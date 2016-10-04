import argparse
#import os
#import sys
#from collections import defaultdict
import itertools
from math import log10
from functools import reduce
import operator
import re
from pdb import set_trace
if False:
    set_trace() # to dodge warnings that pdb isn't being used.
import time
import pickle
import copy
#import sys
#from matplotlib import pyplot as plt
desc = 'Use cross-chamber information to call chamber alleles in SISSOR data'

default_input_file = 'test_subsample.txt'
default_output_file = 'output_calls.txt'

###############################################################################
# PARSE STDIN
###############################################################################

def parseargs():

    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-i', '--input_file', nargs='?', type = str, help='input file, pileup with 3 cells x 24 chambers', default=default_input_file)
    parser.add_argument('-o', '--output_file', nargs='?', type = str, help='file to write output to', default=default_output_file)
    #parser.add_argument('-r', '--region', nargs='?', type = str, help='region to process in format {CHR}:{START}-{END} to process (END excluded)', default=None)

    # default to help option. credit to unutbu: http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
    #if len(sys.argv) <= 1:
    #    parser.print_help()
    #    sys.exit(1)

    args = parser.parse_args()
    return args

###############################################################################
# HELPER FUNCTIONS
###############################################################################

# next function that returns 0 instead of raising StopIteration
# this is convenient for iterating over file 2 at a time
def safenext(iterator):
    try:
        nextval = next(iterator)
        return nextval
    except StopIteration:
        return 0

short_chrom_names = set([str(x) for x in range(1,23)]+['X','Y'])
chrom_names = set(['chr'+str(x) for x in range(1,23)]+['chrX','chrY'])

def format_chrom(chrom_str):
    if chrom_str in short_chrom_names:
        chrom_str = 'chr' + chrom_str
    assert(chrom_str in chrom_names)
    return chrom_str

# product of list
def prod(iterable): # credit to: http://stackoverflow.com/questions/7948291/is-there-a-built-in-product-in-python
    return reduce(operator.mul, iterable, 1)

# sum of two probabilities that are in log form
def addlogs(a,b):
    if a > b:
        return a + log10(1 + pow(10, b - a))
    else:
        return b + log10(1 + pow(10, a - b))

# difference of two probabilities that are in log form
def subtractlogs(a,b):
    if a > b:
        return a + log10(1 - pow(10, b - a))
    else:
        return b + log10(1 - pow(10, a - b))

def remove_multiple_strings(cur_string, replace_list):  # credit to http://stackoverflow.com/questions/30606124/most-efficient-way-to-remove-multiple-substrings-from-string
  for cur_word in replace_list:
    cur_string = cur_string.replace(cur_word, '')
  return cur_string

###############################################################################
# CONSTANTS
###############################################################################

coverage_cut = 3
max_cov = 100
min_nonref = 2
chroms = ['chr{}'.format(i) for i in range(1,23)] + ['chrX','chrY']
cells = ['PGP1_21','PGP1_22','PGP1_A1']
n_chambers = 24
chambers = list(range(1,n_chambers+1))
bases = ['A','T','G','C']
genotypes = list(itertools.combinations_with_replacement(bases,2))
parameters_dir = 'parameters'
n_cells = 3
alleles = ['A','T','G','C']
mixed_alleles = list(set([tuple(sorted(list(set(x)))) for x in itertools.product(alleles,alleles)]))

###############################################################################
# LOAD PARAMETERS
###############################################################################



p_null = pickle.load(open("{}/p_null.p".format(parameters_dir), "rb" )) # estimated p_null, the probability of a strand not being sampled.
ch_priors = pickle.load(open("{}/ch_priors.p".format(parameters_dir), "rb" )) # prior probability of sampling from a chamber
cov_frac_dist = pickle.load(open( "parameters/cov_frac_dist.p", "rb")) # PROBABILITY OF SEEING PARENT 1 ALLELE IN MIXED ALLELE CHAMBER
hom_config_probs = pickle.load(open( "parameters/hom_config_probs.p", "rb")) # STRAND-TO-CHAMBER CONFIGURATIONS
het_config_probs = pickle.load(open( "parameters/het_config_probs.p", "rb"))
genotype_priors = pickle.load(open( "parameters/genotype_priors.p", "rb"))
omega = pickle.load(open( "parameters/omega.p", "rb")) # per-base probability of MDA error
omega_nolog = 10**omega

# INPUT
# G: a tuple that specifies genotype e.g. ('A','T')
# nonzero_chambers: a list containing the indices of chambers 0..23 (or 24, unsampled) where strands may be (have read coverage)
# OUTPUT
# a list of configurations represented by tuples.
# the tuples contain (sl,prob) where sl is a 4-tuple of strand chambers (0..23 or 24) and prob is log probability of config occuring
def singlecell_config(het,nonzero_chambers):

    nz = nonzero_chambers + [n_chambers]
    nzs = set(nonzero_chambers)

    configs = []
    done    = set()
    for c1,c2,c3,c4 in itertools.product(nz,nz,nz,nz):
        if set.intersection({c1,c2,c3,c4},nzs) != nzs:
            continue

        # c1..c4 represent the chamber placements of strands 1..4
        sl = tuple(sorted([c1,c2,c3,c4]))
        if (not het and sl not in hom_config_probs) or (het and sl not in het_config_probs):
            raise Exception('Strand configuration not in config_probs dictionary')

        if sl in done:
            continue

        done.add(sl)

        if het:
            prob = het_config_probs[sl]
        else:
            prob = hom_config_probs[sl]

        configs.append(((c1,c2,c3,c4),prob))

    total = configs[0][1]

    for i in range(1,len(configs)):
        total = addlogs(total,configs[i][1])

    for i in range(len(configs)):
        configs[i] = (configs[i][0], configs[i][1] - total)

    return configs

# INPUT
# G: a tuple that specifies genotype e.g. ('A','T')
# nonzero_chambers: a list containing one list per cell.
# the inner list should contain the indices of chambers 0..23 (or 24, unsampled) where strands are found
# OUTPUT
# a list of configurations represented by tuples.
# the tuple contains one inner tuple per cell.
# these inner tuples contain (sl,prob) where sl is a 4-tuple of strand chambers (0..23 or 24) and prob is log probability of config occuring
def multicell_config(het,nonzero_chambers):

    config_sets = [singlecell_config(het,nz) for nz in nonzero_chambers]

    return itertools.product(*config_sets)

###############################################################################
# LIKELIHOOD CALCULATIONS
###############################################################################

# PRIOR PROBABILITY OF SEEING ALLELES MIXED IN A GIVEN PROPORTION IN A CHAMBER
# accessed as mixed_allele_priors[x][y]
# where x is the number of chambers with reads
# and y is the allele mixture as an alphabetically sorted tuple
def compute_mixed_allele_priors():

    mixed_allele_priors = dict()

    for ref_allele in alleles:

        mixed_allele_priors[ref_allele] = {1:dict(),2:dict(),3:dict(),4:dict()}

        for i in range(1,5):
            for allele in mixed_alleles:
                mixed_allele_priors[ref_allele][i][allele] = -1e10
            for G in genotypes:
                nz = list(range(i))
                het = (G[0] != G[1])
                cfgs = singlecell_config(het,nz)

                for (c1,c2,c3,c4), p in cfgs:
                    # probability of configuration
                    p += genotype_priors[ref_allele][G] - log10(len({c1,c2,c3,c4}))

                    for j in {c1,c2,c3,c4}:

                        alleles_present = []
                        if c1 == j:
                            alleles_present.append(G[0])
                        if c2 == j:
                            alleles_present.append(G[0])
                        if c3 == j:
                            alleles_present.append(G[1])
                        if c4 == j:
                            alleles_present.append(G[1])

                        alleles_present = tuple(sorted(list(set(alleles_present))))

                        if alleles_present not in mixed_allele_priors[ref_allele][i]:
                            mixed_allele_priors[ref_allele][i][alleles_present] = p
                        else:
                            mixed_allele_priors[ref_allele][i][alleles_present] = addlogs(mixed_allele_priors[ref_allele][i][alleles_present], p)

    return mixed_allele_priors

one_allele_cache = dict()
two_allele_cache = dict()
def compute_caches():
    for i in range(0,31):
        for j in range(i+1,31):
            for q in range(33,130):
                for a1_match in [False,True]:
                    for a2_match in [False,True]:
                        qual = 10**((q - 33) * -0.1)
                        frac = i / j if i > 1 else 1e-10
                    
                        x1 = (1.0-qual)*(1.0-omega_nolog) + omega_nolog*qual
                        x2 = omega_nolog*(1.0-qual) + (1-omega_nolog)*qual
                        p1 = x1 if a1_match else x2
                        p2 = x1 if a2_match else x2
                        two_allele_cache[((i / j),qual,a1_match,a2_match)] = log10(frac*p1 + (1-frac)*p2)

                        
    for q in range(33,130):
        qual = 10**((q - 33) * -0.1)

        one_allele_cache[(qual,True)] = log10((1.0-qual)*(1.0-omega_nolog) + omega_nolog*qual)
        one_allele_cache[(qual,False)] = log10(omega_nolog*(1.0-qual) + (1-omega_nolog)*qual)

compute_caches()

def pr_one_chamber_data(alleles_present, base_data, qual_data):

    p = 0
    n = len(base_data)
    assert(n == len(qual_data))

    if len(alleles_present) == 1:

        for base,qual in zip(base_data, qual_data):
            p += one_allele_cache[(qual, (alleles_present[0] == base))]

    elif len(alleles_present) == 2:

        L = len(cov_frac_dist[n])
        for i, p0 in enumerate(cov_frac_dist[n]):
            frac = i / L
            for base,qual in zip(base_data, qual_data):
                p0 += two_allele_cache[(frac, qual, (alleles_present[0] == base), (alleles_present[1] == base))]

            p = addlogs(p, p0)

    return p

def precompute_pr_one_chamber(base_data, qual_data, nonzero_chambers):

    pr_one_ch = dict()
    
    for cell in range(n_cells):
        for chamber in nonzero_chambers[cell]:
            
            # for alleles that are not present at all in the data, the probability of seeing
            # them should be exactly the same.
            present = dict()
            for base in bases:
                present[(base,)] = 0
            for base in base_data[cell][chamber]:
                present[(base,)] = 1
            for allele in mixed_alleles:
                if len(allele) == 2:
                    if (not present[(allele[0],)]) and (not present[(allele[1],)]):
                        present[allele] = 0
                    else:
                        present[allele] = 1

            not_present_val_one_allele = None
            not_present_val_two_allele = None
            
            for allele in mixed_alleles:
                
                if present[allele]:
                    pr_one_ch[(cell, chamber, allele)] = pr_one_chamber_data(allele, base_data[cell][chamber], qual_data[cell][chamber])
                
                elif len(allele) == 1:
                    if not_present_val_one_allele == None:
                        not_present_val_one_allele = pr_one_chamber_data(allele, base_data[cell][chamber], qual_data[cell][chamber])
                    
                    #assert(not_present_val_one_allele == pr_one_chamber_data(allele, base_data[cell][chamber], qual_data[cell][chamber]))
                    pr_one_ch[(cell, chamber, allele)] = not_present_val_one_allele
                elif len(allele) == 2:
                    if not_present_val_two_allele == None:
                        not_present_val_two_allele = pr_one_chamber_data(allele, base_data[cell][chamber], qual_data[cell][chamber])
                    
                    #assert(not_present_val_two_allele == pr_one_chamber_data(allele, base_data[cell][chamber], qual_data[cell][chamber]))
                    pr_one_ch[(cell, chamber, allele)] = not_present_val_two_allele

    return pr_one_ch

# probability that *possible mixed* allele is present to chamber
# allele:    mixed-allele in tuple form that we are testing in chamber
# chamber:   chamber that we are testing for presence of allele in (0-23)
# base_data: list length 24, containing 1 list per chamber. inner list has base pairs called in chamber.
# qual_data: list length 24, containing 1 list per chamber. inner list has q values (p(base call error)) callin in chamber.
# configs:   list of configurations and their probabilities, see earlier

def pr_all_chamber_data(allele, ref_allele, cell, chamber, pr_one_ch, nonzero_chambers, hom_het_configs):

    p_total = None
    p_total = -1e50
    for G in genotypes:

        configs = copy.copy(hom_het_configs[0]) if G[0] == G[1] else copy.copy(hom_het_configs[1])

        for config in configs:
            p = genotype_priors[ref_allele][G]

            (c1,c2,c3,c4), p_cell_cfg = config[cell]

            alleles_present = []

            if c1 == chamber:
                alleles_present.append(G[0])
            if c2 == chamber:
                alleles_present.append(G[0])
            if c3 == chamber:
                alleles_present.append(G[1])
            if c4 == chamber:
                alleles_present.append(G[1])

            alleles_present = tuple(sorted(list(set(alleles_present))))

            if alleles_present != allele:
                continue

            for i in range(0,n_cells):
                if len(nonzero_chambers[i]) == 0:
                    continue

                (c1,c2,c3,c4), p_cell_cfg = config[i]
                p += p_cell_cfg

                for j in nonzero_chambers[i]:

                    alleles_present = []
                    if c1 == j:
                        alleles_present.append(G[0])
                    if c2 == j:
                        alleles_present.append(G[0])
                    if c3 == j:
                        alleles_present.append(G[1])
                    if c4 == j:
                        alleles_present.append(G[1])

                    alleles_present = tuple(sorted(list(set(alleles_present))))

                    p += pr_one_ch[(i,j,alleles_present)]

            p_total = addlogs(p_total,p) if p_total != None else p

    return p_total

# probability that allele is present in chamber
def pr_allele(cell, chamber, pr_one_ch, nonzero_chambers, mixed_allele_priors, ref_allele):

    res = []
    probs = []

    num_nonzero = len(nonzero_chambers[cell])
    hom_het_configs = [multicell_config(False,nonzero_chambers),multicell_config(True,nonzero_chambers)]

    for allele in mixed_alleles:

        probs.append(pr_all_chamber_data(allele, ref_allele, cell, chamber, pr_one_ch, nonzero_chambers, hom_het_configs) + mixed_allele_priors[ref_allele][num_nonzero][allele])

    # denominator for bayes rule posterior calculation
    total = None
    for p in probs:
        total = addlogs(total,p) if total != None else p

    for a,p in zip(mixed_alleles,probs):

        posterior = p - total
        res.append((a,10**posterior))

    return res

###############################################################################
# MAIN FUNCTION AND PARSING
###############################################################################

def call_chamber_alleles(input_file, output_file, SNPs_only=True):

    mixed_allele_priors = compute_mixed_allele_priors()

    with open(input_file,'r') as ipf, open(output_file,'w') as opf:
        for line in ipf:

            el = line.strip().split('\t')

            tags = []
            nonzero_chambers = [[] for i in range(n_cells)]
            nonzero_chamber_count = 0
            base_data = [[[] for i in range(n_chambers)] for j in range(n_cells)]
            qual_data = [[[] for i in range(n_chambers)] for j in range(n_cells)]

            chrom = el[0]
            pos   = int(el[1]) - 1
            ref_base = str.upper(el[2])

            if ref_base == 'N':
                continue
            assert(ref_base in bases)

            total_nonref = 0
            
            for cell_num in range(n_cells):
                for ch_num in range(n_chambers):

                    flat_ix = cell_num * n_chambers + ch_num
                    col_ix = 3 + 4 * flat_ix

                    depth = int(el[col_ix])

                    if depth < coverage_cut or el[col_ix + 1] == '*':
                        continue

                    nonzero_chambers[cell_num].append(ch_num)
                    nonzero_chamber_count += 1

                    bd = str.upper(re.sub(r'\^.|\$|\+[0-9]+[ACGTNacgtn]+|-[0-9]+[ACGTNacgtn]+', '', el[col_ix + 1]))
                    bd = re.sub(r'\.|\,', ref_base, bd)
                    qd = [10**((ord(q) - 33) * -0.1) for q in el[col_ix + 2]]
                    assert(len(bd) == len(qd))
                    assert(len(bd) == depth)

                    paired_bd_qd = [(b,q) for b,q in zip(bd,qd) if b not in ['>','<','*']]
                    if len(paired_bd_qd) < coverage_cut:
                        continue
    
                    bd, qd = zip(*paired_bd_qd)


                    base_data[cell_num][ch_num] = bd[:max_cov]
                    qual_data[cell_num][ch_num] = qd[:max_cov]

                    for b in bd:
                        assert(b in bases)
                        if b != ref_base:
                            total_nonref += 1
            
            if SNPs_only and total_nonref < min_nonref:
                continue
            
            continue
            
            outline_el = ['*']*(n_chambers*n_cells)

            too_many_chambers = False
            for nz in nonzero_chambers:
                if len(nz) > 4:
                    too_many_chambers = True
                    tags.append('TOO_MANY_CHAMBERS')
                    break

            if nonzero_chamber_count > 0 and not too_many_chambers:

                pr_one_ch = precompute_pr_one_chamber(base_data, qual_data, nonzero_chambers)

                for cell in range(0,n_cells):  # for each cell
                    for chamber in nonzero_chambers[cell]:

                        res = pr_allele(cell, chamber, pr_one_ch, nonzero_chambers, mixed_allele_priors, ref_base)

                        out_str = ';'.join(['{}:{}'.format(''.join(a),p) for a,p in res])

                        outline_el[cell*n_chambers + chamber] = out_str

            outline_el = [chrom, str(pos+1), ref_base] + outline_el
            tag_info = ';'.join(tags) if tags != [] else 'N/A'
            outline_el.append(tag_info)
            outline = '\t'.join(outline_el)

            print(outline, file=opf)

if __name__ == '__main__':
    t1 = time.time()
    args = parseargs()
    call_chamber_alleles(args.input_file, args.output_file)
    t2 = time.time()

    print("TOTAL TIME: {} s".format(int(t2-t1)))
