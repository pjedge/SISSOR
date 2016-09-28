import argparse
import numpy as np
#import os
#import sys
#from collections import defaultdict
import itertools
from math import log10, factorial
from functools import reduce
import operator
import re
from pdb import set_trace
if False:
    set_trace() # to dodge warnings that pdb isn't being used.
from copy import deepcopy
from collections import defaultdict
import time
import math
import pickle
from scipy.stats import chisquare
#import sys
#from matplotlib import pyplot as plt
desc = 'Use cross-chamber information to call chamber alleles in SISSOR data'

default_input_dir = 'pileups_subsample'
default_output_dir = 'output_calls.txt'

###############################################################################
# PARSE STDIN
###############################################################################

def parseargs():

    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-i', '--input_dir', nargs='?', type = str, help='input data dir, format {DIR}/PGP1_*/ch*.pileup', default=default_input_dir)
    parser.add_argument('-o', '--output', nargs='?', type = str, help='file to write output to', default=default_output_dir)
    parser.add_argument('-r', '--region', nargs='?', type = str, help='region to process in format {CHR}:{START}-{END} to process (END excluded)', default=None)

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

chroms = ['chr{}'.format(i) for i in range(1,23)] + ['chrX','chrY']
cells = ['PGP1_21','PGP1_22','PGP1_A1']
n_chambers = 24
chambers = list(range(1,n_chambers+1))
bases = ['A','T','G','C']
genotypes = list(itertools.combinations_with_replacement(bases,2))
parameters_dir = 'parameters'
het_snp_rate = 0.0005
hom_snp_rate = 0.001
n_cells = 3
alleles = ['A','T','G','C']

###############################################################################
# ESTIMATE VARIOUS PARAMETERS
###############################################################################

grams_DNA_before_MDA  = 0.25e-12  # 0.25 pg
grams_DNA_after_MDA = 6e-9        # 6 ng
GAIN = grams_DNA_after_MDA / grams_DNA_before_MDA
omega = log10(3.2e-6 / 2 * math.log2(GAIN)) # probability of MDA error.
# formula from this paper: http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0105585

mixed_alleles = list(set([tuple(sorted(list(set(x)))) for x in itertools.product(alleles,alleles)]))
#azip = list(zip(alleles, allele_priors))
transition = {'A':'G','G':'A','T':'C','C':'T'}

total_sampled_cell_positions = pickle.load(open("{}/total_sampled_cell_positions.p".format(parameters_dir), "rb"))
# estimate proportions of strands in chambers
chamber_position_counts = pickle.load(open("{}/chamber_position_counts.p".format(parameters_dir), "rb" ))
total_position_counts = sum([chamber_position_counts[chamber] for chamber in range(0,24)])
chamber_proportions = [chamber_position_counts[chamber]/total_position_counts for chamber in range(0,24)]

chamber_proportions = [x if x != 0 else 1e-10 for x in chamber_proportions]
                       
# estimate p_null, the probability of a strand not being sampled.
# we'll iterate over many possible values and choose the one that minimizes p-value for chi-square goodness of fit
strand_coverage_counts = pickle.load(open("{}/strand_coverage_counts.p".format(parameters_dir), "rb" ))
strand_coverage_proportions = []
tot = sum(list(strand_coverage_counts.values()))
for i in range(0,5):
    strand_coverage_proportions.append(strand_coverage_counts[i] / tot)
print("STRAND COVERAGE PROPORTIONS")
print(strand_coverage_proportions)
p_null = None #0.75
ch_priors = None #[log10((1-p_null)/24)]*24 + [log10(p_null)]

max_likelihood = -float("Inf")

p_null_samples = 3
print("EXPECTED PROPORTIONS")
for i in range(1,p_null_samples):
    
    putative_p_null = i/p_null_samples
    putative_ch_priors = [log10((1-putative_p_null) * x) for x in chamber_proportions] + [log10(putative_p_null)]   # prior probability of sampling a given chamber

    # given these proportions, compute the expected proportions of positions
    # with a given number of strands present.
    
    poss_ch = list(range(0,25)) # the last element is "no chamber"
    coverage_dict = dict()
    for j in range(0,5):
        coverage_dict[j] = -1e10
    seen_cfgs = set()
    plist = []
    for c1,c2,c3,c4 in itertools.product(poss_ch,poss_ch,poss_ch,poss_ch):
        # because of symmetry we ignore cases where we can flip the strands for the same result
        # e.g. we ignore cases where strand 1 is in a later chamber than strand 2 and say that
        # the cases where it is in an earlier chamber are twice as likely. Same for strands 3,4
        if c1 > c2:
            temp = c1
            c1 = c2
            c2 = temp
    
        if c3 > c4:
            temp = c3
            c3 = c4
            c4 = temp  
    
        if (c1,c2,c3,c4) in seen_cfgs:
            continue
        
        if c1 < c2 and c3 < c4:
            perm = log10(4)
        elif (c1 < c2 and c3 == c4) or (c1 == c2 and c3 < c4):
            perm = log10(2)
        elif c1 == c2 and c3 == c4:
            perm = log10(1)
        else:
            continue
        seen_cfgs.add((c1,c2,c3,c4))
        strands_present = 4
        for c in [c1,c2,c3,c4]:
            if c == 24:
                strands_present -= 1
                
        cfg_prob = perm+putative_ch_priors[c1]+putative_ch_priors[c2]+putative_ch_priors[c3]+putative_ch_priors[c4]
        coverage_dict[strands_present] = addlogs(coverage_dict[strands_present],cfg_prob) # probability of configuration

    likelihood = 0
    for j in range(0,5):
        # the likelihood of the observed data given this value for p_null
        # is the product of likelihoods of each observed coverage value
        likelihood += strand_coverage_counts[j] * coverage_dict[j]


    expected_proportions = [10**x for x in coverage_dict.values()]

    print("i={} ".format(i),end='')
    print(expected_proportions,end='')
    print(likelihood)
    
    if likelihood > max_likelihood:
        max_likelihood = likelihood
        p_null = putative_p_null
        ch_priors = putative_ch_priors

print("P(strand not sampled) = {}".format(p_null))
# estimate prior probability of genotypes using strategy described here:
# http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2694485/
# "prior probability of each genotype"

genotype_priors = dict()
for allele in alleles:
    
    # priors on haploid alleles
    p_a = dict()
    p_a[allele] = 1 - hom_snp_rate
    p_a[transition[allele]] = hom_snp_rate / 6 * 4
    for transversion in alleles:
        if transversion in p_a:
            continue
        p_a[transversion] =  hom_snp_rate / 6
    
    genotype_priors[allele] = dict()
    for G in genotypes:
        g1,g2 = G
        # probability of homozygous reference is the probability of neither het or hom SNP
        if g1 == g2 and g1 == allele:
            genotype_priors[allele][G] = 1.0 - het_snp_rate - hom_snp_rate
        elif g1 == g2 and g1 != allele:
            # transitions are 4 times as likely as transversions
            if g1 == transition[allele]:
                genotype_priors[allele][G] = het_snp_rate / 6 * 4 
            else:
                genotype_priors[allele][G] = het_snp_rate / 6
        else: # else it's the product of the haploid priors
                genotype_priors[allele][G] = p_a[g1] * p_a[g2]

    # convert to log
    for G in genotype_priors[allele].keys():
        genotype_priors[allele][G] = log10(genotype_priors[allele][G])


# PROBABILITY OF SEEING PARENT 1 ALLELE IN MIXED ALLELE CHAMBER

#P_parent1_lst = [(0.001,log10(0.495)),(0.5,log10(0.01)),(0.999,log10(0.495))]
cov_frac_dist_raw = pickle.load(open( "parameters/cov_frac_dist.p", "rb"))
cov_frac_dist = defaultdict(list)

lim = 30
for i in range(1,lim+1):
    cov_frac_dist[i] = cov_frac_dist_raw[i]

def chunks(l, n): # credit to http://stackoverflow.com/questions/312443/how-do-you-split-a-list-into-evenly-sized-chunks
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]

# given a list with length N*binsize, sum consectutive binsize-sized chunks into a length N list
def condense(l, binsize):
    return [sum(x) for x in chunks(l,binsize)]
    

last_lst = []
for i in range(lim,2000,lim):
    binsize = int(i / lim)
    lst = condense(cov_frac_dist_raw[i], binsize)
    if lst == []:
        lst = last_lst
    for j in range(i,i+lim):
        cov_frac_dist[j] = lst
    last_lst = lst

for i in range(lim,2000):
    assert(len(cov_frac_dist[i]) == lim)

for i in range(1,200,10):
    plt.figure()
    y = cov_frac_dist[i]
    x = list(range(len(cov_frac_dist[i])))
    plt.plot(x,y)
    plt.title(str(i))

assert(False)
        
###############################################################################
# STRAND-TO-CHAMBER CONFIGURATIONS
###############################################################################

# helper data structures that pre-compute probabilities of strand configurations
poss_ch = list(range(0,25)) # the last element is "no chamber"
het_config_probs = dict()

for c1,c2,c3,c4 in itertools.product(poss_ch,poss_ch,poss_ch,poss_ch):
    # because of symmetry we ignore cases where we can flip the strands for the same result
    # e.g. we ignore cases where strand 1 is in a later chamber than strand 2 and say that
    # the cases where it is in an earlier chamber are twice as likely. Same for strands 3,4
    if c1 > c2:
        temp = c1
        c1 = c2
        c2 = temp

    if c3 > c4:
        temp = c3
        c3 = c4
        c4 = temp  

    if (c1,c2,c3,c4) in het_config_probs:
        continue
    
    if c1 < c2 and c3 < c4:
        perm = log10(4)
    elif (c1 < c2 and c3 == c4) or (c1 == c2 and c3 < c4):
        perm = log10(2)
    elif c1 == c2 and c3 == c4:
        perm = log10(1)
    else:
        continue

    het_config_probs[(c1,c2,c3,c4)] = perm+ch_priors[c1]+ch_priors[c2]+ch_priors[c3]+ch_priors[c4] # probability of configuration
                    
hom_config_probs = dict()

for c1,c2,c3,c4 in itertools.product(poss_ch,poss_ch,poss_ch,poss_ch):
    
    lst  = sorted([c1,c2,c3,c4])
    if tuple(lst) in hom_config_probs:
        continue
    uniq = list(set(lst))
    counts  = [lst.count(u) for u in uniq]
    perm = log10(factorial(len(lst))/prod([factorial(c) for c in counts]))
    
    hom_config_probs[tuple(lst)] = perm+ch_priors[c1]+ch_priors[c2]+ch_priors[c3]+ch_priors[c4] # probability of configuration 

                     
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

def pr_one_chamber_data(alleles_present, base_data, qual_data):

    p = 0
    n = len(base_data)
    assert(n == len(qual_data))  
    
    if len(alleles_present) == 1:
        
        for base,qual in zip(base_data, qual_data):
            if alleles_present[0] == base:
                p += addlogs(subtractlogs(0,qual)+subtractlogs(0,omega), qual+omega)
            else:
                p += addlogs(omega+subtractlogs(0,qual), qual+subtractlogs(0,omega))   
                
    elif len(alleles_present) == 2:
        

        a1 = alleles[0]
        a2 = alleles[1]
        cov = len(base_data)
        for i, p0 in enumerate(cov_frac_dist[cov]):
            frac = i / len(cov_frac_dist[cov]) if i > 1 else 1e-10
            for base,qual in zip(base_data, qual_data):
                
                x1 = addlogs(subtractlogs(0,qual)+subtractlogs(0,omega), qual+omega)
                x2 = addlogs(omega+subtractlogs(0,qual), qual+subtractlogs(0,omega))
                
                p1 = x1 if a1 == base else x2
                p2 = x1 if a2 == base else x2
                
                p1 += log10(frac)
                p2 += log10(1-frac)
                
                p0 += addlogs(p1, p2)
                        
            p = addlogs(p, p0)
            
    return p

def precompute_pr_one_chamber(base_data, qual_data, nonzero_chambers):
    
    pr_one_ch = dict()
    
    for cell in range(n_cells):
        for chamber in nonzero_chambers[cell]:
            for allele in mixed_alleles:
                
                pr_one_ch[(cell, chamber, allele)] = pr_one_chamber_data(allele, base_data[cell][chamber], qual_data[cell][chamber])
        
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

        configs = deepcopy(hom_het_configs[0]) if G[0] == G[1] else deepcopy(hom_het_configs[1])

        for config in configs:
            p = genotype_priors[G]
                
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
    
def call_chamber_alleles(input_dir, output, region):
    
    if region != None:
        region_chrom, bounds = region.split(':')
        region_start, region_end = bounds.split('-')
        region_start = int(region_start)
        region_end = int(region_end)

    mixed_allele_priors = compute_mixed_allele_priors()
    # open all 72 files (3 cells x 24 chambers) into a list of handles
    # input files will be read and processed in parallel for equivalent indices
    # output files will contain the same information but will be corrected based on information from other chambers
    
    input_files = []

    for cell in cells:
            
        for chamber in chambers:
            
            infile_name = '{}/{}/ch{}.pileup'.format(input_dir,cell,chamber)

            input_files.append(open(infile_name,'r'))
        
    # step over files in parallel
    chrom = chroms.pop(0)
    pos = 0
    
    lines = []
    
    for file in input_files:
        
        lines.append(safenext(file))
    
        
    N = len(input_files)

    # print header lines to output files 
    for ix in range(N):
        while 1:
            
            if lines[ix] and lines[ix][0] == '#':
                lines[ix] = safenext(input_files[ix])
            else:
                break
            
    seen_region_chrom = False
    
    coverage_cut = 3
    PROCESSED = 0
    # until all files have reached stopiteration
    with open(output,'w') as opf:
        while sum([not l for l in lines]) < N:

            # "element list"
            # el_lst[x] contains a list with the current line elements for file x (unique cell and chamber)
            el_lst = []
            for l in lines:
                if not l:
                    el_lst.append(0)
                else:
                    el_lst.append(l.strip().split('\t'))
    
            # fix inconsistent chromosome names
            for i in range(len(el_lst)):
                if el_lst[i]:
                    el_lst[i][0] = format_chrom(el_lst[i][0])
                
            on_chrom = [el and el[0] == chrom for el in el_lst]
            on_chrom_pos = [el and el[0] == chrom and int(el[1])-1 == pos for el in el_lst]
    
            # if there are no more files on our chromosome, we move onto the next chromosome
            if sum(on_chrom) == 0:
                if chroms == []:
                    break
                else:
                    chrom = chroms.pop(0)
                    pos = 0
                    continue
            
            # now we know that some subset of the N files are processing our chromosome
            # if none of the files are on our genomic index we also need to iterate forward
            # until some are
            
            if sum(on_chrom_pos) == 0:
                pos = float('Inf')  # positions for each file
                for el in el_lst:
                    if el and el[0] == chrom and int(el[1])-1 < pos:
                        pos = int(el[1])-1  # select minimum position index of all currently considered lines
            
            on_chrom_pos = [el and el[0] == chrom and int(el[1])-1 == pos for el in el_lst]
            assert(sum(on_chrom_pos) > 0)
            
            # now we have a valid chromosome, being considered by some file handles
            # we also have a valid position index, of which some file handles are on it.
            
            # process each infile with a 1 in on_chrom_pos
            # then iterate each of those files
    
            # do something with line element that is on the current genomic index
            if region != None:
                if chrom == region_chrom:
                    seen_region_chrom = True
                if (seen_region_chrom and chrom != region_chrom) or (chrom == region_chrom and pos >= region_end): # past region, we're done
                    for handle in input_files:
                        handle.close()
                    return
                    
                elif chrom != region_chrom or pos < region_start: # before region, skip ahead
                    # iterate to the next lines for each line on the current index        
                    for i in np.where(on_chrom_pos)[0]:
            
                        lines[i] = safenext(input_files[i])
                    continue
                                    
            print("PROCESSED {}".format(PROCESSED))
            tags = []
            nonzero_chambers = [[] for i in range(n_cells)]
            nonzero_chambers_flatix = []
            base_data = [[[] for i in range(n_chambers)] for j in range(n_cells)]
            qual_data = [[[] for i in range(n_chambers)] for j in range(n_cells)]
            ref_base = None
            
            for i in np.where(on_chrom_pos)[0]:
                
                cell_num = int(i / n_chambers)
                ch_num   = i % 24
                
                el = el_lst[i]
                if not el:
                    continue
                
                if int(el[3]) < coverage_cut:
                    continue
                if ref_base == None:
                    ref_base = str.upper(el[2])
                else:
                    try:
                        assert(ref_base == str.upper(el[2]))
                    except:
                        set_trace()
                if ref_base == 'N':
                    continue
                
                assert(ref_base in bases)
                
                nonzero_chambers_flatix.append(i)
                nonzero_chambers[cell_num].append(ch_num)        
        
                bd = str.upper(re.sub(r'\^.|\$', '', el[4]))
                
                bd = re.sub(r'\.|\,', ref_base, bd)
                qd = [((ord(q) - 33) * -0.1) for q in el[5]]

                bd, qd = zip(*[(b,q) for b,q in zip(bd,qd) if b not in ['>','<','*']])

                
                assert(len(bd) == len(qd))

                
                base_data[cell_num][ch_num] = bd
                qual_data[cell_num][ch_num] = qd
            
            outline_el = ['-']*(n_chambers*n_cells)
            
            too_many_chambers = False
            
            for nz in nonzero_chambers:
                if len(nz) > 4:
                    too_many_chambers = True
                    tags.append('TOO_MANY_CHAMBERS')
                    break
               
            if nonzero_chambers_flatix != [] and not too_many_chambers:
                
                pr_one_ch = precompute_pr_one_chamber(base_data, qual_data, nonzero_chambers)
                
                for cell in range(0,n_cells):  # for each cell
                    for chamber in nonzero_chambers[cell]:
                        
                        res = pr_allele(cell, chamber, pr_one_ch, nonzero_chambers, mixed_allele_priors, ref_base)
                        
                        out_str = ';'.join(['{}:{}'.format(''.join(a),p) for a,p in res])
                        
                        outline_el[cell*n_chambers + chamber] = out_str

            
            outline_el = [chrom, str(pos+1)] + outline_el
            tag_info = ';'.join(tags) if tags != [] else 'N/A'
            outline_el.append(tag_info)
            outline = '\t'.join(outline_el)

            print(outline, file=opf)
            
            PROCESSED += 1
            if PROCESSED > 500:
                return
            # iterate to the next lines for each line on the current index        
            for i in np.where(on_chrom_pos)[0]:
    
                lines[i] = safenext(input_files[i])

    # close all chambers
            
    for handle in input_files:
        handle.close()
        

if __name__ == '__main__':
    t1 = time.time()
    args = parseargs()
    call_chamber_alleles(args.input_dir, args.output, args.region)
    t2 = time.time()
    
    print("TOTAL TIME: {} s".format(int(t2-t1)))
        
