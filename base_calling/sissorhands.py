import argparse
import itertools

from math import log10
from pdb import set_trace
if False:
    set_trace() # to dodge warnings that pdb isn't being used.
import time
import re
import pickle
from collections import defaultdict

#from numpy import logaddexp as addlogs

def addlogs(a,b):
    if a > b:
        return a + log10(1 + pow(10, b - a))
    else:
        return b + log10(1 + pow(10, a - b))


from collections import defaultdict
desc = 'Use cross-chamber information to call chamber alleles in SISSOR data'

###############################################################################
# CONSTANTS
###############################################################################

cells = ['PGP1_21','PGP1_22','PGP1_A1']
coverage_cut = 2 # less than this many reads in a chamber are ignored
max_cov = 100
min_nonref = 3
chroms = ['chr{}'.format(i) for i in range(1,23)] + ['chrX','chrY']
n_chambers = 24
chambers = list(range(1,n_chambers+1))
bases = ['A','T','G','C']
genotypes = list(itertools.combinations_with_replacement(bases,2))
parameters_dir = 'parameters'
n_cells = 3
numbers = '0123456789'
base_letters = 'ACGTacgt'
mixed_alleles = [('A',),('T',),('G',),('C',), ('A', 'T'),('A', 'G'),('A', 'C'),('T', 'G'),('T', 'C'),('G', 'C')]
tinylog = -1e5
INDEL_COUNT_LIMIT = 3

chr_num = dict()
for i,chrom in enumerate(chroms):
    chr_num[chrom] = i

###############################################################################
# DEFAULTS...
###############################################################################

default_input_file = 'pileup_example.txt'#'test_subsample.txt'
default_output_file =  'output.txt'#'test_subsample.output'
BASIC = True

default_fragment_boundaries = None
#default_fragment_boundaries = []
#for cell in cells:
#    for chamber in chambers:
#        filename = 'fragment_boundary_beds/{}/ch{}.bed'.format(cell,chamber)
#        default_fragment_boundaries.append(filename)


###############################################################################
# PARSE STDIN
###############################################################################


def parseargs():

    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-i', '--input_file', nargs='?', type = str, help='input file, pileup with n cells x 24 chambers', default=default_input_file)
    parser.add_argument('-o', '--output_file', nargs='?', type = str, help='file to write output to', default=default_output_file)
    parser.add_argument('-b', '--fragment_boundaries', nargs='+', type = str, help='', default=None)
    parser.add_argument('-n', '--num_cells', nargs='?', type = int, help='number of cells being analyzed, default: 3', default=3)

    args = parser.parse_args()
    return args


###############################################################################
# PARAMETERS
###############################################################################

## INITIALIZED AS None TO ENSURE GLOBAL ACCESS
## READ THEM IN LATER USING initialize_parameters() FUNCTION

p_null = pickle.load(open("{}/p_null.p".format(parameters_dir), "rb" )) # estimated p_null, the probability of a strand not being sampled.
ch_priors = pickle.load(open("{}/ch_priors.p".format(parameters_dir), "rb" )) # prior probability of sampling from a chamber
cov_frac_dist = pickle.load(open( "parameters/cov_frac_dist.p", "rb")) # PROBABILITY OF SEEING PARENT 1 ALLELE IN MIXED ALLELE CHAMBER
hom_config_probs = pickle.load(open( "parameters/hom_config_probs.p", "rb")) # STRAND-TO-CHAMBER CONFIGURATIONS
het_config_probs = pickle.load(open( "parameters/het_config_probs.p", "rb"))
genotype_priors = pickle.load(open( "parameters/genotype_priors.p", "rb"))
omega = pickle.load(open( "parameters/omega.p", "rb")) # per-base probability of MDA error
omega_nolog = 10**omega

MDA_ERR = 1e-5 # the probability of consensus error due to MDA
MDA_COR = log10(1 - MDA_ERR)
MDA_ERR = log10(MDA_ERR) - log10(3) # the log probability of an MDA error, divided by 3

INTRACELL_HAP_PENALTY = log10(1e-10) # intra-cell penalty for assigning chambers of different haps to same haps
INTERCELL_HAP_PENALTY = log10(1e-5) # inter-cell penalty for assigning chambers of different haps to same haps

# INPUT
# G: a tuple that specifies genotype e.g. ('A','T')
# nonzero_chambers: a list containing the indices of chambers 0..23 (or 24, unsampled) where strands may be (have read coverage)
# OUTPUT
# a list of configurations represented by tuples.
# the tuples contain (sl,prob) where sl is a 4-tuple of strand chambers (0..23 or 24) and prob is log probability of config occuring
def singlecell_config(het,nonzero_chambers):

    nz = nonzero_chambers + [n_chambers]
    nzs = set(nonzero_chambers)

    configs = defaultdict(lambda: tinylog)
    done    = set()
    for c1,c2,c3,c4 in itertools.product(nz,nz,nz,nz):
        if set.intersection({c1,c2,c3,c4},nzs) != nzs:
            continue

        # c1..c4 represent the chamber placements of strands 1..4
        if het:
            if c1 > c2:
                temp = c1
                c1 = c2
                c2 = temp

            if c3 > c4:
                temp = c3
                c3 = c4
                c4 = temp
            sl = (c1,c2,c3,c4)
        else:
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

        if (c1 == c2 and c1 != 24):
            c2 = 24
        if (c3 == c4 and c3 != 24):
            c4 = 24

        configs[(c1,c2,c3,c4)] = addlogs(configs[(c1,c2,c3,c4)], prob)

        #configs.append(((c1,c2,c3,c4),prob))

    configs = list(configs.items())

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

    return list(itertools.product(*config_sets))

###############################################################################
# LIKELIHOOD CALCULATIONS
###############################################################################

# PRIOR PROBABILITY OF SEEING ALLELES MIXED IN A GIVEN PROPORTION IN A CHAMBER
# accessed as mixed_allele_priors[x][y]
# where x is the number of chambers with reads
# and y is the allele mixture as an alphabetically sorted tuple
def compute_mixed_allele_priors():

    mixed_allele_priors = dict()

    for ref_allele in bases:

        mixed_allele_priors[ref_allele] = {1:dict(),2:dict(),3:dict(),4:dict()}

        for i in range(1,5):
            for allele in mixed_alleles:
                mixed_allele_priors[ref_allele][i][allele] = tinylog
            for G in genotypes:
                nz = list(range(i))
                het = (G[0] != G[1])
                cfgs = singlecell_config(het,nz)

                for (c1,c2,c3,c4), p in cfgs:
                    # probability of configuration
                    p += genotype_priors[ref_allele][G] - log10(len({c1,c2,c3,c4}))

                    for j in {c1,c2,c3,c4}:

                        G1 = G[0]
                        G2 = G[1]
                        G1_present = (c1 == j or c2 == j)
                        G2_present = (c3 == j or c4 == j)

                        if G1_present and G2_present:
                            if G1 != G2:
                                alleles_present = G
                            else:
                                alleles_present = (G1,)
                        elif G1_present and not G2_present:
                            alleles_present = (G1,)
                        elif G2_present and not G1_present:
                            alleles_present = (G2,)

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

def pr_one_chamber_data(alleles_present, base_data, qual_data, fast_mode):


    n = len(base_data)
    assert(n == len(qual_data))
    p = 0

    if len(alleles_present) == 1:

        for base,qual in zip(base_data, qual_data):

            p += one_allele_cache[(qual, (alleles_present[0] == base))]

    elif len(alleles_present) == 2:

        if fast_mode:
            for base,qual in zip(base_data, qual_data):
                p += two_allele_cache[(0.5, qual, (alleles_present[0] == base), (alleles_present[1] == base))]
        else:
            p = -1e10
            L = len(cov_frac_dist[n])
            for i, p0 in enumerate(cov_frac_dist[n]):
                p0 = log10(p0) if p0 > 0 else tinylog
                frac = i / L
                for base,qual in zip(base_data, qual_data):
                    p0 += two_allele_cache[(frac, qual, (alleles_present[0] == base), (alleles_present[1] == base))]

                p = addlogs(p, p0)

    return p

def precompute_pr_one_chamber(base_data, qual_data, nonzero_chambers, fast_mode):

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
                    pr_one_ch[(cell, chamber, allele)] = pr_one_chamber_data(allele, base_data[cell][chamber], qual_data[cell][chamber], fast_mode)

                elif len(allele) == 1:
                    if not_present_val_one_allele == None:
                        not_present_val_one_allele = pr_one_chamber_data(allele, base_data[cell][chamber], qual_data[cell][chamber], fast_mode)

                    #assert(not_present_val_one_allele == pr_one_chamber_data(allele, base_data[cell][chamber], qual_data[cell][chamber]))
                    pr_one_ch[(cell, chamber, allele)] = not_present_val_one_allele
                elif len(allele) == 2:
                    if not_present_val_two_allele == None:
                        not_present_val_two_allele = pr_one_chamber_data(allele, base_data[cell][chamber], qual_data[cell][chamber], fast_mode)

                    #assert(not_present_val_two_allele == pr_one_chamber_data(allele, base_data[cell][chamber], qual_data[cell][chamber]))
                    pr_one_ch[(cell, chamber, allele)] = not_present_val_two_allele

    return pr_one_ch

# probability that *possible mixed* allele is present to chamber
# allele:    mixed-allele in tuple form that we are testing in chamber
# chamber:   chamber that we are testing for presence of allele in (0-23)
# base_data: list length 24, containing 1 list per chamber. inner list has base pairs called in chamber.
# qual_data: list length 24, containing 1 list per chamber. inner list has q values (p(base call error)) callin in chamber.
# configs:   list of configurations and their probabilities, see earlier

def pr_genotype(pr_one_ch, nonzero_chambers, ref_allele, condensed_genotype_set):

    probs = dict()
    outline_el = ['*'] + ['*'] * (n_cells * n_chambers)

    hom_het_configs = [multicell_config(False,nonzero_chambers),multicell_config(True,nonzero_chambers)]
    p_assignment = defaultdict(lambda: tinylog)
    genotype_set = genotypes if condensed_genotype_set == None else condensed_genotype_set

    for G in genotype_set:

        configs = hom_het_configs[0] if G[0] == G[1] else hom_het_configs[1]
        p_total = tinylog

        for config in configs:

            p = 0
            assignments = set()

            for i in range(0,n_cells):
                if len(nonzero_chambers[i]) == 0:
                    continue

                (c1,c2,c3,c4), p_cell_cfg = config[i]
                p += p_cell_cfg

                for j in nonzero_chambers[i]:

                    G1 = G[0]
                    G2 = G[1]
                    G1_present = (c1 == j or c2 == j)
                    G2_present = (c3 == j or c4 == j)

                    if G1_present and G2_present:
                        if G1 != G2:
                            alleles_present = G
                        else:
                            alleles_present = (G1,)
                    elif G1_present and not G2_present:
                        alleles_present = (G1,)
                    elif G2_present and not G1_present:
                        alleles_present = (G2,)
                    else:
                        print("ERROR")
                        exit(1)

                    if len(alleles_present) == 1:

                        # we must account for the fact that the entire consensus may
                        # be wrong due to MDA error
                        p0 = tinylog
                        for base in bases:
                            if base == alleles_present[0]:
                                p1 = pr_one_ch[(i,j,alleles_present)] + MDA_COR
                            else:
                                p1 = pr_one_ch[(i,j,(base,))] + MDA_ERR

                            p0 = addlogs(p0, p1)

                        p += p0
                    else:

                        allele1 = alleles_present[0]
                        allele2 = alleles_present[1]

                        p0 = tinylog
                        #prob_tot = tinylog
                        for base1, base2 in genotypes:

                            x1 = MDA_COR if base1 == allele1 else MDA_ERR
                            x2 = MDA_COR if base2 == allele2 else MDA_ERR

                            if (base1,base2) == alleles_present:
                                p1 = pr_one_ch[(i,j,(base1,base2))] + MDA_COR + MDA_COR
                                p2 = pr_one_ch[(i,j,(base1,base2))] + MDA_ERR + MDA_ERR
                                p1 = addlogs(p1,p2)
                                #prob_tot = addlogs(prob_tot, x1+x2)
                                #prob_tot = addlogs(prob_tot, MDA_ERR)
                            elif base1 == base2:
                                p1 = pr_one_ch[(i,j,(base1,))] + x1 + x2
                                #prob_tot = addlogs(prob_tot, x1+x2)
                            else:
                                p1 = pr_one_ch[(i,j,(base1,base2))] + x1 + x2 + log10(2)
                                #prob_tot = addlogs(prob_tot, x1+x2+log10(2))

                            p0 = addlogs(p0, p1)

                        #assert(abs(1 - 10**prob_tot) < 0.001)

                        p += p0

                    assignments.add((i,j,alleles_present))


            ########### NEW, 10/13/2016 #######################################
            # IF WE KNOW TWO CHAMBERS c_a, c_b BELONG TO SAME HAPLOTYPE
            # incur a penalty for assigning those to non-complementary strands
            # IF WE KNOW TWO CHAMBERS c_a, c_b BELONG TO DIFFERENT HAPLOTYPE
            # incur a penalty for assigning those to complementary strands
            
            # this penalty should be large within-cell and slightly smaller between cell.
            # also, a chamber should get a special tag denoting its status:
            # SAME CELL SAME HAP
            # SAME CELL DIFFERENT HAP
            # DIFFERENT CELL SAME HAP
            # DIFFERENT CELL DIFFERENT HAP
            same_haplotype = defaultdict(lambda: None) ## need to build this
            same_haplotype[(0,1),(0,2)] = True
            same_haplotype[(0,2),(0,1)] = True
            done = set()
            
            # TODO: strand-matched chambers with the same data get different posteriors, figure out why            
            
            # ADDRESS THIS ISSUE, bundled/multicounted configurations
            
            #if c1 < c2 and c3 < c4:
            #    perm = log10(4)
            #elif (c1 < c2 and c3 == c4) or (c1 == c2 and c3 < c4):
            #    perm = log10(2)
            #elif c1 == c2 and c3 == c4:
            #    perm = log10(1)    
            
            
            for CELL1 in range(0,n_cells):
                cell1_chambers, p_cell1_cfg = config[CELL1]
                
                for CELL2 in range(0,n_cells):
                    cell2_chambers, p_cell2_cfg = config[CELL2]

                    for x1,y1 in enumerate(cell1_chambers):
                        for x2, y2 in enumerate(cell2_chambers):
                            
                            if x1 == n_chambers or x2 == n_chambers: # assignments to null chamber (not sampled)
                                continue
                            
                            if ((CELL1,x1),(CELL2,x2)) in done or ((CELL2,x2),(CELL1,x1)) in done:
                                continue
                            
                            if CELL1 == CELL2 and y1 == y2:
                                continue

                            if same_haplotype[(CELL1,x1),(CELL2,x2)] == None:
                                continue # no haplotype information between these chambers   
                            
                            if ((same_haplotype[(CELL1,x1),(CELL2,x2)]
                             and not ((y1 < 2 and y2 < 2) or (y1 >= 2 and y2 >= 2)))
                             or ((not same_haplotype[(CELL1,x1),(CELL2,x2)])
                                and ((y1 < 2 and y2 < 2) or (y1 >= 2 and y2 >= 2)))):   # 0,1 -> parent1 , 2,3 -> parent2
                                
                                # same haplotype but this configuration assigns them to different haplotypes
                                # or diff haplotype but this configuration assigns them to same haplotypes

                                if CELL1 == CELL2:
                                    p += INTRACELL_HAP_PENALTY
                                else: 
                                    p += INTERCELL_HAP_PENALTY
                                    
                            done.add(((CELL1,x1),(CELL2,x2)))
                
            ###################################################################


            p_total = addlogs(p_total, p) # update genotype likelihood sum

            for assignment in assignments:
                p_assignment[assignment] = addlogs(p_assignment[assignment], p + genotype_priors[ref_allele][G])

        #if p_total == tinylog:
        #    set_trace()

        probs[G] = p_total + genotype_priors[ref_allele][G]

    # denominator for bayes rule posterior calculation
    total = tinylog
    for G,p in probs.items():
        total = addlogs(total,p)

    res = []
    SNP = False
    max_posterior = -1
    max_G = ('N','N')
    for G in genotypes:
        if G in probs:
            posterior = 10**(probs[G] - total)
            res.append((G,posterior))
            if posterior > max_posterior:
                max_posterior = posterior
                max_G = G
        else:
            res.append((G,0.0))

    if max_G != (ref_allele,ref_allele):
        SNP = True

    out_str = ';'.join(['{}:{}'.format(''.join(g),p) for g,p in res])
    outline_el[0] = out_str

    for cell in range(0,n_cells):  # for each cell
        for chamber in nonzero_chambers[cell]:

            #num_nonzero = len(nonzero_chambers[cell])
            #pr_all_chamber_data(allele, ref_allele, cell, chamber, pr_one_ch, nonzero_chambers, hom_het_configs)
            probs = []
            for allele in mixed_alleles:
                prob = p_assignment[(cell,chamber,allele)] #+ mixed_allele_priors[ref_allele][num_nonzero][allele]
                probs.append((allele,prob))

            # denominator for bayes rule posterior calculation
            total = tinylog
            for a,p in probs:
                total = addlogs(total,p)

            res = []
            max_allele = ('N',)
            max_posterior = -1
            for a,p in probs:
                posterior = 10**(p - total)
                res.append((a,posterior))

                if posterior > max_posterior:
                    max_posterior = posterior
                    max_allele = a

            if max_allele != (ref_allele,):
                SNP = True

            out_str = ';'.join(['{}:{}'.format(''.join(a),p) for a,p in res])

            outline_el[1 + cell*n_chambers + chamber] = out_str

    return outline_el, SNP


def pr_chambers_BASIC(pr_one_ch, nonzero_chambers, mixed_allele_priors, ref_allele, condensed_genotype_set):

    probs = dict()
    outline_el = ['*'] + ['*'] * (n_cells * n_chambers)

    for cell in range(0,n_cells):  # for each cell
        for chamber in nonzero_chambers[cell]:

            num_nonzero = len(nonzero_chambers[cell])
            #pr_all_chamber_data(allele, ref_allele, cell, chamber, pr_one_ch, nonzero_chambers, hom_het_configs)
            probs = []
            for allele in mixed_alleles:
                prob = pr_one_ch[(cell,chamber,allele)] + mixed_allele_priors[ref_allele][num_nonzero][allele]
                probs.append((allele,prob))

            # denominator for bayes rule posterior calculation
            total = tinylog
            for a,p in probs:
                total = addlogs(total,p)

            res = []
            max_allele = ('N',)
            max_posterior = -1
            for a,p in probs:
                posterior = 10**(p - total)
                res.append((a,posterior))

                if posterior > max_posterior:
                    max_posterior = posterior
                    max_allele = a

            if max_allele != (ref_allele,):
                SNP = True

            out_str = ';'.join(['{}:{}'.format(''.join(a),p) for a,p in res])

            outline_el[1 + cell*n_chambers + chamber] = out_str

    return outline_el, SNP

###############################################################################
# MAIN FUNCTION AND PARSING
###############################################################################

caret_money = re.compile(r'\^.|\$')
comma_or_period = re.compile(r'\.|\,')
plus_or_minus = re.compile(r'\+|-')

def parse_mpileup_base_qual(raw_bd, raw_qd, ref_base):

    bd = [] # base data
    qd = [] # qual data

    indel_count = len(re.findall(caret_money, raw_bd))

    if not re.search(plus_or_minus,raw_bd):

        bd = str.upper(re.sub(caret_money,'', raw_bd))
        bd = re.sub(comma_or_period, ref_base, bd)
        qd = [10**((ord(q) - 33) * -0.1) for q in raw_qd]
        paired_bd_qd = [(b,q) for b,q in zip(bd,qd) if b not in ['>','<','*','n','N']]
        if len(paired_bd_qd) < coverage_cut:
            return [],[], 0

        bd, qd = zip(*paired_bd_qd)
        bd = list(bd)
        qd = list(qd)

    else:

        i = 0   # index for base data
        j = 0   # index for qual data

        while i < len(raw_bd):
            if raw_bd[i] == '.' or raw_bd[i] == ',':   # reference base call
                bd.append(ref_base)
                qd.append(10**((ord(raw_qd[j]) - 33) * -0.1))
                i += 1
                j += 1
            elif raw_bd[i] in base_letters:            # variant call
                bd.append(str.upper(raw_bd[i]))
                qd.append(10**((ord(raw_qd[j]) - 33) * -0.1))
                i += 1
                j += 1
            elif raw_bd[i] == '+' or raw_bd[i] == '-': # indel
                indel_count += 1
                num = int(raw_bd[i+1])
                i += 2
                while(raw_bd[i] in numbers):
                    num *= 10
                    num += int(raw_bd[i])
                    i += 1
                i += num
            elif raw_bd[i] == '^':
                i += 2
            elif raw_bd[i] == '$':
                i += 1
            elif raw_bd[i] in '><*nN':                 # reference skip or deletion
                i += 1
                j += 1

        assert(i == len(raw_bd))
        assert(j == len(raw_qd))

    assert(len(bd) == len(qd))

    for b in bd:
        assert(b in bases)

    return bd, qd, indel_count


def parse_bedfile(input_file):

    boundaries = []
    with open(input_file,'r') as inf:
        for line in inf:

            if len(line) < 3:
                continue

            el = line.strip().split('\t')

            chrom = el[0]
            start = int(el[1])
            stop  = int(el[2])

            boundaries.append((chrom, start, stop))

    return boundaries

# boundary_files[cell*chamber+chamber] should hold name of bed file with fragment boundaries
# use no_cross_chamber_info to ignore cross-chamber information
def call_chamber_alleles(input_file, output_file, boundary_files=None, SNPs_only=False, no_cross_chamber_info=False, minimum_coverage=5):

    if boundary_files != None:
        fragment_boundaries = [[] for i in range(n_cells)]
        for cell in range(0,n_cells):  # for each cell
            for chamber in range(0,n_chambers):
                bfile = boundary_files[cell*n_chambers + chamber]
                fragment_boundaries[cell].append(parse_bedfile(bfile))

    mixed_allele_priors = compute_mixed_allele_priors()
    processed = 0
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
            base_count = {'A':0,'G':0,'T':0,'C':0}
            indel_count = 0
            for cell_num in range(n_cells):
                for ch_num in range(n_chambers):

                    # ensure that position falls inside a called fragment
                    if boundary_files != None:
                        if fragment_boundaries[cell_num][ch_num] == []:
                            continue
                    
                        f_chrom, f_start, f_end = fragment_boundaries[cell_num][ch_num][0]

                        # if we're behind fragment start, skip this spot
                        if chr_num[chrom] < chr_num[f_chrom] or (chrom == f_chrom and pos < f_start):
                            continue

                        # if we're ahead of fragment start, skip to later fragment boundaries
                        while 1:
                            if fragment_boundaries[cell_num][ch_num] == []:
                                break
                            f_chrom, f_start, f_end = fragment_boundaries[cell_num][ch_num][0]
                            if chr_num[chrom] > chr_num[f_chrom] or (chrom == f_chrom and pos >= f_end):
                                fragment_boundaries[cell_num][ch_num].pop(0)
                            else:
                                break
    
                        # if we're not inside fragment, continue
                        if not(chrom == f_chrom and pos >= f_start and pos < f_end):
                            continue

                    flat_ix = cell_num * n_chambers + ch_num
                    col_ix = 3 + 4 * flat_ix

                    depth = int(el[col_ix])

                    if depth < coverage_cut or el[col_ix + 1] == '*':
                        continue

                    raw_bd = el[col_ix + 1]
                    raw_qd = el[col_ix + 2]

                    bd, qd, ic = parse_mpileup_base_qual(raw_bd, raw_qd, ref_base)

                    indel_count += ic

                    if len(bd) < coverage_cut:
                        continue

                    nonzero_chambers[cell_num].append(ch_num)
                    nonzero_chamber_count += 1

                    for b in bd:
                        base_count[b] += 1
                        if b != ref_base:
                            total_nonref += 1

                    base_data[cell_num][ch_num] = bd[:max_cov]
                    qual_data[cell_num][ch_num] = qd[:max_cov]


            base_count = sorted(list(base_count.items()),key=lambda x: x[1],reverse=True)
            major_allele = base_count[0][0]

            fast_mode = False
            condensed_genotype_set = None

            if indel_count > INDEL_COUNT_LIMIT:
                tags.append('ADJACENT_INDEL_OR_CLIP')

            if total_nonref < min_nonref:
                fast_mode = True # skip heavy computation for likely non-SNV
                condensed_genotype_set = [G for G in genotypes if major_allele in G]

            if SNPs_only and total_nonref < min_nonref:
                continue

            outline_el = ['*'] + ['*']*(n_chambers*n_cells)

            too_many_chambers = False
            for nz in nonzero_chambers:
                if len(nz) > 4:
                    too_many_chambers = True
                    tags.append('TOO_MANY_CHAMBERS')
                    break

            maj_alleles = set()
            for cell in range(0,n_cells):  # for each cell
                for chamber in nonzero_chambers[cell]:

                    basecounts = defaultdict(int)
                    for base in base_data[cell][chamber]:
                        basecounts[base] += 1

                    max_k = 'N'
                    max_v = -1
                    for k,v in basecounts.items():
                        if v > max_v:
                            max_v = v
                            max_k = k

                    if(max_k == 'N'):
                        set_trace()
                    maj_alleles.add(max_k)

            if len(maj_alleles) >= 3:
                tags.append('TOO_MANY_ALLELES')

            if nonzero_chamber_count > 0 and not too_many_chambers:

                pr_one_ch = precompute_pr_one_chamber(base_data, qual_data, nonzero_chambers, fast_mode)

                if no_cross_chamber_info:
                    outline_el, SNP = pr_chambers_BASIC(pr_one_ch, nonzero_chambers, mixed_allele_priors, ref_base, condensed_genotype_set)
                    tags.append('NO_CROSS_CHAMBER_INFO')
                else:
                    outline_el, SNP = pr_genotype(pr_one_ch, nonzero_chambers, ref_base, condensed_genotype_set)

                if SNP:
                    tags.append('SNP')

            for cell in range(n_cells):
                for chamber in range(n_chambers):
                    if len(base_data[cell][chamber]) < minimum_coverage:
                        outline_el[1 + cell*n_chambers + chamber] = '*'

            outline_el = [chrom, str(pos+1), ref_base] + outline_el
            tag_info = ';'.join(tags) if tags != [] else 'N/A'
            outline_el.append(tag_info)
            outline_el.append('PILEUP:')
            outline_el = outline_el + el[3:]
            outline = '\t'.join(outline_el)

            print(outline, file=opf)

            processed += 1


if __name__ == '__main__':
    t1 = time.time()
    args = parseargs()
    n_cells = args.num_cells

    call_chamber_alleles(args.input_file, args.output_file, args.fragment_boundaries, SNPs_only=False, no_cross_chamber_info=False)
    t2 = time.time()

    print("TOTAL TIME: {} s".format(int(t2-t1)))
