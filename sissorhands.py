import argparse
import itertools
from math import log10
import time
import pickle
from collections import defaultdict
from file_processing import parse_bedfile, parse_mpileup_base_qual

def addlogs(a,b):
    if a > b:
        return a + log10(1 + pow(10, b - a))
    else:
        return b + log10(1 + pow(10, a - b))

desc = 'Use cross-chamber information to call chamber alleles in SISSOR data'

###############################################################################
# CONSTANTS
###############################################################################

cells = ['PGP1_21','PGP1_22','PGP1_A1']
coverage_cut = 3 # less than this many reads in a chamber are ignored
max_cov = 100
min_nonref = 3
chroms = ['chr{}'.format(i) for i in range(1,23)] + ['chrX','chrY']
n_chambers = 24
chambers = list(range(1,n_chambers+1))
bases = ['A','T','G','C']
diploid_genotypes = list(itertools.combinations_with_replacement(bases,2))
parameters_dir = 'parameters'
n_cells = 3
numbers = '0123456789'
base_letters = 'ACGTacgt'
mixed_alleles = [('A',),('T',),('G',),('C',), ('A', 'T'),('A', 'G'),('A', 'C'),('T', 'G'),('T', 'C'),('G', 'C')]
tinylog = -1e5
INDEL_COUNT_LIMIT = 3
hemizygous_chroms = {'chrX','X','chrY','Y'}
haploid_genotypes = [('A',),('T',),('G',),('C')]

chr_num = dict()
for i,chrom in enumerate(chroms):
    chr_num[chrom] = i

###############################################################################
# PARSE STDIN
###############################################################################

def parseargs():

    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-i', '--input_file', nargs='?', type = str, help='input file, pileup with n cells x 24 chambers', default=None)
    parser.add_argument('-o', '--output_file', nargs='?', type = str, help='file to write output to', default=None)
    parser.add_argument('-b', '--fragment_boundaries', nargs='+', type = str, help='', default=None)
    parser.add_argument('-n', '--num_cells', nargs='?', type = int, help='number of cells being analyzed, default: 3', default=3)

    args = parser.parse_args()
    return args

###############################################################################
# PARAMETERS
###############################################################################

p_null = pickle.load(open("{}/p_null.p".format(parameters_dir), "rb" )) # estimated p_null, the probability of a strand not being sampled.
MDA_dist = pickle.load(open("{}/MDA_dist.p".format(parameters_dir), "rb" )) # distribution of observed MDA error
ch_priors = pickle.load(open("{}/ch_priors.p".format(parameters_dir), "rb" )) # prior probability of sampling from a chamber
cov_frac_dist = pickle.load(open( "parameters/cov_frac_dist.p", "rb")) # PROBABILITY OF SEEING PARENT 1 ALLELE IN MIXED ALLELE CHAMBER
hom_config_probs = pickle.load(open( "parameters/hom_config_probs.p", "rb")) # STRAND-TO-CHAMBER CONFIGURATIONS
het_config_probs = pickle.load(open( "parameters/het_config_probs.p", "rb"))
diploid_genotype_priors = pickle.load(open( "parameters/diploid_genotype_priors.p", "rb"))
haploid_genotype_priors = pickle.load(open( "parameters/haploid_genotype_priors.p", "rb"))

MDA_ERR = 1e-5 # the probability of consensus error due to MDA
MDA_COR = log10(1 - MDA_ERR)
MDA_ERR = log10(MDA_ERR) - log10(3) # the log probability of an MDA error, divided by 3

NUM_BINS = 20
COV_INTERVAL = 10

# a heuristic for finding filtering mismapped positions, CNV, other issues.
# if the same mixture of alleles (e.g. 'A'+'T') is observed in two separate libraries
# above a threshold then there is likely a systemic issue like this
# MDA error should be random and strand-overlap shouldn't happen in multiple chambers often.
def same_mixture_occurs_twice(base_data,nonzero_chambers):

    mixtures_seen = set()

    for cell in range(0,n_cells):  # for each cell
        for chamber in nonzero_chambers[cell]:

            pileup = base_data[cell][chamber]
            counter = defaultdict(int)
            for base in pileup:
                counter[base] += 1

            bases_seen = set()
            for base, v in counter.items():
                if v >= 3:
                    bases_seen.add(base)

            if len(bases_seen) >= 3:
                return True

            elif len(bases_seen) == 2:

                mixture_str = ''.join(sorted(list(bases_seen)))

                if mixture_str in mixtures_seen:
                    return True

                mixtures_seen.add(mixture_str)

    return False

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
# accessed as mixed_allele_priors[r][x][y]
# r is the reference allele
# where x is the number of chambers with reads
# and y is the allele mixture as an alphabetically sorted tuple
def compute_mixed_allele_priors():

    mixed_allele_priors = dict()

    for ref_allele in bases:

        mixed_allele_priors[ref_allele] = {1:dict(),2:dict(),3:dict(),4:dict()}

        for i in range(1,5):
            for allele in mixed_alleles:
                mixed_allele_priors[ref_allele][i][allele] = tinylog
            for G in diploid_genotypes:
                nz = list(range(i))
                het = (G[0] != G[1])
                cfgs = singlecell_config(het,nz)

                for (c1,c2,c3,c4), p in cfgs:
                    # probability of configuration
                    p += diploid_genotype_priors[ref_allele][G] - log10(len({c1,c2,c3,c4}))

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

# create dictionaries that pre-compute some probabilities
# to speed up likelihood calculations for read base data observed in a chamber
one_allele_cache = dict()
two_allele_cache = dict()
def compute_caches():
    for i in range(0,NUM_BINS+1):
        for j in range(i+1,NUM_BINS+1):
            for q in range(33,130):
                for a1_match in [False,True]:
                    for a2_match in [False,True]:
                        qual = 10**((q - 33) * -0.1)
                        frac = i / j if i > 1 else 1e-10

                        x1 = 1.0-qual
                        x2 = qual
                        p1 = x1 if a1_match else x2
                        p2 = x1 if a2_match else x2
                        p3 = frac*p1 + (1-frac)*p2
                        two_allele_cache[((i / j),qual,a1_match,a2_match)] = log10(p3) if p3 > 0 else tinylog

    for q in range(33,130):
        qual = 10**((q - 33) * -0.1)

        one_allele_cache[(qual,True)]  = log10(1.0-qual) if 1.0-qual > 0 else tinylog   #log10((1.0-qual)*(1.0-omega_nolog) + omega_nolog*qual)
        one_allele_cache[(qual,False)] = log10(qual)       #log10(omega_nolog*(1.0-qual) + (1-omega_nolog)*qual)

compute_caches()

# compute the probability of the read data in a single chamber,
# given that the alleles actually present in the chamber are alleles_present.
# base_data[i][j] contains the base pileup for cell i chamber j
# qual_data[i][j] contains the quality scores for cell i chamber j
# fast_mode is a parameter to speed up on positions that have only a couple variant bases observed,
#    so they are almost certainly reference calls. It considers simple 50% odds for
#    observing a given allele when there are overlapped strands in the same chamber,
#   rather than summing over a complete estimated distrbution.
def pr_one_chamber_data(alleles_present, base_data, qual_data, fast_mode):


    n = len(base_data)
    assert(n == len(qual_data))
    p = 0

    # note that to save computation we don't allow both MDA error "mixture" and
    # strand mixture at the same time.
    # this would be really slow and messy without too much added benefit

    # if we are considering that a single strand has fallen in the chamber
    if len(alleles_present) == 1:

        p = tinylog

        L = len(MDA_dist[n])
        for i, p_frac in enumerate(MDA_dist[n]): # MDA_dist[0] has MDA error distribution for 1..10, MDA_error_dist[1] has MDA error distribution for 11..20
            if p_frac == 0:
                continue
            p_frac = log10(p_frac / 3) # divide by 3 because there are 3 possible bases to error to

            frac = i / L

            for MDA_allele in bases:
                if MDA_allele == alleles_present[0]: # can't have main and secondary allele be the same...
                    continue
                p0 = p_frac # 3 possible bases to have an MDA error to
                for base,qual in zip(base_data, qual_data): # sum over all base/quals in pileup
                    p0 += two_allele_cache[(frac, qual, (MDA_allele == base), (alleles_present[0] == base))]

                p = addlogs(p, p0)

    # if we are considering that two strands have fallen in the chamber
    elif len(alleles_present) == 2:

        # if fast mode (mostly reference bases) we use faster computation, assuming basically binomial
        if fast_mode:
            for base,qual in zip(base_data, qual_data):
                p += two_allele_cache[(0.5, qual, (alleles_present[0] == base), (alleles_present[1] == base))]
        # we've observed a lot of non-reference bases. slower computation which sums over entire strand mixture probability distribution
        else:
            p = tinylog
            L = len(cov_frac_dist[n])
            for i, p0 in enumerate(cov_frac_dist[n]):
                p0 = log10(p0) if p0 > 0 else None
                if p0 == None:
                    continue
                frac = i / L
                for base,qual in zip(base_data, qual_data):
                    p0 += two_allele_cache[(frac, qual, (alleles_present[0] == base), (alleles_present[1] == base))]

                p = addlogs(p, p0)

    return p

# pr_one_ch[(i,j,(a,b))] is a precomputed dictionary. It returns the probability that cell i, chamber j, has allele (a,b) present (may be a single base, or a mixture from contamination)
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

# core of the SISSORhands approach
# compute the probability that each allele is present, by computing the probability of
# each genotype with a bayesian calculation.
# to calculate the likelihood of data, sum over all possible events of DNA strands
# being distributed to SISSOR reaction chambers for amplification.

# pr_one_ch[(i,j,(a,b))] is a precomputed dictionary. It returns the probability that cell i, chamber j, has allele (a,b) present (may be a single base, or a mixture from contamination)
# nonzero_chambers[i] is a list of chambers for cell i that have data present
# mixed_allele_priors is a dictionary containing prior probability for a given mixture of alleles in a chamber. This is basically not used at all since the solo-chamber calls are not being used.
# ref_allele is the reference allele at the position we are calling
# condensed_genotype_set is a restricted set of genotypes. This parameter is currently not used at all, but can be used to limit the set of genotypes being considered, for a performance speedup.
# haploid is True if this is chrX or chrY and false otherwise.
def pr_genotype(pr_one_ch, nonzero_chambers, mixed_allele_priors, ref_allele, condensed_genotype_set, haploid):

    probs = dict()
    outline_el = ['*','*'] + sum([(['CELL{}'.format(i)]+['*'] * n_chambers) for i in range(1,n_cells+1)],[])

    hom_het_configs = [multicell_config(False,nonzero_chambers),multicell_config(True,nonzero_chambers)]
    p_assignment = defaultdict(lambda: tinylog)
    genotype_set = diploid_genotypes if condensed_genotype_set == None else condensed_genotype_set

    for G in genotype_set:

        configs = hom_het_configs[0] if G[0] == G[1] else hom_het_configs[1]
        p_total = tinylog

        if haploid:
            configs = hom_het_configs[1] # spoofed as heterozygous
            if G[0] != G[1]: # only consider
                continue


        for config in configs:

            p = 0
            assignments = set()

            if haploid:
                skip = False
                for i in range(0,n_cells):
                    (c1,c2,c3,c4), p_cell_cfg = config[i]
                    if not (c3 == n_chambers and c4 == n_chambers):
                        skip = True
                        break

                if skip:
                    continue

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

                        # we must account for the fact that the entire consensus may
                        # be wrong due to MDA error
                        for base1, base2 in diploid_genotypes:

                            x1 = MDA_COR if base1 == allele1 else MDA_ERR
                            x2 = MDA_COR if base2 == allele2 else MDA_ERR

                            if (base1,base2) == alleles_present:
                                p1 = pr_one_ch[(i,j,(base1,base2))] + MDA_COR + MDA_COR
                                p2 = pr_one_ch[(i,j,(base1,base2))] + MDA_ERR + MDA_ERR
                                p1 = addlogs(p1,p2)
                            elif base1 == base2:
                                p1 = pr_one_ch[(i,j,(base1,))] + x1 + x2
                            else:
                                p1 = pr_one_ch[(i,j,(base1,base2))] + x1 + x2 + log10(2)

                            p0 = addlogs(p0, p1)

                        p += p0

                    assignments.add((i,j,alleles_present))

            p_total = addlogs(p_total, p) # update genotype likelihood sum

            for assignment in assignments:

                if not haploid:
                    p_assignment[assignment] = addlogs(p_assignment[assignment], p + diploid_genotype_priors[ref_allele][G])
                else:
                    p_assignment[assignment] = addlogs(p_assignment[assignment], p + haploid_genotype_priors[ref_allele][G[0]])


        if not haploid:
            probs[G] = p_total + diploid_genotype_priors[ref_allele][G]
        else:
            probs[(G[0],)] = p_total + haploid_genotype_priors[ref_allele][G[0]]
    # denominator for bayes rule posterior calculation
    total = tinylog
    for G,p in probs.items():
        total = addlogs(total,p)

    res = []
    SNP = False
    max_posterior = -1
    max_G = ('N','N')
    allele_count = defaultdict(lambda: tinylog)

    current_genotypes = diploid_genotypes if not haploid else haploid_genotypes
    for G in current_genotypes:
        if G in probs:
            for allele in bases:
                if allele in G:
                    allele_count[allele] = addlogs(allele_count[allele],probs[G])
            posterior = 10**(probs[G] - total)
            res.append((G,posterior))
            if posterior > max_posterior:
                max_posterior = posterior
                max_G = G
        else:
            res.append((G,0.0))

    if max_G != (ref_allele,ref_allele) and max_G != (ref_allele,):
        SNP = True

    res2 = [(allele,10**(allele_count[allele] - total)) for allele in bases]

    allele_str = ';'.join(['{}:{}'.format(a,p) for a,p in res2])
    gen_str = ';'.join(['{}:{}'.format(''.join(g),p) for g,p in res])
    outline_el[0] = allele_str
    outline_el[1] = gen_str


    # determine most likely allele in each chamber, in light of information from other chambers

    for cell in range(0,n_cells):  # for each cell
        for chamber in nonzero_chambers[cell]:

            #num_nonzero = len(nonzero_chambers[cell])
            #pr_all_chamber_data(allele, ref_allele, cell, chamber, pr_one_ch, nonzero_chambers, hom_het_configs)
            probs = []
            for allele in mixed_alleles:
                if haploid and len(allele) > 1:
                    continue
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

            outline_el[3 + cell + cell*n_chambers + chamber] = out_str


    # determine most likely allele in each chamber, IGNORING information from other chambers
    # this is not used in any way in the SISSOR manuscript

    for cell in range(0,n_cells):  # for each cell
        for chamber in nonzero_chambers[cell]:

            num_nonzero = len(nonzero_chambers[cell])
            #pr_all_chamber_data(allele, ref_allele, cell, chamber, pr_one_ch, nonzero_chambers, hom_het_configs)
            probs = []
            for allele in mixed_alleles:
                if not haploid:
                    prob = pr_one_ch[(cell,chamber,allele)] + mixed_allele_priors[ref_allele][num_nonzero][allele]
                else:
                    if len(allele) > 1:
                        continue
                    prob = pr_one_ch[(cell,chamber,allele)] + haploid_genotype_priors[ref_allele][G[0]]
                probs.append((allele,prob))

            # denominator for bayes rule posterior calculation
            total = tinylog
            for a,p in probs:
                total = addlogs(total,p)

            res = []
            max_posterior = -1
            for a,p in probs:
                posterior = 10**(p - total)
                res.append((a,posterior))

                if posterior > max_posterior:
                    max_posterior = posterior

            out_str = ';'.join(['{}:{}'.format(''.join(a),p) for a,p in res])

            outline_el[3 + cell + cell*n_chambers + chamber] += '|'+out_str


    return outline_el, SNP

###############################################################################
# MAIN FUNCTION AND PARSING
###############################################################################

# given an mpileup file over N cells x 24 SISSOR libraries, calls the most likely
# consensus allele as well as the most likely allele in each chamber
# input file is a samtools mpileup file
# output is a special format ("chamber call file") that has a consensus allele probabilities,
#  probabilities for each allele in each chamber at each genomic position,
# as well as a cleaned up version of the original pileups and various tags.
def call_chamber_alleles(input_file, output_file, boundary_files=None, SNPs_only=False):

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

            haploid = (chrom in hemizygous_chroms)

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
            #major_allele = base_count[0][0]

            fast_mode = False
            condensed_genotype_set = None

            if indel_count > INDEL_COUNT_LIMIT:
                tags.append('ADJACENT_INDEL_OR_CLIP')

            if total_nonref < min_nonref:
                fast_mode = True # skip heavy computation for likely non-SNV
                #condensed_genotype_set = [G for G in genotypes if major_allele in G]

            if SNPs_only and total_nonref < min_nonref:
                continue

            outline_el = ['*','*'] + ['*']*((n_chambers+1)*n_cells)

            too_many_chambers = False
            for nz in nonzero_chambers:
                if (not haploid and len(nz) > 4) or (haploid and len(nz) > 2):
                    too_many_chambers = True
                    tags.append('TOO_MANY_CHAMBERS')
                    break

            maj_alleles = set()
            for cell in range(0,n_cells):  # for each cell
                for chamber in nonzero_chambers[cell]:

                    basecounts = defaultdict(int)

                    if len(base_data[cell][chamber]) < 3: # we can expect some weirdness at low coverage
                        continue

                    for base in base_data[cell][chamber]:
                        basecounts[base] += 1

                    max_k = 'N'
                    max_v = -1
                    for k,v in basecounts.items():
                        if v > max_v:
                            max_v = v
                            max_k = k

                    maj_alleles.add(max_k)

            if len(maj_alleles) >= 5:
                tags.append('TOO_MANY_ALLELES')

            if same_mixture_occurs_twice(base_data, nonzero_chambers):
                tags.append('SAME_MIXTURE_OCCURS_TWICE')

            if nonzero_chamber_count > 0 and not too_many_chambers:

                pr_one_ch = precompute_pr_one_chamber(base_data, qual_data, nonzero_chambers, fast_mode)

                outline_el, SNP = pr_genotype(pr_one_ch, nonzero_chambers,mixed_allele_priors, ref_base, condensed_genotype_set, haploid)

                if SNP:
                    tags.append('SNP')

            for cell in range(n_cells):
                for chamber in range(n_chambers):
                    if outline_el[3 + cell + cell*n_chambers + chamber] != '*':
                        outline_el[3 + cell + cell*n_chambers + chamber] = outline_el[3 + cell + cell*n_chambers + chamber] + '|'+''.join(base_data[cell][chamber])

            outline_el = [chrom, str(pos+1), ref_base] + outline_el
            tag_info = ';'.join(tags) if tags != [] else 'N/A'
            outline_el.append(tag_info)
            outline_el = outline_el
            outline = '\t'.join(outline_el)

            print(outline, file=opf)

            processed += 1

# main function
# parse input and run function to call alleles
if __name__ == '__main__':
    t1 = time.time()
    args = parseargs()
    n_cells = args.num_cells

    call_chamber_alleles(args.input_file, args.output_file, args.fragment_boundaries, SNPs_only=False)
    t2 = time.time()

    print("TOTAL TIME: {} s".format(int(t2-t1)))
