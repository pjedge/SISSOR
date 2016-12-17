#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 14:51:49 2016

@author: peter
"""

import pickle
import os
import itertools
from math import log10, factorial
from functools import reduce
import operator
from pdb import set_trace
from collections import defaultdict
import random
import bisect
import math
if False:
    set_trace() # to dodge warnings that pdb isn't being used.

desc = 'Use cross-chamber information to call chamber alleles in SISSOR data'

default_input_dir = 'pileups_subsample'
default_output_dir = 'parameters'
###############################################################################
# CONSTANTS
###############################################################################

chroms = ['chr{}'.format(i) for i in range(1,23)] + ['chrX','chrY']
cells = ['PGP1_21','PGP1_22','PGP1_A1']
n_chambers = 24
chambers = list(range(1,n_chambers+1))
bases = ['A','T','G','C']
genotypes = list(itertools.combinations_with_replacement(bases,2))
het_snp_rate = 0.0005
hom_snp_rate = 0.001
n_cells = 3
chr_num = dict()
for i,chrom in enumerate(chroms):
    chr_num[chrom] = i

###############################################################################
# HELPER FUNCTIONS
###############################################################################

def addlogs(a,b):
    if a > b:
        return a + log10(1 + pow(10, b - a))
    else:
        return b + log10(1 + pow(10, a - b))

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

# product of list
def prod(iterable): # credit to: http://stackoverflow.com/questions/7948291/is-there-a-built-in-product-in-python
    return reduce(operator.mul, iterable, 1)

# credit to David for this function: http://stackoverflow.com/questions/526255/probability-distribution-in-python
def weighted_choice_bisect_compile(items):
    """Returns a function that makes a weighted random choice from items."""
    added_weights = []
    last_sum = 0

    for item, weight in items:
        last_sum += weight
        added_weights.append(last_sum)

    def choice(rnd=random.random, bis=bisect.bisect):
        return items[bis(added_weights, rnd() * last_sum)][0]
    return choice

###############################################################################
# MAIN FUNCTIONS AND PARSING
###############################################################################

def estimate_parameters(suffixes, grams_DNA_before_MDA, grams_DNA_after_MDA):

    chrX_covs = defaultdict(int)
    chamber_position_counts = defaultdict(int)
    strand_coverage_counts = defaultdict(int)
    total_sampled_cell_positions = 0

    for suffix in suffixes:

        # LOAD PICKLE FILES
        partial_chrX_covs = pickle.load(open("parameters/split/chrX_covs.{}.p".format(suffix), "rb" ))
        partial_chamber_position_counts = pickle.load(open("parameters/split/chamber_position_counts.{}.p".format(suffix), "rb" ))
        partial_strand_coverage_counts = pickle.load(open("parameters/split/strand_coverage_counts.{}.p".format(suffix), "rb" ))
        partial_total_sampled_cell_positions = pickle.load(open("parameters/split/total_sampled_cell_positions.{}.p".format(suffix), "rb" ))

        for k,v in partial_chrX_covs.items():
            chrX_covs[k] += v

        for k,v in partial_chamber_position_counts.items():
            chamber_position_counts[k] += v

        for k,v in partial_strand_coverage_counts.items():
            strand_coverage_counts[k] += v

        total_sampled_cell_positions += partial_total_sampled_cell_positions


    chrX_covs_tuple = list(chrX_covs.items())
    total = sum([b for a,b in chrX_covs_tuple])
    chrX_covs_probs = [(a,b/total) for a,b in chrX_covs_tuple]
    chrX_covs_func  = weighted_choice_bisect_compile(chrX_covs_probs)
    # cov_frac_dist_raw is used in the following way:
    # cov_frac_dist_raw[cov][i] where i=1..cov tells the probability of allele p1 being present in fraction
    num_samples = 100000000
    cov_frac_dist_raw = defaultdict(list)
    cov_frac_dist_raw_counts = defaultdict(list)
    cov_frac_dist = defaultdict(list)

    for i in range(num_samples):
        cov1 = chrX_covs_func()
        cov2 = chrX_covs_func()
        tcov = cov1 + cov2
        if cov_frac_dist_raw_counts[tcov] == []:
            cov_frac_dist_raw_counts[tcov] = [0]*tcov
        cov_frac_dist_raw_counts[tcov][cov1] += 1

    for cov in cov_frac_dist_raw_counts.keys():
        for i, count in enumerate(cov_frac_dist_raw_counts[cov]):
            cov_frac_dist_raw[cov].append(count/sum(cov_frac_dist_raw_counts[cov]))
        assert(len(cov_frac_dist_raw[cov]) == cov)

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

    total_position_counts = sum([chamber_position_counts[chamber] for chamber in range(0,24)])
    chamber_proportions = [chamber_position_counts[chamber]/total_position_counts for chamber in range(0,24)]
    chamber_proportions = [x if x != 0 else 1e-10 for x in chamber_proportions]

    p_null = None
    ch_priors = None
    max_likelihood = -float("Inf")
    p_null_samples = 1000

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

    # estimate prior probability of genotypes using strategy described here:
    # http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2694485/
    # "prior probability of each genotype"
    genotype_priors = dict()
    transition = {'A':'G','G':'A','T':'C','C':'T'}
    alleles = ['A','T','G','C']

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

    # ESTIMATE PER-BASE PROBABILITY OF MDA ERROR
    # formula from this paper: http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0105585
    GAIN = grams_DNA_after_MDA / grams_DNA_before_MDA
    omega = log10(3.2e-6 / 2 * math.log2(GAIN)) # probability of MDA error.

   # WRITE PARAMETERS TO PICKLE FILES
    pickle.dump(chrX_covs, open("parameters/chrX_covs.p", "wb"))
    pickle.dump(chamber_position_counts, open("parameters/chamber_position_counts.p","wb"))
    pickle.dump(strand_coverage_counts, open("parameters/strand_coverage_counts.p","wb"))
    pickle.dump(total_sampled_cell_positions, open("parameters/total_sampled_cell_positions.p","wb"))
    pickle.dump(cov_frac_dist, open("parameters/cov_frac_dist.p","wb"))
    pickle.dump(p_null, open("parameters/p_null.p", "wb"))
    pickle.dump(ch_priors, open("parameters/ch_priors.p","wb"))
    pickle.dump(hom_config_probs, open("parameters/hom_config_probs.p","wb"))
    pickle.dump(het_config_probs, open("parameters/het_config_probs.p","wb"))
    pickle.dump(genotype_priors, open("parameters/genotype_priors.p","wb"))
    pickle.dump(omega, open("parameters/omega.p","wb"))

def obtain_counts_parallel(input_file, boundary_files=None, suffix=''):

    coverage_cut = 3
    
    if boundary_files != None:
        fragment_boundaries = [[] for i in range(n_cells)]
        for cell in range(0,n_cells):  # for each cell
            for chamber in range(0,n_chambers):
                bfile = boundary_files[cell*n_chambers + chamber]
                fragment_boundaries[cell].append(parse_bedfile(bfile))

    chrX_covs = defaultdict(int)
    chamber_position_counts = defaultdict(int)
    strand_coverage_counts = defaultdict(int)
    total_sampled_cell_positions = 0

    with open(input_file,'r') as ipf:
        for line in ipf:

            el = line.strip().split('\t')

            chrom = el[0]
            pos   = int(el[1]) - 1
            ref_base = str.upper(el[2])

            if ref_base == 'N':
                continue
            assert(ref_base in bases)

            strand_counts = defaultdict(int)

            for cell_num in range(n_cells):
                for ch_num in range(n_chambers):
                    
                    if boundary_files != None:

                        # ensure that position falls inside a called fragment
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

                    if chrom == 'chrX':
                        chrX_covs[depth] += 1
                    elif chrom != 'chrY':
                        strand_counts[cell_num] += 1
                        chamber_position_counts[ch_num] += 1

            if chrom != 'chrX' and chrom != 'chrY':
                for cell_num in range(n_cells):
                    total_sampled_cell_positions += n_cells
                    strand_coverage_counts[strand_counts[cell_num]] += 1

    output_dir = 'parameters/split'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # WRITE PARAMETERS TO PICKLE FILES
    pickle.dump(chrX_covs, open( "parameters/split/chrX_covs.{}.p".format(suffix), "wb" ))
    pickle.dump(chamber_position_counts, open( "parameters/split/chamber_position_counts.{}.p".format(suffix), "wb" ))
    pickle.dump(strand_coverage_counts, open( "parameters/split/strand_coverage_counts.{}.p".format(suffix), "wb" ))
    pickle.dump(total_sampled_cell_positions, open( "parameters/split/total_sampled_cell_positions.{}.p".format(suffix), "wb" ))
