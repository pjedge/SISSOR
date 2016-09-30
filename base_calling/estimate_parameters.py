#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 14:51:49 2016

@author: peter
"""

import argparse
import pickle
import os
#import sys
#from collections import defaultdict
import itertools
from math import log10
from functools import reduce
import operator
from pdb import set_trace
from collections import defaultdict
import random
import bisect
if False:
    set_trace() # to dodge warnings that pdb isn't being used.
#from copy import deepcopy
#import sys

desc = 'Use cross-chamber information to call chamber alleles in SISSOR data'

default_input_dir = 'pileups_subsample'
default_output_dir = 'parameters'

###############################################################################
# PARSE STDIN
###############################################################################

def parseargs():

    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-i', '--input_file', nargs='?', type = str, help='input file with 3 cells x 24 chambers pileups', default=default_input_dir)
    parser.add_argument('-o', '--output_dir', nargs='?', type = str, help='file to write output to', default=default_output_dir)

    args = parser.parse_args()
    return args

###############################################################################
# CONSTANTS
###############################################################################

sample_rate = 0.001
chroms = ['chr{}'.format(i) for i in range(1,23)] + ['chrX','chrY']
cells = ['PGP1_21','PGP1_22','PGP1_A1']
n_chambers = 24
chambers = list(range(1,n_chambers+1))
bases = ['A','T','G','C']
genotypes = list(itertools.combinations_with_replacement(bases,2))

short_chrom_names = set([str(x) for x in range(1,23)]+['X','Y'])
chrom_names = set(['chr'+str(x) for x in range(1,23)]+['chrX','chrY'])    
    
###############################################################################
# TEMPORARILY HARDCODED (need to be estimated form data in final version)
###############################################################################

p_null = 0.5          # probability that strand is not sampled in chamber. hardcoded until we can estimate it
ch_priors = [log10((1-p_null)/n_chambers)]*n_chambers + [log10(p_null)]   # prior probability of sampling a given chamber

n_cells = 3

# PROBABILITY OF SEEING PARENT 1 ALLELE IN MIXED ALLELE CHAMBER

P_parent1_lst = [(0.001,log10(0.495)),(0.5,log10(0.01)),(0.999,log10(0.495))]
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

def combine_parameters(suffixes):

    chrX_covs = defaultdict(int)
    chamber_position_counts = defaultdict(int)
    strand_coverage_counts = defaultdict(int)
    total_sampled_cell_positions = 0
    
    for suffix in suffixes:
        
        # LOAD PICKLE FILES
        partial_chrX_covs = pickle.load(open("parameters/split/chrX_covs{}.p".format(suffix), "rb" ))
        partial_chamber_position_counts = pickle.load(open("parameters/split/chamber_position_counts{}.p".format(suffix), "rb" ))
        partial_strand_coverage_counts = pickle.load(open("parameters/split/strand_coverage_counts{}.p".format(suffix), "rb" ))
        partial_total_sampled_cell_positions = pickle.load(open("parameters/split/total_sampled_cell_positions{}.p".format(suffix), "rb" ))
        
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
    # cov_frac_dist is used in the following way:
    # cov_frac_dist[cov][i] where i=1..cov tells the probability of allele p1 being present in fraction 
    num_samples = 100000000
    cov_frac_dist = defaultdict(list)
    cov_frac_dist_counts = defaultdict(list)
    
    for i in range(num_samples):
        cov1 = chrX_covs_func()
        cov2 = chrX_covs_func()
        tcov = cov1 + cov2
        if cov_frac_dist_counts[tcov] == []:
            cov_frac_dist_counts[tcov] = [0]*tcov
        cov_frac_dist_counts[tcov][cov1] += 1
    
    for cov in cov_frac_dist_counts.keys():
        for i, count in enumerate(cov_frac_dist_counts[cov]):
            cov_frac_dist[cov].append(count/sum(cov_frac_dist_counts[cov]))
        assert(len(cov_frac_dist[cov]) == cov)
        
   # WRITE PARAMETERS TO PICKLE FILES
    pickle.dump(chrX_covs, open("parameters/chrX_covs.p", "wb"))
    pickle.dump(chamber_position_counts, open("parameters/chamber_position_counts.p","wb"))
    pickle.dump(strand_coverage_counts, open("parameters/strand_coverage_counts.p","wb"))
    pickle.dump(total_sampled_cell_positions, open("parameters/total_sampled_cell_positions.p","wb"))
    pickle.dump(cov_frac_dist, open("parameters/cov_frac_dist.p","wb"))

def estimate_parameters(input_file, suffix=''):

    coverage_cut = 3

    chrX_covs = defaultdict(int)
    chamber_position_counts = defaultdict(int)
    strand_coverage_counts = defaultdict(int)
    total_sampled_cell_positions = 0
    
    with open(input_file,'r') as ipf:
        for line in ipf:
            
            el = line.strip().split('\t')

            chrom = el[0]
            ref_base = str.upper(el[2])

            if ref_base == 'N':
                continue
            assert(ref_base in bases)
                
            strand_counts = defaultdict(int)             

            for cell_num in range(n_cells):
                for ch_num in range(n_chambers):
                
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
    pickle.dump(chrX_covs, open( "parameters/split/chrX_covs{}.p".format(suffix), "wb" ))
    pickle.dump(chamber_position_counts, open( "parameters/split/chamber_position_counts{}.p".format(suffix), "wb" ))
    pickle.dump(strand_coverage_counts, open( "parameters/split/strand_coverage_counts{}.p".format(suffix), "wb" ))
    pickle.dump(total_sampled_cell_positions, open( "parameters/split/total_sampled_cell_positions{}.p".format(suffix), "wb" ))

if __name__ == '__main__':
    args = parseargs()
    estimate_parameters(args.input_file, args.output_dir)
