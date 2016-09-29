#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 14:51:49 2016

@author: peter
"""

import argparse
import numpy as np
import pickle
import os
#import sys
#from collections import defaultdict
import itertools
from math import log10
from functools import reduce
import operator
import re
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
    parser.add_argument('-i', '--input_dir', nargs='?', type = str, help='input data dir, format {DIR}/PGP1_*/ch*.pileup', default=default_input_dir)
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
# MAIN FUNCTION AND PARSING
###############################################################################
    
def estimate_parameters(input_dir, output_dir):

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
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
                
    coverage_cut = 2

    chrX_covs = defaultdict(int)
    
    #MDA_sum   = -1e100
    #MDA_total = 0
    chamber_position_counts = defaultdict(int)
    strand_coverage_counts = defaultdict(int)
    total_sampled_cell_positions = 0
    
    # until all files have reached stopiteration
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
        
        if chrom == 'chrX':
      
            for i in np.where(on_chrom_pos)[0]:
                
                cell_num = int(i / n_chambers)
                ch_num   = i % 24
                
                el = el_lst[i]
                if not el:
                    continue

                d = int(el[3])
                
                if d < coverage_cut:
                    continue

                chrX_covs[d] += 1

        elif chrom != 'chrX' and chrom != 'chrY':  #elif random.random() < sample_rate and chrom != 'chrX' and chrom != 'chrY':
                
            total_sampled_cell_positions += n_cells
            strand_counts = defaultdict(int)         
            for i in np.where(on_chrom_pos)[0]:
                
                cell_num = int(i / n_chambers)
                ch_num   = i % 24
                
                el = el_lst[i]
                if not el:
                    continue
                
                if int(el[3]) < coverage_cut:
                    continue

                strand_counts[cell_num] += 1
                chamber_position_counts[ch_num] += 1
                
            for cell_num in range(n_cells):
                strand_coverage_counts[strand_counts[cell_num]] += 1
                
        # iterate to the next lines for each line on the current index        
        for i in np.where(on_chrom_pos)[0]:

            lines[i] = safenext(input_files[i])
            
    
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
    
    # estimate rate of MDA error
    #MDA_error_rate = MDA_sum - log10(MDA_total)
    
    # WRITE PARAMETERS TO PICKLE FILES
    pickle.dump(chrX_covs_probs, open( "{}/chrX_covs_probs.p".format(output_dir), "wb" ))
    pickle.dump(cov_frac_dist, open( "{}/cov_frac_dist.p".format(output_dir), "wb" ))
    #pickle.dump(MDA_error_rate, open( "{}/MDA_error_rate.p".format(output_dir), "wb" ))
    pickle.dump(chamber_position_counts, open( "{}/chamber_position_counts.p".format(output_dir), "wb" ))
    pickle.dump(strand_coverage_counts, open( "{}/strand_coverage_counts.p".format(output_dir), "wb" ))
    pickle.dump(total_sampled_cell_positions, open( "{}/total_sampled_cell_positions.p".format(output_dir), "wb" ))
    # close all chambers

    for handle in input_files:
        handle.close()
        
if __name__ == '__main__':
    args = parseargs()
    estimate_parameters(args.input_dir, args.output_dir)
