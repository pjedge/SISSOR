# -*- coding: utf-8 -*-
"""
Created on Thu May 19 21:44:48 2016

@author: peter
"""

import argparse
import sys

def parse_args():

    parser = argparse.ArgumentParser(description='Split SISSOR fragments with switch errors')
    # paths to samfiles mapping the same ordered set of RNA reads to different genomes
    parser.add_argument('-f', '--fragments', nargs='?', type = str, help='fragment file to process')
    parser.add_argument('-o', '--output', nargs='?', type = str, help='output fragments')
    parser.add_argument('-t', '--threshold', nargs='?', type = str, help='number of consecutive mismatches with respect to another fragment that count as contamination')


    # default to help option. credit to unutbu: http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
    if len(sys.argv) < 3:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    return args

def overlap(f1, f2):
    
    assert(f1[0][0] <= f2[0][0])
    
    return not (f2[0][0] > f1[-1][0]) 

# determine if two fragments are consistent with each other
def consistent(f1, f2, t):

    assert(f1[0][0] <= f2[0][0])    

    # for allele in fragment1            
    f2iter = iter(f2) # iterator for fragment2
    a2 = f2iter.next()    

    count1 = 0
    count2 = 0
    
    for a1 in f1:

        if a1[0] < a2[0]:
            continue
        
        assert(a1[0] == a2[0])
        
        if a1[1] == a2[1]:
            count1 += 1
        else:
            count2 += 1
        
        a2 = f2iter.next()
    
    same_hap = (count1 > count2)
    
    f2iter = iter(f2) # iterator for fragment2
    a2 = f2iter.next()    

    mismatches = 0

    for a1 in f1:

        if a1[0] < a2[0]:
            continue
        
        assert(a1[0] == a2[0])
        
        if (same_hap and a1[1] == a2[1]) or (not same_hap and a1[1] != a2[1]):
            mismatches = 0
        else:
            mismatches += 1
            if mismatches > t:
                return False
    
        a2 = f2iter.next()
        
    return True
        
        
def read_fragment_matrix(frag_matrix):

    flist = []
    
    with open(frag_matrix,"r") as fm:
        for line in fm:
            if len(line) < 2:
                continue

            el = line.strip().split()

            num_blks      = int(el[0])

            call_list  = el[2:(2+2*num_blks)]              # extract base call part of line
            call_list  = zip(*[iter(call_list)]*2)             # list -> tuple list conversion: credit to http://stackoverflow.com/questions/23286254/convert-list-to-a-list-of-tuples-python
            call_list  = [(int(a)-1, b) for a,b in call_list]  # convert index to 0-based integer
            call_list2 = []

            for ix, blk in call_list:
                curr_ix = ix
                for a in blk:
                    call_list2.append((curr_ix, a))
                    curr_ix += 1

            qlist = el[-1]
            qlist = [10**((ord(q) - 33) * -0.1) for q in qlist]

            alist= [(a,b,c) for ((a,b),c) in zip(call_list2,qlist)]
            
            flist.append(alist)

def fix_chamber_contamination(fragments,outfile):
    
    # read fragments into fragment data structure
    flist = read_fragment_matrix(fragments)
    
    # create a bipartite graph representing consistency of reads
    
    
if __name__ == '__main__':
    args = parse_args()
    fix_chamber_contamination(args.fragments,args.output)