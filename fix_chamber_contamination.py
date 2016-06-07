# -*- coding: utf-8 -*-
"""
Created on Thu May 19 21:44:48 2016

@author: peter
"""

import argparse
import sys
import numpy as np
from collections import defaultdict
import itertools

DEBUG = False

def parse_args():

    parser = argparse.ArgumentParser(description='Split SISSOR fragments with switch errors')
    # paths to samfiles mapping the same ordered set of RNA reads to different genomes
    parser.add_argument('-f', '--fragments', nargs='?', type = str, help='fragment file to process')
    parser.add_argument('-o', '--output', nargs='?', type = str, help='output fragments')
    parser.add_argument('-t', '--threshold', nargs='?', type = float, help='number of consecutive mismatches with respect to another fragment that count as contamination', default = 4)
    parser.add_argument('-c', '--filter_cov1', action='store_true', help='filter out positions with coverage 1', default = False)


    # default to help option. credit to unutbu: http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
    if len(sys.argv) < 3:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    return args

def overlap(f1, f2):

    assert(f1[0][0] <= f2[0][0])

    return not (f2[0][0] > f1[-1][0])

# next function that returns 0 instead of raising StopIteration
# this is convenient for iterating over file 2 at a time
def safenext(iterator):
    try:
        nextval = next(iterator)
        return nextval
    except StopIteration:
        return 0

# determine if two fragments are consistent with each other
def consistent(f1, f2, i,j, contam_bounds, t):

    answer = True
    assert(f1[0][0] <= f2[0][0])


    f2iter = iter(f2) # iterator for fragment2
    a2 = safenext(f2iter)
    if not a2:
        return answer, contam_bounds

    mismatches = 0
    done  = False
    first_pos = True
    switched = False
    possible_switch_pos = None

    for a1 in f1:

        if done:
            break

        while a1[0] > a2[0]:
            a2 = safenext(f2iter)
            if not a2:
                done = True
                break

        if done or a1[0] < a2[0]:
            continue

        assert a1[0] == a2[0]
        if first_pos:
            switched = (a1[1] != a2[1])
            first_pos = False

        if (a1[1] == a2[1] and not switched) or (a1[1] != a2[1] and switched):
            mismatches = 0
            possible_switch_pos = None
        else:
            if possible_switch_pos == None:
                possible_switch_pos = a1[0]
            mismatches += 1
            if mismatches >= t:
                contam_bounds[(i,j)].add(possible_switch_pos)
                answer = False
                switched = not switched
                mismatches = 0
                possible_switch_pos = None

        a2 = safenext(f2iter)
        if not a2:
            break

    contam_bounds[(j,i)] = contam_bounds[(i,j)]

    return answer, contam_bounds


def read_fragment_matrix(frag_matrix):

    names = []
    flist = []

    with open(frag_matrix,"r") as fm:
        for line in fm:
            if len(line) < 2:
                continue

            el = line.strip().split()

            num_blks      = int(el[0])
            name = el[1]

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
            #qlist = [10**((ord(q) - 33) * -0.1) for q in qlist]

            alist= [(a,b,c) for ((a,b),c) in zip(call_list2,qlist)]

            flist.append(alist)
            names.append(name)

    zipped = zip(names,flist)
    sorted_zipped = sorted(zipped,key=lambda x: x[1][0][0])
    sorted_names, sorted_flist = [list(z) for z in zip(*sorted_zipped)]

    return sorted_flist, sorted_names

def matrixify_flist(flist,names, outfile):

    max_ix = 0

    for f in flist:

        for a in f:

            if a[0] > max_ix:

                max_ix = a[0]

    max_name = 0

    for name in names:
        if len(name) > max_name:
            max_name = len(name)

    with open(outfile,'w') as o:
        for f,name in zip(flist,names):

            line = [name.ljust(max_name+1)]+['-'] * max_ix
            for a in f:
                line[a[0]] = a[1]

            pline = ''.join(line)

            print(pline,file=o)


def fix_chamber_contamination(fragments,outfile, threshold=4, filter_cov1=False):

    contam_bounds = defaultdict(set)

    # read fragments into fragment data structure
    flist, names = read_fragment_matrix(fragments)

    if DEBUG:
        matrixify_flist(flist,names, 'sim_data/pretty_fragmatrix')

    N = len(flist)
    assert(N == len(names))

    cov_counts = defaultdict(int)

    for f in flist:

        for a in f:

            cov_counts[a[0]] += 1

    # filter out positions with coverage 1
    if filter_cov1:
        
        filtered_flist = []
        filtered_names = []        
        
        # for each fragment and name
        for f,n in zip(flist,names):
            
            new_f = []
            
            for a in f:
                
                # if coverage greater than 1
                if cov_counts[a[0]] > 1:
                    
                    new_f.append(a)
                    
            # if the fragment has at least 2 alleles after filtering add it to new list
            if len(new_f) >= 2:
                
                filtered_flist.append(new_f)
                filtered_names.append(n)
                
        # overwrite the old fragments and names
        flist = filtered_flist
        names = filtered_names

    #for k in sorted(list(cov_counts.keys())):
    #    print("{}\t{}".format(k,cov_counts[k]))

    # create a graph representing consistency of reads
    #  0: reads do not overlap
    # -1: reads are inconsistent
    #  1: reads are consistent
    C = np.zeros((N,N),dtype='int')

    for i in range(N):
        for j in range(i+1,N):

            f1 = flist[i]
            f2 = flist[j]

            if not overlap(f1,f2):
                continue

            cons,contam_bounds = consistent(f1,f2,i,j,contam_bounds, threshold)

            if cons:
                C[i,j] = 1
                C[j,i] = 1

            else:
                C[i,j] = -1
                C[j,i] = -1

    break_list = [set() for x in range(N)]

    bad = np.sum((C == -1),1) # find how many inconsistents each fragment has

    for i in range(N):

        if bad[i] == 0:
            continue

        f1 = flist[i]
        name = names[i]

        # consider repairing f1.
        # if f1 is inconsistent with multiple reads then we just split f1 where inconsistent.
        # if f1 is inconsistent with one read, then we split both reads at the inconsistent location.


        blist = []
        for j in range(N):
            if C[i,j] != -1:
                continue
            blist.append(j)

        for b1, b2 in itertools.combinations(blist,2):

            common_error = set.intersection(contam_bounds[(i,b1)],contam_bounds[(i,b2)])
            if common_error != set():
                break_list[i] = break_list[i].union(common_error)
                bad[i] -= len(common_error)
                for bx in blist:
                    bad[bx] -= len(common_error)
                    contam_bounds[(i,bx)] -= common_error

                if bad[i] == 0:
                    break

    new_flist = []
    new_names = []

    for i in range(N):

        f1 = flist[i]
        name = names[i]

        # consider repairing f1.
        if bad[i] > 0:
            for j in range(N):

                if C[i,j] != -1:
                    continue
                if contam_bounds[(i,j)] != set(): # inconsistencies weren't previously fixed
                    break_list[i] = contam_bounds[(i,j)]

        if break_list[i] == set():
            # add unedited fragment to new fragment list
            new_flist.append(f1)
            new_names.append(name)
        else:
            # chop up fragment where necessary and add pieces to fragment list
            new_f = []
            name_ctr = 1
            for allele in f1:

                if allele[0] in break_list[i]:
                    if DEBUG:
                        print("Breaking {} at pos {}, with coverage {}".format(name,allele[0], cov_counts[allele[0]]))
                    else:
                        print("Breaking {} at pos {}".format(name,allele[0]))
                    if len(new_f) > 1:
                        new_flist.append(new_f)
                        new_names.append("{}:S{}".format(name,name_ctr))
                        name_ctr += 1

                    new_f = []

                new_f.append(allele)

            if len(new_f) > 1:
                new_flist.append(new_f)
                new_names.append("{}:S{}".format(name,name_ctr))

    # WRITE TO FILE

    lines = []

    for name, fs in zip(new_names, new_flist):
        if len(fs) < 2:
            continue

        fragstr = ''
        num_pairs = 0
        prev_snp_ix = -2
        qual = ' '
        firstpos = fs[0][0]
        for snp_ix, allele, q_char in fs:

            diff = snp_ix - prev_snp_ix

            if diff == 1:
                fragstr += allele
            else:
                num_pairs += 1
                fragstr += ' {} {}'.format(snp_ix+1, allele)

            prev_snp_ix = snp_ix
            qual += q_char

        fragstr += qual

        prefix = '{} {}'.format(num_pairs,name)
        fragstr = prefix + fragstr

        lines.append((firstpos, fragstr))

    lines.sort()

    with open(outfile, 'w') as opf:
        for firstpos, line in lines:
            print(line, file=opf)

    if DEBUG:
        matrixify_flist(new_flist,new_names, 'sim_data/pretty_fragmatrix_fixed')


if __name__ == '__main__':
    args = parse_args()
    fix_chamber_contamination(args.fragments,args.output, args.threshold, args.filter_cov1)
