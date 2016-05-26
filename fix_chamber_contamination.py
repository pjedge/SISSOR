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

def parse_args():

    parser = argparse.ArgumentParser(description='Split SISSOR fragments with switch errors')
    # paths to samfiles mapping the same ordered set of RNA reads to different genomes
    parser.add_argument('-f', '--fragments', nargs='?', type = str, help='fragment file to process')
    parser.add_argument('-o', '--output', nargs='?', type = str, help='output fragments')
    parser.add_argument('-t', '--threshold', nargs='?', type = float, help='number of consecutive mismatches with respect to another fragment that count as contamination', default = 4)

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
    
def check_same_hap(f1,f2):
    count1 = 0
    count2 = 0
    # for allele in fragment1
    f2iter = iter(f2) # iterator for fragment2
    a2 = safenext(f2iter)
    if not a2:
        return (count1 > count2)
    
    for a1 in f1:

        while a1[0] > a2[0]:
            a2 = safenext(f2iter)
            if not a2:
                return (count1 > count2)
                
        if a1[0] < a2[0]:
            continue
        assert a1[0] == a2[0]
        
        if a1[1] == a2[1]:
            count1 += 1
        else:
            count2 += 1
        
        a2 = safenext(f2iter)
        if not a2:
            return (count1 > count2)

    return (count1 > count2)    
    
# determine if two fragments are consistent with each other
def consistent(f1, f2, i,j, contam_bounds, t):
    
    answer = True
    assert(f1[0][0] <= f2[0][0])    
    
    same_hap = check_same_hap(f1,f2)     
    
    f2iter = iter(f2) # iterator for fragment2
    a2 = safenext(f2iter)    
    if not a2:
        return answer, contam_bounds

    mismatches = 0
    start = -1
    end   = -1
    done  = False
    
    for a1 in f1:
        
        if done:
            break
        
        while a1[0] > a2[0]:
            a2 = safenext(f2iter)
            if not a2:
                done = True
                break
        
        if a1[0] < a2[0] or done:
            continue
        assert a1[0] == a2[0]
        
        if (same_hap and a1[1] == a2[1]) or (not same_hap and a1[1] != a2[1]):
            mismatches = 0
            if start != -1 and mismatches > t:
                end = a1[0]
                contam_bounds[(i,j)].add(start)
                contam_bounds[(i,j)].add(end)
                contam_bounds[(j,i)] = contam_bounds[(i,j)]

                start = -1
                end   = -1
        else:
            if start == -1:
                start = a1[0] # need to break fragment before "start"
            mismatches += 1
            if mismatches > t:
                answer = False
        
        a2 = safenext(f2iter)
        if not a2:
            break
            
    if start != -1 and mismatches > t:
        contam_bounds[(i,j)].add(start)
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
    
    return flist, names

def print_coverage_counts(flist):
    
    counts = defaultdict(int)
    
    for f in flist:
        
        for a in f:
            
            counts[a[0]] += 1
            
    print("POS\tCOVERAGE")
    
    for k in sorted(list(counts.keys())):
        
        print("{}\t{}".format(k,counts[k]))
        
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
            
            
def fix_chamber_contamination(fragments,outfile, threshold=4):
    
    contam_bounds = defaultdict(set)
    
    # read fragments into fragment data structure
    flist, names = read_fragment_matrix(fragments)

    #print_coverage_counts(flist)

    #matrixify_flist(flist,names, 'sim_data/pretty_fragmatrix')

    assert(len(flist) == len(names))
    N = len(flist)
    
    # create a graph representing consistency of reads
    #  0: reads do not overlap
    # -1: reads are inconsistent
    #  1: reads are consistent
    C = np.zeros((N,N),dtype='int')
    
    for i in range(len(flist)):
        for j in range(i+1,len(flist)):
            
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

    new_flist = []
    new_names = []
    
    bad = np.sum((C == -1),1) # find how many inconsistents each fragment has
    bad_ix = np.argsort(bad)[::-1]  # sort fragments by number of inconsistenciess
    
    mod_dict = defaultdict(list)
    
    while(sum(bad != 0)):
        for u in range(len(flist)):
            
            i = bad_ix[u]
            f1 = flist[i]
            name = names[i]
            # consider repairing f1.
            # if f1 is inconsistent with multiple reads then we just split f1 where inconsistent.
            # if f1 is inconsistent with one read, then we split both reads at the inconsistent location.
            
            break_pos = set()
            if bad[i] == 1:
                                
                for j in range(len(flist)):
    
                    if C[i,j] != -1:
                        continue
                    if contam_bounds[(i,j)] == set(): # inconsistencies were previously fixed
                        bad[i] = 0
                    else:
                        break_pos = contam_bounds[(i,j)]
                                                  
                    break
                    
            elif bad[i] > 1:
                
                blist = []
                for j in range(len(flist)):
                    if C[i,j] != -1:
                        continue
                    blist.append(j)
                   
                for b1, b2 in itertools.combinations(blist,2):
                    
                    common_error = set.intersection(contam_bounds[(i,b1)],contam_bounds[(i,b2)])
                    if common_error != set():
                        break_pos = break_pos.union(common_error)
                        for bx in blist:
                            contam_bounds[(i,bx)] -= common_error
                        
            if bad[i] == 0:
                # add unedited fragment to new fragment list
                pass
                #new_flist.append(f1)
                #new_names.append(name)
            else:
                # chop up fragment where necessary and add pieces to fragment list
                new_f = []
                name_ctr = 1
                for allele in f1:
                    if allele[0] in break_pos:
                        
                        print("Breaking {} at pos {}".format(name,allele[0]))
    
                        if len(new_f) > 1:
                            mod_dict[name].append(("{}:PART{}".format(name,name_ctr),new_f))
                            #new_flist.append(new_f)
                            #new_names.append("{}:PART{}".format(name,name_ctr))
                            name_ctr += 1
                            
                        new_f = []
                    else:
                        new_f.append(allele)
                        
                    #new_flist.append(new_f)
                    #new_names.append("{}:PART{}".format(name,name_ctr))
                    mod_dict[name].append(("{}:PART{}".format(name,name_ctr),new_f))
    
    
    new_names = []
    new_flist = []
    for name,f in zip(names,flist):
        if name in mod_dict:
            n2, f2 = mod_dict[name]
            new_names += n2
            new_flist += f2
        else:
            new_names.append(name)
            new_flist.append(f)
    
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
            
        
if __name__ == '__main__':
    args = parse_args()
    fix_chamber_contamination(args.fragments,args.output, args.threshold)