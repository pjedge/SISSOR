# -*- coding: utf-8 -*-
"""
Created on Thu May 19 21:44:48 2016

@author: peter
"""

import fragment
import argparse
import sys
from collections import defaultdict
import itertools

VERBOSE = False
DEBUG = False
ID_LEN_NO_MODS = 5

def parse_args():

    parser = argparse.ArgumentParser(description='Split SISSOR fragments with switch errors')
    # paths to samfiles mapping the same ordered set of RNA reads to different genomes
    parser.add_argument('-f', '--fragments', nargs='?', type = str, help='fragment file to process')
    parser.add_argument('-o', '--output', nargs='?', type = str, help='output fragments')
    parser.add_argument('-t', '--threshold', nargs='?', type = int, help='number of consecutive mismatches with respect to another fragment that count as contamination', default = 4)
    parser.add_argument('-f', '--filter', nargs='?', type = int, help='filter for only positions with coverage >= this amount', default = 0)


    # default to help option. credit to unutbu: http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
    if len(sys.argv) < 3:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    return args

def overlap(f1, f2, amt=1):

    assert(f1.seq[0][0] <= f2.seq[0][0])

    s1 = set([a[0] for a in f1.seq])
    s2 = set([a[0] for a in f2.seq])
    inter = set.intersection(s1,s2)

    return len(inter) >= amt


def total_SER(f1, f2):

    assert(f1.seq[0][0] <= f2.seq[0][0])

    f2iter = iter(f2.seq) # iterator for fragment2
    a2 = safenext(f2iter)
    if not a2:
        return 0

    switches = 0
    total_pos = 0
    done  = False
    first_pos = True
    switched = False

    for a1 in f1.seq:

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
            switched = (a1[2] != a2[2])
            first_pos = False

        total_pos += 1
        if (a1[2] != a2[2] and not switched) or (a1[2] == a2[2] and switched):
            switches += 1
            switched = not switched

        a2 = safenext(f2iter)
        if not a2:
            break

    return switches, total_pos

def get_het_breakpoints(frag,min_het_count=3,max_het_frac=0.25):

    breaks = []

    for i in range(len(frag.seq)):
        if frag.seq[i][2] != 'M':
            continue

        hom_count = 0
        het_count = 0

        for j in range(i,len(frag.seq)):
            if frag.seq[j][2] == 'M':
                het_count += 1

                het_frac = het_count/(hom_count+het_count)

                if het_count >= min_het_count and het_frac >= max_het_frac:
                    breaks.append((i,j))

            else:
                hom_count += 1

    breakpos = set()
    for i,j in breaks:
        for k in range(i,j+1):
            breakpos.add(frag.seq[k][0])

    return breakpos

def split_fragment_hets(frag,min_het_count=3,max_het_frac=0.25):

    breakpos = get_het_breakpoints(frag,min_het_count,max_het_frac)

    new_fragment_pieces = []
    new_piece = []
    split_occured = False
    for pos, gpos, call, qual in frag.seq:
        if pos in breakpos:
            split_occured = True
            if len(new_piece) >= 2:
                new_fragment_pieces.append(new_piece)
            new_piece = []
        else:
            if call != 'M':
                new_piece.append((pos,gpos,call,qual))

    if len(new_piece) >= 2:
        new_fragment_pieces.append(new_piece)

    new_piece = []
    if new_fragment_pieces == []:
        return []

    new_fragment_piece_objects = []
    if split_occured:
        for i, piece in enumerate(new_fragment_pieces):
            new_name = '{}:H{}'.format(frag.name,i+1)
            new_fragment_piece_objects.append(fragment.fragment(piece,new_name))
    else:
        new_fragment_piece_objects.append(fragment.fragment(new_fragment_pieces[0],frag.name))
    return new_fragment_piece_objects

def split_flist_hets(flist,min_het_count=3,max_het_frac=0.25):

    new_flist = []

    for f in flist:
        new_flist += split_fragment_hets(f,min_het_count,max_het_frac)

    return new_flist


def filter_flist_heterozygousity(flist,max_het=0.1):

    new_flist = []
    for f in flist:
        hom_count = 0
        het_count = 0

        for pos, gpos, call, qual in f.seq:
            if call == 'M':
                het_count += 1
            elif call == '1' or call == '0':
                hom_count += 1

        het_frac = het_count / (hom_count + het_count)

        if het_frac < max_het:
            new_flist.append(f)

    return new_flist

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
    assert(f1.seq[0][0] <= f2.seq[0][0])

    f2iter = iter(f2.seq) # iterator for fragment2
    a2 = safenext(f2iter)
    if not a2:
        return answer, contam_bounds

    mismatches = 0
    done  = False
    first_pos = True
    switched = False
    possible_switch_pos = None
    last_pos = None

    for a1 in f1.seq:

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
            switched = (a1[2] != a2[2])
            first_pos = False

        if (a1[2] == a2[2] and not switched) or (a1[2] != a2[2] and switched):
            mismatches = 0
            possible_switch_pos = None
        else:
            if possible_switch_pos == None:
                possible_switch_pos = set(range(last_pos+1,a1[0]+1))
            mismatches += 1
            if mismatches >= t:
                contam_bounds[(i,j)] = set.union(contam_bounds[(i,j)], possible_switch_pos)
                answer = False
                switched = not switched
                mismatches = 0
                possible_switch_pos = None

        last_pos = a1[0]

        a2 = safenext(f2iter)
        if not a2:
            break

    contam_bounds[(j,i)] = contam_bounds[(i,j)]

    return answer, contam_bounds

def fragment_comparison_split(flist, threshold=3):

    flist = sorted(flist,key=lambda x: x.seq[0][0])

    contam_bounds = defaultdict(set)

    N = len(flist)

    #for k in sorted(list(cov_counts.keys())):
    #    print("{}\t{}".format(k,cov_counts[k]))

    for i in range(N):
        for j in range(i+1,N):

            f1 = flist[i]
            f2 = flist[j]

            if not overlap(f1,f2):
                continue

            cons,contam_bounds = consistent(f1,f2,i,j,contam_bounds, threshold)

    break_list = [set() for x in range(N)]

    for i in range(N):

        f1 = flist[i]

        # consider repairing f1.
        # if f1 is inconsistent with multiple reads then we just split f1 where inconsistent.
        # if f1 is inconsistent with one read, then we split both reads at the inconsistent location.


        blist = []
        for j in range(N):
            if (i,j) not in contam_bounds or contam_bounds[(i,j)] == set():
                continue
            blist.append(j)

        for b1, b2 in itertools.combinations(blist,2):

            common_error = set.intersection(contam_bounds[(i,b1)], contam_bounds[(i,b2)])

            if common_error != set():
                break_list[i] = break_list[i].union(common_error)
                for bx in blist:
                    inter = set.intersection(common_error, contam_bounds[(i,bx)])
                    if inter != set():
                        rightward = next(iter(inter))
                        leftward = rightward-1
                        while rightward in contam_bounds[(i,bx)]:
                            contam_bounds[(i,bx)].remove(rightward)
                            rightward += 1
                        while leftward in contam_bounds[(i,bx)]:
                            contam_bounds[(i,bx)].remove(leftward)
                            leftward -= 1

    new_flist = []

    for i in range(N):

        f1 = flist[i]

        # consider repairing f1.
        for j in range(N):

            if (i,j) in contam_bounds and contam_bounds[(i,j)] != set(): # inconsistencies weren't previously fixed
                break_list[i] = set.union(break_list[i],contam_bounds[(i,j)])

        if break_list[i] == set():
            # add unedited fragment to new fragment list
            new_flist.append(f1)
        else:
            # chop up fragment where necessary and add pieces to fragment list
            new_seq = []
            name_ctr = 1
            for allele in f1.seq:

                if allele[0] in break_list[i]:
                    if VERBOSE:
                        print("Breaking {} at pos {}".format(f1.name,allele[0]))
                    if len(new_seq) > 1:
                        new_name = "{}:S{}".format(f1.name,name_ctr)
                        new_flist.append(fragment.fragment(new_seq,new_name))
                        name_ctr += 1

                    new_seq = []

                new_seq.append(allele)

            if len(new_seq) > 1:
                new_name = "{}:S{}".format(f1.name,name_ctr)
                new_flist.append(fragment.fragment(new_seq,new_name))

    return new_flist

def filter_discordant_fragments_SER(flist, SER_threshold):

    flist = sorted(flist,key=lambda x: x.seq[0][0])

    N = len(flist)
    #SER_total = np.zeros((N,N),dtype='float')
    #SER_err = np.zeros((N,N),dtype='float')
    new_flist = []

    for i in range(N):
        f1 = flist[i]
#        total = 0
#        err   = 0
        errs = []
        e_total = 0
        t_total = 0
        for j in range(i+1,N):

            f2 = flist[j]

            if not overlap(f1,f2,2):
                #SER_total[i,j] = -1
                #SER_total[j,i] = -1
                #SER_err[i,j] = -1
                #SER_err[j,i] = -1
                continue

            e, t = total_SER(f1,f2)
            e_total += e
            t_total += t
            #SER_total[i,j] = t
            #SER_total[j,i] = t
            #SER_err[i,j] = e
            #SER_err[j,i] = e
#            total += t
#            err   += e
            errs.append(e/t)

        if t_total <= 10 or e_total / t_total < SER_threshold: #len(errs) <= 1 or min(errs) < SER_threshold:
            new_flist.append(f1)

    return new_flist

def filter_fragment_coverage(flist,coverage):

    cov_counts = defaultdict(int)

    for f in flist:

        for a in f.seq:

            cov_counts[a[0]] += 1


    filtered_flist = []

    # for each fragment and name
    for f in flist:

        new_seq = []

        for a in f.seq:

            # if coverage greater than 1
            if cov_counts[a[0]] >= coverage:

                new_seq.append(a)

        # if the fragment has at least 2 alleles after filtering add it to new list
        if len(new_seq) >= 2:

            filtered_flist.append(fragment.fragment(new_seq,f.name))

    return filtered_flist

def add_new_boundaries(flist):

    for i in range(len(flist)):
        flist[i].name = '{}:{}-{}'.format(flist[i].name,flist[i].seq[0][1]+1,flist[i].seq[-1][1]+1)

        for a in flist[i].seq:
            assert(a[2] != 'M') # one last check that we didn't leave in a heterozygous call

    return flist

def filter_het_positions(flist):

    filtered_flist = []

    # for each fragment and name
    for f in flist:

        new_seq = []

        for a in f.seq:

            # if coverage greater than 1
            if a[2] != 'M':

                new_seq.append(a)

        # if the fragment has at least 2 alleles after filtering add it to new list
        if len(new_seq) >= 2:

            filtered_flist.append(fragment.fragment(new_seq,f.name))

    return filtered_flist

def fix_chamber_contamination(fragmentfile, vcf_file, outfile, threshold=2, min_coverage=0, mode='strict'):

    # READ FRAGMENT MATRIX
    flist = fragment.read_fragment_matrix(fragmentfile,vcf_file)

    if mode == 'none':

        # SPLIT FRAGMENTS BASED ON FRAGMENT-FRAGMENT COMPARISON
        flist = filter_het_positions(flist)

    elif mode == 'basic':

        # SPLIT FRAGMENTS BASED ON FRAGMENT-FRAGMENT COMPARISON
        flist = filter_het_positions(flist)
        flist = fragment_comparison_split(flist,threshold)

    elif mode == 'strict':

        # FILTER OUT HIGHLY HETEROZYGOUS FRAGMENTS
        flist = filter_flist_heterozygousity(flist, max_het=0.05)

        # SPLIT FRAGMENTS AT HETEROZYGOUS SPOTS
        flist = split_flist_hets(flist, min_het_count=3,max_het_frac=0.25)

        # FILTER OUT FRAGMENTS THAT ARE HIGHLY DISCORDANT
        flist = filter_discordant_fragments_SER(flist, 0.3)

        # SPLIT FRAGMENTS BASED ON FRAGMENT-FRAGMENT COMPARISON
        flist = fragment_comparison_split(flist,threshold)
    else:
        print("INVALID MODE")
        sys.exit(1)

    # FILTER OUT FRAGMENTS BELOW MIN COVERAGE
    if min_coverage > 1:
        flist = filter_fragment_coverage(flist, min_coverage)

    flist = add_new_boundaries(flist)

    # WRITE TO FILE
    fragment.write_fragment_matrix(flist, outfile)

if __name__ == '__main__':
    args = parse_args()
    fix_chamber_contamination(args.fragments,args.output, args.threshold, args.filter)
