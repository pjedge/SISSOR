#! /bin/python3.5
# -*- coding: utf-8 -*-

# Author : Peter Edge
# Email  : pedge@eng.ucsd.edu
import argparse
import sys

desc = '''computes the mismatch and switch error rates between fragments in a
HapCUT format fragment file and a HapCUT hapblock file'''

# parse arguments
def parse_args():

    parser = argparse.ArgumentParser(description='')
    parser.add_argument('frag_file', nargs='?', type = str, help='path to fragment file')
    parser.add_argument('hapcut_block_file', nargs='?', type = str, help='path to hapcut block file')
    parser.add_argument('output_file', nargs='?', type = str, help='path to output_file')
    print(len(sys.argv))
    # default to help option. credit to unutbu: http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
    if len(sys.argv) < 4:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    return args

# hapcut_blks is a list of blocks, represented as lists of (pos, allele1, allele2)
# frags is a list of lists of (pos, allele)
def single_fragment_error_rate(hapcut_blks, frag, switch_threshold=10):

    switched       = False
    last_base_was_switch = False
    comparisons = 0
    blocknum = 0

    # build a dictionary for rapid indexing into the hapcut block list
    ix_dict = {}
    for i, blk in enumerate(hapcut_blks):
        for j, (pos, a1, a2) in enumerate(blk):
            ix_dict[pos] = (i,j) # for a given SNP index, save the block and index where we can find it

    # iterate over SNPs in the true and assembled haplotypes in parallel
    # i is the index of the current base. x is the current base in the true haplotype. y is the current base in the assembled haplotype.

    # figure out how many blocks are visited
    visited_blocks = set()
    for frag_pos, frag_allele in frag:
        if frag_pos not in ix_dict:
            continue
        i,j = ix_dict[frag_pos]
        if hapcut_blks[i][j][1] in ['0','1']:
            visited_blocks.add(i)

    blocknum = len(visited_blocks)

    # keep separate error counts for the haplotype and its complement
    # they can differ by 1 and we want to minimize switches
    switches   = [0,0]
    mismatches = [0,0]

    for a in [0,1]: # choose which allele to score

        prev_block_ix  = -1

        for frag_pos, frag_allele in frag:

            if frag_pos not in ix_dict:
                continue

            blk_ix, snpblk_ix = ix_dict[frag_pos]

            y = frag_allele
            x = hapcut_blks[blk_ix][snpblk_ix][1+a]

            if x == '-' or y == '-':
                continue # skip cases where fragment or known haplotype are missing data

            if a == 0:
                comparisons += 1

            if blk_ix != prev_block_ix: # this is the first SNP in a new block
                switched = (x != y)
                last_base_was_switch = switched
                prev_block_ix = blk_ix
                continue

            # if there is a mismatch against the true haplotype and we are in a normal state,
            # or if there is a "match" that isn't a match because we are in a switched state,
            # then we need to flip the state again and iterate the count
            if (x != y and not switched) or (x == y and switched): # current base is mismatched, implying a switch

                switched = not switched  # flip the "switched" status

                if last_base_was_switch: # if last base was a switch then this is actually a single-base mismatch
                    # count the 2 switches as a single-base mismatch instead
                    mismatches[a] += 1
                    switches[a] -= 1      # undo count from last base switch
                    print("Flip at {}".format(frag_pos+1))

                    if (switches[a] < 0):
                        switches[a] = 0
                    last_base_was_switch = False

                else:

                    switches[a] += 1
                    print("Switch at {}".format(frag_pos+1))

                    last_base_was_switch = True

            else: # current base is not mismatched
                last_base_was_switch = False

            prev_block_ix = blk_ix

        # special case for switch on last base of previous block; should count as a mismatch
        if last_base_was_switch:
            # count the switch as a single-base mismatch instead
            mismatches[a] += 1
            switches[a] -= 1

            print("Flip at {}".format(frag_pos+1))

            if (switches[a] < 0):
                switches[a] = 0

        print("*****************")
        
    if switches[0] < switches[1]:
        return (switches[0], mismatches[0], comparisons, blocknum)
    else:
        return (switches[1], mismatches[1], comparisons, blocknum)

# parses a hapcut block file into a useful data structure
# list of lists of tuples
# each inner list represents
def parse_hapblock_file(hapblock_file):

    blocklist = []

    try:
        with open(hapblock_file, 'r') as hbf:

            for line in hbf:
                if 'BLOCK' in line:
                    blocklist.append([])
                    continue

                elements = line.strip().split('\t')
                if len(elements) < 3:
                    continue

                pos = int(elements[0])-1
                allele1 = elements[1]
                allele2 = elements[2]

                blocklist[-1].append((pos, allele1, allele2))

    except FileNotFoundError:
        # most of the time, this should mean that the program timed out and therefore didn't produce a phase.
        pass

    return blocklist

# takes a split line from a fragment file and converts it to a list of (pos,allele) tuples
def fragdata2list(fragdata):

    fraglist = []
    frag_iter = iter(fragdata[2:-1])

    for ix_str in frag_iter:
        seq_str = next(frag_iter)
        pos = int(ix_str) - 1

        for i, allele in enumerate(seq_str):
            fraglist.append((pos+i, allele))

    return fraglist

# main function
def fragment_error_rate(frag_file, hapblock_file, output_file):

    hapcut_blks = parse_hapblock_file(hapblock_file)

    with open(frag_file, 'r') as infile:
        with open(output_file, 'w') as outfile:

            # output header
            print("FRAG_ID\tSWITCH_COUNT\tMISMATCH_COUNT\tCOMPARISONS\tBLOCKS_OVERLAPPED\n",file=outfile)

            # processes each fragment
            for line in infile:
                if len(line) < 2:
                    continue  # skip empty lines

                frag_data = line.strip().split()   # split fragment data into list
                frag = fragdata2list(frag_data)    # convert to more useful representation
                frag_name = frag_data[1]           # fragment id

                # compute error rate and print
                switch_count, mismatch_count, comparisons, block_num = single_fragment_error_rate(hapcut_blks, frag)
                print("{}\t{}\t{}\t{}\t{}".format(frag_name, switch_count, mismatch_count,comparisons, block_num),file=outfile)

# parse args and execute main function
if __name__ == '__main__':
    args = parse_args()
    fragment_error_rate(args.frag_file, args.hapcut_block_file, args.output_file)
