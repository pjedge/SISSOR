#! /bin/python3.5
# -*- coding: utf-8 -*-

# Author : Peter Edge
# Email  : pedge@eng.ucsd.edu

# hapcut_blks is a list of blocks, represented as lists of (pos, allele1, allele2)
# frags is a list of lists of (pos, allele)
def single_fragment_error_rate(hapcut_blks, frag):

    switched = False
    last_base_was_switch = False
    comparisons = 0
    blocknum = 0

    # build a dictionary for rapid indexing into the hapcut block list
    ix_dict = {}
    for i, blk in enumerate(hapcut_blks):
        for j, (chrom, ix, pos, a1, a2) in enumerate(blk):
            ix_dict[(chrom,ix)] = (i,j) # for a given SNP index, save the block number, and block index where we can find it

    # iterate over SNPs in the true and assembled haplotypes in parallel
    # i is the index of the current base. x is the current base in the true haplotype. y is the current base in the assembled haplotype.

    # figure out how many blocks are visited
    visited_blocks = set()
    for ix, pos, allele in frag.seq:
        if (frag.chrom,ix) not in ix_dict:
            continue
        i,j = ix_dict[(frag.chrom,ix)]
        if hapcut_blks[i][j][3] in ['0','1']:
            visited_blocks.add(i)

    blocknum = len(visited_blocks)

    # keep separate error counts for the haplotype and its complement
    # they can differ by 1 and we want to minimize switches
    switches   = [0,0]
    mismatches = [0,0]

    for a in [0,1]: # choose which allele to score

        prev_block_ix  = -1

        for frag_ix, frag_pos, frag_allele in frag.seq:

            if (frag.chrom,frag_ix) not in ix_dict:
                continue

            i, j = ix_dict[(frag.chrom,frag_ix)]

            y = frag_allele
            x = hapcut_blks[i][j][3+a]

            if x == '-' or y == '-':
                continue # skip cases where fragment or known haplotype are missing data

            if a == 0:
                comparisons += 1

            if i != prev_block_ix: # this is the first SNP in a new block
                switched = (x != y)
                last_base_was_switch = switched
                prev_block_ix = i
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

                    if (switches[a] < 0):
                        switches[a] = 0
                    last_base_was_switch = False

                else:

                    switches[a] += 1

                    last_base_was_switch = True

            else: # current base is not mismatched
                last_base_was_switch = False

            prev_block_ix = i

        # special case for switch on last base of previous block; should count as a mismatch
        if last_base_was_switch:
            # count the switch as a single-base mismatch instead
            mismatches[a] += 1
            switches[a] -= 1

            if (switches[a] < 0):
                switches[a] = 0

    if switches[0] < switches[1]:
        return (switches[0], mismatches[0], comparisons, blocknum)
    else:
        return (switches[1], mismatches[1], comparisons, blocknum)

# main function
def fragment_error_rate(frag_list, hapblock_list, output_file):

    with open(output_file, 'w') as outfile:

        # output header
        print("FRAG_ID\tSWITCH_COUNT\tMISMATCH_COUNT\tCOMPARISONS\tBLOCKS_OVERLAPPED\n",file=outfile)

        for frag in frag_list:

            # compute error rate and print
            switch_count, mismatch_count, comparisons, block_num = single_fragment_error_rate(hapblock_list, frag)

            print("{}\t{}\t{}\t{}\t{}".format(frag.id, switch_count, mismatch_count,comparisons, block_num),file=outfile)
