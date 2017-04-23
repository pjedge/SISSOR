# -*- coding: utf-8 -*-
"""
Created on Thu May 19 21:44:48 2016

@author: peter
"""
import sys
#import fileIO
import fragment
import sys
import fileIO
from math import log10
from collections import defaultdict
from itertools import combinations
from copy import copy

# OBJECTIVE:
# pair haplotype strands and determine three classes
# SPSC = SAME PARENT SAME CELL
# DPSC = DIFFERENT PARENT SAME CELL
# SPDC = SAME PARENT DIFFERENT CELL
# DPDC = DIFFERENT PARENT DIFFERENT CELL

# we want to read in a haplotype fragment file in HAPCUT2 format, find overlaps,
# and filter based on number of het SNVs overlapped and the consistency
# perhaps find a set of SNVs of largest span that has no switch errors?

# output a bed file with chromosoe, start, stop of region (based on het SNV span)
# and the cells and chambers that contain the paired fragments
# might as well also give the fragment IDs too.

# This bed file will allow us to parse a CCF file and mark pairs of chambers
# based on the verification they've undergone, i.e.
# SPSC.1.2.15 = same parent same cell, match, to cell 2, chamber 15
# SPSC.0.2.15 = same parent same cell, MISMATCH, to cell 2, chamber 15


VERBOSE = False
DEBUG = False
ID_LEN_NO_MODS = 5

tinylog = -1e5
n_chambers = 24
n_cells = 3

THRESHOLD = 0.8 # fraction of match to say haplotype is the same
filter_inconsistent_haplotypes = False
OVERLAP2 = 1


chroms = ['chr{}'.format(i) for i in range(1,23)] + ['chrX','chrY']

chr_num = dict()
for i,chrom in enumerate(chroms):
    chr_num[chrom] = i

# next function that returns 0 instead of raising StopIteration
# this is convenient for iterating over file 2 at a time
def safenext(iterator):
    try:
        nextval = next(iterator)
        return nextval
    except StopIteration:
        return 0

def addlogs(a,b):
    if a > b:
        return a + log10(1 + pow(10, b - a))
    else:
        return b + log10(1 + pow(10, a - b))

def overlap(f1, f2, haplotype_dict=None):

    assert(f1.seq[0][0] <= f2.seq[0][0])
    assert(f1.seq[0][1] <= f2.seq[0][1])


    el1 = f1.name.split(':')
    start1, end1 = [int(x) for x in el1[-1].split('-')]

    el2 = f2.name.split(':')
    start2, end2 = [int(x) for x in el2[-1].split('-')]

    if end1 - start2 < 3: #start2 >= end1:
        return None, None, None

    end = min([end1,end2])
    start = start2

    #s1 = [a for a in f1.seq if a[1] > start and a[1] < end]
    #s2 = [a for a in f2.seq if a[1] > start and a[1] < end]

    match1 = 0
    total1 = 0
    match2 = 0
    total2 = 0
    for (snp_ix, gen_ix, a, q) in f1.seq: #s1:
        if haplotype_dict[gen_ix] not in ['0','1'] or a not in ['0','1']:
            continue
        if haplotype_dict[gen_ix] == a:
            match1 += 1
        total1 += 1

    for (snp_ix, gen_ix, a, q) in f2.seq: #s2:
        if haplotype_dict[gen_ix] not in ['0','1'] or a not in ['0','1']:
            continue
        if haplotype_dict[gen_ix] == a:
            match2 += 1
        total2 += 1

    #if total1 != 0 and total2 != 0:
    #    print(total1,end=' ')
    #    print(total2)

    if total1 < OVERLAP2 or total2 < OVERLAP2:
        return None, None, None

    if (match1/total1 > THRESHOLD and match2/total2 > THRESHOLD):

        return start, end, 'P1'

    elif (1-(match1/total1) > THRESHOLD and 1-(match2/total2) > THRESHOLD):

        return start, end, 'P2'

    else:

        return None, None, None



cell_map = {'PGP1_21':0,'PGP1_22':1,'PGP1_A1':2}

def assign_fragments(flist, outputfile, haplotype_file=None):#hapblocks):

    total = 0
    paired_fragments = []

    t_blocklist = fileIO.parse_hapblock_file(haplotype_file,use_SNP_index=False)

    for t_block in t_blocklist:
         # convert t_block to a dict for convenience

        haplotype_dict = defaultdict(lambda: '-')
        for i, a1, a2 in t_block:
            haplotype_dict[i] = a1 # index by genomic pos (1-indexed), return one haplotype

        for f1,f2 in combinations(flist, 2):


            if f1.seq[0][0] > f2.seq[0][0]:

                temp = f1
                f1 = f2
                f2 = temp

            assert(f1.seq[0][0] <= f2.seq[0][0])
            assert(f1.seq[0][1] <= f2.seq[0][1])

            start_SNP, end_SNP, parent = overlap(f1,f2, haplotype_dict)

            if start_SNP == None:
                continue

            el1 = f1.name.split(':')
            chrom = el1[0]
            start1, end1 = [int(x) for x in el1[-1].split('-')]
            cell1 = cell_map[el1[2]]
            chamber1 = el1[3]
            assert(chamber1[0:2] == 'CH')
            chamber1 = int(chamber1[2:]) - 1

            el2 = f2.name.split(':')
            start2, end2 = [int(x) for x in el2[-1].split('-')]
            cell2 = cell_map[el2[2]]
            chamber2 = el2[3]
            assert(chamber2[0:2] == 'CH')
            chamber2 = int(chamber2[2:]) - 1

            end = min([end1,end2])
            start = start2

            total += end - start

                # order the pairs in increasing cell, chamber
            if not (cell1 < cell2 or (cell1 == cell2 and chamber1 < chamber2)):
                temp = (cell1,chamber1)
                (cell1, chamber1) = (cell2, chamber2)
                (cell2, chamber2) = temp

            paired_fragments.append((chrom, start, end, cell1, chamber1, cell2, chamber2, parent))

    paired_fragments.sort(key=lambda x: (chr_num[x[0]],x[1]))

    with open(outputfile, 'w') as opf:

        for (chrom, start, end, cell1, chamber1, cell2, chamber2, parent) in paired_fragments:
            print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(chrom, start, end, cell1, chamber1, cell2, chamber2, parent), file=opf)

    print("TOTAL          : {}".format(total))


def annotate_paired_strands(chamber_call_file, fragment_assignment_file, output_file):

    # read in file with spans of paired strands with matching haplotypes
    paired_fragments = [] # list of (start, end, ix1, ix2) tuples
    with open(fragment_assignment_file,'r') as infile:
        for line in infile:
            el = line.strip().split()
            chrom = el[0]
            start = int(el[1])
            end   = int(el[2])
            cell1 = int(el[3])
            chamber1 = int(el[4])
            cell2 = int(el[5])
            chamber2 = int(el[6])
            parent = el[7]

            paired_fragments.append((chrom, start, end, cell1, chamber1, cell2, chamber2, parent))

    paired_fragments.sort(key=lambda x: (chr_num[x[0]],x[1]))

    current_paired_fragments = []

    with open(chamber_call_file,'r') as ccf, open(output_file,'w') as opf:
        #print('chr\tpos\tsissor_call\tCGI_allele\tref',file=mof)

        for line in ccf:
            ccf_line = line.strip().split('\t')
            ccf_chrom = ccf_line[0]
            ccf_pos   = int(ccf_line[1])

            # pop info for haplotype-paired regions into the "current" list
            while paired_fragments != [] and (chr_num[paired_fragments[0][0]] < chr_num[ccf_chrom] or (paired_fragments[0][0] == ccf_chrom and paired_fragments[0][1] <= ccf_pos)):
                current_paired_fragments.append(paired_fragments.pop(0))

            # filter out paired fragments that we're past
            criteria = lambda pf: (pf[0] == ccf_chrom and pf[1] <= ccf_pos and ccf_pos <= pf[2])
            current_paired_fragments = list(filter(criteria,current_paired_fragments))

            tags = ccf_line[80].split(';')
            new_tags = []

            for chrom, start, end, cell1, chamber1, cell2, chamber2, parent in current_paired_fragments:

                assert(cell1 < cell2 or (cell1 == cell2 and chamber1 < chamber2))

                new_tag = 'HP:{}:{},{}:{},{}'.format(parent,cell1,chamber1,cell2,chamber2)
                new_tags.append(new_tag)

            if tags == ['N/A'] and new_tags != []:
                tags = new_tags
            else:
                tags += new_tags

            ccf_line[80] = ';'.join(tags)
            new_line = '\t'.join(ccf_line)
            print(new_line, file=opf)

def pair_strands(fragmentfile, vcf_file, outputfile, haplotype_file):

    # READ HAPLOTYPE BLOCKS

    # READ FRAGMENT MATRIX
    flist = fragment.read_fragment_matrix(fragmentfile,vcf_file)

    # ASSIGN FRAGMENTS TO HAPLOTYPES
    assign_fragments(flist, outputfile, haplotype_file)

def count_not_matchable(fragmentfile, vcf_file, haplotype_file):

    flist = fragment.read_fragment_matrix(fragmentfile, vcf_file)
     # convert t_block to a dict for convenience
    haplotype_dict = None
    if haplotype_file != None:
        haplotype_dict = defaultdict(lambda: '-')
        t_blocklist = fileIO.parse_hapblock_file(haplotype_file,use_SNP_index=False)
        for t_block in t_blocklist:
            for i, a1, a2 in t_block:
                haplotype_dict[i] = a1 # index by genomic pos (1-indexed), return one haplotype

    unmatchable = defaultdict(int)
    for f1 in flist:
        el1 = f1.name.split(':')
        start1, end1 = [int(x) for x in el1[-1].split('-')]
        cell = el1[2]
        match1 = 0
        total1 = 0

        for (snp_ix, gen_ix, a, q) in f1.seq:  # s1:
            if haplotype_dict[gen_ix] not in ['0', '1'] or a not in ['0', '1']:
                continue
            if haplotype_dict[gen_ix] == a:
                match1 += 1
            total1 += 1

        if total1 < 1:

            unmatchable[(cell,'too_short')] += end1 - start1


        elif not (match1/total1 > THRESHOLD or 1-match1/total1 > THRESHOLD):

            unmatchable[(cell,'mismatch')] += end1 - start1

    return unmatchable


def test_pair_strands():
    #chroms = ['chr{}'.format(x) for x in range(1,23)]
    chroms = ['chr20']

    TOTAL = 0

    for chrom in chroms:
        fragmentfile = '../haplotyping/sissor_project/data/PGP1_ALL/fragmat/cov1_basic/{}'.format(chrom)
        vcf_file = '../haplotyping/sissor_project/data/PGP1_VCFs_BACindex/{}.vcf'.format(chrom)
        outputfile = 'strand_pairs_new/{}'.format(chrom)
        haplotype_file = '../haplotyping/sissor_project/experiments/hapcut2_PGP1_ALL/cov1_basic/{}.output'.format(chrom)
        total = pair_strands(fragmentfile, vcf_file, outputfile, haplotype_file)
        TOTAL += total

    print(TOTAL)

if __name__ == '__main__':
    test_pair_strands()
    #pair_strands()
