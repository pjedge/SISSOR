# -*- coding: utf-8 -*-
"""
Created on Thu May 19 21:44:48 2016

@author: peter
"""
import sys
sys.path.append('../../HapTools')
sys.path.append('../haplotyping')
#import fileIO
import fragment
import sys
from math import log10
from itertools import combinations

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

CUTOFF = 0.99

def overlap(f1, f2, amt=3):

    if f1.seq[0][0] > f2.seq[0][0]:

        temp = f1
        f1 = f2
        f2 = temp

    assert(f1.seq[0][0] <= f2.seq[0][0])

    s1 = set([a[0] for a in f1.seq])
    s2 = set([a[0] for a in f2.seq])

    inter = set.intersection(s1,s2)

    if len(inter) < amt:
        return None, None

    inter_seq1 = [x for x in f1.seq if x[0] in inter]
    inter_seq2 = [x for x in f2.seq if x[0] in inter]

    gpos = [x[1] for x in inter_seq1]
    start = min(gpos)
    end   = max(gpos)

    p_samehap = 0
    p_diffhap = 0
    total = 0
    matches = 0
    for (snp_ix1, gen_ix1, a1, q1), (snp_ix2, gen_ix2, a2, q2) in zip(inter_seq1, inter_seq2):


        q1 = 10**((ord(q1)-33)/-10)
        q2 = 10**((ord(q2)-33)/-10)

        total += 1
        if a1 == a2:
            matches += 1
            p_samehap += log10((1-q1)*(1-q2)+q1*q2)
            p_diffhap += log10((1-q1)*q2+(1-q2)*q1)
        else:
            p_diffhap += log10((1-q1)*(1-q2)+q1*q2)
            p_samehap += log10((1-q1)*q2+(1-q2)*q1)

    #post_samehap = 10**(p_samehap - addlogs(p_samehap, p_diffhap))
    #print("{} {} {}".format(matches, total, post_samehap))

    if matches/total > 0.8:#post_samehap > CUTOFF:

        return start, end

    else:
        return None, None

cell_map = {'PGP1_21':0,'PGP1_22':1,'PGP1_A1':2}

def assign_fragments(flist, outputfile):#hapblocks):

    total = 0

    with open(outputfile, 'w') as opf:
        for f1,f2 in combinations(flist, 2):

            start, end = overlap(f1,f2,4)

            if start == None:
                continue

            total += end - start

            el1 = f1.name.split(':')
            chrom = el1[0]
            cell1 = cell_map[el1[2]]
            chamber1 = el1[3]
            assert(chamber1[0:2] == 'CH')
            chamber1 = int(chamber1[2:])

            el2 = f2.name.split(':')
            cell2 = cell_map[el2[2]]
            chamber2 = el2[3]
            assert(chamber2[0:2] == 'CH')
            chamber2 = int(chamber2[2:])

            print('{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(chrom, start, end, cell1, chamber1, cell2, chamber2), file=opf)

    print(total)


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

            paired_fragments.append((chrom, start, end, cell1, chamber1, cell2, chamber2))

    paired_fragments.sort(key=lambda x: (chr_num[x[0]],x[1]))

    current_paired_fragments = []

    with open(chamber_call_file,'r') as ccf, open(output_file,'w') as opf:
        #print('chr\tpos\tsissor_call\tCGI_allele\tref',file=mof)

        snp_ix = -1
        prev_ccf_chrom = None
        for line in ccf:
            ccf_line = line.strip().split('\t')
            ccf_chrom = ccf_line[0]

            if ccf_chrom != prev_ccf_chrom:
                snp_ix = -1
            prev_ccf_chrom = ccf_chrom

            snp_ix += 1

            ccf_pos   = int(ccf_line[1])

            if paired_fragments != []:

                while(chr_num[paired_fragments[0][0]] <= chr_num[ccf_chrom] and paired_fragments[0][1] <= ccf_pos):
                    current_paired_fragments.append(paired_fragments.pop(0))

            # filter out paired fragments that we're past

            current_paired_fragments = [
            (chrom, start, end, cell1, chamber1, cell2, chamber2)
            for (chrom, start, end, cell1, chamber1, cell2, chamber2) in current_paired_fragments
            if (chrom == ccf_chrom and start <= ccf_pos and ccf_pos <= end)
            ]

            if ccf_chrom in ['chrX','chrY']:
                continue

            tags = ccf_line[80].split(';')
            new_tags = []

            for chrom, start, end, cell1, chamber1, cell2, chamber2 in current_paired_fragments:
                    # order the pairs in increasing cell, chamber
                if not (cell1 < cell2 or (cell1 == cell2 and chamber1 < chamber2)):
                    temp = (cell1,chamber1)
                    (cell1, chamber1) = (cell2, chamber2)
                    (cell2, chamber2) = temp

                new_tag = 'HP:{},{}:{},{}'.format(cell1,chamber1,cell2,chamber2)
                new_tags.append(new_tag)

            if tags == ['N/A'] and new_tags != []:
                tags = new_tags
            else:
                tags += new_tags

            ccf_line[80] = ';'.join(tags)
            new_line = '\t'.join(ccf_line)
            print(new_line, file=opf)

def pair_strands(fragmentfile, vcf_file, outputfile):

    # READ HAPLOTYPE BLOCKS

    # READ FRAGMENT MATRIX
    flist = fragment.read_fragment_matrix(fragmentfile,vcf_file)

    # ASSIGN FRAGMENTS TO HAPLOTYPES
    flist = assign_fragments(flist, outputfile)

if __name__ == '__main__':
    pair_strands()
