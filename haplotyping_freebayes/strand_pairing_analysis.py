# -*- coding: utf-8 -*-
"""
Created on Thu May 19 21:44:48 2016

@author: peter
"""
import sys
sys.path.append('/home/pedge/git/HapTools')
sys.path.append('/home/peter/git/HapTools')

#import fileIO
import fragment
import sys
import fileIO
from math import log10
from collections import defaultdict
from itertools import combinations
from copy import copy
from DJSF import Union, Find, Node #MakeSet
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

THRESHOLD = 0.9 # fraction of match to say haplotype is the same
filter_inconsistent_haplotypes = False
OVERLAP1 = 3
OVERLAP2 = 3


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

    if haplotype_dict == None:
        s1 = set([a[0] for a in f1.seq])
        s2 = set([a[0] for a in f2.seq])

        inter = set.intersection(s1,s2)

        if len(inter) < OVERLAP1:
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

        if matches/total > THRESHOLD:#post_samehap > CUTOFF:

            return start, end

        else:
            return None, None

    else:

        el1 = f1.name.split(':')
        start1, end1 = [int(x) for x in el1[-1].split('-')]


        el2 = f2.name.split(':')
        start2, end2 = [int(x) for x in el2[-1].split('-')]

        if end1 - start2 < 3: #start2 >= end1:
            return None, None

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
            return None, None

        if (match1/total1 > THRESHOLD and match2/total2 > THRESHOLD) or (1-(match1/total1) > THRESHOLD and 1-(match2/total2) > THRESHOLD):#post_samehap > CUTOFF:

            return start, end

        else:
            return None, None



cell_map = {'PGP1_21':0,'PGP1_22':1,'PGP1_A1':2}

def assign_fragments(flist, outputfile, haplotype_file=None):#hapblocks):

    total = 0
    paired_fragments = []

     # convert t_block to a dict for convenience
    haplotype_dict = None
    if haplotype_file != None:
        haplotype_dict = defaultdict(lambda: '-')
        t_blocklist = fileIO.parse_hapblock_file(haplotype_file,use_SNP_index=False)
        for t_block in t_blocklist:
            for i, a1, a2 in t_block:
                haplotype_dict[i] = a1 # index by genomic pos (1-indexed), return one haplotype


    for f1,f2 in combinations(flist, 2):


        if f1.seq[0][0] > f2.seq[0][0]:

            temp = f1
            f1 = f2
            f2 = temp

        assert(f1.seq[0][0] <= f2.seq[0][0])
        assert(f1.seq[0][1] <= f2.seq[0][1])

        start_SNP, end_SNP = overlap(f1,f2, haplotype_dict)

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

        paired_fragments.append((chrom, start, end, cell1, chamber1, cell2, chamber2))

    paired_fragments.sort(key=lambda x: (chr_num[x[0]],x[1]))


    if filter_inconsistent_haplotypes:

        paired_fragments_copy = copy(paired_fragments)
        current_paired_fragments = []
        firstpos = min([x[1] for x in paired_fragments])
        lastpos = max([x[2] for x in paired_fragments])
        chrom = paired_fragments[0][0]
        bad_fragments = set()
        failed = False
        fail_positions = []
        for pos in range(firstpos, lastpos+1):
            if pos % 1000000 == 0:
                print('{} Mb...'.format(int(pos/1000000)))
            num_popped = 0
            while paired_fragments != [] and (chr_num[paired_fragments[0][0]] < chr_num[chrom] or (paired_fragments[0][0] == chrom and paired_fragments[0][1] <= pos)):
                current_paired_fragments.append(paired_fragments.pop(0))
                num_popped += 1

            l1 = len(current_paired_fragments)
            # filter out paired fragments that we're past
            criteria = lambda pf: (pf[0] == chrom and pf[1] <= pos and pos <= pf[2])
            current_paired_fragments = list(filter(criteria,current_paired_fragments))
            l2 = len(current_paired_fragments)

            num_removed = l1 - l2

            if num_popped == 0 and num_removed == 0:
                if failed == True:
                    fail_positions.append(pos)
                continue

            # we'll determine if haplotype configuration is possible using a disjoint-set-forest
            # make a dictionary to be able to find our disjoint-set-forest nodes
            nodedict = dict()
            elements = [(x[3],x[4]) for x in current_paired_fragments] + [(x[5],x[6]) for x in current_paired_fragments]
            for element in elements:
                nodedict[element] = Node(element)

            # union together paired fragments, saying "we know these are same haplotype"
            for (chrom, start, end, cell1, chamber1, cell2, chamber2) in current_paired_fragments:
                Union(nodedict[(cell1,chamber1)],nodedict[(cell2,chamber2)])

            # parentdict has an arbitrary key, and the value is a set of same-haplotype elements
            parentdict = defaultdict(set)
            for v in nodedict.values():
                parentdict[Find(v).label].add(v.label)

            # traverse our haplotype sets and determine if any of them are fishy
            # e.g. three strands from the same haplotype in the same cell
            failed = False
            fail_cell = None
            for hapset in parentdict.values():

                cellcounter = defaultdict(int)
                for cell,chamber in hapset:
                    cellcounter[cell] += 1

                if cellcounter[cell] > 2:
                    fail_cell = cell
                    failed = True
                    break

            if failed:
                print("FAILURE at {} {}: cell {} has >=3 strands same haplotype.".format(chrom, pos, fail_cell))
                print("haplotype sets:")
                for hapset in parentdict.values():
                    print(hapset)
                print("bad_fragment pairs:")
                for frag in current_paired_fragments:
                    print(frag)
                print("----------------------------------------")
                bad_fragments = bad_fragments.union(set(current_paired_fragments))

        filtered_paired_fragments = list(set(paired_fragments_copy) - bad_fragments)
        filtered_paired_fragments.sort(key=lambda x: (chr_num[x[0]],x[1]))
        import pickle
        pickle.dump(fail_positions,open('fail_positions.p','wb'))
    else:

        filtered_paired_fragments = paired_fragments

    total_filtered = 0

    with open(outputfile, 'w') as opf:

        for (chrom, start, end, cell1, chamber1, cell2, chamber2) in filtered_paired_fragments:
            total_filtered += end - start
            print('{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(chrom, start, end, cell1, chamber1, cell2, chamber2), file=opf)

    print("TOTAL          : {}".format(total))
    if filter_inconsistent_haplotypes:
        print("TOTAL, FILTERED: {}".format(total_filtered))

    return total


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

            for chrom, start, end, cell1, chamber1, cell2, chamber2 in current_paired_fragments:

                assert(cell1 < cell2 or (cell1 == cell2 and chamber1 < chamber2))

                new_tag = 'HP:{},{}:{},{}'.format(cell1,chamber1,cell2,chamber2)
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
    total = assign_fragments(flist, outputfile, haplotype_file)

    return total

def parse_bedfile(input_file):

    boundaries = []
    with open(input_file,'r') as inf:
        for line in inf:

            if len(line) < 3:
                continue

            el = line.strip().split('\t')

            chrom = el[0]
            start = int(el[1])
            stop  = int(el[2])

            boundaries.append((chrom, start, stop))

    return boundaries

def pair_strands_XY(bedfile_lst,outfile):

    def filter_XY(bedlst):
        return [(chrom, start, stop) for (chrom, start, stop) in bedlst if chrom in ['chrX','chrY','X','Y']]

    bed_data = [filter_XY(parse_bedfile(bf)) for bf in bedfile_lst]
    with open(outfile,'w') as of:
        for i in range(len(bed_data)):
            for j in range(i+1, len(bed_data)):

                bd1 = bed_data[i]
                bd2 = bed_data[j]

                for (chrom1,start1,end1) in bd1:
                    for (chrom2, start2, end2) in bd2:

                        if chrom1 != chrom2:
                            continue

                        if start1 > start2:
                            temp1 = start1
                            temp2 = end1
                            start1 = start2
                            end1   = end2
                            start2 = temp1
                            end2   = temp2

                        cell1 = int(i / n_chambers)
                        chamber1 = int(i % n_chambers)
                        cell2 = int(j / n_chambers)
                        chamber2 = int(j % n_chambers)

                        if start2 < end1:
                            #overlap
                            end_overlap = min([end1, end2])
                            start_overlap = start2
                            el = [chrom1,start_overlap,end_overlap,cell1,chamber1,cell2,chamber2]
                            line = '\t'.join([str(x) for x in el])
                            print(line,file=of)


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
    pass
    #test_pair_strands()
    #pair_strands()
