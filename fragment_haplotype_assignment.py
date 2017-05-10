# -*- coding: utf-8 -*-
"""
Created on Thu May 19 21:44:48 2016

@author: peter
"""
import sys
import fragment
import fileIO
from collections import defaultdict
from itertools import combinations
from chamber_allele_call_accuracy import parse_bedfile

ID_LEN_NO_MODS = 5

tinylog = -1e5
n_chambers = 24
n_cells = 3

THRESHOLD = 0.8 # fraction of match to say haplotype is the same

chroms = ['chr{}'.format(i) for i in range(1,23)] + ['chrX','chrY']

chr_num = dict()
for i,chrom in enumerate(chroms):
    chr_num[chrom] = i

# parse a haplotype block file from HapCUT2 and return a list of blocks.
# each block is a list of
def parse_hapblock_file(hapblock_file,use_SNP_index=True):

    blocklist = [] # data will be a list of blocks, each with a list tying SNP indexes to haplotypes

    try:
        with open(hapblock_file, 'r') as hbf:

            for line in hbf:
                if len(line) < 3: # empty line
                    continue
                if 'BLOCK' in line:
                    blocklist.append([])
                    continue

                elements = line.strip().split('\t')
                if len(elements) < 3: # not enough elements to have a haplotype
                    continue

                pos = int(elements[4])-1 if not use_SNP_index else int(elements[0])-1
                allele1 = elements[1]
                allele2 = elements[2]

                blocklist[-1].append((pos, allele1, allele2))

    except FileNotFoundError:
        # most of the time, this should mean that the program timed out and therefore didn't produce a phase.
        print("File {} was missing. Error calculation will continue without it.".format(hapblock_file), file=sys.stderr)
        pass

    return blocklist

# assign a single fragment object to a haplotype, given a haplotype dictionary containing phase
def assign_fragment(f, haplotype_dict):

    match = 0
    total = 0
    el = f.name.split(':')
    start, end = [int(x) for x in el[-1].split('-')]

    for (snp_ix, gen_ix, a, q) in f.seq: #s1:
        if haplotype_dict[gen_ix] not in ['0','1'] or a not in ['0','1']:
            continue
        if haplotype_dict[gen_ix] == a:
            match += 1
        total += 1

    if total == 0:

        return None, None, None

    elif match/total > THRESHOLD:

        return start, end, 'P1'

    elif 1 - match/total > THRESHOLD:

        return start, end, 'P2'

    else:

        return None, None, None

cell_map = {'PGP1_21':0,'PGP1_22':1,'PGP1_A1':2}


# assign fragments in a fragment block file (processed SISSOR fragments) to P1 or P2 haplotype
# print out the chromosome and boundaries of the fragment, the cell and chamber in which it appears, and the P1/P2 assignment.
def assign_fragment_haplotypes(fragmentfile, vcf_file, outputfile, haplotype_file):

    flist = fragment.read_fragment_matrix(fragmentfile, vcf_file)

    assigned_fragments = []
    t_blocklist = parse_hapblock_file(haplotype_file,use_SNP_index=False)

    for t_block in t_blocklist:

        # complexity of this approach is poor
        # but it's straightforward
        # and necessary because haplotype blocks may be split
        # after assembly, resulting in fragments that span two haplotype blocks.
         # convert t_block to a dict for convenience
        haplotype_dict = defaultdict(lambda: '-')
        for i, a1, a2 in t_block:
            haplotype_dict[i] = a1 # index by genomic pos (1-indexed), return one haplotype

        for f in flist:

            start_SNP, end_SNP, parent = assign_fragment(f, haplotype_dict)

            if start_SNP == None:
                continue

            el = f.name.split(':')
            chrom = el[0]
            start, end = [int(x) for x in el[-1].split('-')]
            cell = cell_map[el[2]]
            chamber = el[3]
            assert(chamber[0:2] == 'CH')
            chamber = int(chamber[2:]) - 1

            assigned_fragments.append((chrom, start, end, cell, chamber, parent))

    assigned_fragments.sort(key=lambda x: (chr_num[x[0]],x[1]))

    with open(outputfile, 'w') as opf:

        for (chrom, start, end, cell, chamber, parent) in assigned_fragments:
            print('{}\t{}\t{}\t{}\t{}\t{}'.format(chrom, start+1, end+1, cell, chamber, parent), file=opf)


# edit the "chamber call files" which have have allele probabilities for each chamber,
# with tags "P1" and "P2" that list which cells and chambers have strands assigned to P1 haplotype
# or P2 haplotype at that position.
# this information will be used in later steps to compare variant calls between haplotype-matched strands
def annotate_assigned_fragments(chamber_call_file, fragment_assignment_file, output_file):

    # read in file with spans of assigned strands with matching haplotypes
    assigned_fragments = [] # list of (start, end, ix1, ix2) tuples
    with open(fragment_assignment_file,'r') as infile:
        for line in infile:
            el = line.strip().split()
            chrom = el[0]
            start = int(el[1])
            end   = int(el[2])
            cell = int(el[3])
            chamber = int(el[4])
            parent = el[5]

            assigned_fragments.append((chrom, start, end, cell, chamber, parent))

    assigned_fragments.sort(key=lambda x: (chr_num[x[0]],x[1]))

    current_assigned_fragments = []

    with open(chamber_call_file,'r') as ccf, open(output_file,'w') as opf:
        #print('chr\tpos\tsissor_call\tCGI_allele\tref',file=mof)

        for line in ccf:
            ccf_line = line.strip().split('\t')
            ccf_chrom = ccf_line[0]
            ccf_pos   = int(ccf_line[1])

            # pop info for haplotype-assigned regions into the "current" list
            while assigned_fragments != [] and (chr_num[assigned_fragments[0][0]] < chr_num[ccf_chrom] or (assigned_fragments[0][0] == ccf_chrom and assigned_fragments[0][1] <= ccf_pos)):
                current_assigned_fragments.append(assigned_fragments.pop(0))

            # filter out assigned fragments that we're past
            criteria = lambda pf: (pf[0] == ccf_chrom and pf[1] <= ccf_pos and ccf_pos <= pf[2])
            current_assigned_fragments = list(filter(criteria,current_assigned_fragments))

            tags = ccf_line[80].split(';')
            new_tags = []

            P1_chambers = []
            P2_chambers = []
            for chrom, start, end, cell, chamber, parent in current_assigned_fragments:
                if parent == 'P1':
                    P1_chambers.append((cell, chamber))
                elif parent == 'P2':
                    P2_chambers.append((cell, chamber))

            P1_chambers.sort()
            P2_chambers.sort()
            new_tag_P1 = ':'.join(['P1'] + ['{},{}'.format(cell, chamber) for (cell, chamber) in P1_chambers])
            new_tag_P2 = ':'.join(['P2'] + ['{},{}'.format(cell, chamber) for (cell, chamber) in P2_chambers])
            new_tags += [new_tag_P1, new_tag_P2]

            if tags == ['N/A'] and new_tags != []:
                tags = new_tags
            else:
                tags += new_tags

            ccf_line[80] = ';'.join(tags)
            new_line = '\t'.join(ccf_line)
            print(new_line, file=opf)

# counts the total length of fragments that cannot be assigned to haplotypes
def count_not_matchable(fragmentfile, vcf_file, haplotype_file):

    flist = fragment.read_fragment_matrix(fragmentfile, vcf_file)
     # convert t_block to a dict for convenience
    haplotype_dict = None
    if haplotype_file != None:
        haplotype_dict = defaultdict(lambda: '-')
        t_blocklist = parse_hapblock_file(haplotype_file,use_SNP_index=False)
        for t_block in t_blocklist:
            for i, a1, a2 in t_block:
                haplotype_dict[i] = a1 # index by genomic pos (1-indexed), return one haplotype

    unmatchable = defaultdict(int)
    for f in flist:
        el = f.name.split(':')
        start, end = [int(x) for x in el[-1].split('-')]
        cell = el[2]
        match = 0
        total = 0

        for (snp_ix, gen_ix, a, q) in f1.seq:  # s1:
            if haplotype_dict[gen_ix] not in ['0', '1'] or a not in ['0', '1']:
                continue
            if haplotype_dict[gen_ix] == a:
                match += 1
            total += 1

        if total < 1 or not (match/total > THRESHOLD or 1-match/total1 > THRESHOLD):

            unmatchable[cell] += end - start

    return unmatchable
