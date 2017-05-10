#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  9 16:35:02 2016

@author: peter
"""

import sys
import fileIO
import pickle


# a class representing a HapCUT2 format haplotype fragment
# this is mostly agnostic to the SISSOR method, except for two special fields: first_piece and last_piece
# these are present so when we are done processing SISSOR fragments to remove contamination and detectable switch errors,
# we can extend the ends of the fragments to the original boundaries rather than the endpieces being truncated
# at the first and last heterozygous SNV from our SNV set used for phasing.
class fragment:

    def __init__(self, seq, name, first_piece = True, last_piece = True):
        self.seq = seq                         # list of (snp index, genomic index, allele call, quality score) tuples
        self.name = name                       # fragment ID / name
        self.first_piece = first_piece
        self.last_piece  = last_piece

    def __str__(self):
        fragstr = ''
        num_pairs = 0
        prev_snp_ix = -2
        qual = ' '
        for snp_ix, genome_ix, allele, q_char in self.seq:

            diff = snp_ix - prev_snp_ix

            if diff == 1:
                fragstr += allele
            else:
                num_pairs += 1
                fragstr += ' {} {}'.format(snp_ix+1, allele)

            prev_snp_ix = snp_ix
            qual += q_char

        fragstr += qual

        prefix = '{} {}'.format(num_pairs,self.name)
        fragstr = prefix + fragstr
        return fragstr

# read in a HapCUT2 format fragment file into a list of fragment objects
def read_fragment_matrix(frag_matrix, vcf_file):

    snp_ix = 0
    vcf_dict = dict()
    with open(vcf_file,'r') as infile:
        for line in infile:
            if line[:1] == '#':
                continue
            el = line.strip().split('\t')
            if len(el) < 5:
                continue

            genomic_pos = int(el[1])-1
            vcf_dict[snp_ix] = genomic_pos
            snp_ix += 1

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

            alist= [(a,vcf_dict[a],b,c) for ((a,b),c) in zip(call_list2,qlist)]

            frag = fragment(alist,name)
            flist.append(frag)

    sorted_flist = sorted(flist,key=lambda x: x.seq[0][0])

    return sorted_flist

# print out a list of fragment objects to a HapCUT2 format fragment file
def write_fragment_matrix(flist,outfile):
    lines = []

    for f in flist:
        if len(f.seq) < 2:
            continue
        firstpos = f.seq[0][0]
        lines.append((firstpos, str(f)))

    lines.sort()

    with open(outfile, 'w') as opf:
        for firstpos, line in lines:
            print(line, file=opf)

# a visualization function -- visualize fragments in full-blown aligned matrix format.
# very space inefficient but useful for analyzing toy examples and specific aligned haplotype fragment cases.
def matrixify_flist(flist, o=None):

    max_ix = 0
    max_name = 0
    min_ix = min([f.seq[0][0] if len(f.seq) > 0 else int(1e9) for f in flist])

    for f in flist:

        if len(f.name) > max_name:
            max_name = len(f.name)

        for a in f.seq:

            if a[0] > max_ix:

                max_ix = a[0]

    if o != None:
        for f in flist:

            line = ['-'] * (max_ix+1-min_ix)
            for a in f.seq:
                line[a[0] - min_ix] = a[2]

            line = [f.name.ljust(max_name+1)] + line
            pline = ''.join(line)

            print(pline,file=o)
    else:
        for f in flist:

            line = ['-'] * (max_ix+1 - min_ix)
            for a in f.seq:
                line[a[0] - min_ix] = a[2]

            line = [f.name.ljust(max_name+1)] + line
            pline = ''.join(line)

            print(pline)

# determine if two fragments overlap
def overlap(f1, f2, amt=1):

    s1 = set([a[0] for a in f1.seq])
    s2 = set([a[0] for a in f2.seq])
    inter = set.intersection(s1,s2)

    return len(inter) >= amt
