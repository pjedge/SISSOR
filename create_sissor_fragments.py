#! /bin/python3.5
# -*- coding: utf-8 -*-

# Author : Peter Edge
# Email  : pedge@eng.ucsd.edu

class fragment:

    def __init__(self):
        self.id    = None
        self.chrom = None
        self.seq   = [] # a list of (ix, pos, allele) tuples

def create_sissor_fragments(sissor_vcf_file, bed_file, hapcut_block_file, pileup_file, covered_file):

    # read bed file boundaries into a list structure
    bed_data = []
    with open(bed_file, 'r') as infile:
        for line in infile:
            if len(line) < 2:
                continue
            bed_data.append(line.strip().split())

    # create a list of (chromosome, pos) tuples specifying each allele present in hapcut haplotype
    # want a set too for fast lookup
    hapblock_file_snps = {}

    with open(hapcut_block_file, 'r') as infile:

        for line in infile:
            if ('BLOCK' in line or '********' in line or len(line) < 2):
                continue
            el = line.strip().split()
            last_el  = el[-1]
            el2 = last_el.split(':')
            last_el2 = el2[-1]

            # we want to filter out lines that end in FV or whose last column's
            # last number is > 2.0
            if last_el2 == 'FV' or float(last_el2) < 2.0:
                continue

            chrom = el[3]
            snp_ix = int(el[0])-1
            pos   = int(el[4])-1

            hapblock_file_snps[(chrom,pos)] = snp_ix

    covered = set()
    if pileup_file != None:
        with open(pileup_file, 'r') as infile:

            for line in infile:

                # read line elements
                el    = line.strip().split()
                chrom = el[0]
                pos   = int(el[1])-1
                qual  = int(el[4])
                depth = int(el[7])

                # filter on quality and depth
                # we only care about this position if it's also in the HapCUT block file
                if qual >= 30 and depth >= 5 and (chrom, pos) in hapblock_file_snps:
                    covered.add((chrom,pos))
    else:
        with open(covered_file, 'r') as infile:
            for line in infile:
                if len(line) < 2:
                    continue
                el = line.strip().split()
                chrom = el[0]
                pos   = int(el[1])
                covered.add((chrom,pos))

    # create a set of indices from the haploid SISSOR VCF where the non-reference alleles are found
    nonrefs = set()
    with open(sissor_vcf_file, 'r') as infile:

        for line in infile:
            if line[0] == '#' or len(line) < 2:
                continue

            # read line elements
            el    = line.strip().split()

            chrom = el[0]
            pos   = int(el[1])-1
            nonrefs.add((chrom,pos))

    # create a list where each element is a list A corresponding to a bed region
    # each A contains tuples for each snp found in the hapcut haplotype that falls in the region
    frag_list = [fragment() for i in bed_data]

    for i, (chrom, p1, p2, l, chamber) in enumerate(bed_data):
        p1 = int(p1)
        p2 = int(p2)

        frag_list[i].id = '{}:{}-{}'.format(chrom[3:],p1,p2)
        frag_list[i].chrom = chrom

        # add SNPs inside current bed region to fs
        for pos in range(p1,p2+1):
            if (chrom, pos) in covered:
                snp_ix = hapblock_file_snps[(chrom,pos)]
                if (chrom, pos) in nonrefs:
                    frag_list[i].seq.append((snp_ix, pos, '1')) # mark this as a nonref allele
                else:
                    frag_list[i].seq.append((snp_ix, pos, '0')) # mark this as a reference allele

    return frag_list
