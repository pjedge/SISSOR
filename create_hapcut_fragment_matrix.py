#! /bin/python3.5
# -*- coding: utf-8 -*-

# Author : Peter Edge
# Email  : pedge@eng.ucsd.edu
from collections import defaultdict
import os
import re
'''
Converts SISSOR haploid vcf files to hapcut fragment matrices.

Designed for work on ONE chromosome across many chambers:
    - The variant_vcf_file and pileup file should be for the same, single chromosome.
    - The sissor_vcf_files and bed_files should be in same order by chamber.
    - Recommended to include all sissor_vcf_files and bed_files for one SISSOR experiment.
    - The output_fragment_file should specify the chromosome label.
'''

def create_hapcut_fragment_matrices_freebayes(sissor_vcf_files, bed_files, variant_vcf_files, output_dir):

    dp_pat = re.compile("DP=(\d+)")
    qr_pat = re.compile(";QR=(\d+)")

    lines = defaultdict(list) # lines[chrom] contains a list of (first_index, line) tuples. these will be sorted by first pos and printed to fragment outfile

    snp_indices = dict() # chrom -> list of snp indices
    ix_counter = defaultdict(int)

    # need to define the columns of the fragment matrix (indexes and positions of SNPs)
    # use a vcf file as a template to know what SNP indices are
    for variant_vcf_file in variant_vcf_files:
        with open(variant_vcf_file, 'r') as infile:

            for line in infile:
                if (line[0] == '#' or len(line) < 2):
                    continue
                el = line.strip().split()
                chrom = el[0]
                pos   = int(el[1])-1
                snp_indices[(chrom,pos)] = ix_counter[chrom] # given a chromosome and genomic index, return the SNP index within that chroms vcf
                ix_counter[chrom] += 1

    # step through triples of a haploid vcf file, the bedfile that defines
    # where fragment boundaries are, and the pileup file.
    # pileup file allows us to know which SNPs have acceptable depth/quality
    for sissor_vcf_file, bed_file in zip(sissor_vcf_files, bed_files):

        call_dict = dict()
        with open(sissor_vcf_file, 'r') as infile:

            for line in infile:
                if len(line) < 3 or line[0] == '#':
                    continue

                # read line elements
                el    = line.strip().split()
                if len(el) < 3:
                    continue

                chrom = el[0]
                pos   = int(el[1])-1
                fields1 = el[7]
                labels = el[8]
                assert(labels[:2] == 'GT')
                fields2 = el[9]
                call = fields2[0]
                if call not in ['0','1']:
                    continue
                qual = int(float(el[5])) if call == '1' else int(float(re.findall(qr_pat,fields1)[0]))

                depth = int(float(re.findall(dp_pat,fields1)[0]))
                q_char = '~' if qual>=93 else chr(33 + qual)

                # filter on quality and depth
                # we only care about this position if it's also in the HapCUT block file
                if qual >= 30 and depth >= 5 and (chrom, pos) in snp_indices:
                    call_dict[(chrom,pos)] = (call, q_char, snp_indices[(chrom,pos)])

        # read bed file boundaries into a list structure
        bed_data = []
        with open(bed_file, 'r') as infile:
            for line in infile:
                if len(line) < 2:
                    continue
                bed_data.append(line.strip().split()[:3])

        # create a list where each element is a list corresponding to a bed region
        # each A contains tuples for each snp found in the hapcut haplotype that falls in the region
        frag_snps = [[] for i in bed_data]
        names     = []

        for i, (chrom, p1, p2) in enumerate(bed_data):
            p1 = int(p1)
            p2 = int(p2)
            curr_name = '{}:{}-{}'.format(chrom,p1,p2)
            names.append(curr_name)
            # consume SNPs leading up to bed region
            # add SNPs inside current bed region to fs
            for snp_pos in range(p1,p2+1):
                if (chrom, snp_pos) not in call_dict:
                    continue

                call, Q, snp_ix = call_dict[(chrom, snp_pos)] # get the quality char from the genotype call

                frag_snps[i].append((snp_ix, chrom, snp_pos, call, Q)) # mark this as a nonref allele

        # now take this information and make it into a fragment matrix file line
        for name, fs in zip(names, frag_snps):
            if len(fs) < 2:
                continue

            fragstr = ''
            num_pairs = 0
            prev_snp_ix = -2
            qual = ''
            chrom    = fs[0][1]
            firstpos = fs[0][2]
            for snp_ix, chrom, pos, allele, q_char in fs:

                diff = snp_ix - prev_snp_ix

                if diff == 1:
                    fragstr += allele
                else:
                    num_pairs += 1
                    fragstr += ' {} {}'.format(snp_ix+1, allele)

                prev_snp_ix = snp_ix
                qual += q_char

            fragstr += ' ' + ''.join(qual)

            prefix = '{} {}'.format(num_pairs,name)
            fragstr = prefix + fragstr

            lines[chrom].append((firstpos, fragstr))

    # go through the output lines we've accrued for each chromosome.
    # sort them and print them to a fragment file for each chromosome.
    for chrom in lines.keys():
        lines[chrom].sort()

        output_fragment_file = os.path.join(output_dir, chrom)

        with open(output_fragment_file, 'w') as opf:
            for firstpos, line in lines[chrom]:
                print(line, file=opf)

def create_hapcut_fragment_matrices(sissor_vcf_files, bed_files, pileup_files, variant_vcf_files, output_dir):

    lines = defaultdict(list) # lines[chrom] contains a list of (first_index, line) tuples. these will be sorted by first pos and printed to fragment outfile

    snp_indices = dict() # chrom -> list of snp indices
    ix_counter = defaultdict(int)

    # need to define the columns of the fragment matrix (indexes and positions of SNPs)
    # use a vcf file as a template to know what SNP indices are
    for variant_vcf_file in variant_vcf_files:
        with open(variant_vcf_file, 'r') as infile:

            for line in infile:
                if (line[0] == '#' or len(line) < 2):
                    continue
                el = line.strip().split()
                chrom = el[0]
                pos   = int(el[1])-1
                snp_indices[(chrom,pos)] = ix_counter[chrom] # given a chromosome and genomic index, return the SNP index within that chroms vcf
                ix_counter[chrom] += 1

    # step through triples of a haploid vcf file, the bedfile that defines
    # where fragment boundaries are, and the pileup file.
    # pileup file allows us to know which SNPs have acceptable depth/quality
    for sissor_vcf_file, bed_file, pileup_file in zip(sissor_vcf_files, bed_files, pileup_files):

        Qdict = dict()
        with open(pileup_file, 'r') as infile:

            for line in infile:

                # read line elements
                el    = line.strip().split()
                chrom = el[0]
                pos   = int(el[1])-1
                qual  = int(el[4])
                depth = int(el[7])
                q_char = '~' if qual>=93 else chr(33 + qual)

                # filter on quality and depth
                # we only care about this position if it's also in the HapCUT block file
                if qual >= 30 and depth >= 5 and (chrom, pos) in snp_indices:
                    Qdict[(chrom,pos)] = (q_char, snp_indices[(chrom,pos)])

        # read bed file boundaries into a list structure
        bed_data = []
        with open(bed_file, 'r') as infile:
            for line in infile:
                if len(line) < 2:
                    continue
                bed_data.append(line.strip().split()[:3])

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

        # create a list where each element is a list corresponding to a bed region
        # each A contains tuples for each snp found in the hapcut haplotype that falls in the region
        frag_snps = [[] for i in bed_data]
        names     = []

        for i, (chrom, p1, p2) in enumerate(bed_data):
            p1 = int(p1)
            p2 = int(p2)
            curr_name = '{}:{}-{}'.format(chrom,p1,p2)
            names.append(curr_name)
            # consume SNPs leading up to bed region
            # add SNPs inside current bed region to fs
            for snp_pos in range(p1,p2+1):
                if (chrom, snp_pos) not in Qdict:
                    continue

                Q, snp_ix = Qdict[(chrom, snp_pos)] # get the quality char from the genotype call

                if (chrom, snp_pos) in nonrefs:
                    frag_snps[i].append((snp_ix, chrom, snp_pos, '1', Q)) # mark this as a nonref allele
                else:
                    frag_snps[i].append((snp_ix, chrom, snp_pos, '0', Q)) # mark this as a ref allele

        # now take this information and make it into a fragment matrix file line
        for name, fs in zip(names, frag_snps):
            if len(fs) < 2:
                continue

            fragstr = ''
            num_pairs = 0
            prev_snp_ix = -2
            qual = ''
            chrom    = fs[0][1]
            firstpos = fs[0][2]
            for snp_ix, chrom, pos, allele, q_char in fs:

                diff = snp_ix - prev_snp_ix

                if diff == 1:
                    fragstr += allele
                else:
                    num_pairs += 1
                    fragstr += ' {} {}'.format(snp_ix+1, allele)

                prev_snp_ix = snp_ix
                qual += q_char

            fragstr += ' ' + ''.join(qual)

            prefix = '{} {}'.format(num_pairs,name)
            fragstr = prefix + fragstr

            lines[chrom].append((firstpos, fragstr))

    # go through the output lines we've accrued for each chromosome.
    # sort them and print them to a fragment file for each chromosome.
    for chrom in lines.keys():
        lines[chrom].sort()

        output_fragment_file = os.path.join(output_dir, chrom)

        with open(output_fragment_file, 'w') as opf:
            for firstpos, line in lines[chrom]:
                print(line, file=opf)
