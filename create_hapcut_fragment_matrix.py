#! /bin/python3.5
# -*- coding: utf-8 -*-

# Author : Peter Edge
# Email  : pedge@eng.ucsd.edu
import argparse
import sys

desc = '''converts SISSOR haploid vcf files to hapcut fragment matrices'''

# parse arguments
def parse_args():

    parser = argparse.ArgumentParser(description='')
    parser.add_argument('sissor_vcf_file', nargs='?', type = str, help='path to SISSOR haploid vcf file containing only non-reference calls')
    parser.add_argument('bed_file', nargs='?', type = str, help='bed file where each entry specifies the boundary of a SISSOR fragment')
    parser.add_argument('variant_vcf_file', nargs='?', type = str, help='vcf file containing reliable variant calls to be used with HapCUT', default = None)
    parser.add_argument('pileup_file', nargs='?', type = str, help='oldschool samtools pileup file from SISSOR raw reads')
    parser.add_argument('output_fragment_file', nargs='?', type = str, help='path to result fragment file')
    # default to help option. credit to unutbu: http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
    if len(sys.argv) < 5:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    return args

def create_hapcut_fragment_matrix(sissor_vcf_file, bed_file, variant_vcf_file, pileup_file, output_fragment_file):

    # read bed file boundaries into a list structure
    bed_data = []
    with open(bed_file, 'r') as infile:
        for line in infile:
            if len(line) < 2:
                continue
            bed_data.append(line.strip().split())

    # create a list of (chromosome, pos) tuples specifying each allele present in hapcut haplotype
    snp_indices_unsorted = []
    # want a set too for fast lookup
    snp_indices_set      = set()

    # need to define the columns of the fragment matrix (indexes and positions of SNPs)
    # use a vcf file as a template to know where SNPs and their indices are
    with open(variant_vcf_file, 'r') as infile:

        ix = 1
        for line in infile:
            if (line[0] == '#' or len(line) < 2):
                continue
            el = line.strip().split()
            chrom = el[0]
            pos   = int(el[1])-1
            snp_indices_unsorted.append((ix,chrom,pos))
            snp_indices_set.add((chrom,pos))
            ix += 1

    snp_indices = sorted(snp_indices_unsorted, key=lambda tup: tup[0])

    covered = set()
    if pileup_file != None:
        with open(pileup_file, 'r') as infile:

            for line in infile:

                # read line elements
                el    = line.strip().split()
                chrom = el[0]
                pos   = int(el[1])-1
                qual  = float(el[4])
                depth = int(el[7])
                Q = 10**(qual*-0.1)

                # filter on quality and depth
                # we only care about this position if it's also in the HapCUT block file
                if qual >= 30 and depth >= 5 and (chrom, pos) in snp_indices_set:
                    covered.add((chrom,pos,Q))

    # create a list where each element is a list corresponding to a bed region
    # each A contains tuples for each snp found in the hapcut haplotype that falls in the region
    frag_snps = [[] for i in bed_data]
    names     = []

    for i, (chrom, p1, p2, l, chamber) in enumerate(bed_data):
        p1 = int(p1)
        p2 = int(p2)
        curr_name = '{}:{}-{}'.format(chrom[3:],p1,p2)
        names.append(curr_name)
        if snp_indices == []:
            break
        snp_ix, snp_chrom, snp_pos = snp_indices.pop(0)
        # consume SNPs leading up to bed region
        while (snp_chrom != chrom or snp_pos < p1) and snp_indices != []:
            snp_ix, snp_chrom, snp_pos = snp_indices.pop(0)
        # add SNPs inside current bed region to fs
        while snp_chrom == chrom and p1 <= snp_pos and snp_pos <= p2 and snp_indices != []:
            if (snp_chrom, snp_pos) in covered:
                frag_snps[i].append((snp_ix, snp_chrom, snp_pos, '0', Q))
            snp_ix, snp_chrom, snp_pos = snp_indices.pop(0)

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

    # go through our list of fragment SNP indices and use the set of nonref alleles
    # to mark each nonreference position as such
    for i in range(0,len(frag_snps)):

        for j in range(0,len(frag_snps[i])):

            snp_ix, chrom, pos, allele, Q = frag_snps[i][j]

            if (chrom, pos) in nonrefs:
                frag_snps[i][j] = (snp_ix, chrom, pos, '1', Q) # mark this as a nonref allele


    # now take this information and print it to a fragment matrix file
    with open(output_fragment_file, 'w') as opf:

        for name, fs in zip(names, frag_snps):
            fragstr = ''
            num_pairs = 0
            prev_snp_ix = -2
            qual = ''

            for snp_ix, chrom, pos, allele, Q in fs:

                diff = snp_ix - prev_snp_ix

                if diff == 1:
                    fragstr += allele
                else:
                    num_pairs += 1
                    fragstr += ' {} {}'.format(snp_ix+1, allele)

                prev_snp_ix = snp_ix
                qual += Q

            fragstr += ' ' + ''.join(qual)

            prefix = '{} {} '.format(num_pairs,name)
            fragstr = prefix + fragstr
            print(fragstr, file=opf)

# parse args and execute main function
if __name__ == '__main__':
    args = parse_args()
    create_hapcut_fragment_matrix(args.sissor_vcf_file, args.bed_file, args.variant_vcf_file, args.pileup_file, args.output_fragment_file)
