#! /bin/python3.5
# -*- coding: utf-8 -*-

# Author : Peter Edge
# Email  : pedge@eng.ucsd.edu
import argparse
import os
import sys
import math
from collections import defaultdict

desc = '''
Converts SISSOR haploid vcf files to hapcut fragment matrices.

Designed for work on ONE chromosome across many chambers:
    - The variant_vcf_file and pileup file should be for the same, single chromosome.
    - The sissor_vcf_files and bed_files should be in same order by chamber.
    - Recommended to include all sissor_vcf_files and bed_files for one SISSOR experiment. 
    - The output_fragment_file should specify the chromosome label.
'''

# parse arguments
def parse_args():

    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-s','--sissor_vcf_files', nargs='+', type = str, help='path to SISSOR haploid vcf file containing only non-reference calls',default=[])
    parser.add_argument('-b','--bed_files', nargs='+', type = str, help='bed file where each entry specifies the boundary of a SISSOR fragment',default=[])
    parser.add_argument('-p','--pileup_files', nargs='+', type = str, help='oldschool samtools pileup file from SISSOR raw reads',default=[])
    parser.add_argument('-v','--variant_vcf_files', nargs='+', type = str, help='vcf file containing reliable variant calls to be used with HapCUT',default=[])
    parser.add_argument('-o','--output_dir', nargs='?', type = str, help='directory to output fragment files')
    parser.add_argument('-c','--condensed_pileup', action='store_true', help='set this flag if using special 3 column condensed pileup')
    # default to help option. credit to unutbu: http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
    if len(sys.argv) < 5:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    return args

def create_hapcut_fragment_matrices(sissor_vcf_files, bed_files, variant_vcf_files, pileup_files, output_dir, condensed_pileup=False):

    if not os.path.exists(output_dir):
        os.path.makedirs(output_dir)

    lines = defaultdict(list) # lines[chrom] contains a list of (first_index, line) tuples. these will be sorted by first pos and printed to fragment outfile
    snp_indices = dict() # chrom -> list of snp indices
    ix_counter = defaultdict(lambda: 1)
    
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
        
        # read bed file boundaries into a list structure
        bed_data = []
        with open(bed_file, 'r') as infile:
            for line in infile:
                if len(line) < 2:
                    continue
                bed_data.append(line.strip().split())

        # make a "Q_dict" dictionary to indicate if a position is sufficiently
        # covered by reads and what its quality score is
        Q_dict = dict()
        if condensed_pileup:
            
            with open(pileup_file, 'r') as infile:
    
                for line in infile:
    
                    # read line elements
                    el    = line.strip().split()
                    chrom = el[0]
                    pos   = int(el[1])
                    qual  = float(el[2])
    
                    # already filtered
                    
                    q_char = '~' if qual < 5.011872336272714e-10 else chr(int(33-10*math.log10(qual)))

                    Q_dict[(chrom,pos)] = q_char
        
        else:
            with open(pileup_file, 'r') as infile:
    
                for line in infile:
    
                    # read line elements
                    el    = line.strip().split()
                    chrom = el[0]
                    pos   = int(el[1])-1
                    qual  = int(el[4])
                    depth = int(el[7])
                    Q = 10**(qual*-0.1)
                    q_char = '~' if qual>=126 else chr(33 + qual)
                    
                    # filter on quality and depth
                    # we only care about this position if it's also in the HapCUT block file
                    if qual >= 30 and depth >= 5 and (chrom, pos) in snp_indices:
                        Q_dict[(chrom,pos)] = q_char

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
    
        for i, (chrom, p1, p2, l, chamber) in enumerate(bed_data):
            p1 = int(p1)
            p2 = int(p2)
            curr_name = '{}:{}-{}'.format(chrom,p1,p2)
            names.append(curr_name)
            # consume SNPs leading up to bed region
            # add SNPs inside current bed region to fs
            for snp_pos in range(p1,p2+1):
                if (chrom, snp_pos) not in Q_dict or (chrom, snp_pos) not in snp_indices:
                    continue
                
                Q = Q_dict[(chrom, snp_pos)] # get the quality char from the genotype call
                snp_ix = snp_indices[(chrom, snp_pos)]
                
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

            
# parse args and execute main function
if __name__ == '__main__':
    args = parse_args()
    create_hapcut_fragment_matrices(args.sissor_vcf_files, args.bed_files, args.variant_vcf_files, args.pileup_files, args.output_dir, args.condensed_pileup)
