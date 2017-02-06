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

HET_SPLIT_MIN_SIZE = 2
FIXED_Q_CHAR = '5'

short_chrom_names = set([str(x) for x in range(1,23)]+['X','Y'])
chrom_names = set(['chr'+str(x) for x in range(1,23)]+['chrX','chrY'])

def format_chrom(chrom_str):
    if chrom_str in short_chrom_names:
        chrom_str = 'chr' + chrom_str
    if chrom_str not in chrom_names:
        return False
    return chrom_str

def create_hapcut_fragment_matrices_freebayes(ploidy1_files, ploidy2_files, bed_files, cell_id, chamber_names, variant_vcf_files, output_dir,hets_in_seq=False):

    dp_pat = re.compile("DP=(\d+)")
    qr_pat = re.compile(";QR=(\d+)")

    lines = defaultdict(list) # lines[chrom] contains a list of (first_index, line) tuples. these will be sorted by first pos and printed to fragment outfile

    snp_indices = dict() # chrom -> list of snp indices

    # need to define the columns of the fragment matrix (indexes and positions of SNPs)
    # use a vcf file as a template to know what SNP indices are
    for variant_vcf_file in variant_vcf_files:
        ix_counter = 0
        with open(variant_vcf_file, 'r') as infile:

            for line in infile:
                if (line[0] == '#' or len(line) < 2):
                    continue
                el = line.strip().split()
                chrom = format_chrom(el[0])
                if not chrom:
                    continue
                pos   = int(el[1])-1
                snp_indices[(chrom,pos)] = ix_counter # given a chromosome and genomic index, return the SNP index within that chroms vcf
                ix_counter += 1

    for ploidy1_file, ploidy2_file, bed_file, chamber_name in zip(ploidy1_files, ploidy2_files, bed_files, chamber_names):

        call_dict = dict()
        with open(ploidy1_file, 'r') as infile:

            for line in infile:
                if len(line) < 3 or line[0] == '#':
                    continue

                # read line elements
                el    = line.strip().split()
                if len(el) < 3:
                    continue

                chrom = format_chrom(el[0])
                if not chrom:
                    continue
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
                #q_char = '~' if qual>=93 else chr(33 + qual)

                # filter on quality and depth
                # we only care about this position if it's also in the HapCUT block file
                if qual >= 30 and depth >= 5 and (chrom, pos) in snp_indices:
                    call_dict[(chrom,pos)] = (call,snp_indices[(chrom,pos)])

        # read bed file boundaries into a list structure
        bed_data = []
        with open(bed_file, 'r') as infile:
            for line in infile:
                if len(line) < 2:
                    continue
                curr_data = line.strip().split()[:3]
                curr_data[0] = format_chrom(curr_data[0])
                if not curr_data[0]:
                    continue
                bed_data.append(curr_data)

        # use diploid FreeBayes file to find positions reliably called as heterozygous and split on them
        het_calls = dict()

        if ploidy2_file != None:

            with open(ploidy2_file,'r') as infile:
                for line in infile:
                    if len(line) < 3 or line[0] == '#':
                        continue

                    # read line elements
                    el    = line.strip().split()
                    if len(el) < 10:
                        continue
                    chrom = format_chrom(el[0])
                    if not chrom:
                        continue
                    pos   = int(el[1])-1
                    fields1 = el[7]
                    labels = el[8]
                    assert(labels[:2] == 'GT')
                    fields2 = el[9]
                    genotype = fields2[0:3]

                    #qual = int(float(el[5]))
                    depth = int(float(re.findall(dp_pat,fields1)[0]))

                    #if qual > 30 and depth >= 5 and genotype in ['0/1','1/0']:
                    if genotype in ['0/1','1/0'] and depth >= 5 and (chrom, pos) in snp_indices:

                        het_calls[(chrom,pos)] = snp_indices[(chrom,pos)]

        # create a list where each element is a list corresponding to a bed region
        # each A contains tuples for each snp found in the hapcut haplotype that falls in the region
        frag_snps = []
        names     = []

        for i, (chrom, p1, p2) in enumerate(bed_data):
            p1 = int(p1)
            p2 = int(p2)
            curr_name = '{}:{}-{}:{}:CH{}'.format(chrom,p1,p2,cell_id,chamber_name)
            frag = []
            # consume SNPs leading up to bed region
            # add SNPs inside current bed region to fs

            #part_counter = 0
            #seen_het = False
            snps_spanned = 0
            hets_seen = 0
            hom_seen = 0

            for snp_pos in range(p1,p2+1):

                if (chrom,snp_pos) in call_dict or (chrom,snp_pos) in het_calls:
                    snps_spanned += 1

                if (chrom,snp_pos) in het_calls:
                    hets_seen += 1
                    if hets_in_seq:
                        snp_ix = het_calls[(chrom,snp_pos)]
                        frag.append((snp_ix, chrom, snp_pos, '2', FIXED_Q_CHAR))
                    continue

                    #if len(frag) >= HET_SPLIT_MIN_SIZE: # and not seen_het:
                    #    print("Splitting {} at {} due to heterozygous call".format(curr_name, snp_pos))
                    #    frag_snps.append(frag)
                    #    part_name = '{}:H{}'.format(curr_name,part_counter)
                    #    part_counter += 1
                    #    names.append(part_name)
                    #frag = []
                    #seen_het = True
                    #continue


                if (chrom, snp_pos) not in call_dict:
                    continue
                hom_seen += 1
                call, snp_ix = call_dict[(chrom, snp_pos)] # get the quality char from the genotype call

                frag.append((snp_ix, chrom, snp_pos, call, FIXED_Q_CHAR)) # mark this as a nonref allele

            if hom_seen >= 2:
                frag_snps.append(frag)
                #if part_counter == 0:
                names.append(curr_name)
                #else:
                #    part_name = '{}:H{}'.format(curr_name,part_counter)
                #    names.append(part_name)

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
                #qual += q_char
                # currently just fixing quality char at 1% error since we can't really justify using genotype Qs...
                qual += FIXED_Q_CHAR

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
