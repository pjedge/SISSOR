#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 18:29:04 2016

@author: peter
"""

# chamber call file has results from SISSOR cross-chamber allele calls
# GFF file has set of known alleles for individual (for comparison)
# VCF file has another set of heterozygous SNVs for individual
# cutoff is the minimum probability of an allele call to use it
# result outfile simply contains the counts of match vs mismatched alleles
# mismatch_outfile prints locations of mismatched allele calls
def chamber_allele_call_accuracy(chamber_call_file, GFF_file, VCF_file, cutoff, result_outfile, mismatch_outfile):

    # add heterozygous alleles observed in one CGI dataset to a dictionary
    cgi_dict = dict()
    with open(GFF_file,'r') as gff:
        for line in gff:
            if line[0] == '#' or len(line) < 3:
                continue
            gff_line = line.strip().split('\t')
            if gff_line[2] != 'SNP':
                continue

            gff_chrom = gff_line[0]
            gff_pos   = int(gff_line[3])

            assert(gff_pos == int(gff_line[4]))

            info = gff_line[8].strip().split(';')
            a_info = info[0]
            assert('alleles' in a_info)
            a_info = a_info[8:]

            #r_info = info[2]
            #assert('ref_allele' in r_info)
            ref_allele = 'N'#r_info[10:]

            if len(a_info) == 1:
                #alleles = {a_info}
                continue
            elif len(a_info) == 3:
                alleles = {a_info[0],a_info[2]}
            else:
                print("Unexpected number of alleles")
                assert(False)

            assert(len(alleles) == 2)

            cgi_dict[(gff_chrom,gff_pos)] = alleles

    # add heterozygous SNVs observed in our second CGI dataset to the dictionary
    with open(VCF_file,'r') as vcf:
        for line in vcf:
            if line[0] == '#' or len(line) < 3:
                continue

            vcf_line = line.strip().split('\t')

            vcf_chrom = vcf_line[0]
            vcf_pos   = int(vcf_line[1])
            alleles = {vcf_line[3],vcf_line[4]}

            if (vcf_chrom,vcf_pos) in cgi_dict:
                if cgi_dict[(vcf_chrom,vcf_pos)] != alleles:
                    continue # vcf data doesn't match GFF CGI data
            else:
                cgi_dict[(vcf_chrom,vcf_pos)] = alleles

    match = 0
    mismatch = 0
    snv_match = 0
    snv_mismatch = 0

    with open(chamber_call_file,'r') as ccf, open(mismatch_outfile,'w') as mof:
        print('chr\tpos\tsissor_call\tCGI_allele\tref',file=mof)
        for line in ccf:
            ccf_line = line.strip().split('\t')
            ccf_chrom = ccf_line[0]
            ccf_pos   = int(ccf_line[1])
            ref_allele = ccf_line[2]

            tags = ccf_line[-1].split(';')

            if ccf_chrom == 'chrX' or ccf_chrom == 'chrY':
                continue

            if (ccf_chrom,ccf_pos) not in cgi_dict:
                continue

            if 'TOO_MANY_ALLELES' in tags or 'TOO_MANY_CHAMBERS' in tags:
                continue

            call = ccf_line[3]

            if call == '*':
                continue

            el2 = call.split(';')

            genotype_prob = -1
            #max_genotype = 'NN'
            for entry in el2:

                genotype, prob = entry.split(':')
                prob = float(prob)

                if prob > genotype_prob:
                    genotype_prob = prob
                    #max_genotype = genotype

            ccf_alleles       = []

            for call in ccf_line[4:-1]:
                if call == '*':
                    continue

                el2 = call.split(';')

                max_prob = -1
                max_allele = 'N'
                for entry in el2:

                    a_info, prob = entry.split(':')
                    prob = float(prob)

                    if len(a_info) == 1:
                        allele = {a_info}
                    elif len(a_info) == 2:
                        allele = {a_info[0],a_info[1]}

                    if prob > max_prob:
                        max_prob = prob
                        max_allele = allele

                if len(max_allele) == 1:
                    ccf_alleles.append(max_allele)

            #has_allele_mix = False
            #for allele in ccf_alleles:
            #    if len(allele) == 2:
            #        has_allele_mix = True

            if genotype_prob < cutoff:
                continue

            for allele in ccf_alleles:

                if len(allele) > 1:
                    continue

                if allele <= cgi_dict[(ccf_chrom,ccf_pos)]:

                    match += 1
                    if allele != {ref_allele}:
                        snv_match += 1

                else:

                    mismatch += 1
                    if allele != {ref_allele}:
                        snv_mismatch += 1

                    print("{}\t{}\t{}\t{}\t{}".format(ccf_chrom,ccf_pos,''.join(allele),''.join(cgi_dict[(ccf_chrom,ccf_pos)]), ref_allele),file=mof)

    with open(result_outfile,'w') as rof:
        print("MATCH:        {}".format(match),file=rof)
        print("MISMATCH:     {}".format(mismatch),file=rof)
        print("ERR RATE:     {}".format(mismatch / (match + mismatch)),file=rof)
        print("SNV MATCH:    {}".format(snv_match),file=rof)
        print("SNV MISMATCH: {}".format(snv_mismatch),file=rof)
        print("SNV FDR:      {}".format(snv_mismatch / (snv_match + snv_mismatch)),file=rof)
