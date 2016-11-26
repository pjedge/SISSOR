#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 18:29:04 2016

@author: peter
"""
import pickle
from pyliftover import LiftOver
import re
lo = LiftOver('hg38ToHg19.over.chain')
GMS_LIMIT = 0.5
chroms = set(['chr{}'.format(i) for i in range(1,23)])
bases = {'A','T','G','C'}

ref_allele_re = re.compile('ref_allele ([ACGT])')


# chamber call file has results from SISSOR cross-chamber allele calls
# GFF file has set of known alleles for individual (for comparison)
# VCF file has another set of heterozygous SNVs for individual
# cutoff is the minimum probability of an allele call to use it
# result outfile simply contains the counts of match vs mismatched alleles
# mismatch_outfile prints locations of mismatched allele calls
def accuracy_count(chamber_call_file, GFF_file, WGS_VCF_file, GMS_file, cutoff, counts_pickle_file, mismatch_outfile):

    dp_pat = re.compile("DP=(\d+)")
    qr_pat = re.compile(";QR=(\d+)")
    type_pat = re.compile("TYPE=(\w+);")

    truth_dict = dict()
    # add SNVs seen in WGS dataset to dictionary
    with open(WGS_VCF_file,'r') as vcf, open(lifted_WGS_VCF_file,'w') as lifted_vcf:
        for line in vcf:
            if line[0] == '#' or len(line) < 3:
                print(line,end='',file=lifted_vcf)
                continue

            pos_type = re.findall(type_pat,line)
            depth = int(float(re.findall(dp_pat,fields1)[0]))
            if depth < 20:
                continue

            # we only care about reference and SNPs, no indels
            if not(pos_type == [] or pos_type[0] == 'snp'):
                continue

            vcf_line = line.strip().split('\t')
            fields = vcf_line[9]
            if vcf_line[8][0:2] != 'GT':
                continue
            genotype = fields[0:3]
            qual = int(float(re.findall(qr_pat,fields1)[0])) if genotype == '0/0' else float(el[5])

            if qual < 100:
                continue

            vcf_chrom = vcf_line[0]
            vcf_pos = int(vcf_line[1])
            ref_allele = vcf_line[3]
            variant_allele = vcf_line[4][0]
            third_allele = 'N'
            if len(vcf_line[4]) == 3:
                third_allele = vcf_line[4][2]

            alleles = set()

            if '0' in genotype:
                assert(ref_allele in bases)
                alleles.add(ref_allele)
            if '1' in genotype:
                assert(variant_allele in bases)
                alleles.add(variant_allele)
            if '2' in genotype:
                assert(third_allele in bases)
                alleles.add(third_allele)

            truth_dict[(vcf_chrom,vcf_pos)] = alleles
            ref[(vcf_chrom,vcf_pos)] = ref_allele

    '''
    # add heterozygous SNVs observed in our second CGI dataset to the dictionary
    with open(CGI_VCF_file,'r') as vcf:
        for line in vcf:
            if line[0] == '#' or len(line) < 3:
                continue

            vcf_line = line.strip().split('\t')

            vcf_chrom = vcf_line[0]
            vcf_pos   = int(vcf_line[1])
            alleles = {vcf_line[3],vcf_line[4]}

            if (vcf_chrom,vcf_pos) in truth_dict:
                truth_dict[(vcf_chrom,vcf_pos)] = truth_dict[(vcf_chrom,vcf_pos)].union(alleles)
            else:
                truth_dict[(vcf_chrom,vcf_pos)] = alleles
    '''

    # add heterozygous alleles observed in one CGI dataset to a dictionary
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

            infoline = gff_line[8]
            info = infoline.strip().split(';')
            a_info = info[0]
            assert('alleles' in a_info)
            a_info = a_info[8:]

            ref_allele = None
            found = re.findall(ref_allele_re, infoline)
            if len(found) > 0:
                ref_allele = found[0]

            if len(a_info) == 1:
                alleles = {a_info}
            elif len(a_info) == 3:
                alleles = {a_info[0],a_info[2]}
            else:
                print("Unexpected number of alleles")
                assert(False)

            if (gff_chrom,gff_pos) in truth_dict:
                truth_dict[(gff_chrom,gff_pos)] = truth_dict[(gff_chrom,gff_pos)].union(alleles)
            else:
                truth_dict[(gff_chrom,gff_pos)] = alleles
                # we observed only 1 SNP allele and hasn't been seen in other datasets
                # add reference base because of this CGI dataset's base false negative rate
                if len(alleles) == 1 and ref_allele != None:
                    truth_dict[(gff_chrom,gff_pos)] = truth_dict[(gff_chrom,gff_pos)].union({ref_allele})

    with open(GMS_file,'r') as GMS:
        for line in GMS:
            if line[0] == '#' or len(line) < 3:
                continue
            el = line.strip().split('\t')

            chrom     = el[0]
            pos       = int(el[1])
            gms_score = float(el[5])

            if (chrom, pos) in truth_dict and gms_score < GMS_LIMIT:
                del truth_dict[(chrom, pos)]

    mismatch     = 0
    total_called = 0
    total_known  = 0
    snv_mismatch = 0
    snv_match    = 0

    with open(chamber_call_file,'r') as ccf, open(mismatch_outfile,'w') as mof:
        #print('chr\tpos\tsissor_call\tCGI_allele\tref',file=mof)
        for line in ccf:
            ccf_line = line.strip().split('\t')
            ccf_chrom = ccf_line[0]
            ccf_pos   = int(ccf_line[1])
            ref_allele = ccf_line[2]

            tags = ccf_line[76].split(';')

            if ccf_chrom == 'chrX' or ccf_chrom == 'chrY':
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

            if genotype_prob < cutoff:
                continue

            total_called += 1

            if (ccf_chrom,ccf_pos) not in truth_dict:
                continue

            total_known += 1

            ccf_alleles       = []

            for call in ccf_line[4:76]:
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
            is_mismatch = False
            is_SNV_match = False
            is_SNV_mismatch = False
            mismatch_allele = {'N'}
            for allele in ccf_alleles:

                # allele should be a subset of observed alleles
                # else we call it a mismatch
                if not (allele <= truth_dict[(ccf_chrom,ccf_pos)]):

                    if allele != {ref_allele}:
                        is_SNV_mismatch = True

                    is_mismatch = True
                    mismatch_allele = allele

                elif allele != {ref_allele} and allele <= truth_dict[(ccf_chrom,ccf_pos)]:

                    is_SNV_match = True

            if is_mismatch:
                mismatch += 1
                cgi_set = truth_dict[(ccf_chrom,ccf_pos)]
                line2 = "{}\tOTHERS: {}\tSISSOR: {}".format(line.strip(),''.join(cgi_set),''.join(mismatch_allele))
                print(line2,file=mof)

            if is_SNV_mismatch:
                snv_mismatch += 1
            elif is_SNV_match:
                snv_match += 1

    counts = (mismatch, total_known, snv_mismatch, snv_match, total_called, cutoff)
    pickle.dump(counts,open(counts_pickle_file,'wb'))

def accuracy_aggregate(counts_pickle_files, result_outfile):

    mismatch     = 0
    total_known  = 0
    total_called = 0
    snv_mismatch = 0
    snv_match    = 0
    cutoff = None
    for pfile in counts_pickle_files:
        (mismatch0, total_known0, snv_mismatch0, snv_match0, total_called0, cutoff0) = pickle.load(open(pfile,'rb'))
        mismatch     += mismatch0
        total_known  += total_known0
        total_called += total_called0
        snv_mismatch += snv_mismatch0
        snv_match    += snv_match0

        if cutoff == None:
            cutoff = cutoff0

    description = '''POSTERIOR CUTOFF is the cutoff used for the posterior probability
    of a positions genotype as computed by sissorhands tool.
    CALLS ABOVE CUTOFF refers to the total number of
    positions (ref or variant) called using sissorhands base caller above the cutoff.
    HAVE TRUTH is the number of positions (including ref calls) above the cutoff that overlap the truth dataset
    we are comparing against (combination of PGP1 WGS data and CGI data)
    MISMATCH refers to positions that had known SNPs in our truth dataset and
    sissorhands called an allele that did not match a known allele (i.e. one allele of
    a heterozygous or the single allele of a homozygous. For positions seen as a single
    allele variant in CGI data the reference was also considered valid). The ERROR RATE is
    the fraction of these two values, MISMATCH/(HAVE TRUTH)

    The SNV MATCH counts positions where sissorhands called an snv in some chamber
    and it is a known SNV (in truth dataset). SNV MISMATCH counts those that are not known SNVs.
    The SNV FDR is the SNP false discovery rate, snv_mismatch / (snv_mismatch + snv_match)
    '''

    with open(result_outfile,'w') as rof:
        print("POSTERIOR CUTOFF:   {}".format(cutoff),file=rof)
        print("CALLS ABOVE CUTOFF: {}".format(total_called),file=rof)
        print("---------------------------",file=rof)
        print("HAVE TRUTH:         {}".format(total_known),file=rof)
        print("MISMATCH:           {}".format(mismatch),file=rof)
        print("ERROR RATE:         {}".format(mismatch / total_known),file=rof)
        print("---------------------------",file=rof)
        print("SNV MATCH:          {}".format(snv_match),file=rof)
        print("SNV MISMATCH:       {}".format(snv_mismatch),file=rof)
        print("SNV FDR:            {}".format(snv_mismatch / (snv_mismatch + snv_match)),file=rof)
        print("---------------------------",file=rof)
        print(description)
