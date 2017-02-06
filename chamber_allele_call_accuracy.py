#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 18:29:04 2016

@author: peter
"""
import pickle
import re
import pysam
import itertools
from collections import defaultdict, namedtuple

count_key = namedtuple('count_key', 'is_SNP matches_CGI CGI_genotype')

GMS_LIMIT = 0.5
chroms = ['chr{}'.format(i) for i in range(1,23)]
bases = {'A','T','G','C'}
n_chambers = 24

FILTER_LOW_MAPPABILITY = False
STRAND_MATCH_MIN_COV = 3

chr_num = dict()
for i,chrom in enumerate(chroms):
    chr_num[chrom] = i

ref_allele_re = re.compile('ref_allele ([ACGT])')

def split_vcf(input_vcf, chunklist, output_vcfs):

    regions_output = list(zip(chunklist, output_vcfs))
    (chrom, start, stop), outputfile = regions_output.pop(0)
    output = open(outputfile, 'w')
    with open(input_vcf,'r') as vcf:
        for line in vcf:
            if line[0] == '#' or len(line) < 3:
                continue

            vcf_line = line.strip().split('\t')

            vcf_chrom = vcf_line[0]
            vcf_pos = int(vcf_line[1])

            if vcf_chrom != chrom or vcf_pos > stop:

                done = False
                while not (vcf_chrom == chrom and vcf_pos >= start and vcf_pos <= stop):
                    output.close()
                    if len(regions_output) == 0:
                        done = True
                        break
                    (chrom, start, stop), outputfile = regions_output.pop(0)
                    output = open(outputfile, 'w')
                if done:
                    break

            # write to lifted over vcf
            print(line,end='',file=output)

    if not output.closed:
        output.close()

# this takes the genome mappability score (GMS) files and splits them into
# convenient 5 Mb chunks, corresponding to regions called by sissorhands
def split_gms(infiles, chunklist, outputlst):

    chunks_output = list(zip(chunklist, outputlst))

    for infile in infiles:
        (chrom, start, stop), outputfile = chunks_output.pop(0)
        assert(start == 1)              # should be beginning of chrom
        assert(chrom+'.gms' in infile)  # file should have chrom name
        output = open(outputfile,'w')
        with open(infile,'r') as gms:
            for line in gms:
                if line[0] == '#' or len(line) < 3:
                    continue
                el = line.strip().split('\t')

                gms_chrom     = el[0]
                gms_pos       = int(el[1])

                if gms_pos > stop:
                    output.close()
                    if chunks_output == []:
                        break
                    (chrom, start, stop), outputfile = chunks_output.pop(0)
                    output = open(outputfile,'w')

                # make sure we're really on the current region
                assert(gms_chrom == chrom)
                assert(gms_pos >= start)
                assert(gms_pos <= stop)

                print(line,end='',file=output)

    if not output.closed:
        output.close()

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


# chamber call file has results from SISSOR cross-chamber allele calls
# GFF file has set of known alleles for individual (for comparison)
# VCF file has another set of heterozygous SNVs for individual
# cutoff is the minimum probability of an allele call to use it
# result outfile simply contains the counts of match vs mismatched alleles
# mismatch_outfile prints locations of mismatched allele calls
def accuracy_count(chamber_call_file, GFF_file, WGS_VCF_file, CGI_VCF_file, GMS_file, HG19, region, cutoff, counts_pickle_file, mismatch_outfile, bedfile_filter=None): #WGS_VCF_file,

    dp_pat = re.compile("DP=(\d+)")
    qr_pat = re.compile(";QR=(\d+)")
    type_pat = re.compile("TYPE=([^;]+);")
    region_chrom, region_start, region_end = region
    truth_dict = dict()
    gms_set    = set()
    hg19_fasta = pysam.FastaFile(HG19)

    # add heterozygous alleles observed in one CGI dataset to a dictionary
    with open(GFF_file,'r') as gff:
        for line in gff:
            if line[0] == '#' or len(line) < 3:
                continue
            gff_line = line.strip().split('\t')

            if gff_line[2] == 'SNP':

                gff_chrom = gff_line[0]
                gff_pos   = int(gff_line[3])
                if not (gff_chrom == region_chrom and gff_pos >= region_start and gff_pos <= region_end):
                    continue

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
                    #if len(alleles) == 1 and ref_allele != None:
                    #    truth_dict[(gff_chrom,gff_pos)] = truth_dict[(gff_chrom,gff_pos)].union({ref_allele})

            elif gff_line[2] == 'REF':

                gff_chrom   = gff_line[0]
                gff_start   = int(gff_line[3])
                gff_end     = int(gff_line[4])
                if not ((gff_chrom == region_chrom and gff_start >= region_start and gff_start <= region_end) or
                        (gff_chrom == region_chrom and gff_end >= region_start and gff_end <= region_end)):
                    continue

                for gff_pos in range(gff_start, gff_end+1):

                    ref_lookup = str.upper(hg19_fasta.fetch(gff_chrom, gff_pos-1, gff_pos))  # pysam uses 0-index

                    if ref_lookup not in bases:  # e.g. 'N' alleles
                        continue

                    if (gff_chrom,gff_pos) in truth_dict:
                        truth_dict[(gff_chrom,gff_pos)] = truth_dict[(gff_chrom,gff_pos)].union({ref_lookup})
                    else:
                        truth_dict[(gff_chrom,gff_pos)] = {ref_lookup}


    # add heterozygous SNVs observed in our second CGI dataset to the dictionary
    with open(CGI_VCF_file,'r') as vcf:
        for line in vcf:
            if line[0] == '#' or len(line) < 3:
                continue

            vcf_line = line.strip().split('\t')

            vcf_chrom = vcf_line[0]
            vcf_pos   = int(vcf_line[1])

            if not (vcf_chrom == region_chrom and vcf_pos >= region_start and vcf_pos <= region_end):
                continue

            alleles = {vcf_line[3],vcf_line[4]}

            if (vcf_chrom,vcf_pos) in truth_dict:
                truth_dict[(vcf_chrom,vcf_pos)] = truth_dict[(vcf_chrom,vcf_pos)].union(alleles)
            else:
                truth_dict[(vcf_chrom,vcf_pos)] = alleles


    # add SNVs seen in WGS dataset to dictionary
    with open(WGS_VCF_file,'r') as vcf:
        for line in vcf:
            if line[0] == '#' or len(line) < 3:
                continue

            vcf_line = line.strip().split('\t')
            fields = vcf_line[7]
            pos_type = re.findall(type_pat,line)
            depth = int(float(re.findall(dp_pat,fields)[0]))
            if depth < 20:
                continue

            # we only care about reference and SNPs, no indels or MNPs
            if not(pos_type == [] or pos_type[0] == 'snp'):
                continue
            if not (len(vcf_line[3]) == 1 and len(vcf_line[4]) == 1):
                continue

            if vcf_line[8][0:2] != 'GT':
                continue
            genotype = vcf_line[9][0:3]
            qual = int(float(re.findall(qr_pat,fields)[0])) if genotype == '0/0' else float(vcf_line[5])

            if qual < 100:
                continue


            vcf_chrom = vcf_line[0]
            vcf_pos = int(vcf_line[1])
            ref_allele = str.upper(vcf_line[3])
            variant_allele = str.upper(vcf_line[4][0])
            third_allele = 'N'
            if len(vcf_line[4]) == 3:
                third_allele = str.upper(vcf_line[4][2])

            alleles = set()

            if '0' in genotype:
                alleles.add(ref_allele)
            if '1' in genotype:
                alleles.add(variant_allele)
            if '2' in genotype:
                alleles.add(third_allele)

            # skip cases with N alleles, or with any other unusual cases
            if len(alleles.intersection(bases)) != len(alleles):
                continue

            if (vcf_chrom,vcf_pos) in truth_dict:
                truth_dict[(vcf_chrom,vcf_pos)] = truth_dict[(vcf_chrom,vcf_pos)].union(alleles)
            #else:
            #    truth_dict[(vcf_chrom,vcf_pos)] = alleles


    with open(GMS_file,'r') as GMS:
        for line in GMS:
            if line[0] == '#' or len(line) < 3:
                continue
            el = line.strip().split('\t')

            chrom     = el[0]
            pos       = int(el[1])
            gms_score = float(el[5])

            if not (chrom == region_chrom and pos >= region_start and pos <= region_end):
                continue

            if gms_score > GMS_LIMIT:
                gms_set.add((chrom, pos))

    counts = defaultdict(int)

    fragment_boundaries = parse_bedfile(bedfile_filter) if bedfile_filter != None else None

    with open(chamber_call_file,'r') as ccf, open(mismatch_outfile,'w') as mof:
        #print('chr\tpos\tsissor_call\tCGI_allele\tref',file=mof)
        for line in ccf:
            ccf_line = line.strip().split('\t')
            ccf_chrom = ccf_line[0]
            ccf_pos   = int(ccf_line[1])

            if not (ccf_chrom == region_chrom and ccf_pos >= region_start and ccf_pos <= region_end):
                continue

            if bedfile_filter != None:
                #######################################################################
                #######################################################################
                # filter out positions not covered by any fragment. Want to know if accuracy within fragments is better
                if fragment_boundaries == []:
                    continue

                f_chrom, f_start, f_end = fragment_boundaries[0]

                # if we're behind fragment start, skip this spot
                if chr_num[ccf_chrom] < chr_num[f_chrom] or (ccf_chrom == f_chrom and ccf_pos < f_start):
                    continue

                # if we're ahead of fragment start, skip to later fragment boundaries
                while 1:
                    if fragment_boundaries == []:
                        break
                    f_chrom, f_start, f_end = fragment_boundaries[0]
                    if chr_num[ccf_chrom] > chr_num[f_chrom] or (ccf_chrom == f_chrom and ccf_pos >= f_end):
                        fragment_boundaries.pop(0)
                    else:
                        break

                # if we're not inside fragment, continue
                if not(ccf_chrom == f_chrom and ccf_pos >= f_start and ccf_pos < f_end):
                    continue
                #######################################################################
                #######################################################################

            if FILTER_LOW_MAPPABILITY and (ccf_chrom, ccf_pos) not in gms_set:  # poor mappability
                continue

            ref_allele = ccf_line[2]

            tags = ccf_line[80].split(';')

            if 'TOO_MANY_ALLELES' in tags or 'TOO_MANY_CHAMBERS' in tags or 'SAME_MIXTURE_OCCURS_TWICE' in tags: #or 'ADJACENT_INDEL_OR_CLIP' in tags: # 'TOO_MANY_ALLELES' in tags or
                continue

            call = ccf_line[3]

            if call == '*':
                continue

            el2 = call.split(';')

            alleles = set()

            for entry in el2:

                allele, prob = entry.split(':')
                prob = float(prob)

                if prob > cutoff:
                    alleles.add(allele)

            if len(alleles) == 0:
                continue

            base_call_list = [x for x in ccf_line[5:80] if 'CELL' not in x]

            for call in base_call_list:
                if call == '*':
                    continue

                xchamber_calls, basic_calls, pileup = call.split('|')



            '''
            ccf_alleles       = []
            has_allele_mix = False

            base_call_list = [x for x in ccf_line[5:80] if 'CELL' not in x]

            for call in base_call_list:
                if call == '*':
                    continue

                xchamber_calls, basic_calls, pileup = call.split('|')

                el2 = xchamber_calls.split(';')

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

                if len(max_allele) == 2:
                    has_allele_mix = True
                    break
                elif len(max_allele) == 1:
                    ccf_alleles.append(max_allele)

            if has_allele_mix:
                continue
            '''

            is_SNP = (alleles != {ref_allele})
            matches_CGI  = None
            CGI_genotype = None
            if (ccf_chrom,ccf_pos) in truth_dict:
                matches_CGI = (alleles <= truth_dict[(ccf_chrom,ccf_pos)])
                if truth_dict[(ccf_chrom,ccf_pos)] == {ref_allele}:
                    CGI_genotype = '0/0'
                elif len(truth_dict[(ccf_chrom,ccf_pos)]) == 1:
                    CGI_genotype = '1/1'
                elif len(truth_dict[(ccf_chrom,ccf_pos)]) == 2:
                    CGI_genotype = '0/1'
                else:
                    CGI_genotype = '!'

            counts[count_key(is_SNP,matches_CGI,CGI_genotype)] += 1

            # alleles should be a subset of observed alleles
            # else we call it a mismatch
            if (ccf_chrom,ccf_pos) in truth_dict and not (alleles <= truth_dict[(ccf_chrom,ccf_pos)]):

                cgi_set = truth_dict[(ccf_chrom,ccf_pos)]
                tag = "OTHERS:{};SISSOR:{}".format(''.join(cgi_set),''.join(alleles))
                el = line.strip().split('\t')
                if el[-1] == 'N/A':
                    el[-1] = tag
                else:
                    el[-1] = '{};{}'.format(el[-1], tag)

                print('\t'.join(el),file=mof)

    #counts = (mismatch, total_known, snv_mismatch, snv_match, total_called, cutoff)
    pickle.dump(counts,open(counts_pickle_file,'wb'))
'''
# chamber call file has results from SISSOR cross-chamber allele calls
# GFF file has set of known alleles for individual (for comparison)
# VCF file has another set of heterozygous SNVs for individual
# cutoff is the minimum probability of an allele call to use it
# result outfile simply contains the counts of match vs mismatched alleles
# mismatch_outfile prints locations of mismatched allele calls
def accuracy_count_strand_pairing(chamber_call_file, GFF_file, WGS_VCF_file, CGI_VCF_file, GMS_file, HG19, region, cutoff, counts_pickle_file, mismatch_outfile, same_cell_only=False): #WGS_VCF_file,

    dp_pat = re.compile("DP=(\d+)")
    qr_pat = re.compile(";QR=(\d+)")
    type_pat = re.compile("TYPE=([^;]+);")
    region_chrom, region_start, region_end = region
    truth_dict = dict()
    gms_set    = set()
    hg19_fasta = pysam.FastaFile(HG19)

    # add heterozygous alleles observed in one CGI dataset to a dictionary
    with open(GFF_file,'r') as gff:
        for line in gff:
            if line[0] == '#' or len(line) < 3:
                continue
            gff_line = line.strip().split('\t')

            if gff_line[2] == 'SNP':

                gff_chrom = gff_line[0]
                gff_pos   = int(gff_line[3])
                if not (gff_chrom == region_chrom and gff_pos >= region_start and gff_pos <= region_end):
                    continue

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

            elif gff_line[2] == 'REF':

                gff_chrom   = gff_line[0]
                gff_start   = int(gff_line[3])
                gff_end     = int(gff_line[4])
                if not ((gff_chrom == region_chrom and gff_start >= region_start and gff_start <= region_end) or
                        (gff_chrom == region_chrom and gff_end >= region_start and gff_end <= region_end)):
                    continue

                for gff_pos in range(gff_start, gff_end+1):

                    ref_lookup = str.upper(hg19_fasta.fetch(gff_chrom, gff_pos-1, gff_pos))  # pysam uses 0-index

                    if ref_lookup not in bases:  # e.g. 'N' alleles
                        continue

                    if (gff_chrom,gff_pos) in truth_dict:
                        truth_dict[(gff_chrom,gff_pos)] = truth_dict[(gff_chrom,gff_pos)].union({ref_lookup})
                    else:
                        truth_dict[(gff_chrom,gff_pos)] = {ref_lookup}


    # add heterozygous SNVs observed in our second CGI dataset to the dictionary
    with open(CGI_VCF_file,'r') as vcf:
        for line in vcf:
            if line[0] == '#' or len(line) < 3:
                continue

            vcf_line = line.strip().split('\t')

            vcf_chrom = vcf_line[0]
            vcf_pos   = int(vcf_line[1])

            if not (vcf_chrom == region_chrom and vcf_pos >= region_start and vcf_pos <= region_end):
                continue

            alleles = {vcf_line[3],vcf_line[4]}

            if (vcf_chrom,vcf_pos) in truth_dict:
                truth_dict[(vcf_chrom,vcf_pos)] = truth_dict[(vcf_chrom,vcf_pos)].union(alleles)
            else:
                truth_dict[(vcf_chrom,vcf_pos)] = alleles


    # add SNVs seen in WGS dataset to dictionary
    with open(WGS_VCF_file,'r') as vcf:
        for line in vcf:
            if line[0] == '#' or len(line) < 3:
                continue

            vcf_line = line.strip().split('\t')
            fields = vcf_line[7]
            pos_type = re.findall(type_pat,line)
            depth = int(float(re.findall(dp_pat,fields)[0]))
            if depth < 20:
                continue

            # we only care about reference and SNPs, no indels or MNPs
            if not(pos_type == [] or pos_type[0] == 'snp'):
                continue
            if not (len(vcf_line[3]) == 1 and len(vcf_line[4]) == 1):
                continue

            if vcf_line[8][0:2] != 'GT':
                continue
            genotype = vcf_line[9][0:3]
            qual = int(float(re.findall(qr_pat,fields)[0])) if genotype == '0/0' else float(vcf_line[5])

            if qual < 100:
                continue


            vcf_chrom = vcf_line[0]
            vcf_pos = int(vcf_line[1])
            ref_allele = str.upper(vcf_line[3])
            variant_allele = str.upper(vcf_line[4][0])
            third_allele = 'N'
            if len(vcf_line[4]) == 3:
                third_allele = str.upper(vcf_line[4][2])

            alleles = set()

            if '0' in genotype:
                alleles.add(ref_allele)
            if '1' in genotype:
                alleles.add(variant_allele)
            if '2' in genotype:
                alleles.add(third_allele)

            # skip cases with N alleles, or with any other unusual cases
            if len(alleles.intersection(bases)) != len(alleles):
                continue

            if (vcf_chrom,vcf_pos) in truth_dict:
                truth_dict[(vcf_chrom,vcf_pos)] = truth_dict[(vcf_chrom,vcf_pos)].union(alleles)
            #else:
            #    truth_dict[(vcf_chrom,vcf_pos)] = alleles


    with open(GMS_file,'r') as GMS:
        for line in GMS:
            if line[0] == '#' or len(line) < 3:
                continue
            el = line.strip().split('\t')

            chrom     = el[0]
            pos       = int(el[1])
            gms_score = float(el[5])

            if not (chrom == region_chrom and pos >= region_start and pos <= region_end):
                continue

            if gms_score > GMS_LIMIT:
                gms_set.add((chrom, pos))

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

            if not (ccf_chrom == region_chrom and ccf_pos >= region_start and ccf_pos <= region_end):
                continue

            if FILTER_LOW_MAPPABILITY and (ccf_chrom, ccf_pos) not in gms_set:  # poor mappability
                continue

            ref_allele = ccf_line[2]

            tags = ccf_line[80].split(';')

            if 'TOO_MANY_ALLELES' in tags or 'SAME_MIXTURE_OCCURS_TWICE' in tags or 'TOO_MANY_CHAMBERS' in tags: #or 'ADJACENT_INDEL_OR_CLIP' in tags:
                continue

            call = ccf_line[3]

            if call == '*':
                continue

            el2 = call.split(';')

            alleles = []
            #for entry in el2:

            #    allele, prob = entry.split(':')
            #    prob = float(prob)

            #    if prob > cutoff:
            #        alleles.add(allele)
                    #max_genotype = genotype

            base_call_list = [x for x in ccf_line[5:80] if 'CELL' not in x]

            for call in base_call_list:
                if call == '*':
                    continue

                xchamber_calls, basic_calls, pileup = call.split('|')

            call_dict = dict()

            for i, call in enumerate(base_call_list):

                cell = int(i / n_chambers)
                chamber = int(i % n_chambers)

                if call == '*':
                    continue

                xchamber_calls, basic_calls, pileup = call.split('|')

                el2 = xchamber_calls.split(';')

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

                if len(max_allele) == 1 and max_prob > cutoff:
                    call_dict[(cell,chamber)] = max_allele


            haplotype_pairs = []
            for tag in tags:
                if tag[0:3] == 'HP:':
                    foo, strand1, strand2 = tag.split(':')
                    cell1,chamber1 = strand1.split(',')
                    cell2,chamber2 = strand2.split(',')
                    cell1 = int(cell1)
                    chamber1 = int(chamber1)
                    cell2 = int(cell2)
                    chamber2 = int(chamber2)

                    if cell1 != cell2 and same_cell_only:
                        continue

                    haplotype_pairs.append(((cell1,chamber1),(cell2,chamber2)))

            for (cell1,chamber1),(cell2,chamber2) in haplotype_pairs:
                if (cell1,chamber1) not in call_dict or (cell2,chamber2) not in call_dict:
                    continue

                if call_dict[(cell1,chamber1)] == call_dict[(cell2,chamber2)]:

                    alleles.append(call_dict[(cell1,chamber1)])


            if len(alleles) == 0:
                continue

            # alleles should be a subset of observed alleles
            # else we call it a mismatch
            for allele in alleles:

                total_called += 1

                if (ccf_chrom,ccf_pos) not in truth_dict:
                    continue

                total_known += 1

                if not (allele <= truth_dict[(ccf_chrom,ccf_pos)]):

                    if allele != {ref_allele}:
                        snv_mismatch += 1

                    mismatch += 1
                    cgi_set = truth_dict[(ccf_chrom,ccf_pos)]
                    tag = "OTHERS:{};SISSOR:{}".format(''.join(cgi_set),''.join(allele))
                    el = line.strip().split('\t')
                    if el[-1] == 'N/A':
                        el[-1] = tag
                    else:
                        el[-1] = '{};{}'.format(el[-1], tag)

                    print('\t'.join(el),file=mof)

                elif allele != {ref_allele} and allele <= truth_dict[(ccf_chrom,ccf_pos)]:

                    snv_match += 1

    counts = (mismatch, total_known, snv_mismatch, snv_match, total_called, cutoff)
    pickle.dump(counts,open(counts_pickle_file,'wb'))

'''

# chamber call file has results from SISSOR cross-chamber allele calls
# GFF file has set of known alleles for individual (for comparison)
# VCF file has another set of heterozygous SNVs for individual
# cutoff is the minimum probability of an allele call to use it
# result outfile simply contains the counts of match vs mismatched alleles
# mismatch_outfile prints locations of mismatched allele calls
def accuracy_count_strand_pairing_genomic(chamber_call_file, GFF_file, WGS_VCF_file, CGI_VCF_file, GMS_file, HG19, region, cutoff, counts_pickle_file, mismatch_outfile, same_cell_only=False): #WGS_VCF_file,

    dp_pat = re.compile("DP=(\d+)")
    qr_pat = re.compile(";QR=(\d+)")
    type_pat = re.compile("TYPE=([^;]+);")
    region_chrom, region_start, region_end = region
    truth_dict = dict()
    gms_set    = set()
    hg19_fasta = pysam.FastaFile(HG19)

    # add heterozygous alleles observed in one CGI dataset to a dictionary
    with open(GFF_file,'r') as gff:
        for line in gff:
            if line[0] == '#' or len(line) < 3:
                continue
            gff_line = line.strip().split('\t')

            if gff_line[2] == 'SNP':

                gff_chrom = gff_line[0]
                gff_pos   = int(gff_line[3])
                if not (gff_chrom == region_chrom and gff_pos >= region_start and gff_pos <= region_end):
                    continue

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

            elif gff_line[2] == 'REF':

                gff_chrom   = gff_line[0]
                gff_start   = int(gff_line[3])
                gff_end     = int(gff_line[4])
                if not ((gff_chrom == region_chrom and gff_start >= region_start and gff_start <= region_end) or
                        (gff_chrom == region_chrom and gff_end >= region_start and gff_end <= region_end)):
                    continue

                for gff_pos in range(gff_start, gff_end+1):

                    ref_lookup = str.upper(hg19_fasta.fetch(gff_chrom, gff_pos-1, gff_pos))  # pysam uses 0-index

                    if ref_lookup not in bases:  # e.g. 'N' alleles
                        continue

                    if (gff_chrom,gff_pos) in truth_dict:
                        truth_dict[(gff_chrom,gff_pos)] = truth_dict[(gff_chrom,gff_pos)].union({ref_lookup})
                    else:
                        truth_dict[(gff_chrom,gff_pos)] = {ref_lookup}


    # add heterozygous SNVs observed in our second CGI dataset to the dictionary
    with open(CGI_VCF_file,'r') as vcf:
        for line in vcf:
            if line[0] == '#' or len(line) < 3:
                continue

            vcf_line = line.strip().split('\t')

            vcf_chrom = vcf_line[0]
            vcf_pos   = int(vcf_line[1])

            if not (vcf_chrom == region_chrom and vcf_pos >= region_start and vcf_pos <= region_end):
                continue

            alleles = {vcf_line[3],vcf_line[4]}

            if (vcf_chrom,vcf_pos) in truth_dict:
                truth_dict[(vcf_chrom,vcf_pos)] = truth_dict[(vcf_chrom,vcf_pos)].union(alleles)
            else:
                truth_dict[(vcf_chrom,vcf_pos)] = alleles


    # add SNVs seen in WGS dataset to dictionary
    with open(WGS_VCF_file,'r') as vcf:
        for line in vcf:
            if line[0] == '#' or len(line) < 3:
                continue

            vcf_line = line.strip().split('\t')
            fields = vcf_line[7]
            pos_type = re.findall(type_pat,line)
            depth = int(float(re.findall(dp_pat,fields)[0]))
            if depth < 20:
                continue

            # we only care about reference and SNPs, no indels or MNPs
            if not(pos_type == [] or pos_type[0] == 'snp'):
                continue
            if not (len(vcf_line[3]) == 1 and len(vcf_line[4]) == 1):
                continue

            if vcf_line[8][0:2] != 'GT':
                continue
            genotype = vcf_line[9][0:3]
            qual = int(float(re.findall(qr_pat,fields)[0])) if genotype == '0/0' else float(vcf_line[5])

            if qual < 100:
                continue


            vcf_chrom = vcf_line[0]
            vcf_pos = int(vcf_line[1])
            ref_allele = str.upper(vcf_line[3])
            variant_allele = str.upper(vcf_line[4][0])
            third_allele = 'N'
            if len(vcf_line[4]) == 3:
                third_allele = str.upper(vcf_line[4][2])

            alleles = set()

            if '0' in genotype:
                alleles.add(ref_allele)
            if '1' in genotype:
                alleles.add(variant_allele)
            if '2' in genotype:
                alleles.add(third_allele)

            # skip cases with N alleles, or with any other unusual cases
            if len(alleles.intersection(bases)) != len(alleles):
                continue

            if (vcf_chrom,vcf_pos) in truth_dict:
                truth_dict[(vcf_chrom,vcf_pos)] = truth_dict[(vcf_chrom,vcf_pos)].union(alleles)
            #else:
            #    truth_dict[(vcf_chrom,vcf_pos)] = alleles


    with open(GMS_file,'r') as GMS:
        for line in GMS:
            if line[0] == '#' or len(line) < 3:
                continue
            el = line.strip().split('\t')

            chrom     = el[0]
            pos       = int(el[1])
            gms_score = float(el[5])

            if not (chrom == region_chrom and pos >= region_start and pos <= region_end):
                continue

            if gms_score > GMS_LIMIT:
                gms_set.add((chrom, pos))

    counts = defaultdict(int)

    with open(chamber_call_file,'r') as ccf, open(mismatch_outfile,'w') as mof:
        #print('chr\tpos\tsissor_call\tCGI_allele\tref',file=mof)
        for line in ccf:
            ccf_line = line.strip().split('\t')
            ccf_chrom = ccf_line[0]
            ccf_pos   = int(ccf_line[1])

            if not (ccf_chrom == region_chrom and ccf_pos >= region_start and ccf_pos <= region_end):
                continue

            if FILTER_LOW_MAPPABILITY and (ccf_chrom, ccf_pos) not in gms_set:  # poor mappability
                continue

            ref_allele = ccf_line[2]

            tags = ccf_line[80].split(';')

            if 'TOO_MANY_ALLELES' in tags or 'SAME_MIXTURE_OCCURS_TWICE' in tags or 'TOO_MANY_CHAMBERS' in tags:# or 'ADJACENT_INDEL_OR_CLIP' in tags:
                continue

            call = ccf_line[3]

            if call == '*':
                continue

            el2 = call.split(';')

            alleles = set()
            #for entry in el2:

            #    allele, prob = entry.split(':')
            #    prob = float(prob)

            #    if prob > cutoff:
            #        alleles.add(allele)
                    #max_genotype = genotype

            base_call_list = [x for x in ccf_line[5:80] if 'CELL' not in x]

            for call in base_call_list:
                if call == '*':
                    continue

                xchamber_calls, basic_calls, pileup = call.split('|')


            call_dict = dict()

            for i, call in enumerate(base_call_list):

                cell = int(i / n_chambers)
                chamber = int(i % n_chambers)

                if call == '*':
                    continue

                xchamber_calls, basic_calls, pileup = call.split('|')

                el2 = xchamber_calls.split(';')

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

                if max_prob > cutoff and len(pileup) >= STRAND_MATCH_MIN_COV:
                    call_dict[(cell,chamber)] = max_allele
                else:

                    # try to call position with multiple bases instead.
                    max_prob = -1
                    max_allele = 'N'
                    for b1,b2 in itertools.combinations(bases,2):

                        totalprob = 0
                        for entry in el2:

                            a_info, prob = entry.split(':')
                            prob = float(prob)

                            if len(a_info) == 1:
                                allele = {a_info}
                            elif len(a_info) == 2:
                                allele = {a_info[0],a_info[1]}

                            if allele <= {b1,b2}:
                                totalprob += prob

                        if totalprob > max_prob:
                            max_prob = totalprob
                            max_allele = {b1,b2}

                    if max_prob > cutoff and len(pileup) >= STRAND_MATCH_MIN_COV:
                        call_dict[(cell,chamber)] = max_allele

            haplotype_pairs = []
            for tag in tags:
                if tag[0:3] == 'HP:':
                    foo, strand1, strand2 = tag.split(':')
                    cell1,chamber1 = strand1.split(',')
                    cell2,chamber2 = strand2.split(',')
                    cell1 = int(cell1)
                    chamber1 = int(chamber1)
                    cell2 = int(cell2)
                    chamber2 = int(chamber2)

                    if cell1 != cell2 and same_cell_only:
                        continue

                    haplotype_pairs.append(((cell1,chamber1),(cell2,chamber2)))

            STRAND_MISMATCH = False
            HALF_MATCH = False
            has_a_pair = False
            mm_list = []
            for (cell1,chamber1),(cell2,chamber2) in haplotype_pairs:
                if (cell1,chamber1) not in call_dict or (cell2,chamber2) not in call_dict:
                    continue

                alleles = alleles.union(call_dict[(cell1,chamber1)])
                alleles = alleles.union(call_dict[(cell2,chamber2)])

                has_a_pair = True
                L1 = len(call_dict[(cell1,chamber1)])
                L2 = len(call_dict[(cell2,chamber2)])

                hm = (call_dict[(cell1,chamber1)] <= call_dict[(cell2,chamber2)] or call_dict[(cell2,chamber2)] <= call_dict[(cell1,chamber1)])

                if ((L1 == 1 and L2 == 2) or (L1 == 2 and L2 == 1)) and hm:
                    HALF_MATCH = True
                    str_allele1 = ''.join(call_dict[(cell1,chamber1)])
                    str_allele2 = ''.join(call_dict[(cell2,chamber2)])
                    mm_list.append(('half_match',cell1,chamber1,str_allele1,cell2,chamber2,str_allele2))

                elif not hm:

                    STRAND_MISMATCH = True
                    str_allele1 = ''.join(call_dict[(cell1,chamber1)])
                    str_allele2 = ''.join(call_dict[(cell2,chamber2)])
                    mm_list.append(('mismatch',cell1,chamber1,str_allele1,cell2,chamber2,str_allele2))

            # require at least one pair of haplotype-matched calls at position
            if not has_a_pair:
                continue

            if len(alleles) == 0:
                continue

            if STRAND_MISMATCH or HALF_MATCH:
                tag = "STRAND_MISMATCHES:{}".format(mm_list)
                el = line.strip().split('\t')
                if el[-1] == 'N/A':
                    el[-1] = tag
                else:
                    el[-1] = '{};{}'.format(el[-1], tag)
                print('\t'.join(el),file=mof)

            if (ccf_chrom,ccf_pos) not in truth_dict:
                CGI_MATCH = None
            elif alleles <= truth_dict[(ccf_chrom,ccf_pos)]:
                CGI_MATCH = 1
            else:
                CGI_MATCH = 0

            SNV = alleles != {ref_allele}

            counts[(SNV,CGI_MATCH,STRAND_MISMATCH,HALF_MATCH)] += 1

    pickle.dump(counts,open(counts_pickle_file,'wb'))


# chamber call file has results from SISSOR cross-chamber allele calls
# GFF file has set of known alleles for individual (for comparison)
# VCF file has another set of heterozygous SNVs for individual
# cutoff is the minimum probability of an allele call to use it
# result outfile simply contains the counts of match vs mismatched alleles
# mismatch_outfile prints locations of mismatched allele calls
def accuracy_count_single_strand_for_ROC(chamber_call_file, GFF_file, WGS_VCF_file, CGI_VCF_file, GMS_file, HG19, region, cutoff, counts_pickle_file, mismatch_outfile, bedfile_filter=None): #WGS_VCF_file,

    dp_pat = re.compile("DP=(\d+)")
    qr_pat = re.compile(";QR=(\d+)")
    type_pat = re.compile("TYPE=([^;]+);")
    region_chrom, region_start, region_end = region
    truth_dict = dict()
    gms_set    = set()
    hg19_fasta = pysam.FastaFile(HG19)

    # add heterozygous alleles observed in one CGI dataset to a dictionary
    with open(GFF_file,'r') as gff:
        for line in gff:
            if line[0] == '#' or len(line) < 3:
                continue
            gff_line = line.strip().split('\t')

            if gff_line[2] == 'SNP':

                gff_chrom = gff_line[0]
                gff_pos   = int(gff_line[3])
                if not (gff_chrom == region_chrom and gff_pos >= region_start and gff_pos <= region_end):
                    continue

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


            elif gff_line[2] == 'REF':

                gff_chrom   = gff_line[0]
                gff_start   = int(gff_line[3])
                gff_end     = int(gff_line[4])
                if not ((gff_chrom == region_chrom and gff_start >= region_start and gff_start <= region_end) or
                        (gff_chrom == region_chrom and gff_end >= region_start and gff_end <= region_end)):
                    continue

                for gff_pos in range(gff_start, gff_end+1):

                    ref_lookup = str.upper(hg19_fasta.fetch(gff_chrom, gff_pos-1, gff_pos))  # pysam uses 0-index

                    if ref_lookup not in bases:  # e.g. 'N' alleles
                        continue

                    if (gff_chrom,gff_pos) in truth_dict:
                        truth_dict[(gff_chrom,gff_pos)] = truth_dict[(gff_chrom,gff_pos)].union({ref_lookup})
                    else:
                        truth_dict[(gff_chrom,gff_pos)] = {ref_lookup}


    # add heterozygous SNVs observed in our second CGI dataset to the dictionary
    with open(CGI_VCF_file,'r') as vcf:
        for line in vcf:
            if line[0] == '#' or len(line) < 3:
                continue

            vcf_line = line.strip().split('\t')

            vcf_chrom = vcf_line[0]
            vcf_pos   = int(vcf_line[1])

            if not (vcf_chrom == region_chrom and vcf_pos >= region_start and vcf_pos <= region_end):
                continue

            alleles = {vcf_line[3],vcf_line[4]}

            if (vcf_chrom,vcf_pos) in truth_dict:
                truth_dict[(vcf_chrom,vcf_pos)] = truth_dict[(vcf_chrom,vcf_pos)].union(alleles)
            else:
                truth_dict[(vcf_chrom,vcf_pos)] = alleles

    '''
    # add SNVs seen in WGS dataset to dictionary
    with open(WGS_VCF_file,'r') as vcf:
        for line in vcf:
            if line[0] == '#' or len(line) < 3:
                continue

            vcf_line = line.strip().split('\t')
            fields = vcf_line[7]
            pos_type = re.findall(type_pat,line)
            depth = int(float(re.findall(dp_pat,fields)[0]))
            if depth < 20:
                continue

            # we only care about reference and SNPs, no indels or MNPs
            if not(pos_type == [] or pos_type[0] == 'snp'):
                continue
            if not (len(vcf_line[3]) == 1 and len(vcf_line[4]) == 1):
                continue

            if vcf_line[8][0:2] != 'GT':
                continue
            genotype = vcf_line[9][0:3]
            qual = int(float(re.findall(qr_pat,fields)[0])) if genotype == '0/0' else float(vcf_line[5])

            if qual < 100:
                continue


            vcf_chrom = vcf_line[0]
            vcf_pos = int(vcf_line[1])
            ref_allele = str.upper(vcf_line[3])
            variant_allele = str.upper(vcf_line[4][0])
            third_allele = 'N'
            if len(vcf_line[4]) == 3:
                third_allele = str.upper(vcf_line[4][2])

            alleles = set()

            if '0' in genotype:
                alleles.add(ref_allele)
            if '1' in genotype:
                alleles.add(variant_allele)
            if '2' in genotype:
                alleles.add(third_allele)

            # skip cases with N alleles, or with any other unusual cases
            if len(alleles.intersection(bases)) != len(alleles):
                continue

            if (vcf_chrom,vcf_pos) in truth_dict:
                truth_dict[(vcf_chrom,vcf_pos)] = truth_dict[(vcf_chrom,vcf_pos)].union(alleles)
            #else:
            #    truth_dict[(vcf_chrom,vcf_pos)] = alleles


    with open(GMS_file,'r') as GMS:
        for line in GMS:
            if line[0] == '#' or len(line) < 3:
                continue
            el = line.strip().split('\t')

            chrom     = el[0]
            pos       = int(el[1])
            gms_score = float(el[5])

            if not (chrom == region_chrom and pos >= region_start and pos <= region_end):
                continue

            if gms_score > GMS_LIMIT:
                gms_set.add((chrom, pos))
    '''

    counts = defaultdict(int)

    fragment_boundaries = parse_bedfile(bedfile_filter) if bedfile_filter != None else None

    with open(chamber_call_file,'r') as ccf, open(mismatch_outfile,'w') as mof:
        #print('chr\tpos\tsissor_call\tCGI_allele\tref',file=mof)
        for line in ccf:
            ccf_line = line.strip().split('\t')
            ccf_chrom = ccf_line[0]
            ccf_pos   = int(ccf_line[1])

            if not (ccf_chrom == region_chrom and ccf_pos >= region_start and ccf_pos <= region_end):
                continue

            if bedfile_filter != None:
                #######################################################################
                #######################################################################
                # filter out positions not covered by any fragment. Want to know if accuracy within fragments is better
                if fragment_boundaries == []:
                    continue

                f_chrom, f_start, f_end = fragment_boundaries[0]

                # if we're behind fragment start, skip this spot
                if chr_num[ccf_chrom] < chr_num[f_chrom] or (ccf_chrom == f_chrom and ccf_pos < f_start):
                    continue

                # if we're ahead of fragment start, skip to later fragment boundaries
                while 1:
                    if fragment_boundaries == []:
                        break
                    f_chrom, f_start, f_end = fragment_boundaries[0]
                    if chr_num[ccf_chrom] > chr_num[f_chrom] or (ccf_chrom == f_chrom and ccf_pos >= f_end):
                        fragment_boundaries.pop(0)
                    else:
                        break

                # if we're not inside fragment, continue
                if not(ccf_chrom == f_chrom and ccf_pos >= f_start and ccf_pos < f_end):
                    continue
                #######################################################################
                #######################################################################

            if FILTER_LOW_MAPPABILITY and (ccf_chrom, ccf_pos) not in gms_set:  # poor mappability
                continue

            ref_allele = ccf_line[2]

            tags = ccf_line[80].split(';')

            if 'TOO_MANY_ALLELES' in tags or 'TOO_MANY_CHAMBERS' in tags or 'SAME_MIXTURE_OCCURS_TWICE' in tags: #or 'ADJACENT_INDEL_OR_CLIP' in tags: # 'TOO_MANY_ALLELES' in tags or
                continue

            call = ccf_line[3]

            if call == '*':
                continue

            ccf_alleles       = []
            has_allele_mix = False

            base_call_list = [x for x in ccf_line[5:80] if 'CELL' not in x]

            for call in base_call_list:
                if call == '*':
                    continue

                xchamber_calls, basic_calls, pileup = call.split('|')

                el2 = xchamber_calls.split(';')

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

            # alleles should be a subset of observed alleles
            # else we call it a mismatch

            for alleles in ccf_alleles:

                is_SNP = (alleles != {ref_allele})
                matches_CGI  = None
                CGI_genotype = None
                if (ccf_chrom,ccf_pos) in truth_dict:
                    matches_CGI = (alleles <= truth_dict[(ccf_chrom,ccf_pos)])
                    if truth_dict[(ccf_chrom,ccf_pos)] == {ref_allele}:
                        CGI_genotype = '0/0'
                    elif len(truth_dict[(ccf_chrom,ccf_pos)]) == 1:
                        CGI_genotype = '1/1'
                    elif len(truth_dict[(ccf_chrom,ccf_pos)]) == 2:
                        CGI_genotype = '0/1'
                    else:
                        CGI_genotype = '!'

                counts[count_key(is_SNP,matches_CGI,CGI_genotype)] += 1

                if (ccf_chrom,ccf_pos) in truth_dict and not (alleles <= truth_dict[(ccf_chrom,ccf_pos)]):

                    cgi_set = truth_dict[(ccf_chrom,ccf_pos)]
                    tag = "CGI:{};SISSOR:{}".format(''.join(cgi_set),''.join(alleles))
                    el = line.strip().split('\t')
                    if el[-1] == 'N/A':
                        el[-1] = tag
                    else:
                        el[-1] = '{};{}'.format(el[-1], tag)

                    print('\t'.join(el),file=mof)

    pickle.dump(counts,open(counts_pickle_file,'wb'))


def accuracy_aggregate(counts_pickle_files, pickle_outfile):

    counts = defaultdict(int)

    for pfile in counts_pickle_files:
        temp_dict = pickle.load(open(pfile,'rb'))
        for k,v in temp_dict.items():
            counts[k] += v

    pickle.dump(counts,open(pickle_outfile,'wb'))
    #print("Strand Mismatch, not in CGI, divided by total:")
    #print(sum([x[1] for x in counts.items() if x[0][1] != None and x[0][2] == True]) / sum([x[1] for x in counts.items() if x[0][1] != None]))
'''
def generate_table(pickle_list):

    countlist = []
    # for each new file, we add its cutoff
    for pickle_file in pickle_list:

        countlist.append(pickle.load(open(pickle_file,'rb')))

    table  = [list(countlist[0]._fields)]

    uniq = set()

    for count in countlist:
        if uniq = uniq.union(set(count.keys()))
'''


'''
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


    error_rate = mismatch / total_known if total_known > 0 else 'N/A'
    snv_fdr    = snv_mismatch / (snv_mismatch + snv_match) if (snv_mismatch + snv_match) > 0 else 'N/A'

    with open(result_outfile,'w') as rof:
        print("POSTERIOR CUTOFF:   {}".format(cutoff),file=rof)
        print("CALLS ABOVE CUTOFF: {}".format(total_called),file=rof)
        print("---------------------------",file=rof)
        print("HAVE TRUTH:         {}".format(total_known),file=rof)
        print("MISMATCH:           {}".format(mismatch),file=rof)
        print("ERROR RATE:         {}".format(error_rate),file=rof)
        print("---------------------------",file=rof)
        print("SNV MATCH:          {}".format(snv_match),file=rof)
        print("SNV MISMATCH:       {}".format(snv_mismatch),file=rof)
        print("SNV FDR:            {}".format(snv_fdr),file=rof)
        print("---------------------------",file=rof)
        #print(description)
'''