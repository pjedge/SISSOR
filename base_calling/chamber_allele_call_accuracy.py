#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 18:29:04 2016

@author: peter
"""
import pickle
import re
import pysam
from collections import defaultdict

GMS_LIMIT = 0.5
chroms = ['chr{}'.format(i) for i in range(1,23)]
bases = {'A','T','G','C'}
n_chambers = 24

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


def same_mixture_occurs_twice(pileups):

    mixtures_seen = set()
    
    for pileup in pileups:
        
        counter = defaultdict(int)
        for base in pileup:
            counter[base] += 1
            
        bases_seen = set()
        for base, v in counter.items():
            if v >= 3:
                bases_seen.add(base)
                
        if len(bases_seen) >= 3:
            return True
            
        elif len(bases_seen) == 2:
            
            mixture_str = ''.join(sorted(list(bases_seen)))
            
            if mixture_str in mixtures_seen:
                return True
            
            mixtures_seen.add(mixture_str)
            
    return False
    

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
def accuracy_count(chamber_call_file, GFF_file, WGS_VCF_file, CGI_VCF_file, GMS_file, HG19, region, cutoff, counts_pickle_file, mismatch_outfile, bedfile_filter): #WGS_VCF_file,

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

    fragment_boundaries = parse_bedfile(bedfile_filter)

    with open(chamber_call_file,'r') as ccf, open(mismatch_outfile,'w') as mof:
        #print('chr\tpos\tsissor_call\tCGI_allele\tref',file=mof)
        for line in ccf:
            ccf_line = line.strip().split('\t')
            ccf_chrom = ccf_line[0]
            ccf_pos   = int(ccf_line[1])

            if not (ccf_chrom == region_chrom and ccf_pos >= region_start and ccf_pos <= region_end):
                continue

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
                        
                        
                        
            if (ccf_chrom, ccf_pos) not in gms_set:  # poor mappability
                continue

            ref_allele = ccf_line[2]

            tags = ccf_line[80].split(';')

            if ccf_chrom == 'chrX' or ccf_chrom == 'chrY':
                continue

            if 'TOO_MANY_CHAMBERS' in tags or 'ADJACENT_INDEL_OR_CLIP' in tags: # 'TOO_MANY_ALLELES' in tags or 
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
            
            pileups = []
            for call in base_call_list:
                if call == '*':
                    continue

                xchamber_calls, basic_calls, pileup = call.split('|')

                pileups.append(pileup)

            
            if same_mixture_occurs_twice(pileups): # likely mapping issue
                continue
            
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
            total_called += 1

            if (ccf_chrom,ccf_pos) not in truth_dict:
                continue

            total_known += 1

            # alleles should be a subset of observed alleles
            # else we call it a mismatch
            if not (alleles <= truth_dict[(ccf_chrom,ccf_pos)]):

                if alleles != {ref_allele}:
                    snv_mismatch += 1

                mismatch += 1
                cgi_set = truth_dict[(ccf_chrom,ccf_pos)]
                tag = "OTHERS:{};SISSOR:{}".format(''.join(cgi_set),''.join(alleles))
                el = line.strip().split('\t')
                if el[-1] == 'N/A':
                    el[-1] = tag
                else:
                    el[-1] = '{};{}'.format(el[-1], tag)

                print('\t'.join(el),file=mof)

            elif alleles != {ref_allele} and alleles <= truth_dict[(ccf_chrom,ccf_pos)]:

                snv_match += 1

    counts = (mismatch, total_known, snv_mismatch, snv_match, total_called, cutoff)
    pickle.dump(counts,open(counts_pickle_file,'wb'))

# chamber call file has results from SISSOR cross-chamber allele calls
# GFF file has set of known alleles for individual (for comparison)
# VCF file has another set of heterozygous SNVs for individual
# cutoff is the minimum probability of an allele call to use it
# result outfile simply contains the counts of match vs mismatched alleles
# mismatch_outfile prints locations of mismatched allele calls
def accuracy_count_strand_pairing(chamber_call_file, GFF_file, WGS_VCF_file, CGI_VCF_file, GMS_file, HG19, region, cutoff, counts_pickle_file, mismatch_outfile): #WGS_VCF_file,

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

            if (ccf_chrom, ccf_pos) not in gms_set:  # poor mappability
                continue

            ref_allele = ccf_line[2]

            tags = ccf_line[80].split(';')

            if ccf_chrom == 'chrX' or ccf_chrom == 'chrY':
                continue

            if 'TOO_MANY_ALLELES' in tags or 'TOO_MANY_CHAMBERS' in tags or 'ADJACENT_INDEL_OR_CLIP' in tags:
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
            
            call_dict = defaultdict(lambda: None)
            
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
                
                if len(max_allele) == 1:
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
                    haplotype_pairs.append(((cell1,chamber1),(cell2,chamber2)))
                    
            for (cell1,chamber1),(cell2,chamber2) in haplotype_pairs:
                
                if call_dict[(cell1,chamber1)] == call_dict[(cell2,chamber2)]:
                    
                    alleles.append(call_dict[(cell1,chamber1)])
                

            if len(alleles) == 0:
                continue

            total_called += 1

            if (ccf_chrom,ccf_pos) not in truth_dict:
                continue

            total_known += 1

            # alleles should be a subset of observed alleles
            # else we call it a mismatch
            if not (alleles <= truth_dict[(ccf_chrom,ccf_pos)]):

                if alleles != {ref_allele}:
                    snv_mismatch += 1

                mismatch += 1
                cgi_set = truth_dict[(ccf_chrom,ccf_pos)]
                tag = "OTHERS:{};SISSOR:{}".format(''.join(cgi_set),''.join(alleles))
                el = line.strip().split('\t')
                if el[-1] == 'N/A':
                    el[-1] = tag
                else:
                    el[-1] = '{};{}'.format(el[-1], tag)

                print('\t'.join(el),file=mof)

            elif allele != {ref_allele} and alleles <= truth_dict[(ccf_chrom,ccf_pos)]:

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

    #description = '''POSTERIOR CUTOFF is the cutoff used for the posterior probability
    #of a positions genotype as computed by sissorhands tool.
    #CALLS ABOVE CUTOFF refers to the total number of
    #positions (ref or variant) called using sissorhands base caller above the cutoff.
    #HAVE TRUTH is the number of positions (including ref calls) above the cutoff that overlap the truth dataset
    #we are comparing against (combination of PGP1 WGS data and CGI data)
    #MISMATCH refers to positions that had known SNPs in our truth dataset and
    #sissorhands called an allele that did not match a known allele (i.e. one allele of
    #a heterozygous or the single allele of a homozygous. For positions seen as a single
    #allele variant in CGI data the reference was also considered valid). The ERROR RATE is
    #the fraction of these two values, MISMATCH/(HAVE TRUTH)

    #The SNV MATCH counts positions where sissorhands called an snv in some chamber
    #and it is a known SNV (in truth dataset). SNV MISMATCH counts those that are not known SNVs.
    #The SNV FDR is the SNP false discovery rate, snv_mismatch / (snv_mismatch + snv_match)
    #'''

    description = '''POSTERIOR CUTOFF is the cutoff used for the posterior probability
    that an allele is present as computed by sissorhands tool.
    CALLS ABOVE CUTOFF refers to the total number of
    positions (ref or variant) called using sissorhands base caller above the cutoff.
    HAVE TRUTH is the number of positions (including ref calls) above the cutoff that overlap the truth dataset
    we are comparing against (PGP1 CGI data)
    MISMATCH refers to positions that had known SNPs in our truth dataset and
    sissorhands called an allele that did not match a known allele (i.e. one allele of
    a heterozygous or the single allele of a homozygous. For positions seen as a single
    allele variant in CGI data the reference was also considered valid). The ERROR RATE is
    the fraction of these two values, MISMATCH/(HAVE TRUTH)

    The SNV MATCH counts positions where sissorhands called an snv in some chamber
    and it is a known SNV (in truth dataset). SNV MISMATCH counts those that are not known SNVs.
    The SNV FDR is the SNP false discovery rate, snv_mismatch / (snv_mismatch + snv_match)
    '''

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
        print(description)
