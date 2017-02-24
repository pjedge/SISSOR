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
import random
from DJSF import Union, Find, Node #MakeSet

count_key = namedtuple('count_key', 'is_SNP matches_CGI CGI_genotype')
count_key_sp = namedtuple('count_key_sp', 'is_SNP matches_CGI third_chamber_match')
count_key_ind = namedtuple('count_key_ind', 'cell is_SNP matches_CGI third_chamber_match')

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

# the original implementation of the haplotype strand-pairing code found it more
# convenient to list off all pairs of chambers that were paired,
# since the criteria for pairing did not result in disjoint sets necessarily.
# in the last implemenation they are always disjoint sets, but the format stayed the same.
# so this function takes a list of ((cell1,chamber1),(cell2,chamber2))
# strand pairs, and returns two disjoint sets for which chambers belong to the two haplotypes.
# note that P1,P2 for different positions are totally unrelated!
# they're arbitrary names unique to the position.
def hap_pairs_to_sets(haplotype_pairs):

    # make a dictionary to be able to find our disjoint-set-forest nodes
    nodedict = dict()
    elements = [x[0] for x in haplotype_pairs] + [x[1] for x in haplotype_pairs]
    for element in elements:
        nodedict[element] = Node(element)

    # union together paired fragments, saying "we know these are same haplotype"
    for ((cell1, chamber1), (cell2, chamber2)) in haplotype_pairs:
        Union(nodedict[(cell1,chamber1)],nodedict[(cell2,chamber2)])

    # parentdict has an arbitrary key, and the value is a set of same-haplotype elements
    parentdict = defaultdict(set)
    for v in nodedict.values():
        parentdict[Find(v).label].add(v.label)

    hapsets = list(parentdict.values())
    assert (len(hapsets) <= 2)
    P1 = set()
    P2 = set()
    if len(hapsets) > 0:
        P1 = hapsets[0]
    if len(hapsets) == 2:
        P2 = hapsets[1]

    return P1,P2

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


def generate_truth_dict(GFF_file, WGS_VCF_file, CGI_VCF_file, HG19, region, truth_dict_pickle):

    dp_pat = re.compile("DP=(\d+)")
    qr_pat = re.compile(";QR=(\d+)")
    type_pat = re.compile("TYPE=([^;]+);")
    region_chrom, region_start, region_end = region
    truth_dict = dict()
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

    pickle.dump(truth_dict,open(truth_dict_pickle,'wb'))

# chamber call file has results from SISSOR cross-chamber allele calls
# GFF file has set of known alleles for individual (for comparison)
# VCF file has another set of heterozygous SNVs for individual
# cutoff is the minimum probability of an allele call to use it
# result outfile simply contains the counts of match vs mismatched alleles
# mismatch_outfile prints locations of mismatched allele calls
def accuracy_count(chamber_call_file, truth_dict_pickle, GMS_file, cutoff, counts_pickle_file, mismatch_outfile, bedfile_filter=None): #WGS_VCF_file,

    truth_dict = pickle.load(open(truth_dict_pickle,'rb'))
    gms_set    = set()

    if FILTER_LOW_MAPPABILITY:
        with open(GMS_file,'r') as GMS:
            for line in GMS:
                if line[0] == '#' or len(line) < 3:
                    continue
                el = line.strip().split('\t')

                chrom     = el[0]
                pos       = int(el[1])
                gms_score = float(el[5])

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

# chamber call file has results from SISSOR cross-chamber allele calls
# GFF file has set of known alleles for individual (for comparison)
# VCF file has another set of heterozygous SNVs for individual
# cutoff is the minimum probability of an allele call to use it
# result outfile simply contains the counts of match vs mismatched alleles
# mismatch_outfile prints locations of mismatched allele calls
def accuracy_count_strand_pairing(chamber_call_file, truth_dict_pickle, GMS_file, cutoff, counts_pickle_file, mismatch_outfile, strand_mismatch_pickle, separate_haplotypes = True, mode='all_cell'): #WGS_VCF_file,

    if mode not in ['same_cell','all_cell','cross_cell','ind_same_cell']:
        raise ValueError('Invalid accuracy count mode.')

    same_cell  = False
    cross_cell = False
    if mode == 'same_cell' or mode == 'ind_same_cell':
        same_cell = True
    if mode == 'cross_cell':
        cross_cell = True

    truth_dict = pickle.load(open(truth_dict_pickle,'rb'))
    gms_set    = set()

    if FILTER_LOW_MAPPABILITY:
        with open(GMS_file,'r') as GMS:
            for line in GMS:
                if line[0] == '#' or len(line) < 3:
                    continue
                el = line.strip().split('\t')

                chrom     = el[0]
                pos       = int(el[1])
                gms_score = float(el[5])

                if gms_score > GMS_LIMIT:
                    gms_set.add((chrom, pos))

    counts = defaultdict(int)
    strand_mismatch_counts = defaultdict(int)

    with open(chamber_call_file,'r') as ccf, open(mismatch_outfile,'w') as mof:
        #print('chr\tpos\tsissor_call\tCGI_allele\tref',file=mof)
        for line in ccf:
            ccf_line = line.strip().split('\t')
            ccf_chrom = ccf_line[0]
            ccf_pos   = int(ccf_line[1])

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

                basecounts = defaultdict(int)
                for base in pileup:
                    basecounts[base] += 1

                max_base = 'N'
                max_v = -1
                for k,v in basecounts.items():
                    if v > max_v:
                        max_v = v
                        max_base = k

                # at very low coverages have to be sort of strict.
                if (len(pileup) == 3 or len(pileup) == 4) and max_v < 3:
                    continue

                el2 = xchamber_calls.split(';')

                max_prob = -1
                max_allele = {'N'}
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

                if len(max_allele) == 1 and max_allele == {max_base} and max_prob > cutoff and len(pileup) >= STRAND_MATCH_MIN_COV:
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

            P1,P2 = hap_pairs_to_sets(haplotype_pairs)

            if mode == 'ind_same_cell':
                hap_sets = []
                for hap_set in [P1,P2]:
                    cell_dict = defaultdict(set)
                    for hcell,hch in hap_set:
                        cell_dict[hcell].add((hcell,hch))

                    hap_sets += list(cell_dict.values())

            else:
                hap_sets = [P1,P2] if separate_haplotypes else [P1.union(P2)]


            random.shuffle(haplotype_pairs)

            for hap_set in hap_sets:

                if len(hap_set) == 0:
                    continue

                matched_chambers = None
                matched_allele = None
                third_chamber_match = False

                for (cell1,chamber1),(cell2,chamber2) in haplotype_pairs:

                    # both chambers should have an allele call.
                    if (cell1,chamber1) not in call_dict or (cell2,chamber2) not in call_dict:
                        continue

                    # skip pair if we should only be comparing within cells
                    if cell1 != cell2 and same_cell:
                        continue

                    # only consider pair if it's in the current haplotype we're considering
                    if not ((cell1,chamber1) in hap_set and (cell2,chamber2) in hap_set):
                        assert((cell1,chamber1) not in hap_set)   # if one not in the haplotype, neither should be
                        assert((cell2,chamber2) not in hap_set)
                        continue

                    # skip pair if we should only be comparing between cells
                    if cell1 == cell2 and cross_cell:
                        continue

                    allele1 = call_dict[(cell1,chamber1)]
                    allele2 = call_dict[(cell2,chamber2)]

                    if allele1 == allele2:

                        strand_mismatch_counts[((cell1,cell2),'match')] += 1

                        # don't replace the old strand pair if this one isn't a SNP
                        if matched_chambers != None and allele1 == {ref_allele}:
                            continue

                        matched_chambers = (cell1,chamber1),(cell2,chamber2)
                        matched_allele = allele1

                    else:
                        strand_mismatch_counts[((cell1,cell2),'mismatch')] += 1

                if matched_chambers == None:
                    continue

                for (cell,chamber),allele in call_dict.items():
                    if (cell,chamber) not in matched_chambers and allele == matched_allele:
                        third_chamber_match = True

                if (ccf_chrom,ccf_pos) not in truth_dict:
                    CGI_MATCH = None
                elif matched_allele <= truth_dict[(ccf_chrom,ccf_pos)]:
                    CGI_MATCH = 1
                else:
                    CGI_MATCH = 0

                cgi_allele = truth_dict[(ccf_chrom,ccf_pos)] if (ccf_chrom,ccf_pos) in truth_dict else {'None'}

                if CGI_MATCH == False:
                    tag = "MISMATCH:{}:{}".format(''.join(list(cgi_allele)),matched_chambers)
                    el = line.strip().split('\t')
                    if el[-1] == 'N/A':
                        el[-1] = tag
                    else:
                        el[-1] = '{};{}'.format(el[-1], tag)
                    print('\t'.join(el),file=mof)

                SNV = matched_allele != {ref_allele}

                if mode == 'ind_same_cell':
                    assert(matched_chambers[0][0] == matched_chambers[1][0]) # have to be same cell...
                    counts[count_key_ind(matched_chambers[0][0],SNV,CGI_MATCH,third_chamber_match)] += 1 # same as count_key_sp but also saves the cell.
                else:
                    counts[count_key_sp(SNV,CGI_MATCH,third_chamber_match)] += 1

    pickle.dump(counts,open(counts_pickle_file,'wb'))
    pickle.dump(strand_mismatch_counts,open(strand_mismatch_pickle,'wb'))

def accuracy_aggregate(counts_pickle_files, pickle_outfile):

    counts = defaultdict(int)

    for pfile in counts_pickle_files:
        temp_dict = pickle.load(open(pfile,'rb'))
        for k,v in temp_dict.items():
            counts[k] += v

    pickle.dump(counts,open(pickle_outfile,'wb'))
    #print("Strand Mismatch, not in CGI, divided by total:")
    #print(sum([x[1] for x in counts.items() if x[0][1] != None and x[0][2] == True]) / sum([x[1] for x in counts.items() if x[0][1] != None]))

# all credit to Alex Martelli on stackoverflow for this function
#http://stackoverflow.com/questions/2166818/python-how-to-check-if-an-object-is-an-instance-of-a-namedtuple
def isnamedtupleinstance(x):
    t = type(x)
    b = t.__bases__
    if len(b) != 1 or b[0] != tuple: return False
    f = getattr(t, '_fields', None)
    if not isinstance(f, tuple): return False
    return all(type(n)==str for n in f)

def generate_table(pickle_file, output_file):
    with open(output_file,'w') as outf:
        counts = pickle.load(open(pickle_file,'rb'))
        sample_key = list(counts.keys())[0]
        if isnamedtupleinstance(sample_key):
            header = list(sample_key._fields) + ['count']
            print('\t'.join(header),file=outf)
        for k,v in counts.items():
            if type(k)==tuple or isnamedtupleinstance(k):
                line = [(x if type(x)!=bool else int(x)) for x in k] + [v]
            else:
                line = [k,v]
            print('\t'.join([str(x) for x in line]), file=outf)
