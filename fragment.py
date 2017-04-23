#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  9 16:35:02 2016

@author: peter
"""

import sys
import fileIO
import error_rates
import pickle

class fragment:

    def __init__(self, seq, name, switch_errors=None, mismatch_errors = None):
        self.seq = seq                         # list of (snp index, genomic index, allele call, quality score) tuples
        self.name = name                       # fragment ID / name

    def __str__(self):
        fragstr = ''
        num_pairs = 0
        prev_snp_ix = -2
        qual = ' '
        for snp_ix, genome_ix, allele, q_char in self.seq:

            diff = snp_ix - prev_snp_ix

            if diff == 1:
                fragstr += allele
            else:
                num_pairs += 1
                fragstr += ' {} {}'.format(snp_ix+1, allele)

            prev_snp_ix = snp_ix
            qual += q_char

        fragstr += qual

        prefix = '{} {}'.format(num_pairs,self.name)
        fragstr = prefix + fragstr
        return fragstr

def read_fragment_matrix(frag_matrix, vcf_file):

    snp_ix = 0
    vcf_dict = dict()
    with open(vcf_file,'r') as infile:
        for line in infile:
            if line[:1] == '#':
                continue
            el = line.strip().split('\t')
            if len(el) < 5:
                continue

            genomic_pos = int(el[1])-1
            vcf_dict[snp_ix] = genomic_pos
            snp_ix += 1

    flist = []

    with open(frag_matrix,"r") as fm:
        for line in fm:
            if len(line) < 2:
                continue

            el = line.strip().split()

            num_blks      = int(el[0])
            name = el[1]

            call_list  = el[2:(2+2*num_blks)]              # extract base call part of line
            call_list  = zip(*[iter(call_list)]*2)             # list -> tuple list conversion: credit to http://stackoverflow.com/questions/23286254/convert-list-to-a-list-of-tuples-python
            call_list  = [(int(a)-1, b) for a,b in call_list]  # convert index to 0-based integer
            call_list2 = []

            for ix, blk in call_list:
                curr_ix = ix
                for a in blk:
                    call_list2.append((curr_ix, a))
                    curr_ix += 1

            qlist = el[-1]
            #qlist = [10**((ord(q) - 33) * -0.1) for q in qlist]

            alist= [(a,vcf_dict[a],b,c) for ((a,b),c) in zip(call_list2,qlist)]

            frag = fragment(alist,name)
            flist.append(frag)

    sorted_flist = sorted(flist,key=lambda x: x.seq[0][0])

    return sorted_flist

def write_fragment_matrix(flist,outfile):
    lines = []

    for f in flist:
        if len(f.seq) < 2:
            continue
        firstpos = f.seq[0][0]
        lines.append((firstpos, str(f)))

    lines.sort()

    with open(outfile, 'w') as opf:
        for firstpos, line in lines:
            print(line, file=opf)

def matrixify_flist(flist, o=None):

    max_ix = 0
    max_name = 0
    min_ix = min([f.seq[0][0] if len(f.seq) > 0 else int(1e9) for f in flist])

    for f in flist:

        if len(f.name) > max_name:
            max_name = len(f.name)

        for a in f.seq:

            if a[0] > max_ix:

                max_ix = a[0]

    if o != None:
        for f in flist:

            line = ['-'] * (max_ix+1-min_ix)
            for a in f.seq:
                line[a[0] - min_ix] = a[2]

            line = [f.name.ljust(max_name+1)] + line
            pline = ''.join(line)

            print(pline,file=o)
    else:
        for f in flist:

            line = ['-'] * (max_ix+1 - min_ix)
            for a in f.seq:
                line[a[0] - min_ix] = a[2]

            line = [f.name.ljust(max_name+1)] + line
            pline = ''.join(line)

            print(pline)

def overlap(f1, f2, amt=1):

    s1 = set([a[0] for a in f1.seq])
    s2 = set([a[0] for a in f2.seq])
    inter = set.intersection(s1,s2)

    return len(inter) >= amt

# compute error rates by using another haplotype block file as ground truth
def fragment_fragment_error_rate(sissor_frags, bac_frags, vcf_file, outfile, vis_outfile, pickle_outfile, chamber_filter=None):

    snp_ix = 0
    vcf_dict = dict()
    with open(vcf_file,'r') as infile:
        for line in infile:
            if line[:1] == '#':
                continue
            el = line.strip().split('\t')
            if len(el) < 5:
                continue

            genomic_pos = int(el[1])-1
            vcf_dict[snp_ix] = genomic_pos
            snp_ix += 1

    # parse and get stuff to compute error rates
    flist_sissor = read_fragment_matrix(sissor_frags,vcf_file)
    flist_bac = read_fragment_matrix(bac_frags,vcf_file)
    ref_name = fileIO.get_ref_name(vcf_file)
    num_snps    = fileIO.count_SNPs(vcf_file)
    num_covered = sum(error_rates.find_covered_positions(sissor_frags, num_snps))

    def flip(allele):
        if allele == '1':
            return '0'
        elif allele == '0':
            return '1'
        elif allele == '-':
            return '-'
        else:
            print("ERROR")
            exit(0)

    total_switch_count   = 0
    total_mismatch_count = 0
    total_poss_mm        = 0

    with open(outfile,'w') as OUTFILE, open(vis_outfile,'w') as VIS_OUTFILE:
        #print("FRAG_ID\tREF_NAME\tSWITCH_CT\tPOSS_SW\tMISMATCH_CT\tPOSS_MM\tSWITCH_LOC\tMISMATCH_LOC",file=OUTFILE)
        total_num_ref = 0
        total_num_alt = 0
        for frag_sissor in flist_sissor:

            if chamber_filter != None and ':CH{0:02d}:'.format(chamber_filter) not in frag_sissor.name:
                continue
            fragment_as_block = [(snp_ix, allele, flip(allele)) for (snp_ix,g_ix,allele,qual) in frag_sissor.seq]

            bac_overlaps      = [f for f in flist_bac if overlap(f, frag_sissor)]
            bac_blocks        = [[(snp_ix, allele, flip(allele)) for (snp_ix,g_ix,allele,qual) in bo.seq] for bo in bac_overlaps]

            #err = error_rates.error_rate_calc(bac_blocks, [fragment_as_block], vcf_file)

            #reslist = [fragment.name,ref_name,err.switch_count, err.poss_sw, err.mismatch_count, err.poss_mm,switchpos,mismatchpos]
            #swloc[f.name] = switchpos
            #mmloc[f.name] = mismatchpos
            #err_rates[f.name] = err.get_switch_mismatch_rate()
            #err_rates[name] = err.get_switch_rate()
            fragment_switch_count = 0
            fragment_mismatch_count = 0
            fragment_poss_mm = 0

            for f, blk in zip(bac_overlaps, bac_blocks):
                err, num_ref, num_alt = error_rates.error_rate_calc([blk], [fragment_as_block], vcf_file)
                total_num_ref += num_ref
                total_num_alt += num_alt
                print("{} {}".format(total_num_ref,total_num_alt))
                current_switch_count = err.get_switch_count()
                current_mismatch_count = err.get_mismatch_count()
                current_poss_mm = err.get_poss_mm()
                fragment_switch_count += current_switch_count
                fragment_mismatch_count += current_mismatch_count
                fragment_poss_mm += current_poss_mm

                errseq = []
                for sw_loc in err.switch_loc[ref_name]:
                    errseq.append((sw_loc,vcf_dict[sw_loc],'S',1))
                for mm_loc in err.mismatch_loc[ref_name]:
                    errseq.append((mm_loc,vcf_dict[mm_loc],'M',1))

                assert(len(errseq) == current_switch_count + current_mismatch_count)

                print("switches:   {}".format(current_switch_count),file=VIS_OUTFILE)
                print("mismatches: {}".format(current_mismatch_count),file=VIS_OUTFILE)
                print("positions:  {}".format(current_poss_mm),file=VIS_OUTFILE)
                print("error rate: {}".format((current_switch_count + current_mismatch_count)/current_poss_mm),file=VIS_OUTFILE) if current_poss_mm > 0 else 0

                errseq.sort()
                error_vis = fragment(errseq,"ERRORS")
                fragment_vis = [frag_sissor,f,error_vis]
                matrixify_flist(fragment_vis,VIS_OUTFILE)
                pos_str = " ".join(["{}:{}".format(ref_name,x[1]+1) for x in errseq])
                print("err pos:    {}".format(pos_str),file=VIS_OUTFILE)
                print("\n***********************************************************************\n",file=VIS_OUTFILE)

            err_rate = (fragment_switch_count + fragment_mismatch_count)/fragment_poss_mm if fragment_poss_mm > 0 else 0
            print("{}\t{}\t{}\t{}\t{}".format(frag_sissor.name,err_rate, fragment_switch_count,fragment_mismatch_count,fragment_poss_mm),file=OUTFILE)
                #print("{}:{}:{}\t".format(f.name,err.get_switch_count(),err.get_mismatch_count()),file=OUTFILE,end='')

            #print("\t".join([str(x) for x in reslist]),file=OUTFILE)
            total_switch_count += fragment_switch_count
            total_mismatch_count += fragment_mismatch_count
            total_poss_mm        += fragment_poss_mm


    print("switches:   {}".format(total_switch_count))
    print("mismatches: {}".format(total_mismatch_count))
    print("positions:  {}".format(total_poss_mm))
    print("error rate: {}".format((total_switch_count + total_mismatch_count)/total_poss_mm)) if total_poss_mm > 0 else 0
    print("num ref mm: {}".format(total_num_ref))
    print("num alt mm: {}".format(total_num_alt))

    res_tuple = (total_switch_count,total_mismatch_count,total_poss_mm)
    pickle.dump(res_tuple,open(pickle_outfile,'wb'))

    #return swloc, mmloc,err_rates

#fragment_fragment_error_rate('haplotyping/sissor_project/data/PGP1_ALL/fragmat/cov1_strict/chr20', 'haplotyping/sissor_project/data/BAC_frags/chr20', 'haplotyping/sissor_project/data/PGP1_VCFs_BACindex/chr20.vcf', 'sissor_bac_fragments_error.txt', 'sissor_bac_visualization.txt','sissor_bac_pickle.p',chamber_filter=None)
