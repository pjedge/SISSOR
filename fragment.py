#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  9 16:35:02 2016

@author: peter
"""



class fragment:

    def __init__(self, seq, name, switch_errors=None, mismatch_errors = None):
        self.seq = seq                         # list of (snp index, genomic index, allele call, quality score) tuples
        self.name = name                       # fragment ID / name
        self.switch_errors = switch_errors     # list of SNP index positions where switch errors were found
        self.mismatch_errors = mismatch_errors # list of SNP index positions where mismatch errors occured
        self.haplotype = None
        
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

def matrixify_flist(flist, outfile=None):

    max_ix = 0
    max_name = 0

    for f in flist:

        if len(f.name) > max_name:
            max_name = len(f.name)

        for a in f.seq:

            if a[0] > max_ix:

                max_ix = a[0]

    if outfile != None:
        with open(outfile,'w') as o:
            for f in flist:

                line = ['-'] * (max_ix+1)
                for a in f.seq:
                    line[a[0]] = a[2]

                line = [f.name.ljust(max_name+1)] + line
                pline = ''.join(line)

                print(pline,file=o)
    else:
        for f in flist:

            line = ['-'] * (max_ix+1)
            for a in f.seq:
                line[a[0]] = a[2]

            line = [f.name.ljust(max_name+1)] + line
            pline = ''.join(line)

            print(pline)

# compute error rates by using another haplotype block file as ground truth
def fragment_hapblock_error_rate(truth_file, frag_file, vcf_file, outfile,chamber_filter=None):

    import sys
    sys.path.append("/home/peter/git/HapTools")
    import fileIO
    import error_rates

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
    flist = read_fragment_matrix(frag_file,vcf_file)
    new_flist = []
    for f in flist:
        new_seq = [(pos,gpos,call,qual) for (pos,gpos,call,qual) in f.seq if call != '2']
        new_flist.append(fragment(new_seq,f.name))
    flist = new_flist

    t_blocklist = fileIO.parse_hapblock_file(truth_file,use_SNP_index=True)
    num_snps    = fileIO.count_SNPs(vcf_file)
    num_covered = sum(error_rates.find_covered_positions(frag_file, num_snps))
    ref_name    = fileIO.get_ref_name(vcf_file)
    # compute error result object, update the runtime and AN50 / completeness
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

    swloc = dict()
    mmloc = dict()
    err_rates = dict()

    with open(outfile,'w') as OUTFILE:
        print("FRAG_ID\tREF_NAME\tSWITCH_CT\tPOSS_SW\tMISMATCH_CT\tPOSS_MM\tSWITCH_LOC\tMISMATCH_LOC",file=OUTFILE)
        for f in flist:

            if chamber_filter != None and ':CH{0:02d}:'.format(chamber_filter) not in f.name:
                continue
            fragment_as_block = [(snp_ix, allele, flip(allele)) for (snp_ix,g_ix,allele,qual) in f.seq]

            err = error_rates.error_rate_calc(t_blocklist, [fragment_as_block], num_snps, num_covered, ref_name, track_error_positions=True)

            switchpos   = [vcf_dict[x]+1 for x in err.switch_loc[ref_name]]
            mismatchpos = [vcf_dict[x]+1 for x in err.mismatch_loc[ref_name]]
            #reslist = [fragment.name,ref_name,err.switch_count, err.poss_sw, err.mismatch_count, err.poss_mm,switchpos,mismatchpos]
            swloc[f.name] = switchpos
            mmloc[f.name] = mismatchpos
            err_rates[f.name] = err.get_switch_mismatch_rate()
            #err_rates[name] = err.get_switch_rate()

            #print("\t".join([str(x) for x in reslist]),file=OUTFILE)

    #return err_rates
    return swloc, mmloc,err_rates
