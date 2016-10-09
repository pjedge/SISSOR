#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 18:29:04 2016

@author: peter
"""

# chamber call file has results from SISSOR cross-chamber allele calls
# GFF file has set of known alleles for individual (for comparison)
# cutoff is the minimum probability of an allele call to use it
# result outfile simply contains the counts of match vs mismatched alleles
# mismatch_outfile prints locations of mismatched allele calls
def chamber_allele_call_accuracy(chamber_call_file, GFF_file, cutoff, result_outfile, mismatch_outfile):

    with open(chamber_call_file,'r') as ccf, open(GFF_file,'r') as gff, open(result_outfile,'w') as rof, open(mismatch_outfile,'w') as mof:
        print('chr\tpos\tsissor_call\tCGI_allele\tref',file=mof)
        chroms = ['chr{}'.format(x) for x in range(1,23)] + ['chrX','chrY']
        def read_ccf_line():
            ccf_line  = ccf.readline()
            if not ccf_line or len(ccf_line) < 3:
                return None,None,None,None
                
            ccf_line = ccf_line.strip().split('\t')
            ccf_chrom = ccf_line[0]
            ccf_pos   = int(ccf_line[1])
            ref_allele = ccf_line[2]

            alleles = []

            for call in ccf_line[3:-1]:
                if call == '*':
                    continue

                el2 = call.split(';')

                for entry in el2:

                    a_info, prob = entry.split(':')
                    prob = float(prob)

                    if prob > cutoff:
                        if len(a_info) == 1:
                            allele = {a_info}
                        elif len(a_info) == 2:
                            allele = {a_info[0],a_info[1]}
                        alleles.append(allele)

            return ccf_chrom, ccf_pos, alleles, ref_allele

        def read_gff_line():
            gff_line  = gff.readline()
            if not gff_line:
                return None,None,None,None
            elif gff_line[0] == '#': # header line
                return read_gff_line()
            gff_line = gff_line.strip().split('\t')
            gff_chrom = gff_line[0]
            gff_pos   = int(gff_line[3])
            if(gff_pos != int(gff_line[4])):
                return read_gff_line()

            info = gff_line[8].strip().split(';')
            a_info = info[0]
            assert('alleles' in a_info)
            a_info = a_info[8:]

            r_info = info[2]
            assert('ref_allele' in r_info)
            ref_allele = r_info[10:]

            if len(a_info) == 1:
                alleles = {a_info}
            elif len(a_info) == 3:
                alleles = {a_info[0],a_info[2]}
            else:
                print("Unexpected number of alleles")
                assert(False)

            return gff_chrom, gff_pos, alleles, ref_allele

        match = 0
        mismatch = 0
            
        i = 0
        ccf_chrom = None
        gff_chrom = None
        ref_allele = None
        while i < len(chroms):
            
            chrom = chroms[i]

            while ccf_chrom != chrom:
                ccf_chrom, ccf_pos, ccf_alleles, ref_allele = read_ccf_line()
                if not ccf_chrom:
                    break
            while gff_chrom != chrom:
                gff_chrom, gff_pos, gff_alleles, ref_allele2 = read_gff_line()
                if not gff_chrom:
                    break                
            # leave these out for now
            if chrom == 'chrX' or chrom == 'chrY':
                i += 1
                continue
            # in this loop, both files are considering positions of same chromosome
            while True:
    
                while ccf_chrom == chrom and ccf_pos < gff_pos:
                    ccf_chrom, ccf_pos, ccf_alleles, ref_allele = read_ccf_line()
                while gff_chrom == chrom and gff_pos < ccf_pos:
                    gff_chrom, gff_pos, gff_alleles, ref_allele2 = read_gff_line()
                                
                if gff_chrom != chrom or ccf_chrom != chrom:
                    i += 1
                    break
                
                assert(ccf_chrom == gff_chrom)
                assert(ccf_pos   == gff_pos)
                
                # do stuff with data lines
                
                for allele in ccf_alleles:
                    
                    if allele <= gff_alleles:
                        
                        match += 1
                    
                    else:
                        
                        mismatch += 1
                        print("{}\t{}\t{}\t{}\t{}".format(ccf_chrom,ccf_pos,''.join(allele),''.join(gff_alleles), ref_allele),file=mof)
                
                
                # we've processed a position so move both files ahead one
                ccf_chrom, ccf_pos, ccf_alleles, ref_allele  = read_ccf_line()
                gff_chrom, gff_pos, gff_alleles, ref_allele2 = read_gff_line()
        
        
        print("MATCH:    {}".format(match),file=rof)
        print("MISMATCH: {}".format(mismatch),file=rof)