#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 18:29:04 2016

@author: peter
"""

# chamber call file has results from SISSOR cross-chamber allele calls
#
def base_calls_to_vcf(chamber_call_files, cutoff, output_vcf, output_ccf):
    
    with open(output_vcf,'w') as ovcf, open(output_ccf,'w') as occf:
        for chamber_call_file in chamber_call_files:
            with open(chamber_call_file,'r') as ccf:
                for raw_line in ccf:
                    
                    ccf_line = raw_line.strip().split('\t')
                    ccf_chrom = ccf_line[0]
                    ccf_pos   = int(ccf_line[1])
        
                    #if (ccf_chrom, ccf_pos) not in gms_set:  # poor mappability
                    #    continue
        
                    ref_allele = ccf_line[2]
        
                    tags = ccf_line[76].split(';')
        
                    if ccf_chrom == 'chrX' or ccf_chrom == 'chrY':
                        continue
        
                    if 'TOO_MANY_ALLELES' in tags or 'TOO_MANY_CHAMBERS' in tags or 'ADJACENT_INDEL_OR_CLIP' in tags:
                        continue
        
                    call = ccf_line[3]
        
                    if call == '*':
                        continue
        
                    el2 = call.split(';')
        
                    genotype_prob = -1
                    max_genotype = 'NN'
                    for entry in el2:
        
                        genotype, prob = entry.split(':')
                        prob = float(prob)
        
                        if prob > genotype_prob:
                            genotype_prob = prob
                            max_genotype = genotype
        
                    if genotype_prob < cutoff:
                        continue
        
                    a1, a2 = max_genotype
                    
                    if a1 == a2:  # not heterozygous
                        continue  
        
                    if a1 == ref_allele:
                        variant_allele = a2
                    elif a2 == ref_allele:
                        variant_allele = a1
                    else:
                        continue
        
                    vcf_line = '{}\t{}\t.\t{}\t{}\t100\tPASS\tSNP:99\tGT\t0/1'.format(ccf_chrom,ccf_pos,ref_allele,variant_allele)
                    
                    print(raw_line,end='',file=occf)
                    print(vcf_line,file=ovcf)  
        
