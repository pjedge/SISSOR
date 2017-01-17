#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 18:05:15 2016

@author: peter
"""
import sys
import pysam
import copy

def fix_sam(insam, outbam):


    new_SQ = [{'LN': 249250621, 'SN': 'chr1'},
    {'LN': 243199373, 'SN': 'chr2'},
    {'LN': 198022430, 'SN': 'chr3'},
    {'LN': 191154276, 'SN': 'chr4'},
    {'LN': 180915260, 'SN': 'chr5'},
    {'LN': 171115067, 'SN': 'chr6'},
    {'LN': 159138663, 'SN': 'chr7'},
    {'LN': 146364022, 'SN': 'chr8'},
    {'LN': 141213431, 'SN': 'chr9'},
    {'LN': 135534747, 'SN': 'chr10'},
    {'LN': 135006516, 'SN': 'chr11'},
    {'LN': 133851895, 'SN': 'chr12'},
    {'LN': 115169878, 'SN': 'chr13'},
    {'LN': 107349540, 'SN': 'chr14'},
    {'LN': 102531392, 'SN': 'chr15'},
    {'LN': 90354753, 'SN': 'chr16'},
    {'LN': 81195210, 'SN': 'chr17'},
    {'LN': 78077248, 'SN': 'chr18'},
    {'LN': 59128983, 'SN': 'chr19'},
    {'LN': 63025520, 'SN': 'chr20'},
    {'LN': 48129895, 'SN': 'chr21'},
    {'LN': 51304566, 'SN': 'chr22'},
    {'LN': 155270560, 'SN': 'chrX'},
    {'LN': 59373566, 'SN': 'chrY'},
    {'LN': 16569, 'SN': 'chrM'}]

    infile = pysam.AlignmentFile(insam,"rb");
    new_header = copy.deepcopy(infile.header)

    convert_name = dict()
    for i in range(1,23):
        convert_name[str(i)] = i-1
    convert_name['X'] = 22
    convert_name['Y'] = 23
    convert_name['MT'] = 24
    new_header['SQ'] = new_SQ
    outfile = pysam.AlignmentFile(outbam,"wb", header=new_header);

    for record in infile:
        if record.reference_id != -1 and record.reference_name in convert_name:
            record.reference_id = convert_name[record.reference_name]
            outfile.write(record)
    infile.close()
    outfile.close()


if __name__ == '__main__':

    insam = sys.argv[1]
    outbam = sys.argv[2]
    fix_sam(insam,outbam)
