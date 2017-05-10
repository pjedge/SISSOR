#! /bin/python3.5
# -*- coding: utf-8 -*-

# Author : Peter Edge
# Email  : pedge@eng.ucsd.edu

import sys

def prune_hapblock_file(hapblock_file, output_file, snp_conf_cutoff, split_conf_cutoff, use_refhap_heuristic):

    with open(hapblock_file,'r') as inf, open(output_file,'w') as of:
        blk_count = 0
        for line in inf:
            if 'BLOCK' in line:
                blk_count = 0
            if len(line) < 3 or 'BLOCK' in line or '****' in line:
                print(line,file=of,end='')
                continue

            el = line.strip().split()
            pruned_refhap_heuristic = int(el[8])
            split_conf = float(el[9]) if el[9] != '.' else 100
            snp_conf   = float(el[10]) if el[10] != '.' else 100

            if split_conf < split_conf_cutoff and blk_count >= 2:
                print('******** ',file=of)
                print('BLOCK: (from split)',file=of)

            if (use_refhap_heuristic and pruned_refhap_heuristic) or (snp_conf < snp_conf_cutoff):
                el[1] = '-'
                el[2] = '-'
                print('\t'.join(el),file=of)
            else:
                print(line,file=of,end='')

            blk_count += 1
