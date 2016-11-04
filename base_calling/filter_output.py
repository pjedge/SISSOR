#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 18:29:04 2016

@author: peter
"""
import sys

infile = sys.argv[1]
CUTOFF = 0.999999

with open(infile,'r') as inf:
    for line in inf:
        el = line.strip().split()
        chrom = el[0]
        pos   = el[1]
        ref_allele = el[2]
        
        for call in el[3:-1]:
            if call == '*':
                continue
            
            el2 = call.split(';')
            
            for entry in el2:
                
                allele, prob = entry.split(':')
                prob = float(prob)
                
                if prob > CUTOFF and allele != ref_allele and len(allele) == 2:
                    print(line,end='')
                    
