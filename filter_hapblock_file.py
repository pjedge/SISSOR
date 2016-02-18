#! /bin/python3.5
# -*- coding: utf-8 -*-

# Author : Peter Edge
# Email  : pedge@eng.ucsd.edu

def filter_hapblock_file(infile, outfile):

    with open(infile, 'r') as inf, open(outfile, 'w') as opf:

        for line in inf:
            # block headers and dividers
            if ('BLOCK' in line or '********' in line or len(line) < 2):
                print(line,file=opf,end='')
                continue
            # data lines
            el = line.strip().split()
            last_el  = el[-1]
            el2 = last_el.split(':')
            last_el2 = el2[-1]

            # we want to filter out lines that end in FV or whose last column's
            # last number is > 2.0
            if last_el2 != 'FV' and float(last_el2) < 2.0:
                print(line,file=opf,end='')
