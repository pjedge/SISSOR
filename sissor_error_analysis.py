#! /bin/python3.5
# -*- coding: utf-8 -*-

# Author : Peter Edge
# Email  : pedge@eng.ucsd.edu
import argparse
import sys
from filter_hapblock_file import filter_hapblock_file
from create_sissor_fragments import create_sissor_fragments
from fragment_error_rate import fragment_error_rate

desc = '''creates a fragment matrix from SISSOR and then funnels this into fragment_error_rate'''


# parse arguments
def parse_args():

    parser = argparse.ArgumentParser(description='')
    parser.add_argument('sissor_vcf_file', nargs='?', type = str, help='path to haploid vcf file. present alleles are variant calls. missing alleles are assumed to be ref')
    parser.add_argument('bed_file', nargs='?', type = str, help='bed file where each entry specifies the boundary of a SISSOR fragment')
    parser.add_argument('hapcut_block_file', nargs='?', type = str, help='path to hapcut block file')
    parser.add_argument('filtered_hapcut_block_file', nargs='?', type = str, help='path to write filtered hapcut block file to')
    parser.add_argument('intermediate_fragment_file', nargs='?', type = str, help='path to write intermediate fragment matrix file')
    parser.add_argument('output_file', nargs='?', type = str, help='path to write error analysis output')

    # default to help option. credit to unutbu: http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
    if len(sys.argv) < 7:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    return args

def sissor_error_analysis(sissor_vcf_file, bed_file, hapcut_block_file, filtered_hapcut_block_file, intermediate_fragment_file, output_file):

    filter_hapblock_file(hapcut_block_file, filtered_hapcut_block_file)
    create_sissor_fragments(sissor_vcf_file, bed_file, filtered_hapcut_block_file, intermediate_fragment_file)
    fragment_error_rate(intermediate_fragment_file, filtered_hapcut_block_file, output_file)

# parse args and execute main function
if __name__ == '__main__':
    args = parse_args()
    sissor_error_analysis(args.sissor_vcf_file, args.bed_file, args.hapcut_block_file, args.filtered_hapcut_block_file, args.intermediate_fragment_file, args.output_file)
