#! /bin/python3.5
# -*- coding: utf-8 -*-

# Author : Peter Edge
# Email  : pedge@eng.ucsd.edu
import argparse
import sys
from create_sissor_fragments import create_sissor_fragments
from fragment_error_rate import fragment_error_rate

desc = '''creates a fragment matrix from SISSOR and then funnels this into fragment_error_rate'''

# parse arguments
def parse_args():

    parser = argparse.ArgumentParser(description='')
    parser.add_argument('sissor_vcf_file', nargs='?', type = str, help='path to haploid vcf file. present alleles are variant calls. missing alleles are assumed to be ref')
    parser.add_argument('bed_file', nargs='?', type = str, help='bed file where each entry specifies the boundary of a SISSOR fragment')
    parser.add_argument('hapcut_block_file', nargs='?', type = str, help='path to hapcut block file')
    parser.add_argument('output_file', nargs='?', type = str, help='path to write error analysis output')
    parser.add_argument('-p','--pileup_file', nargs='?', type = str, help='oldschool samtools pileup file from SISSOR raw reads for determining positions with coverage >= 5 and quality >= 30',default=None)
    parser.add_argument('-c', '--covered_file', nargs='?', type = str, help='small file with intersection of HapCUT positions and pileup positions meeting threshold',default=None)

    # default to help option. credit to unutbu: http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
    if len(sys.argv) < 5:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    return args

def sissor_error_analysis(sissor_vcf_file, bed_file, hapcut_block_file, output_file, pileup_file, covered_file):

    frag_list = create_sissor_fragments(sissor_vcf_file, bed_file, hapcut_block_file, pileup_file, covered_file)
    fragment_error_rate(frag_list, hapcut_block_file, output_file)

# parse args and execute main function
if __name__ == '__main__':
    args = parse_args()
    if args.pileup_file == None and args.covered_file == None:
        print('Must specify either pileup file (-p) or covered file (-c)')
        exit(0)

    sissor_error_analysis(args.sissor_vcf_file, args.bed_file, args.hapcut_block_file, args.output_file, args.pileup_file, args.covered_file)
