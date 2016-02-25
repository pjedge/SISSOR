
import sys
import argparse

# parse arguments
def parse_args():

    parser = argparse.ArgumentParser(description='')
    parser.add_argument('hapcut_block_file', nargs='?', type = str, help='path to hapcut block file')
    parser.add_argument('pileup_file', nargs='?', type = str, help='oldschool samtools pileup file from SISSOR raw reads for determining positions with coverage >= 5 and quality >= 30')
    parser.add_argument('output_file', nargs='?', type = str, help='path to output file')
    # default to help option. credit to unutbu: http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
    if len(sys.argv) < 4:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    return args

def get_covered_positions(hapcut_block_file, pileup_file, output_file):

    # want a set too for fast lookup
    snp_indices_set      = {}

    with open(hapcut_block_file, 'r') as infile:

        for line in infile:
            if ('BLOCK' in line or '********' in line or len(line) < 2):
                continue
            el = line.strip().split()
            chrom = el[3]
            pos   = int(el[4])-1
            snp_indices_set.add((chrom,pos))

    nonrefs = set()
    with open(pileup_file, 'r') as infile:

        for line in infile:

            # read line elements
            el    = line.strip().split()
            chrom = el[0]
            pos   = int(el[1])-1
            qual  = int(el[4])
            depth = int(el[7])

            # filter on quality and depth
            # we only care about this position if it's also in the HapCUT block file
            if qual >= 30 and depth >= 5 and (chrom, pos) in snp_indices_set:
                nonrefs.add((chrom,pos))

    with open(output_file,'w') as opf:

        for chrom, pos in nonrefs:

            print("{}\t{}".format(chrom, pos),file=opf)


if __name__ == '__main__':
    args = parse_args()
    get_covered_positions(args.hapcut_block_file, args.pileup_file, args.output_file)
