import sys
pileup_file = sys.argv[1]
variant_vcf_file = sys.argv[2]

snp_indices_set      = set()
# need to define the columns of the fragment matrix (indexes and positions of SNPs)
# use a vcf file as a template to know where SNPs and their indices are
with open(variant_vcf_file, 'r') as infile:

    for line in infile:
        if (line[0] == '#' or len(line) < 2):
            continue
        el = line.strip().split()
        chrom = el[0]
        pos   = int(el[1])-1
        snp_indices_set.add((chrom,pos))

with open(pileup_file, 'r') as infile:

    for line in infile:

        # read line elements
        el    = line.strip().split()
        chrom = el[0]
        pos   = int(el[1])-1
        qual  = float(el[4])
        depth = int(el[7])

        Q = 10**(qual*-0.1)

        # filter on quality and depth
        # we only care about this position if it's also in the HapCUT block file
        if qual >= 30 and depth >= 5 and (chrom, pos) in snp_indices_set:
            print("{}\t{}\t{}".format(chrom,pos,Q))
