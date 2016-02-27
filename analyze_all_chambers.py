# script to analyze all chambers

from sissor_error_analysis import sissor_error_analysis
import os

hapblock_file = '/home/pedge/git/SISSOR/data/haplotype.all.bac.txt'
if not os.path.exists('sissor_error_output'):
    os.mkdir('sissor_error_output')

for c in range(2,25):
    c_fill = str(c).zfill(2)
    vcf = '/oasis/tscc/scratch/wkchu/SISSOR/PGP1_21_highoutputmem/vcf/PGP1_21_ch{}.sorted.fragment.depth5.vcf'.format(c_fill)
    bed = '/oasis/tscc/scratch/wkchu/SISSOR/PGP1_21_highoutputmem/fragmentboundarybed/PGP1_21_FragmentBoundaryCh{}.bed'.format(c)
    output = 'sissor_error_output/chamber{}_error.txt'.format(c)
    pfile = '/oasis/tscc/scratch/wkchu/SISSOR/PGP1_21_highoutputmem/pileup/PGP1_21_ch{}.sorted.fragment.quality.pileup'.format(c_fill)

    sissor_error_analysis(vcf, bed, hapblock_file, output, pileup_file=pfile)
