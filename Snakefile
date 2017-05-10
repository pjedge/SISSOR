from collections import defaultdict
import pickle
import random
import estimate_parameters
import fix_sam
import chamber_allele_call_accuracy as caca
import base_calls_to_vcf
import fragment_haplotype_assignment
from itertools import product
#sys.path.append('/home/pedge/git/HapTools')
from file_processing import prune_hapblock_file, filter_wgs_vcf, split_vcf
from plot_sissor import plot_sissor
from create_hapcut_fragment_matrix import create_hapcut_fragment_matrix
from fix_chamber_contamination import fix_chamber_contamination
import generate_tables
#import run_tools
import fragment
import fileIO
import os

import filter_vcf
import calculate_haplotype_statistics as chs
localrules: all, simlinks, pileup_test, make_accuracy_table


BCFTOOLS = '/path/to/bcftools' #'/opt/biotools/bcftools/bin/bcftools'
SAMTOOLS = '/path/to/samtools' #'/opt/biotools/samtools/1.3/bin/samtools'
PICARD   = '/path/to/picard.jar' #'/home/pedge/installed/picard.jar'
CROSSMAP = '/path/to/CrossMap.py' #'/home/pedge/installed/opt/python/bin/CrossMap.py'
HAPCUT2  = '/path/to/HAPCUT2' #'/home/pedge/git/hapcut2/build/HAPCUT2'
PYPY     = '/path/to/pypy3.3' #'/home/pedge/installed/pypy3.3-5.5-alpha-20161013-linux_x86_64-portable/bin/pypy3.3'
WGS_VCF_URL = 'https://www.encodeproject.org/files/ENCFF995BBX/@@download/ENCFF995BBX.vcf.gz'
HG19     = '/path/to/hg19.fa' #'/oasis/tscc/scratch/pedge/data/genomes/hg19/hg19.fa' #'/home/wkchu/zhang_lab_oasis/resources_v2.8_b37/human_g1k_v37_decoy.fasta'
CGI_SNPs1 = '/path/to/CGI_reference.gff' # Complete Genomics reference dataset for PGP1 in GFF format #'/oasis/tscc/scratch/wkchu/SISSOR/PGP1_A1/BAM/ns.gff'
WGS_BAM_URL = 'https://www.encodeproject.org/files/ENCFF713HUF/@@download/ENCFF713HUF.bam' # URL for 60x WGS of PGP1f cells from ENCODE
GRCH38   = '/path/to/grch38.fa' #'/oasis/tscc/scratch/pedge/data/genomes/grch38/grch38.fa'
GRCH38_URL = 'https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz'

# a collection of BAC Libraries for PGP1, each in VCF format.
BAC_VCFs = ['/path/to/BAC_library1.vcf','/path/to/BAC_library2.vcf'] #['/home/wkchu/BACPoolPileup/Indx73.bac.pileup.vcf',
#'/home/wkchu/BACPoolPileup/Indx74.bac.pileup.vcf',
#'/home/wkchu/BACPoolPileup/Indx75.bac.pileup.vcf',
#'/home/wkchu/BACPoolPileup/Indx76.bac.pileup.vcf',
#'/home/wkchu/BACPoolPileup/Indx77.bac.pileup.vcf',
#'/home/wkchu/BACPoolPileup/Indx78.bac.pileup.vcf',
#'/home/wkchu/BACPoolPileup/Indx79.bac.pileup.vcf',
#'/home/wkchu/BACPoolPileup/Indx80.bac.pileup.vcf',
#'/home/wkchu/BACPoolPileup/Indx81.bac.pileup.vcf',
#'/home/wkchu/BACPoolPileup/Indx82.bac.pileup.vcf',
#'/home/wkchu/BACPoolPileup/Indx83.bac.pileup.vcf',
#'/home/wkchu/BACPoolPileup/Indx84.bac.pileup.vcf',
#'/home/wkchu/BACPoolPileup/Indx85.2.bac.pileup.vcf',
#'/home/wkchu/BACPoolPileup/Indx85.bac.pileup.vcf',
#'/home/wkchu/BACPoolPileup/Indx86.bac.pileup.vcf',
#'/home/wkchu/BACPoolPileup/Indx87.bac.pileup.vcf',
#'/home/wkchu/BACPoolPileup/Indx88.bac.pileup.vcf',
#'/home/wkchu/BACPoolPileup/Indx89.bac.pileup.vcf',
#'/home/wkchu/BACPoolPileup/Indx90.bac.pileup.vcf',
#'/home/wkchu/BACPoolPileup/Indx91.bac.pileup.vcf',
#'/home/wkchu/BACPoolPileup/Indx92.2.bac.pileup.vcf',
#'/home/wkchu/BACPoolPileup/Indx92.bac.pileup.vcf',
#'/home/wkchu/BACPoolPileup/Indx93.bac.pileup.vcf',
#'/home/wkchu/BACPoolPileup/Indx94.bac.pileup.vcf',
#'/home/wkchu/BACPoolPileup/Indx95.bac.pileup.vcf',
#'/home/wkchu/BACPoolPileup/Indx96.bac.pileup.vcf']

# path to directories with the 3 SISSOR Library BAM files.
#PGP1_21_dir = '/path/to/SISSOR_library1_dir' #'/oasis/tscc/scratch/wkchu/SISSOR/PGP1_21_highoutputmem/BAM'
#PGP1_22_dir = '/path/to/SISSOR_library2_dir' #'/oasis/tscc/scratch/wkchu/SISSOR/PGP1_22/2016OctMergedBAM'
#PGP1_A1_dir = '/path/to/SISSOR_library3_dir' #'/oasis/tscc/scratch/wkchu/SISSOR/PGP1_A1/2016OctMergedBAM'

# cutoff values are phred scaled probability of error
# so 30 => 0.999 accuracy
# and 50 => 0.99999 accuracy
variant_calling_cutoffs = [7,8,9,10,30,50,70,90,110,130,150]

chroms = ['chr{}'.format(i) for i in range(1,23)]
chroms_XY = chroms + ['chrX','chrY']
chroms_set_XY = set(chroms_XY)

chunksize = int(5e6)

# break the genome into 5 Mb chunks for parallel base calling with SISSORhands and Freebayes

grch38_size_list = [('chr1', 248956422),
 ('chr2', 242193529),
 ('chr3', 198295559),
 ('chr4', 190214555),
 ('chr5', 181538259),
 ('chr6', 170805979),
 ('chr7', 159345973),
 ('chr8', 145138636),
 ('chr9', 138394717),
 ('chr10', 133797422),
 ('chr11', 135086622),
 ('chr12', 133275309),
 ('chr13', 114364328),
 ('chr14', 107043718),
 ('chr15', 101991189),
 ('chr16', 90338345),
 ('chr17', 83257441),
 ('chr18', 80373285),
 ('chr19', 58617616),
 ('chr20', 64444167),
 ('chr21', 46709983),
 ('chr22', 50818468),
 ('chrX', 156040895),
 ('chrY', 57227415)]

hg19_size_list = [('chr1', 249250621),
 ('chr2', 243199373),
 ('chr3', 198022430),
 ('chr4', 191154276),
 ('chr5', 180915260),
 ('chr6', 171115067),
 ('chr7', 159138663),
 ('chr8', 146364022),
 ('chr9', 141213431),
 ('chr10', 135534747),
 ('chr11', 135006516),
 ('chr12', 133851895),
 ('chr13', 115169878),
 ('chr14', 107349540),
 ('chr15', 102531392),
 ('chr16', 90354753),
 ('chr17', 81195210),
 ('chr18', 78077248),
 ('chr19', 59128983),
 ('chr20', 63025520),
 ('chr21', 48129895),
 ('chr22', 51304566),
 ('chrX', 155270560),
 ('chrY', 59373566)]


# create chunks of grch38
# return list of (chrom,start,stop) tuples. stop is inclusive
grch38_chunklist = []
for chrom, chrlen in grch38_size_list:
    for start in range(1,chrlen+1,chunksize):
        end = start+chunksize-1 if start+chunksize-1 < chrlen else chrlen
        grch38_chunklist.append((chrom,start,end))

grch38_regions = ['{}.{}.{}'.format(chrom,start,stop) for chrom,start,stop in grch38_chunklist]

grch38_regions_noXY = ['{}.{}.{}'.format(chrom,start,stop) for chrom,start,stop in grch38_chunklist if chrom not in ['chrX','chrY']]
grch38_regions_XYonly = ['{}.{}.{}'.format(chrom,start,stop) for chrom,start,stop in grch38_chunklist if chrom in ['chrX','chrY']]
grch38_chunklist_noXY = [(chrom,start,stop) for chrom,start,stop in grch38_chunklist if chrom not in ['chrX','chrY']]
grch38_chunklist_XYonly = [(chrom,start,stop) for chrom,start,stop in grch38_chunklist if chrom in ['chrX','chrY']]


# create chunks of hg19
# return list of (chrom,start,stop) tuples. stop is inclusive
chunklist = []

for chrom, chrlen in hg19_size_list:
    for start in range(1,chrlen+1,chunksize):
        end = start+chunksize-1 if start+chunksize-1 < chrlen else chrlen
        chunklist.append((chrom,start,end))

regions = ['{}.{}.{}'.format(chrom,start,stop) for chrom,start,stop in chunklist]

#regions = regions[5:8]

cells = ['PGP1_21','PGP1_22','PGP1_A1']
chambers = list(range(1,25))

cells = ['PGP1_21','PGP1_22','PGP1_A1']
chambers = list(range(1,25))
chambers_pad = ['{0:02d}'.format(c) for c in chambers]

HAPLOTYPE_ERROR_RATES_DIR = 'haplotyping/error_rates'
HAPLOTYPE_PLOTS_DIR       = 'haplotyping/plots'
HAPLOTYPE_EXP_DIR         = 'haplotyping/experiments'
haplotyping_VCF_dir      = 'haplotyping/data/haplotyping_VCFs'
BAC_VCF_dir    = 'haplotyping/data/PGP1_VCFs_BACindex'

modes = ['same_cell','all_cell','cross_cell','ind_same_cell']

# generate all tables and results.
rule all:
    input:
        'accuracy_reports/tables/unphased_call_accuracy_master_table.txt',
        'accuracy_reports/tables/phased_call_accuracy_master_table.txt',
        'accuracy_reports/tables/strand_strand_mismatch_rates_table.txt',
        'accuracy_reports/tables/strand_strand_nucleotide_substitutions.png',
        "{}/sissor_haplotype_error_genome.png".format(HAPLOTYPE_PLOTS_DIR),

# make a table showing the accuracy of same-haplotype strand-matched calls from SISSOR.
rule make_phased_accuracy_master_table:
    params: job_name = 'make_phased_accuracy_master_table'
    input:  same_cell = 'accuracy_reports/same_cell/cutoff10.counts.p',
            all_cell = 'accuracy_reports/all_cell/cutoff10.counts.p',
            cross_cell = 'accuracy_reports/cross_cell/cutoff10.counts.p',
            ind_same_cell = 'accuracy_reports/ind_same_cell/cutoff10.counts.p',
    output: table = 'accuracy_reports/tables/phased_call_accuracy_master_table.txt',
    run:
        generate_tables.generate_phased_calls_table(output.table, input.same_cell, input.all_cell, input.cross_cell, input.ind_same_cell)

# make a table showing the accuracy of the SISSORhands consensus variant caller
rule make_unphased_accuracy_master_table:
    params: job_name = 'make_unphased_accuracy_master_table'
    input:  counts_files = expand('accuracy_reports/unphased/cutoff{c}.counts.p',c=variant_calling_cutoffs)
    output: table = 'accuracy_reports/tables/unphased_call_accuracy_master_table.txt'
    run:
        generate_tables.generate_unphased_calls_table(output.table, input.counts_files, variant_calling_cutoffs)

# make a table showing rates of strand-strand mismatches, presumably due to MDA
rule make_strand_strand_mismatch_table:
    params: job_name = 'make_strand_strand_mismatch_table'
    input:  same_cell = 'accuracy_reports/strand_mismatch_same_cell/cutoff10.counts.p',
            all_cell = 'accuracy_reports/strand_mismatch_all_cell/cutoff10.counts.p',
            cross_cell = 'accuracy_reports/strand_mismatch_cross_cell/cutoff10.counts.p'
    output: table = 'accuracy_reports/tables/strand_strand_mismatch_rates_table.txt',
    run:
        generate_tables.generate_strand_strand_mismatch_table(output.table, input.same_cell, input.all_cell, input.cross_cell)

# make a plot showing the rate of nucleotide substitutions between SISSOR strands
rule generate_nucleotide_substitution_plot:
    params: job_name = 'generate_nucleotide_substitution_plot'
    input:  same_cell = 'accuracy_reports/strand_mismatch_same_cell/cutoff10.counts.p',
    output: plot = 'accuracy_reports/tables/strand_strand_nucleotide_substitutions.png',
    run:
        generate_tables.generate_nucleotide_substitution_plot(output.plot, input.same_cell)

# generate a table including everything in the raw accuracy count dictionaries
rule make_accuracy_table:
    params: job_name = 'make_accuracy_table.{report_name}.{cut}'
    input:  counts = 'accuracy_reports/{report_name}/cutoff{cut}.counts.p',
    output: table = 'accuracy_reports/tables/{report_name}.{cut}.table.txt',
    run:
        caca.generate_table(input.counts, output.table)

# aggregate the accuracy counts for phased or unphased SISSOR variant calls.
rule accuracy_aggregate_counts:
    params: job_name = 'accuracy_aggregate_counts.{report_name}.{cut}'
    input: counts = expand('accuracy_reports/{{report_name}}/split/{r}/{{cut}}/counts.p',r=regions),
    output: counts = 'accuracy_reports/{report_name}/cutoff{cut}.counts.p',
    run:
        caca.accuracy_aggregate(input.counts,output.counts)

# aggregate the separate files that print out the positions with mismatches
# this is only for inspecting individual SNVs that differ from PGP1 reference,
# or that have strand-strand MDA-related mismatches, etc.
rule accuracy_aggregate_mismatches:
    params: job_name = 'accuracy_aggregate_mismatches.{cut}'
    input:  mof = expand('accuracy_reports/{{report_name}}/split/{r}/{{cut}}/mismatches',r=regions),
    output: mof = 'accuracy_reports/{report_name}/cutoff{cut}.mismatches'
    shell:  'cat {input.mof} > {output.mof}'

# generate a dictionary with SNV/reference calls for every 5mb chunk of the genome.
# we intersect a CGI dataset for PGP1 with a 60x Illumina WGS dataset,
# to get highly accurate calls for 2.7 Gb Reference positions and 3Mb SNV positions
# that are shared between PGP1 fibroblast and lymphocyte cell lines.

# we then extract SNVs from a set of PGP1 BAC libraries to use as secondary validation for
# SNVs not seen in the other datasets.
rule generate_ref_dict:
    params: job_name = 'generate_ref_dict.{chrom}.{start}.{end}'
    input:  gff = CGI_SNPs1,
            hg19 = HG19,
            wgs_done = 'wgs/freebayes_hg19/{chrom}.done', #wgs_vcf = 'wgs/freebayes_hg19/{chrom}.{start}.{end}.vcf',
            bac_vcfs = BAC_VCFs,
    output: pickle = 'ref_dicts/{chrom}.{start}.{end}.ref_dict.p'
    run:
        reg_chrom = str(wildcards.chrom)
        reg_start = int(wildcards.start)
        reg_end   = int(wildcards.end)
        region = (reg_chrom, reg_start, reg_end)
        wgs_vcf = 'wgs/freebayes_hg19/{}.{}.{}.vcf'.format(reg_chrom, reg_start, reg_end)
        assert(os.path.isfile(wgs_vcf))
        caca.generate_ref_dict(input.gff, wgs_vcf, input.bac_vcfs, input.hg19, region, output.pickle)

# count numbers for accuracy of unphased, consensus variant calls from SISSORhands.
rule accuracy_count_unphased:
    params: job_name = 'accuracy_count_unphased.{r}'
    input:  ccf = 'base_calls/unphased/{r}.out',
            ref_dict = 'ref_dicts/{r}.ref_dict.p',
    output: counts_lst = expand('accuracy_reports/unphased/split/{{r}}/{cut}/counts.p',cut=variant_calling_cutoffs),
            mof_lst = expand('accuracy_reports/unphased/split/{{r}}/{cut}/mismatches',cut=variant_calling_cutoffs),
    run:
        for counts, mof, cut in zip(output.counts_lst, output.mof_lst, variant_calling_cutoffs):
            caca.accuracy_count_unphased(input.ccf, input.ref_dict, (1.0-10**(-0.1*float(cut))), counts, mof)

# count numbers for the accuracy of SISSOR using same-haplotype strand-matching technique
# there are 3 primary ways:
# same cell: get variant calls unique to a single cell
# all cell: use calls from pairs within OR between cells, to maximize coverage/number of calls
# cross-cell: compare between strands in different cells only, to gauge error rate of strand pairing approach without it being inflated by cell-specific mutations
rule accuracy_count_phased:
    params: job_name = 'accuracy_count_{mode,(same_cell|all_cell|cross_cell|ind_same_cell)}.{r}'
    input:  ccf = 'base_calls/phased/{r}.ccf',
            ref_dict = 'ref_dicts/{r}.ref_dict.p',
    output: counts = 'accuracy_reports/{mode,(same_cell|all_cell|cross_cell|ind_same_cell)}/split/{r}/10/counts.p',
            mof = 'accuracy_reports/{mode,(same_cell|all_cell|cross_cell|ind_same_cell)}/split/{r}/10/mismatches',
            smf = 'accuracy_reports/strand_mismatch_{mode,(same_cell|all_cell|cross_cell|ind_same_cell)}/split/{r}/10/counts.p'
    run:
        CUT = 0.9
        caca.accuracy_count_phased(input.ccf, input.ref_dict, CUT, output.counts, output.mof, output.smf,separate_haplotypes=True,mode=wildcards.mode)

###################################################################

# OBTAIN SET OF REFERENCE AND VARIANT CALLS FROM PGP1F DATASET
# procedure is to call freebayes on short chunks of bam,
# combine each individual chrom,
# sort the chromosome separately,
# then split into short regions again.

# re-split the lifted over and sorted VCFs into 5Mb chunks corresponding to hg19.
rule split_freebayes_wgs:
    params: job_name  = 'split_freebayes_wgs_{c}',
    input:  vcf  = 'wgs/lifted_sorted/{c}.vcf'
    output: done  = touch('wgs/freebayes_hg19/{c}.done')
    run:
        outfiles   = ['wgs/freebayes_hg19/{}.{}.{}.vcf'.format(chrom,start,end) for (chrom,start,end) in chunklist if chrom == wildcards.c]
        chr_chunks = [(chrom,start,end) for (chrom,start,end) in chunklist if chrom == wildcards.c]
        split_vcf(input.vcf, chr_chunks, outfiles)

# combine, and sort, the lifted over VCFs for PGP1f WGS SNV/ref calls for every genomic position
rule sort_freebayes_wgs:
    params: job_name  = 'sort_freebayes_wgs.{c}',
    input:  expand('wgs/lifted/{r}.vcf',r=grch38_regions),
    output: 'wgs/lifted_sorted/{c}.vcf'
    run:
        infiles =  [f for (chrom,start,end),f in zip(grch38_chunklist,input) if chrom == wildcards.c]
        shell('''
        grep '^#' {infiles[0]} > {output}
        cat {infiles} |
        grep -v '^#' |
        vcf-sort -c -p 4 >> {output}
        ''')

# liftover individual VCF files to hg19
rule liftover_wgs:
    params: job_name  = 'liftover.{r}',
    input:  vcf = 'wgs/freebayes/{r}.vcf',
            chain = 'hg38ToHg19.over.chain'
    output: vcf = 'wgs/lifted/{r}.vcf',
    shell:
        '''
        {CROSSMAP} \
        vcf \
        {input.chain} \
        {input.vcf} \
        {HG19} \
        {output.vcf}
        '''

# run freebayes to call reference and variants for every genomic position for PGP1f WGS bam file.
# as with everything else, this is done in 5Mb chunks to allow cluster execution
rule run_freebayes_wgs:
    params: job_name  = 'freebayes.{chr}.{start}.{stop}',
    input:  bam   = 'wgs/wgs.bam',
            bai   = 'wgs/wgs.bam.bai',
            ref   = GRCH38
    output: region_vcf = 'wgs/freebayes/{chr}.{start}.{stop}.vcf'
    shell:
        '''
        freebayes -f {GRCH38} \
        --standard-filters \
        --report-monomorphic \
        --region {wildcards.chr}:{wildcards.start}..{wildcards.stop} \
         {input.bam} > {output.region_vcf}
        '''

# download grch38 genome
rule download_grch38:
    params: job_name = 'download_grch38'
    output: fa = GRCH38
    shell:
        '''
        wget {GRCH38_URL} -O {output.fa}.gz
        gunzip {output.fa}.gz
        '''

# download raw bam of 60x Illumina WGS on PGP1f cells
rule download_WGS_bam:
    params: job_name = 'download_WGS_bam'
    output: bam = 'wgs/wgs.bam'
    shell:
        '''
        wget {WGS_BAM_URL} -O {output.bam}
        '''
###################################################################

# add phasing information (P1 and P2 info tags) to unphased chamber call files
# this will allow us to traverse the files and pair up calls in different chambers
# to perform same-haplotype strand-matching.
rule annotate_assigned_fragments:
    params: job_name = 'annotate_assigned_fragments.{r}'
    input:  ccf = 'base_calls/unphased/{r}.out',
            asn = 'fragment_haplotype_assignments/all',
    output: ccf = 'base_calls/phased/{r}.ccf',
    run:
        fragment_haplotype_assignment.annotate_assigned_fragments(input.ccf,input.asn,output.ccf)

# combine separate files for fragment-to-haplotype assignments
rule combine_assigned_fragment_files:
    params: job_name = 'combine_assigned_fragment_files'
    input:  sep = expand('fragment_haplotype_assignments/{chrom}',chrom=chroms),  #+['chrXY']
    output: combined = 'fragment_haplotype_assignments/all',
    shell:  'cat {input.sep} > {output.combined}'

# assign every SISSOR haplotype fragment back to the assembled haplotype it matches best
rule assign_fragment_haplotypes:
    params: job_name = 'assign_fragment_haplotypes.chr{chrom,\d+}'
    input: frag = 'haplotyping/data/PGP1_ALL/fragmat/cov1_strict/chr{chrom,\d+}',
           vcf  = 'haplotyping/data/haplotyping_VCFs/chr{chrom,\d+}.vcf',
           hap  = 'haplotyping/experiments/hapcut2_PGP1_ALL/cov1_strict/chr{chrom,\d+}.output'
    output: sep = 'fragment_haplotype_assignments/chr{chrom,\d+}',
    run:
        fragment_haplotype_assignment.assign_fragment_haplotypes(input.frag,input.vcf,output.sep,input.hap)

# make a bar chart showing the accuracy (switch, mismatch) and completeness (N50, AN50) for SISSOR haplotypes
rule plot_hapcut2_results:
    params:
        job_name = "plot_hapcut2_sissor"
    input:
        stats_file  = "{}/hapcut2.stats.p".format(HAPLOTYPE_ERROR_RATES_DIR),
        labels_file = "{}/hapcut2.labels.p".format(HAPLOTYPE_ERROR_RATES_DIR)
    output:
        plot2 = "%s/sissor_haplotype_error_genome.png" % HAPLOTYPE_PLOTS_DIR
    run:
        data = pickle.load(open(input.stats_file,"rb"))
        labels = pickle.load(open(input.labels_file,"rb"))
        plot_sissor(data,labels,output.plot2)

# calculate switch and mismatch error rates of the assemble SISSOR haplotypes
# also calculates N50 and AN50, phase rate, etc
exp =['cov1_none','cov1_strict']
exp_labels=['Unprocessed Fragments','Processed Fragments']
rule calculate_error_rates:
    params:
        job_name = "hapcut2_error_rates"
    input:
        hapblocks = expand("{E}/hapcut2_PGP1_ALL/{exp}/{c}.output",E=HAPLOTYPE_EXP_DIR,exp=exp,c=chroms),
        fragmats  = expand("haplotyping/data/PGP1_ALL/fragmat/{exp}/{c}",exp=exp,c=chroms),
        var_vcfs  = expand("{v}/{c}.vcf",v=haplotyping_VCF_dir, c=chroms),
        bac_haps  = expand("haplotyping/data/BAC/{c}.filtered",c=chroms),
        contig_size_file = 'hg19.chrom.sizes'
    output:
        stats_file  = "{}/hapcut2.stats.p".format(HAPLOTYPE_ERROR_RATES_DIR),
        labels_file = "{}/hapcut2.labels.p".format(HAPLOTYPE_ERROR_RATES_DIR)
    run:
        # list of lists of error results,
        data   = []
        labels = exp_labels

        for x in exp: # sample

            assembly_files = ["{}/hapcut2_PGP1_ALL/{}/{}.output".format(HAPLOTYPE_EXP_DIR,x,c) for c in chroms]
            frag_files     = ["haplotyping/data/PGP1_ALL/fragmat/{}/{}".format(x,c) for c in chroms]
            vcf_files      = ["{}/{}.vcf".format(haplotyping_VCF_dir,c) for c in chroms]
            truth_files    = ["haplotyping/data/BAC/{}.filtered".format(c) for c in chroms]
            truth_vcf_files      = ["{}/{}.vcf".format(BAC_VCF_dir,c) for c in chroms]

            err = chs.hapblock_hapblock_error_rate_multiple(truth_files, truth_vcf_files, assembly_files, frag_files, vcf_files, input.contig_size_file)

            print("{} results over all chromosomes:".format(x))
            print(err)

            data.append(err)

        pickle.dump(data,open(output.stats_file,"wb"))
        pickle.dump(labels,open(output.labels_file,"wb"))

# Prune SNVs from SISSOR haplotypes and split haplotype blocks at possible switch errors, at moderate stringency
rule prune_haplotype:
    params:
        job_name = "{s}.{x}.prune_haplotype",
    input:
        hapblocks = expand("{E}/hapcut2_{{s}}/{{x}}/{c}.output.uncorrected",E=HAPLOTYPE_EXP_DIR,c=chroms)
    output:
        hapblocks = expand("{E}/hapcut2_{{s}}/{{x}}/{c}.output",E=HAPLOTYPE_EXP_DIR,c=chroms)
    run:
        for i, o in zip(input.hapblocks,output.hapblocks):
            prune_hapblock_file(i, o, snp_conf_cutoff=13.01, split_conf_cutoff=13.01, use_refhap_heuristic=True) #split_conf_cutoff=0.9999

# Run HapCUT2 to assemble haplotypes from SISSOR fragments (raw fragments and post-processed fragments to compare completeness and accuracy)
rule run_hapcut2:
    params:
        job_name = "{s}.{c}.{x}.hapcut",
    input:
        frag_file = "haplotyping/data/{s}/fragmat/{x}/{c}",
        vcf_file  = lambda wildcards: expand("{v}/{c}.vcf", v=haplotyping_VCF_dir,c=wildcards.c)
    output:
        hapblocks = "{E}/hapcut2_{s}/{x}/{c}.output.uncorrected",
    shell:
        '''
        {HAPCUT2} --fragments {input.frag_file} --vcf {input.vcf_file} --output {output.hapblocks} --ea 1
        '''

# we now have fragment files that are nearly in HapCUT2 format but
# with one alteration -- they've been augmented with 'M' symbols for "mixed calls"
# this rule applies the following steps:
# filter out fragments with many mixed calls
# split fragments at clusters of mixed calls
# remove fragments that are highly discordant with the phase of other fragments
# split fragments at detectable switch errors to other fragments (2 or more SNVs switched with respect to another fragment)
# if a fragment has the same switch to 2 or more fragments then only that offending fragment is split
# if it is ambiguous (2 fragments, switched only w.r.t. each other, both are split)
rule fix_fragmat:
    params:
        job_name  = "fix_fragmat.{x}",
    input:
        var_vcfs = expand("{v}/{CHR}.vcf",v=haplotyping_VCF_dir,CHR=chroms),
        P_ALL = expand("haplotyping/data/PGP1_ALL/augmented_fragmat/{c}",c=chroms)
    output:
        fixed = expand("haplotyping/data/PGP1_ALL/fragmat/{{x}}/{c}",c=chroms)
    run:
        if 'none' in wildcards.x:
            mode = 'none'
        elif 'basic' in wildcards.x:
            mode = 'basic'
        elif 'strict' in wildcards.x:
            mode = 'strict'

        for i,v,o in zip(input.P_ALL, input.var_vcfs, output.fixed):
            fix_chamber_contamination(i,v,o,threshold=2, min_coverage=0,mode=mode)

# generate a HapCUT2 format fragment file from our SISSORhands "chamber call file" variant calls.
rule generate_fragmatrix:
    params: job_name = 'generate_fragmatrix'
    input:  ccf      = 'haplotyping/data/PGP1_ALL.het.whole_genome.ccf',
            vcfs     = expand('{v}/{c}.vcf',v=haplotyping_VCF_dir,c=chroms),
            bounds   = expand('eric_fragment_boundary_beds/{P[0]}/ch{P[1]}.bed',P=product(cells,chambers)),
    output: fragmat  = expand('haplotyping/data/PGP1_ALL/augmented_fragmat/{c}',c=chroms),

    run:
        odir = 'haplotyping/data/PGP1_ALL/augmented_fragmat'
        # generate fragment matrix from sissorhands base calls
        create_hapcut_fragment_matrix(chamber_call_file=input.ccf, variant_vcf_files=input.vcfs,fragment_boundary_files=input.bounds, output_dir=odir)

# combine the filtered "chamber call files" into one file.
# the calls in this file will be used to generate HapCUT2 fragment matrix files
rule combine_filtered:
    params: job_name = 'combine_filtered'
    input:  ccf = expand('haplotyping/data/ccf_split/{r}.ccf',r=regions),
            vcf = expand('haplotyping/data/vcf_split/{r}.vcf',r=regions)
    output: ccf = 'haplotyping/data/PGP1_ALL.het.whole_genome.ccf',
            vcf = 'haplotyping/data/PGP1_ALL.het.whole_genome.vcf',
    shell:
        '''cat {input.ccf} > {output.ccf}
           cat {input.vcf} > {output.vcf}'''

# we will use per-chamber variant calls from SISSORhands to assemble haplotypes from 3 SISSOR cells.
# we start by taking the output of SISSORhands and filtering them for only positions in our heterozygous SNV set (from PGP1f Illmina WGS).
# we do this because the output of SISSORhands is massive and scattered in many files,
# as it has calls and probability information for every chamber at every genomic position.
rule filter_CCF:
    params: job_name = 'filter_CCF_{chrom}.{start}.{end}'
    input:  ccf      = 'base_calls/unphased/{chrom}.{start}.{end}.out',
            vcf      = 'haplotyping/data/haplotyping_VCFs/all.vcf'
    output: ccf      = 'haplotyping/data/ccf_split/{chrom}.{start}.{end}.ccf',
            vcf      = 'haplotyping/data/vcf_split/{chrom}.{start}.{end}.vcf',
    run:
        het_pos = dict()
        reg_chrom = str(wildcards.chrom)
        reg_start = int(wildcards.start)
        reg_end   = int(wildcards.end)
        with open(input.vcf,'r') as vcf_file:
            for line in vcf_file:
                if len(line) < 3:
                    continue
                el = line.strip().split('\t')
                chrom = el[0]
                pos   = int(el[1])
                if chrom == reg_chrom and pos <= reg_end and pos >= reg_start:
                    het_pos[(chrom,pos)] = line.strip()

        with open(input.ccf,'r') as infile, open(output.ccf,'w') as ccf_out, open(output.vcf,'w') as vcf_out:
            for line in infile:
                if len(line) < 3:
                    continue
                el = line.strip().split('\t')
                chrom = el[0]
                pos   = int(el[1])
                if (chrom,pos) in het_pos:
                    print(line.strip(),file=ccf_out)
                    print(het_pos[(chrom,pos)],file=vcf_out)

# prune the BAC haplotypes and split blocks to remove switch errors,
# at extremely high stringency.
# this will be our reference set to gauge SISSOR haplotyping error against so we
# have a low tolerance for error
rule prune_BAC:
    params:
        job_name = "BAC.prune_haplotype",
    input:
        hapblocks = expand("haplotyping/data/BAC/{c}",c=chroms)
    output:
        hapblocks = expand("haplotyping/data/BAC/{c}.filtered",c=chroms)
    run:
        for i, o in zip(input.hapblocks,output.hapblocks):
            prune_hapblock_file(i, o, snp_conf_cutoff=40, split_conf_cutoff=40, use_refhap_heuristic=True) #split_conf_cutoff=0.9999

# run HapCUT2 to assemble high-confidence reference BAC haplotype.
rule run_hapcut2_BAC:
    params:
        job_name = "{c}.BAC.hapcut2",
    input:
        frag_file = "haplotyping/data/BAC_frags_fixed/{c}",
        vcf_file  = "%s/{c}.vcf" % BAC_VCF_dir
    output:
        hapblocks = "haplotyping/data/BAC/{c}"
    shell:
        '''
        {HAPCUT2} --fragments {input.frag_file} --vcf {input.vcf_file} --output {output.hapblocks} --ea 1
        '''

# remove switch errors from BAC fragments.
# split each fragment at any position with 2 or more SNVs switched with respect to another fragment
# if a fragment has the same switch to 2 or more fragments then only that offending fragment is split
# if it is ambiguous (2 fragments, switched only w.r.t. each other, both are split)
rule fix_fragmat_BAC:
    params:
        job_name  = "fix_fragmat_BAC_{c}",
    input:
        var_vcfs = "%s/{c}.vcf" % BAC_VCF_dir,
        frags = "haplotyping/data/BAC_frags/{c}",
        WGS_variants = "%s/{c}.vcf" % haplotyping_VCF_dir,
    output:
        fixed = "haplotyping/data/BAC_frags_fixed/{c}"
    run:
        fix_chamber_contamination(input.frags,input.var_vcfs,output.fixed,threshold=2, min_coverage=0,mode='basic',vcf_filter=input.WGS_variants)

# create a combined version of the haplotyping VCF with only heterozygous SNVs
rule combine_haplotyping_VCFs:
    params: job_name = "combine_haplotyping_VCFs",
    input: expand("haplotyping/data/haplotyping_VCFs/{c}.vcf",c=chroms)
    output: 'haplotyping/data/haplotyping_VCFs/all.vcf'
    shell:
        '''cat {input} > {output}'''

# filter the haplotyping VCFs for quality and coverage
# keep only heterozygous SNVs because that's all we care about
# return separate VCFs per chromosome
rule filter_haplotyping_VCFs:
    params: job_name  = "filter_haplotyping_VCF",
    input: vcf = "wgs/wgs_hg19.vcf"
    output: vcfs = expand("haplotyping/data/haplotyping_VCFs/{c}.vcf",c=chroms)
    run:
        filter_wgs_vcf(input.vcf,output.vcfs,chroms) # input.filter_set

# sort the variants that we'll use for phasing.
rule sort_wgs_SNPs:
    params: job_name  = 'liftover_SNPs',
    input:  vcf = 'wgs/wgs_hg19.autosomesXY.unsorted.vcf',
    output: vcf = 'wgs/wgs_hg19.vcf'
    shell:
        '''cat {input} | vcf-sort -c -p 4 > {output}'''

# remove all chromosomes except chr1-chr22,chrX,chrY.
# this is just to make sorting work properly.
rule filter_wgs_autosomes:
    params: job_name  = 'filter_wgs_autosomes',
    input:  vcf = 'wgs/wgs_hg19.unsorted.vcf',
    output: vcf = 'wgs/wgs_hg19.autosomesXY.unsorted.vcf'
    run:
        with open(input.vcf,'r') as inf, open(output.vcf,'w') as outf:
            for line in inf:
                if line[:1] == '#':
                    continue
                el = line.strip().split('\t')
                if len(el) < 5:
                    continue
                if el[0] in chroms_set_XY:
                    print(line.strip(),file=outf)

# our data is aligned to hg19.
# carry over the 60x PGP1 fibroblast varants to hg19 so we can use them with
# HapCUT2 to phase.
rule liftover_wgs_SNPs:
    params: job_name  = 'liftover_SNPs',
    input:  vcf = 'wgs/wgs_grch38.vcf',
            chain = 'hg38ToHg19.over.chain',
            hg19 = HG19
    output: vcf = 'wgs/wgs_hg19.unsorted.vcf',
    shell:
        '''
        {CROSSMAP} \
        vcf \
        {input.chain} \
        {input.vcf} \
        {input.hg19} \
        {output.vcf}
        '''

# download VCF file with variants from a 60x Illumina WGS experiment on PGP1 fibroblast cells
# these variants will be used to assemble haplotypes with HapCUT2
rule download_WGS_vcf:
    params: job_name = 'download_WGS_vcf'
    output: vcf = 'wgs/wgs_grch38.vcf'
    shell:
        '''
        wget {WGS_VCF_URL} -O {output.vcf}.gz
        mv {output.vcf}.gz {output.vcf}
        '''

# run the SISSORhands custom variant caller for SISSOR technology to call variants using
# cross-chamber information.
rule call_alleles:
    params: job_name = 'call_alleles.{r}'
    input:  pileup = 'pileups/split/{r}.pileup',
            MDA_dist = 'parameters/MDA_dist.p',
            cov_frac_dist = 'parameters/cov_frac_dist.p',
            p_null = 'parameters/p_null.p',
            ch_priors = 'parameters/ch_priors.p',
            hom_config_probs = 'parameters/hom_config_probs.p',
            het_config_probs = 'parameters/het_config_probs.p',
            haploid_genotype_priors = 'parameters/haploid_genotype_priors.p',
            diploid_genotype_priors = 'parameters/diploid_genotype_priors.p'
    output: ccf = 'base_calls/unphased/{r}.out'
    shell:  '{PYPY} sissorhands.py -i {input.pileup} -o {output.ccf}' # this is the slowest step and requires massive parallel compute time. PYPY makes this a little faster for each job.
    #run:
        #import sissorhands
        #sissorhands.call_chamber_alleles(input.pileup, output.ccf)

# conglomerate statistics from the previous rule into the statistics needed for base calling
rule estimate_parameters:
    params: job_name = 'estimate_parameters'
    input:  expand('parameters/split/chrX_MDA_fracs.{r}.p',r=regions),
            expand('parameters/split/chrX_covs.{r}.p',r=regions),
            expand('parameters/split/chamber_position_counts.{r}.p',r=regions),
            expand('parameters/split/strand_coverage_counts.{r}.p',r=regions),
            #expand('parameters/split/total_sampled_cell_positions.{r}.p',r=regions),
    output: 'parameters/MDA_dist.p',
            'parameters/chrX_covs.p',
            'parameters/chamber_position_counts.p',
            'parameters/strand_coverage_counts.p',
            'parameters/cov_frac_dist.p',
            'parameters/p_null.p',
            'parameters/ch_priors.p',
            'parameters/hom_config_probs.p',
            'parameters/het_config_probs.p',
            'parameters/diploid_genotype_priors.p',
            'parameters/haploid_genotype_priors.p',
    run:
        estimate_parameters.estimate_parameters(regions)

# parse through bam files and count various statistics from the SISSOR libraries
# for the purpose of estimating things like the probability of strands falling in a chamber,
# or the distribution of coverage for chromosome X (used for estimation of strand overlap error)
rule obtain_counts_parallel:
    params: job_name = 'obtain_counts_parallel.{r}'
    input:  pileup = 'pileups/split/{r}.pileup',
    output: 'parameters/split/chrX_MDA_fracs.{r}.p',
            'parameters/split/chrX_covs.{r}.p',
            'parameters/split/chamber_position_counts.{r}.p',
            'parameters/split/strand_coverage_counts.{r}.p',
    run:
        estimate_parameters.obtain_counts_parallel(input.pileup, boundary_files=None, suffix=wildcards.r)

# perform a massive multi-pileup of the reads for all 72 SISSOR libraries bam files (3 cells x 24 chambers)
# pileup in small chunks of 5Mb to allow parallel processing of chunks on the cluster
rule pileup:
    params: job_name = 'pileup.{r}',
            chrom = lambda wildcards: wildcards.r.split('.')[0],
            start = lambda wildcards: wildcards.r.split('.')[1],
            stop  = lambda wildcards: wildcards.r.split('.')[2]
    input:  bams = expand('bams/{P[0]}/ch{P[1]}.bam',P=product(cells, chambers)),
            bais = expand('bams/{P[0]}/ch{P[1]}.bam.bai',P=product(cells, chambers))
    output: pileup = 'pileups/split/{r}.pileup'
    shell:  '{SAMTOOLS} mpileup --region {params.chrom}:{params.start}-{params.stop} --adjust-MQ 50 --max-depth 100 --output-MQ --min-MQ 30 --min-BQ 20 --fasta-ref {HG19} --output {output.pileup} {input.bams}'

# index bamfile for pileup
rule index_bam:
    params: job_name = lambda wildcards: 'index_bam.{}'.format(str(wildcards.x).replace("/", "."))
    input:  bam = '{x}.bam'
    output: bai = '{x}.bam.bai'
    shell:  '{SAMTOOLS} index {input.bam} {output.bai}'

# PGP1_22 and PGP1_A1 bam files don't have chr before chromosome labels. necessary
# to change the chromosome names for mpileup to include reference call.
# credit to petervangalen and Pierre Lindenbaum for this rule's code: https://www.biostars.org/p/13462/
rule fix_sam:
    params: job_name = 'fix_sam.{ce}.{ch}'
    input:  sam = 'bams/{ce}/old.ch{ch}.sam'
    output: bam = 'bams/{ce}/ch{ch}.bam'
    run:
        fix_sam.fix_sam(input.sam,output.bam)

# convert bam file to sam file
rule bam2sam:
    params: job_name = 'bam2sam.{ce}.{ch}'
    input:  'bams/{ce}/old.ch{ch}.bam'
    output: temp(rules.fix_sam.input)
    shell: '{SAMTOOLS} view -h {input} > {output}'

# simlink data to make path naming scheme consistent between PGP1_21 and PGP1_22
#rule simlinks:
#    run:
#        for ch,chpad in zip(chambers,chambers_pad):
#            shell('''
#            mkdir -p bams/PGP1_21 bams/PGP1_22 bams/PGP1_A1
#            ln -s {PGP1_21_dir}/PGP1_21_ch{chpad}.sorted.chr.bam bams/PGP1_21/ch{ch}.bam
#            ln -s {PGP1_22_dir}/PGP1_22_ch{chpad}.sorted.bam bams/PGP1_22/old.ch{ch}.bam
#            ln -s {PGP1_A1_dir}/PGP1_A1_ch{chpad}.sorted.bam bams/PGP1_A1/old.ch{ch}.bam
#            ''')
