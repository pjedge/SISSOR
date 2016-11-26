from collections import defaultdict
import pickle
import random
import estimate_parameters
import fix_sam
import chamber_allele_call_accuracy as caca
import base_calls_to_vcf

localrules: all, simlinks

chroms = ['chr{}'.format(i) for i in range(1,23)]

# cutoff values are coded as 1.0-10**(-1*x)
# so 3 => 0.999
# and 5 => 0.99999
cutoffs = list(range(1,11))
#cutoffs = [5]

rule all:
    input:
        expand('gms/{chr}.gms',chr=chroms)
        #expand('accuracy_reports/cutoff{cut}.results',cut=cutoffs)

cells = ['PGP1_21','PGP1_22','PGP1_A1']
chambers = list(range(1,25))
chambers_pad = ['{0:02d}'.format(c) for c in chambers]
grams_DNA_before_MDA  = 0.25e-12  # 0.25 pg
grams_DNA_after_MDA = 6e-9        # 6 ng
chunksize = int(5e6)

SAMTOOLS = '/opt/biotools/samtools/1.3/bin/samtools'
HG19     = '/oasis/tscc/scratch/pedge/data/genomes/hg19/hg19.fa' #'/home/wkchu/zhang_lab_oasis/resources_v2.8_b37/human_g1k_v37_decoy.fasta'
CGI_SNPs1 = '/oasis/tscc/scratch/wkchu/SISSOR/PGP1_A1/BAM/ns.gff'
CGI_SNPs2 = '/oasis/tscc/scratch/pedge/sissor_project/data/PGP1_VCFs/all.vcf'
PGP1_21_dir = '/oasis/tscc/scratch/wkchu/SISSOR/PGP1_21_highoutputmem/BAM'
PGP1_22_dir = '/oasis/tscc/scratch/wkchu/SISSOR/PGP1_22/previous/BAM'
PGP1_A1_dir = '/oasis/tscc/scratch/wkchu/SISSOR/PGP1_A1/HiSeqCombinedBAM'

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

# create chunks of hg19
# return list of (chrom,start,stop) tuples. stop is inclusive
chunklist = []
for chrom, chrlen in hg19_size_list:
    for start in range(1,chrlen+1,chunksize):
        end = start+chunksize-1 if start+chunksize-1 < chrlen else chrlen
        chunklist.append((chrom,start,end))

regions = ['{}.{}.{}'.format(chrom,start,stop) for chrom,start,stop in chunklist]

rule download_GMS:
    params: job_name = 'download_GMS.{chr}'
    output: 'gms/{chr}.gms',
    shell:
    '''
        wget http://labshare.cshl.edu/shares/schatzlab/www-data/darkmatter/by_tech/hg19_illumina/{chr}.gms.gz -O {output}.gz
        gunzip {output}.gz
    '''


rule het_vcf_file:
    params: job_name = 'het_vcf_file.{cut}'
    input:  ccf = 'output/chamber_allele_calls.out'
    output: vcf = 'het_vcfs/cutoff{cut}.whole_genome.vcf',
    run:
        base_calls_to_vcf.base_calls_to_vcf(input.ccf,(1.0-10**(-1*int(wildcards.cut))),output.vcf)

rule accuracy_aggregate:
    params: job_name = 'accuracy_aggregate.{cut}'
    input: counts = expand('accuracy_reports/split/{r}/{{cut}}/counts.p',r=regions),
            mof = expand('accuracy_reports/split/{r}/{{cut}}/mismatches',r=regions),
    output: report = 'accuracy_reports/cutoff{cut}.results',
            mof = 'accuracy_reports/cutoff{cut}.mismatches'
    run:
        caca.accuracy_aggregate(input.counts,output.report)
        shell('cat {input.mof} > {output.mof}')

rule accuracy_count:
    params: job_name = 'accuracy_count'
    input:  ccf = 'output/split/{r}.out',
            gff = CGI_SNPs1,
            vcf = CGI_SNPs2,
    output: counts_lst = expand('accuracy_reports/split/{{r}}/{cut}/counts.p',cut=cutoffs),
            mof_lst = expand('accuracy_reports/split/{{r}}/{cut}/mismatches',cut=cutoffs),
    run:
        for counts, mof, cut in zip(output.counts_lst, output.mof_lst, cutoffs):
            caca.accuracy_count(input.ccf,input.gff,input.vcf,(1.0-10**(-1*int(cut))),counts,mof)

'''
rule combine_base_calls:
    params: job_name = 'combine_base_calls'
    input:  expand('output/split/{r}.out',r=regions)
    output: 'output/chamber_allele_calls.out'
    shell: 'cat {input} > {output}'

rule call_alleles:
    params: job_name = 'call_alleles.{r}'
    input:  pileup = 'pileups/split/{r}.pileup',
            boundary_files = expand('fragment_boundary_beds/{ce}/ch{ch}.bed',ch=chambers,ce=cells),
            cov_frac_dist = 'parameters/cov_frac_dist.p',
            p_null = 'parameters/p_null.p',
            ch_priors = 'parameters/ch_priors.p',
            hom_config_probs = 'parameters/hom_config_probs.p',
            het_config_probs = 'parameters/het_config_probs.p',
            genotype_priors = 'parameters/genotype_priors.p',
            omega = 'parameters/omega.p'
    output: 'output/split/{r}.out'
    run:
        import sissorhands
        sissorhands.call_chamber_alleles(input.pileup, output[0], input.boundary_files)

rule estimate_parameters:
    params: job_name = 'estimate_parameters'
    input:  expand('parameters/split/chrX_covs.{r}.p',r=regions),
            expand('parameters/split/chamber_position_counts.{r}.p',r=regions),
            expand('parameters/split/strand_coverage_counts.{r}.p',r=regions),
            expand('parameters/split/total_sampled_cell_positions.{r}.p',r=regions),
    output: 'parameters/chrX_covs.p',
            'parameters/chamber_position_counts.p',
            'parameters/strand_coverage_counts.p',
            'parameters/total_sampled_cell_positions.p',
            'parameters/cov_frac_dist.p',
            'parameters/p_null.p',
            'parameters/ch_priors.p',
            'parameters/hom_config_probs.p',
            'parameters/het_config_probs.p',
            'parameters/genotype_priors.p',
            'parameters/omega.p'
    run:
        estimate_parameters.estimate_parameters(regions,grams_DNA_before_MDA,grams_DNA_after_MDA)

rule obtain_counts_parallel:
    params: job_name = 'obtain_counts_parallel.{r}'
    input:  pileup = 'pileups/split/{r}.pileup',
            boundary_files = expand('fragment_boundary_beds/{ce}/ch{ch}.bed',ch=chambers,ce=cells)
    output: 'parameters/split/chrX_covs.{r}.p',
            'parameters/split/chamber_position_counts.{r}.p',
            'parameters/split/strand_coverage_counts.{r}.p',
            'parameters/split/total_sampled_cell_positions.{r}.p',
    run:
        estimate_parameters.obtain_counts_parallel(input.pileup, input.boundary_files, wildcards.r)

rule pileup:
    params: job_name = 'pileup.{r}',
            chrom = lambda wildcards: wildcards.r.split('.')[0],
            start = lambda wildcards: wildcards.r.split('.')[1],
            stop  = lambda wildcards: wildcards.r.split('.')[2]
    input:  bams = expand('bams/{ce}/ch{ch}.bam',ch=chambers,ce=cells),
            bais = expand('bams/{ce}/ch{ch}.bam.bai',ch=chambers,ce=cells)
    output: pileup = 'pileups/split/{r}.pileup'
    shell:  '{SAMTOOLS} mpileup --region {params.chrom}:{params.start}-{params.stop} --adjust-MQ 50 --max-depth 100 --output-MQ --min-MQ 30 --min-BQ 20 --fasta-ref {HG19} --output {output.pileup} {input.bams}'

rule index_bam:
    params: job_name = 'index_bam{ce}.{ch}'
    input:  bam = 'bams/{ce}/ch{ch}.bam'
    output: bai = 'bams/{ce}/ch{ch}.bam.bai'
    shell:  '{SAMTOOLS} index {input.bam} {output.bai}'

# PGP1_22 and PGP1_A1 don't have chr before chromosome labels. necessary for pileup to include reference
# credit to petervangalen and Pierre Lindenbaum for this rule's code: https://www.biostars.org/p/13462/

rule fix_sam:
    params: job_name = 'fix_sam.{ce}.{ch}'
    input:  sam = 'bams/{ce}/old.ch{ch}.sam'
    output: bam = 'bams/{ce}/ch{ch}.bam'
    run:
        fix_sam.fix_sam(input.sam,output.bam)

rule bam2sam:
    params: job_name = 'bam2sam.{ce}.{ce}'
    input:  'bams/{ce}/old.ch{ce}.bam'
    output: temp(rules.fix_sam.input)
    shell: '{SAMTOOLS} view -h {input} > {output}'
'''
# simlink data to make path naming scheme consistent between PGP1_21 and PGP1_22
rule simlinks:
    run:
        for ch,chpad in zip(chambers,chambers_pad):
            shell('''
            mkdir -p bams/PGP1_21 bams/PGP1_22 bams/PGP1_A1
            ln -s {PGP1_21_dir}/PGP1_21_ch{chpad}.sorted.chr.bam bams/PGP1_21/ch{ch}.bam
            ln -s {PGP1_22_dir}/PGP1_22_ch{chpad}.sorted.bam bams/PGP1_22/old.ch{ch}.bam
            ln -s {PGP1_A1_dir}/PGP1_A1_ch{chpad}.sorted.bam bams/PGP1_A1/old.ch{ch}.bam
            ''')
