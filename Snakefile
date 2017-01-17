from collections import defaultdict
import pickle
import random
import estimate_parameters
import fix_sam
import chamber_allele_call_accuracy as caca
import base_calls_to_vcf
import strand_pairing_analysis
from itertools import product

localrules: all, simlinks, pileup_test

# cutoff values are coded as 1.0-10**(-1*x)
# so 3 => 0.999
# and 5 => 0.99999
cutoffs = [1,3,5,7,9,11,13,15] #list(range(1,17))

paired_cutoffs = [2,3]

het_vcf_cutoffs = list(range(1,6))
chroms = ['chr{}'.format(i) for i in range(1,23)]
chroms_XY = chroms + ['chrX','chrY']
chunksize = int(5e6)

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
#regions = ['chr6.135000001.140000000']
cells = ['PGP1_21','PGP1_22','PGP1_A1']
chambers = list(range(1,25))

cells = ['PGP1_21','PGP1_22','PGP1_A1']
chambers = list(range(1,25))
chambers_pad = ['{0:02d}'.format(c) for c in chambers]
grams_DNA_before_MDA  = 0.25e-12  # 0.25 pg
grams_DNA_after_MDA = 6e-9        # 6 ng

SAMTOOLS = '/opt/biotools/samtools/1.3/bin/samtools'
PICARD   = '/home/pedge/installed/picard.jar'
CROSSMAP = '/home/pedge/installed/opt/python/bin/CrossMap.py'
HAPCUT2  = '/home/pedge/git/hapcut2/build/HAPCUT2'
PYPY     = '/home/pedge/installed/pypy3.3-5.5-alpha-20161013-linux_x86_64-portable/bin/pypy3.3'
WGS_VCF_URL = 'https://www.encodeproject.org/files/ENCFF995BBX/@@download/ENCFF995BBX.vcf.gz'
HG19     = '/oasis/tscc/scratch/pedge/data/genomes/hg19/hg19.fa' #'/home/wkchu/zhang_lab_oasis/resources_v2.8_b37/human_g1k_v37_decoy.fasta'
CGI_SNPs1 = '/oasis/tscc/scratch/wkchu/SISSOR/PGP1_A1/BAM/ns.gff'
CGI_SNPs2 = '/oasis/tscc/scratch/pedge/sissor_project/data/PGP1_VCFs_BACindex/all.vcf'
PGP1_21_dir = '/oasis/tscc/scratch/wkchu/SISSOR/PGP1_21_highoutputmem/BAM'
PGP1_22_dir = '/oasis/tscc/scratch/wkchu/SISSOR/PGP1_22/2016OctMergedBAM'#'/oasis/tscc/scratch/wkchu/SISSOR/PGP1_22/previous/BAM'
PGP1_A1_dir = '/oasis/tscc/scratch/wkchu/SISSOR/PGP1_A1/2016OctMergedBAM'#'/oasis/tscc/scratch/wkchu/SISSOR/PGP1_A1/HiSeqCombinedBAM'

HAPLOTYPE_ERROR_RATES_DIR = 'haplotyping/error_rates'
HAPLOTYPE_PLOTS_DIR       = 'haplotyping/plots'
HAPLOTYPE_EXP_DIR         = 'haplotyping/experiments'
BAC_vcf_dir               = 'haplotyping/data/PGP1_VCFs_BACindex'

rule all:
    input:
        "{}/sissor_haplotype_error_chromosome.png".format(HAPLOTYPE_PLOTS_DIR),
        "{}/sissor_haplotype_error_genome.png".format(HAPLOTYPE_PLOTS_DIR),
        #expand('accuracy_reports/unpaired/cutoff{cut}.results',cut=cutoffs),
        #expand('accuracy_reports/paired_same_cell/cutoff{cut}.results',cut=paired_cutoffs),
        #expand('accuracy_reports/paired_all_cells/cutoff{cut}.results',cut=paired_cutoffs),


rule het_vcf_file:
    params: job_name = 'het_vcf_file.{cut}'
    input:  ccfs    = expand('base_calls/unpaired/{r}.out',r=regions)
    output: vcf     = 'het_vcfs/cutoff{cut}/whole_genome.vcf',
            ccf_het = 'het_vcfs/cutoff{cut}/whole_genome.out'
    run:
        base_calls_to_vcf.base_calls_to_vcf(input.ccfs,(1.0-10**(-1*int(wildcards.cut))),output.vcf, output.ccf_het)

rule accuracy_count_strand_pairing:
    params: job_name = 'accuracy_count_strand_pairing.{mode}.{chrom}.{start}.{end}'
    input:  ccf = 'base_calls/paired/{chrom}.{start}.{end}.out',
            gff = CGI_SNPs1,
            cgi_vcf  = CGI_SNPs2,
            gms = 'gms/{chrom}.{start}.{end}.gms',
            hg19 = HG19,
            wgs_vcf = 'wgs/wgs_hg19.vcf'
    output: counts_lst = temp(expand('accuracy_reports/paired_{{mode}}/split/{{chrom}}.{{start}}.{{end}}/{cut}/counts.p',cut=paired_cutoffs)),
            mof_lst = temp(expand('accuracy_reports/paired_{{mode}}/split/{{chrom}}.{{start}}.{{end}}/{cut}/mismatches',cut=paired_cutoffs)),
    run:
        reg_chrom = str(wildcards.chrom)
        reg_start = int(wildcards.start)
        reg_end   = int(wildcards.end)
        region = (reg_chrom, reg_start, reg_end)

        if wildcards.mode not in ['same_cell','all_cells']:
            exit(1)
        same_cell_only = (wildcards.mode == 'same_cell')


        for counts, mof, cut in zip(output.counts_lst, output.mof_lst, paired_cutoffs):
            CUT = (1.0-10**(-1*int(cut))) if cut != None else None
            caca.accuracy_count_strand_pairing(input.ccf, input.gff, input.wgs_vcf, input.cgi_vcf, input.gms, input.hg19, region, CUT, counts, mof, same_cell_only=same_cell_only)

rule annotate_strand_pairs:
    params: job_name = 'annotate_strand_pairs.{r}'
    input:  ccf = 'base_calls/unpaired/{r}.out',
            pairs = 'strand_pairs/all',
    output: ccf = 'base_calls/paired/{r}.out',
    run:
        strand_pairing_analysis.annotate_paired_strands(input.ccf,input.pairs,output.ccf)

rule combine_paired_strand_files:
    params: job_name = 'combine_paired_strand_files'
    input:  sep = expand('strand_pairs/{chrom}',chrom=chroms),
    output: combined = 'strand_pairs/all',
    shell:  'cat {input.sep} > {output.combined}'

rule pair_strands:
    params: job_name = 'pair_strands'
    input: frag = 'haplotyping/data/PGP1_ALL/fragmat/cov1_basic/{chrom}',
           vcf  = 'haplotyping/data/PGP1_VCFs_BACindex/{chrom}.vcf',
           hap  = 'haplotyping/experiments/hapcut2_PGP1_ALL/cov1_basic/{chrom}.output'
    output: sep = 'strand_pairs/{chrom}',
    run:
        strand_pairing_analysis.pair_strands(input.frag,input.vcf,output.sep,input.hap)


# PLOT RESULTS
rule plot_hapcut2_results:
    params:
        job_name = "plot_hapcut2_sissor"
    input:
        stats_file  = "{}/hapcut2.stats.p".format(HAPLOTYPE_ERROR_RATES_DIR),
        labels_file = "{}/hapcut2.labels.p".format(HAPLOTYPE_ERROR_RATES_DIR)
    output:
        plot1 = "%s/sissor_haplotype_error_chromosome.png" % HAPLOTYPE_PLOTS_DIR,
        plot2 = "%s/sissor_haplotype_error_genome.png" % HAPLOTYPE_PLOTS_DIR
    run:
        data = pickle.load(open(input.stats_file,"rb"))
        labels = pickle.load(open(input.labels_file,"rb"))
        plot_data.plot_experiment_sissor(data,labels,[],output.plot1)
        plot_sissor(data,labels,output.plot2)

exp =['cov1_none','cov1_basic','cov1_strict','cov2_basic']
exp_labels=['No processing','Basic Processing','Strict Processing','Basic Processing, coverage >= 2']
rule calculate_error_rates:
    params:
        job_name = "hapcut2_error_rates"
    input:
        hapblocks = expand("{E}/hapcut2_PGP1_ALL/{exp}/{c}.output",E=experiments_dir,exp=exp,c=chroms),
        runtimes  = expand("{E}/hapcut2_PGP1_ALL/{exp}/{c}.runtime",E=experiments_dir,exp=exp,c=chroms),
        fragmats  = expand("{d}/PGP1_ALL/fragmat/{exp}/{c}",d=data_dir,exp=exp,c=chroms),
        var_vcfs  = expand("{v}/{c}.vcf",v=BAC_vcf_dir, c=chroms),
        bac_haps  = expand("sissor_project/data/BAC/{c}.filtered",c=chroms)
    output:
        stats_file  = "{}/hapcut2.stats.p".format(HAPLOTYPE_ERROR_RATES_DIR),
        labels_file = "{}/hapcut2.labels.p".format(HAPLOTYPE_ERROR_RATES_DIR)
    run:
        # list of lists of error results,
        data   = [] # index of the list is a 'condition' each inner list has 23 error results (each chrom).
        labels = exp_labels

        for x in exp: # sample
            datalist = []
            for c in chroms: # chromosome
                assembly_file = "{}/hapcut2_PGP1_ALL/{}/{}.output".format(HAPLOTYPE_EXP_DIR,x,c)
                runtime_file  = "{}/hapcut2_PGP1_ALL/{}/{}.runtime".format(HAPLOTYPE_EXP_DIR,x,c)
                frag_file     = "{}/PGP1_ALL/fragmat/{}/{}".format(data_dir,x,c)
                vcf_file      = "{}/{}.vcf".format(BAC_vcf_dir,c)
                truth_file    = "haplotyping/data/BAC/{}.filtered".format(c)
                err = error_rates.hapblock_hapblock_error_rate(truth_file, assembly_file, frag_file, vcf_file, runtime_file, use_SNP_index=True)
                datalist.append(err)

            print("{} results over all chromosomes:".format(x))
            print(sum(datalist,error_rates.error_result()))

            data.append(datalist)

        pickle.dump(data,open(output.stats_file,"wb"))
        pickle.dump(labels,open(output.labels_file,"wb"))


# RUN HAPCUT2
rule prune_haplotype:
    params:
        job_name = "{s}.{x}.prune_haplotype",
    input:
        hapblocks = expand("{E}/hapcut2_{{s}}/{{x}}/{c}.output.uncorrected",E=experiments_dir,c=chroms)
    output:
        hapblocks = expand("{E}/hapcut2_{{s}}/{{x}}/{c}.output",E=experiments_dir,c=chroms)
    run:
        for i, o in zip(input.hapblocks,output.hapblocks):
            fileIO.prune_hapblock_file(i, o, snp_conf_cutoff=0.95, split_conf_cutoff=-1, use_refhap_heuristic=True) #split_conf_cutoff=0.9999

# RUN HAPCUT2
rule run_hapcut2:
    params:
        job_name = "{s}.{c}.{x}.hapcut",
    input:
        frag_file = "%s/{s}/fragmat/{x}/{c}" % data_dir,
        vcf_file  = lambda wildcards: expand("{v}/{c}.vcf", v=BAC_vcf_dir,c=wildcards.c)
    output:
        hapblocks = "{E}/hapcut2_{s}/{x}/{c}.output.uncorrected",
        runtime = "{E}/hapcut2_{s}/{x}/{c}.runtime"
    run:
        # run hapcut
        runtime = run_tools.run_hapcut2(HAPCUT2, input.frag_file, input.vcf_file, output.hapblocks, 5, 0.8, '--ea 1')
        with open(output.runtime,'w') as rf:
            print(runtime, file=rf)

# COMBINE FRAGMENT MATRICES
rule fix_fragmat:
    params:
        job_name  = "fix_fragmat.{x}",
    input:
        var_vcfs = expand("{v}/{CHR}.vcf",v=BAC_vcf_dir,CHR=chroms),
        P_ALL = expand("{dat}/PGP1_ALL/augmented_fragmat/{c}",dat=data_dir,c=chroms)
    output:
        fixed = expand("{{data_dir}}/PGP1_ALL/fragmat/{{x}}/{c}",c=chroms)
    run:
        mincov = 0
        if 'cov2' in wildcards.x:
            mincov = 2
        elif 'cov3' in wildcards.x:
            mincov = 3

        if 'none' in wildcards.x:
            mode = 'none'
        elif 'basic' in wildcards.x:
            mode = 'basic'
        elif 'strict' in wildcards.x:
            mode = 'strict'

        for i,v,o in zip(input.P_ALL, input.var_vcfs, output.fixed):
            fix_chamber_contamination(i,v,o,threshold=2, min_coverage=mincov,mode=mode)

rule generate_fragmatrix:
    params: job_name = 'generate_fragmatrix'
    input:  ccf      = 'haplotyping/data/PGP1_ALL.het.whole_genome.ccf',
            vcfs     = expand('{v}/{c}.vcf',v=BAC_vcf_dir,c=chroms),
            bounds   = expand('eric_fragment_boundary_beds/{P[0]}/ch{P[1]}.bed',P=product(cells,chambers)),
    output: fragmat  = expand('haplotyping/data/PGP1_ALL/augmented_fragmat/{c}',c=chroms),

    run:
        odir = 'haplotyping/data/PGP1_ALL/augmented_fragmat'
        # generate fragment matrix from sissorhands base calls
        create_hapcut_fragment_matrix_CCF(chamber_call_file=input.ccf, variant_vcf_files=input.vcfs,fragment_boundary_files=input.bounds, output_dir=odir)

rule combine_filtered:
    params: job_name = 'combine_filtered'
    input:  ccf = expand('haplotyping/data/ccf_split/{r}.ccf',r=regions),
            vcf = expand('haplotyping/data/vcf_split/{r}.vcf',r=regions)
    output: ccf = 'haplotyping/data/PGP1_ALL.het.whole_genome.ccf',
            vcf = 'haplotyping/data/PGP1_ALL.het.whole_genome.vcf',
    shell:
        '''cat {input.ccf} > {output.ccf}
           cat {input.vcf} > {output.vcf}'''

rule filter_CCF:
    params: job_name = 'filter_CCF_{chrom}.{start}.{end}'
    input:  ccf      = 'base_calls/split/{chrom}.{start}.{end}.out',
            vcf      = 'haplotyping/data/PGP1_VCFs_BACindex/all.vcf'
    output: ccf      = temp('haplotyping/data/ccf_split/{chrom}.{start}.{end}.ccf'),
            vcf      = temp('haplotyping/data/vcf_split/{chrom}.{start}.{end}.vcf'),
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

# assemble BAC haplotypes

rule prune_BAC:
    params:
        job_name = "BAC.prune_haplotype",
    input:
        hapblocks = expand("haplotyping/data/BAC/{c}",c=chroms)
    output:
        hapblocks = expand("haplotyping/data/BAC/{c}.filtered",c=chroms)
    run:
        for i, o in zip(input.hapblocks,output.hapblocks):
            fileIO.prune_hapblock_file(i, o, snp_conf_cutoff=0.9999, split_conf_cutoff=0.9999, use_refhap_heuristic=True) #split_conf_cutoff=0.9999

# RUN HAPCUT2

rule run_hapcut2_BAC:
    params:
        job_name = "{c}.BAC.hapcut2",
    input:
        frag_file = "%s/BAC_frags_fixed/{c}" % data_dir,
        vcf_file  = "%s/{c}.vcf" % BAC_vcf_dir
    output:
        hapblocks = "haplotyping/data/BAC/{c}"
    shell:
        '''
        {HAPCUT2} --fragments {input.frag_file} --vcf {input.vcf_file} --output {output.hapblocks} --ea 1
        '''

# REMOVE SWITCH ERRORS FROM BAC FRAGMENTS
rule fix_fragmat_BAC:
    params:
        job_name  = "fix_fragmat_BAC_{c}",
    input:
        var_vcfs = "%s/{c}.vcf" % BAC_vcf_dir,
        frags = "haplotyping/data/BAC_frags/{c}"
    output:
        fixed = "haplotyping/data/BAC_frags_fixed/{c}"
    run:
        fix_chamber_contamination(input.frags,input.var_vcfs,output.fixed,threshold=2, min_coverage=0,mode='basic')



rule accuracy_aggregate:
    params: job_name = 'accuracy_aggregate.{cut}'
    input: counts = expand('accuracy_reports/{{report_name}}/split/{r}/{{cut}}/counts.p',r=regions),
            mof = expand('accuracy_reports/{{report_name}}/split/{r}/{{cut}}/mismatches',r=regions),
    output: report = 'accuracy_reports/{report_name}/cutoff{cut}.results',
            mof = 'accuracy_reports/{report_name}/cutoff{cut}.mismatches'
    run:
        caca.accuracy_aggregate(input.counts,output.report)
        shell('cat {input.mof} > {output.mof}')

rule accuracy_count_unpaired:
    params: job_name = 'accuracy_count_unpaired.{chrom}.{start}.{end}'
    input:  ccf = 'base_calls/unpaired/{chrom}.{start}.{end}.out',
            gff = CGI_SNPs1,
            cgi_vcf  = CGI_SNPs2,
            gms = 'gms/{chrom}.{start}.{end}.gms',
            hg19 = HG19,
            wgs_vcf = 'wgs/wgs_hg19.vcf',
            #bed_filter = 'eric_fragment_boundary_beds/sorted_merged.bed'
    output: counts_lst = temp(expand('accuracy_reports/unpaired/split/{{chrom}}.{{start}}.{{end}}/{cut}/counts.p',cut=cutoffs)),
            mof_lst = temp(expand('accuracy_reports/unpaired/split/{{chrom}}.{{start}}.{{end}}/{cut}/mismatches',cut=cutoffs)),
    run:
        reg_chrom = str(wildcards.chrom)
        reg_start = int(wildcards.start)
        reg_end   = int(wildcards.end)
        region = (reg_chrom, reg_start, reg_end)
        for counts, mof, cut in zip(output.counts_lst, output.mof_lst, cutoffs):
            caca.accuracy_count(input.ccf, input.gff, input.wgs_vcf, input.cgi_vcf, input.gms, input.hg19, region, (1.0-10**(-1*int(cut))), counts, mof)

rule liftover_wgs_SNPs:
    params: job_name  = 'liftover_SNPs',
    input:  vcf = 'wgs/wgs_grch38.vcf',
            chain = 'hg38ToHg19.over.chain',
            hg19 = HG19
    output: vcf = 'wgs/wgs_hg19.vcf',
    shell:
        '''
        {CROSSMAP} \
        vcf \
        {input.chain} \
        {input.vcf} \
        {input.hg19} \
        {output.vcf}
        '''

rule download_WGS_vcf:
    params: job_name = 'download_WGS_vcf'
    output: vcf = 'wgs/wgs_grch38.vcf'
    shell:
        '''
        wget {WGS_VCF_URL} -O {output.vcf}.gz
        mv {output.vcf}.gz {output.vcf}
        '''

rule split_gms:
    params: job_name  = 'split_gms',
    input:  gms = expand('gms/{c}.gms',c=chroms_XY)
    output: gms = expand('gms/{r}.gms',r=regions)
    run:
        caca.split_gms(input.gms, chunklist, output.gms)

rule download_GMS:
    params: job_name = 'download_GMS.{chr}'
    output: 'gms/{chr}.gms',
    shell:
        '''
        wget http://labshare.cshl.edu/shares/schatzlab/www-data/darkmatter/by_tech/hg19_illumina/{wildcards.chr}.gms.gz -O {output}.gz
        gunzip {output}.gz
        '''
'''
rule call_alleles:
    params: job_name = 'call_alleles.{r}'
    input:  pileup = 'pileups/split/{r}.pileup',
            #boundary_files = expand('fragment_boundary_beds/{P[0]}/ch{P[1]}.bed',P=product(cells, chambers)),
            MDA_dist = 'parameters/MDA_dist.p',
            cov_frac_dist = 'parameters/cov_frac_dist.p',
            p_null = 'parameters/p_null.p',
            ch_priors = 'parameters/ch_priors.p',
            hom_config_probs = 'parameters/hom_config_probs.p',
            het_config_probs = 'parameters/het_config_probs.p',
            haploid_genotype_priors = 'parameters/haploid_genotype_priors.p',
            diploid_genotype_priors = 'parameters/diploid_genotype_priors.p'
    output: ccf = 'base_calls/unpaired/{r}.out'
    shell:  '{PYPY} sissorhands.py -i {input.pileup} -o {output.ccf}'
    #run:
        #import sissorhands
        #sissorhands.call_chamber_alleles(input.pileup, output.ccf)

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
            #'parameters/total_sampled_cell_positions.p',
            #'parameters/omega.p'
    run:
        estimate_parameters.estimate_parameters(regions,grams_DNA_before_MDA,grams_DNA_after_MDA)

rule obtain_counts_parallel:
    params: job_name = 'obtain_counts_parallel.{r}'
    input:  pileup = 'pileups/split/{r}.pileup',
            #boundary_files = expand('fragment_boundary_beds/{P[0]}/ch{P[1]}.bed',P=product(cells, chambers))
    output: 'parameters/split/chrX_MDA_fracs.{r}.p',
            'parameters/split/chrX_covs.{r}.p',
            'parameters/split/chamber_position_counts.{r}.p',
            'parameters/split/strand_coverage_counts.{r}.p',
            #'parameters/split/total_sampled_cell_positions.{r}.p',
    run:
        estimate_parameters.obtain_counts_parallel(input.pileup, boundary_files=None, suffix=wildcards.r)

rule pileup:
    params: job_name = 'pileup.{r}',
            chrom = lambda wildcards: wildcards.r.split('.')[0],
            start = lambda wildcards: wildcards.r.split('.')[1],
            stop  = lambda wildcards: wildcards.r.split('.')[2]
    input:  bams = expand('bams/{P[0]}/ch{P[1]}.bam',P=product(cells, chambers)),
            bais = expand('bams/{P[0]}/ch{P[1]}.bam.bai',P=product(cells, chambers))
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
    params: job_name = 'bam2sam.{ce}.{ch}'
    input:  'bams/{ce}/old.ch{ch}.bam'
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
