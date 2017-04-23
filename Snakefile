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
import fileIO
from plot_sissor import plot_sissor
from create_hapcut_fragment_matrix_CCF import create_hapcut_fragment_matrix_CCF
from fix_chamber_contamination import fix_chamber_contamination
#import run_tools
import fragment
import fileIO

import filter_vcf
import calculate_haplotype_statistics as chs
localrules: all, simlinks, pileup_test, make_accuracy_table

# cutoff values are phred scaled probability of error
# so 30 => 0.999 accuracy
# and 50 => 0.99999 accuracy
variant_calling_cutoffs = [7,8,9,10,30,50,70,90,110,130,150]

chroms = ['chr{}'.format(i) for i in range(1,23)]
chroms_set = set(chroms)
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

#regions = regions[5:8]

cells = ['PGP1_21','PGP1_22','PGP1_A1']
chambers = list(range(1,25))

cells = ['PGP1_21','PGP1_22','PGP1_A1']
chambers = list(range(1,25))
chambers_pad = ['{0:02d}'.format(c) for c in chambers]

BCFTOOLS = '/opt/biotools/bcftools/bin/bcftools'
SAMTOOLS = '/opt/biotools/samtools/1.3/bin/samtools'
PICARD   = '/home/pedge/installed/picard.jar'
CROSSMAP = '/home/pedge/installed/opt/python/bin/CrossMap.py'
HAPCUT2  = '/home/pedge/git/hapcut2/build/HAPCUT2'
PYPY     = '/home/pedge/installed/pypy3.3-5.5-alpha-20161013-linux_x86_64-portable/bin/pypy3.3'
WGS_VCF_URL = 'https://www.encodeproject.org/files/ENCFF995BBX/@@download/ENCFF995BBX.vcf.gz'
HG19     = '/oasis/tscc/scratch/pedge/data/genomes/hg19/hg19.fa' #'/home/wkchu/zhang_lab_oasis/resources_v2.8_b37/human_g1k_v37_decoy.fasta'
CGI_SNPs1 = '/oasis/tscc/scratch/wkchu/SISSOR/PGP1_A1/BAM/ns.gff'
CGI_SNPs2 = '/oasis/tscc/scratch/wkchu/SISSOR/NewPGP1REF/PGP-Harvard-var-hu43860C-20160106T060652Z.vcf'
BAC_VCFs = ['/home/wkchu/BACPoolPileup/Indx73.bac.pileup.vcf',
'/home/wkchu/BACPoolPileup/Indx74.bac.pileup.vcf',
'/home/wkchu/BACPoolPileup/Indx75.bac.pileup.vcf',
'/home/wkchu/BACPoolPileup/Indx76.bac.pileup.vcf',
'/home/wkchu/BACPoolPileup/Indx77.bac.pileup.vcf',
'/home/wkchu/BACPoolPileup/Indx78.bac.pileup.vcf',
'/home/wkchu/BACPoolPileup/Indx79.bac.pileup.vcf',
'/home/wkchu/BACPoolPileup/Indx80.bac.pileup.vcf',
'/home/wkchu/BACPoolPileup/Indx81.bac.pileup.vcf',
'/home/wkchu/BACPoolPileup/Indx82.bac.pileup.vcf',
'/home/wkchu/BACPoolPileup/Indx83.bac.pileup.vcf',
'/home/wkchu/BACPoolPileup/Indx84.bac.pileup.vcf',
'/home/wkchu/BACPoolPileup/Indx85.2.bac.pileup.vcf',
'/home/wkchu/BACPoolPileup/Indx85.bac.pileup.vcf',
'/home/wkchu/BACPoolPileup/Indx86.bac.pileup.vcf',
'/home/wkchu/BACPoolPileup/Indx87.bac.pileup.vcf',
'/home/wkchu/BACPoolPileup/Indx88.bac.pileup.vcf',
'/home/wkchu/BACPoolPileup/Indx89.bac.pileup.vcf',
'/home/wkchu/BACPoolPileup/Indx90.bac.pileup.vcf',
'/home/wkchu/BACPoolPileup/Indx91.bac.pileup.vcf',
'/home/wkchu/BACPoolPileup/Indx92.2.bac.pileup.vcf',
'/home/wkchu/BACPoolPileup/Indx92.bac.pileup.vcf',
'/home/wkchu/BACPoolPileup/Indx93.bac.pileup.vcf',
'/home/wkchu/BACPoolPileup/Indx94.bac.pileup.vcf',
'/home/wkchu/BACPoolPileup/Indx95.bac.pileup.vcf',
'/home/wkchu/BACPoolPileup/Indx96.bac.pileup.vcf']

PGP1_21_dir = '/oasis/tscc/scratch/wkchu/SISSOR/PGP1_21_highoutputmem/BAM'
PGP1_22_dir = '/oasis/tscc/scratch/wkchu/SISSOR/PGP1_22/2016OctMergedBAM'#'/oasis/tscc/scratch/wkchu/SISSOR/PGP1_22/previous/BAM'
PGP1_A1_dir = '/oasis/tscc/scratch/wkchu/SISSOR/PGP1_A1/2016OctMergedBAM'#'/oasis/tscc/scratch/wkchu/SISSOR/PGP1_A1/HiSeqCombinedBAM'

HAPLOTYPE_ERROR_RATES_DIR = 'haplotyping/error_rates'
HAPLOTYPE_PLOTS_DIR       = 'haplotyping/plots'
HAPLOTYPE_EXP_DIR         = 'haplotyping/experiments'

haplotyping_VCF_dir      = 'haplotyping/data/haplotyping_VCFs'
BAC_VCF_dir    = 'haplotyping/data/PGP1_VCFs_BACindex'

modes = ['same_cell','all_cell','cross_cell','ind_same_cell']
#modes = ['cross_cell']
rule all:
    input:
        expand('accuracy_reports/{mode}/cutoff10.counts.p',mode=modes),
        expand('accuracy_reports/{mode}/cutoff10.mismatches',mode=modes),
        expand('accuracy_reports/tables/strand_mismatch_{mode}.10.table.txt',mode=modes),
        'accuracy_reports/tables/unphased.10.table.txt',

        #expand('accuracy_reports/unphased/cutoff{cut}.counts.p',cut=variant_calling_cutoffs),
        #"{}/sissor_haplotype_error_genome.png".format(HAPLOTYPE_PLOTS_DIR),
        #expand('ref_dicts/{r}.ref_dict.p',r=regions),
        #expand('accuracy_reports/tables/{mode}.10.table.txt',mode=modes),

rule make_accuracy_table:
    params: job_name = 'make_accuracy_table.{report_name}.{cut}'
    input:  counts = 'accuracy_reports/{report_name}/cutoff{cut}.counts.p',
    output: table = 'accuracy_reports/tables/{report_name}.{cut}.table.txt',
    run:
        caca.generate_table(input.counts, output.table)

rule accuracy_aggregate_counts:
    params: job_name = 'accuracy_aggregate_counts.{report_name}.{cut}'
    input: counts = expand('accuracy_reports/{{report_name}}/split/{r}/{{cut}}/counts.p',r=regions),
    output: counts = 'accuracy_reports/{report_name}/cutoff{cut}.counts.p',
    run:
        caca.accuracy_aggregate(input.counts,output.counts)

rule accuracy_aggregate_mismatches:
    params: job_name = 'accuracy_aggregate_mismatches.{cut}'
    input:  mof = expand('accuracy_reports/{{report_name}}/split/{r}/{{cut}}/mismatches',r=regions),
    output: mof = 'accuracy_reports/{report_name}/cutoff{cut}.mismatches'
    shell:  'cat {input.mof} > {output.mof}'

rule generate_ref_dict:
    params: job_name = 'generate_ref_dict.{chrom}.{start}.{end}'
    input:  gff = CGI_SNPs1,
            cgi_vcf  = CGI_SNPs2,
            hg19 = HG19,
            wgs_vcf = 'wgs/wgs_hg19.vcf',
            bac_vcfs = BAC_VCFs
    output: pickle = 'ref_dicts/{chrom}.{start}.{end}.ref_dict.p'
    run:
        reg_chrom = str(wildcards.chrom)
        reg_start = int(wildcards.start)
        reg_end   = int(wildcards.end)
        region = (reg_chrom, reg_start, reg_end)
        caca.generate_ref_dict(input.gff, input.wgs_vcf, input.cgi_vcf, input.bac_vcfs, input.hg19, region, output.pickle)

rule accuracy_count_unphased:
    params: job_name = 'accuracy_count_unphased.{r}'
    input:  ccf = 'base_calls/unphased/{r}.out',
            ref_dict = 'ref_dicts/{r}.ref_dict.p',
            gms = 'gms/{r}.gms'
    output: counts_lst = expand('accuracy_reports/unphased/split/{{r}}/{cut}/counts.p',cut=variant_calling_cutoffs),
            mof_lst = expand('accuracy_reports/unphased/split/{{r}}/{cut}/mismatches',cut=variant_calling_cutoffs),
    run:
        for counts, mof, cut in zip(output.counts_lst, output.mof_lst, variant_calling_cutoffs):
            caca.accuracy_count_unphased(input.ccf, input.ref_dict, input.gms, (1.0-10**(-0.1*float(cut))), counts, mof)

rule accuracy_count_phased:
    params: job_name = 'accuracy_count_{mode,(same_cell|all_cell|cross_cell|ind_same_cell)}.{r}'
    input:  ccf = 'base_calls/phased/{r}.ccf',
            ref_dict = 'ref_dicts/{r}.ref_dict.p',
            gms = 'gms/{r}.gms'
    output: counts = 'accuracy_reports/{mode,(same_cell|all_cell|cross_cell|ind_same_cell)}/split/{r}/10/counts.p',
            mof = 'accuracy_reports/{mode,(same_cell|all_cell|cross_cell|ind_same_cell)}/split/{r}/10/mismatches',
            smf = 'accuracy_reports/strand_mismatch_{mode,(same_cell|all_cell|cross_cell|ind_same_cell)}/split/{r}/10/counts.p'
    run:
        CUT = 0.9
        caca.accuracy_count_phased(input.ccf, input.ref_dict, input.gms, CUT, output.counts, output.mof, output.smf,separate_haplotypes=True,mode=wildcards.mode)

# add phasing information to unphased chamber call files
rule annotate_assigned_fragments:
    params: job_name = 'annotate_assigned_fragments.{r}'
    input:  ccf = 'base_calls/unphased/{r}.out',
            asn = 'fragment_haplotype_assignments/all',
    output: ccf = 'base_calls/phased/{r}.ccf',
    run:
        fragment_haplotype_assignment.annotate_assigned_fragments(input.ccf,input.asn,output.ccf)

# combine separate files for assembled fragment haplotypes
rule combine_assigned_fragment_files:
    params: job_name = 'combine_assigned_fragment_files'
    input:  sep = expand('fragment_haplotype_assignments/{chrom}',chrom=chroms),
    output: combined = 'fragment_haplotype_assignments/all',
    shell:  'cat {input.sep} > {output.combined}'

# assign fragments to assembled haplotypes
rule assign_fragment_haplotypes:
    params: job_name = 'assign_fragment_haplotypes.{chrom}'
    input: frag = 'haplotyping/data/PGP1_ALL/fragmat/cov1_strict/{chrom}',
           vcf  = 'haplotyping/data/haplotyping_VCFs/{chrom}.vcf',
           hap  = 'haplotyping/experiments/hapcut2_PGP1_ALL/cov1_strict/{chrom}.output'
    output: sep = 'fragment_haplotype_assignments/{chrom}',
    run:
        fragment_haplotype_assignment.assign_fragment_haplotypes(input.frag,input.vcf,output.sep,input.hap)

# PLOT RESULTS
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


# RUN HAPCUT2
rule prune_haplotype:
    params:
        job_name = "{s}.{x}.prune_haplotype",
    input:
        hapblocks = expand("{E}/hapcut2_{{s}}/{{x}}/{c}.output.uncorrected",E=HAPLOTYPE_EXP_DIR,c=chroms)
    output:
        hapblocks = expand("{E}/hapcut2_{{s}}/{{x}}/{c}.output",E=HAPLOTYPE_EXP_DIR,c=chroms)
    run:
        for i, o in zip(input.hapblocks,output.hapblocks):
            fileIO.prune_hapblock_file(i, o, snp_conf_cutoff=0.95, split_conf_cutoff=0.95, use_refhap_heuristic=True) #split_conf_cutoff=0.9999

# RUN HAPCUT2
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

# COMBINE FRAGMENT MATRICES
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

rule generate_fragmatrix:
    params: job_name = 'generate_fragmatrix'
    input:  ccf      = 'haplotyping/data/PGP1_ALL.het.whole_genome.ccf',
            vcfs     = expand('{v}/{c}.vcf',v=haplotyping_VCF_dir,c=chroms),
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

#rule combine_intersect:
#    params: job_name = "combine_intersect",
#    input: expand("haplotyping/data/haplotyping_VCFs/{c}.vcf",c=chroms)
#    output: 'haplotyping/data/haplotyping_VCFs/all.vcf'
#    shell:
#        '''cat {input} > {output}'''

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
        frag_file = "haplotyping/data/BAC_frags_fixed/{c}",
        vcf_file  = "%s/{c}.vcf" % BAC_VCF_dir
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
        var_vcfs = "%s/{c}.vcf" % BAC_VCF_dir,
        frags = "haplotyping/data/BAC_frags/{c}",
        #WGS_variants = "%s/{c}.vcf" % haplotyping_VCF_dir,
    output:
        fixed = "haplotyping/data/BAC_frags_fixed/{c}"
    run:
        fix_chamber_contamination(input.frags,input.var_vcfs,output.fixed,threshold=2, min_coverage=0,mode='basic',vcf_filter=None) #input.WGS_variants)

#rule filter_CGI_VCF:
#    params: job_name  = "filter_CGI_VCF",
#    input: vcfs = expand("{d}/{c}.vcf",d=BAC_VCF_dir,c=chroms),
#           filter_set = "/oasis/tscc/scratch/wkchu/SISSOR/NewPGP1REF/PGP-Harvard-var-hu43860C-20160106T060652Z.vcf"
#    output: vcfs = expand("{d}/{c}.vcf",d=haplotyping_VCF_dir,c=chroms)
#    run:
#        filter_vcf.filter_vcf(input.vcfs,output.vcfs,input.filter_set)

rule combine_haplotyping_VCFs:
    params: job_name = "combine_haplotyping_VCFs",
    input: expand("haplotyping/data/haplotyping_VCFs/{c}.vcf",c=chroms)
    output: 'haplotyping/data/haplotyping_VCFs/all.vcf'
    shell:
        '''cat {input} > {output}'''

rule filter_haplotyping_VCFs:
    params: job_name  = "filter_haplotyping_VCF",
    input: vcf = "wgs/wgs_hg19.vcf"
    output: vcfs = expand("{d}/{c}.vcf",d=haplotyping_VCF_dir,c=chroms)
    run:
        filter_vcf.filter_wgs_vcf(input.vcf,output.vcfs,chroms) # input.filter_set

rule sort_wgs_SNPs:
    params: job_name  = 'liftover_SNPs',
    input:  vcf = 'wgs/wgs_hg19.autosomes.unsorted.vcf',
    output: vcf = 'wgs/wgs_hg19.vcf'
    shell:
        '''cat {input} | vcf-sort -c -p 4 > {output}'''

rule filter_wgs_autosomes:
    params: job_name  = 'filter_wgs_autosomes',
    input:  vcf = 'wgs/wgs_hg19.unsorted.vcf',
    output: vcf = 'wgs/wgs_hg19.autosomes.unsorted.vcf'
    run:
        with open(input.vcf,'r') as inf, open(output.vcf,'w') as outf:
            for line in inf:
                if line[:1] == '#':
                    continue
                el = line.strip().split('\t')
                if len(el) < 5:
                    continue
                if el[0] in chroms_set:
                    print(line.strip(),file=outf)

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
    output: ccf = 'base_calls/unphased/{r}.out'
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
        estimate_parameters.estimate_parameters(regions)

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
    params: job_name = lambda wildcards: 'index_bam.{}'.format(str(wildcards.x).replace("/", "."))
    input:  bam = '{x}.bam'
    output: bai = '{x}.bam.bai'
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
# simlink data to make path naming scheme consistent between PGP1_21 and PGP1_22
'''
rule simlinks:
    run:
        for ch,chpad in zip(chambers,chambers_pad):
            shell('''
            mkdir -p bams/PGP1_21 bams/PGP1_22 bams/PGP1_A1
            ln -s {PGP1_21_dir}/PGP1_21_ch{chpad}.sorted.chr.bam bams/PGP1_21/ch{ch}.bam
            ln -s {PGP1_22_dir}/PGP1_22_ch{chpad}.sorted.bam bams/PGP1_22/old.ch{ch}.bam
            ln -s {PGP1_A1_dir}/PGP1_A1_ch{chpad}.sorted.bam bams/PGP1_A1/old.ch{ch}.bam
            ''')
