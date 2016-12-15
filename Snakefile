# author: Peter Edge
# 3/29/2016
# email: pedge@eng.ucsd.edu

# this is a snakemake Snakefile, written using snakemake 3.5.5
# configure paths and values in cluster.yaml and config.yaml to your system
# example execution for a TORQUE cluster:
# snakemake -j 200 --cluster-config cluster.yaml --cluster "qsub -A {cluster.group} -V -q {cluster.queue} -o {cluster.qsub_stdout_dir}/{params.job_name}.o -e {cluster.qsub_stdout_dir}/{params.job_name}.e -N {params.job_name} -l nodes=1:ppn={cluster.ppn} -l walltime={cluster.walltime} -M {cluster.email} -m e -d {cluster.working_dir}" --local-cores 1
sys.path.append('/home/peter/git/HapTools')
configfile: "config.yaml"
localrules: all, plot_hapcut2_results, simlinks, clean

import sys
sys.path.append(config['haptools_dir'])
import run_tools
import error_rates
import fileIO
import plot_data
import os
from os.path import join
from create_hapcut_fragment_matrix import create_hapcut_fragment_matrices_freebayes
from create_hapcut_fragment_matrix_CCF import create_hapcut_fragment_matrix_CCF
from fix_chamber_contamination import fix_chamber_contamination
import pickle
from collections import defaultdict
from plot_sissor import plot_sissor

# qsub stderr and stdout directories are not automatically created by snakemake!
if not os.path.exists(config['qsub_stdout_dir']):
    os.makedirs(config['qsub_stdout_dir'])

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
 ('chr22', 51304566)]
 #
 #('chrX', 155270560),
 #('chrY', 59373566)]

# create chunks of hg19
# return list of (chrom,start,stop) tuples. stop is inclusive
chunklist = []

for chrom, chrlen in hg19_size_list:
    for start in range(1,chrlen+1,chunksize):
        end = start+chunksize-1 if start+chunksize-1 < chrlen else chrlen
        chunklist.append((chrom,start,end))

regions = ['{}.{}.{}'.format(chrom,start,stop) for chrom,start,stop in chunklist]

# specify chambers
data_dir = config["data_dir"]
experiments_dir = config["experiments_dir"]
plots_dir = config["plots_dir"]
chambers = list(range(1,25))
chambers_pad = ['{0:02d}'.format(c) for c in chambers]
chroms  = ['chr{}'.format(i) for i in range(1,23)]
samples = ['PGP1_ALL','PGP1_21','PGP1_22','PGP1_A1']
cells=samples[1:]

variant_vcf_dir = join(data_dir,'PGP1_VCFs')
variant_vcf_dir_fp = join(data_dir,'PGP1_VCFs_w_false_positives')

rule all:
    input:
        expand("{P}/sissor_haplotype_error_chromosome.png",P=config['plots_dir']),
        expand("{P}/sissor_haplotype_error_genome.png",P=config['plots_dir']),
        #expand("{d}/{s}/fixed_beds/ch{ch}.bed",d=data_dir,s=cells,ch=chambers)

# PLOT RESULTS
rule plot_hapcut2_results:
    params:
        job_name = "plot_hapcut2_sissor"
    input:
        stats_file  = "{}/hapcut2.stats.p".format(config['error_rates_dir']),
        labels_file = "{}/hapcut2.labels.p".format(config['error_rates_dir'])
    output:
        plot1 = "%s/sissor_haplotype_error_chromosome.png" % config['plots_dir'],
        plot2 = "%s/sissor_haplotype_error_genome.png" % config['plots_dir']

    run:
        data = pickle.load(open(input.stats_file,"rb"))
        labels = pickle.load(open(input.labels_file,"rb"))
        plot_data.plot_experiment_sissor(data[1:3],labels[1:3],[],output.plot1)
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
        var_vcfs  = expand("{v}/{c}.vcf",v=variant_vcf_dir,c=chroms),
        bac_haps  = expand("{bac}/{c}.filtered",bac=config['BAC_hapblocks'], c=chroms)
    output:
        stats_file  = "{}/hapcut2.stats.p".format(config['error_rates_dir']),
        labels_file = "{}/hapcut2.labels.p".format(config['error_rates_dir'])
    run:
        # list of lists of error results,
        data   = [] # index of the list is a 'condition' each inner list has 23 error results (each chrom).
        labels = exp_labels

        for x in exp: # sample
            datalist = []
            for c in chroms: # chromosome
                assembly_file = "{}/hapcut2_PGP1_ALL/{}/{}.output".format(config['experiments_dir'],x,c)
                runtime_file  = "{}/hapcut2_PGP1_ALL/{}/{}.runtime".format(config['experiments_dir'],x,c)
                frag_file     = "{}/PGP1_ALL/fragmat/{}/{}".format(data_dir,x,c)
                vcf_file      = "{}/{}.vcf".format(variant_vcf_dir,c)
                truth_file    = "{}/{}.filtered".format(config['BAC_hapblocks'], c)
                err = error_rates.hapblock_hapblock_error_rate(truth_file, assembly_file, frag_file, vcf_file, runtime_file, use_SNP_index=False)
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
            fileIO.prune_hapblock_file(i, o, snp_conf_cutoff=0.95, split_conf_cutoff=0.49, use_refhap_heuristic=True) #split_conf_cutoff=0.9999

# RUN HAPCUT2
rule run_hapcut2:
    params:
        job_name = "{s}.{c}.{x}.hapcut",
    input:
        frag_file = "%s/{s}/fragmat/{x}/{c}" % data_dir,
        vcf_file  = lambda wildcards: expand("{v}/{c}.vcf", v=variant_vcf_dir,c=wildcards.c)
    output:
        hapblocks = "{E}/hapcut2_{s}/{x}/{c}.output.uncorrected",
        runtime = "{E}/hapcut2_{s}/{x}/{c}.runtime"
    run:
        # run hapcut
        runtime = run_tools.run_hapcut2(config['hapcut2'], input.frag_file, input.vcf_file, output.hapblocks, 5, 0.8, '--ea 1')
        with open(output.runtime,'w') as rf:
            print(runtime, file=rf)

# CREATE BED FILES FOR THE NEW FIXED FRAGMENTS
'''
rule make_fixed_beds:
    params:
        job_name  = "make_fixed_beds",
    input:
        expand("{d}/PGP1_ALL/fragmat/cov1_strict/{c}",d=data_dir,c=chroms)
    output:
        expand("{d}/{s}/fixed_beds/ch{ch}.bed",d=data_dir,s=cells,ch=chambers)
    run:
        cell_boundaries = defaultdict(list) # key: (cell,chamber#,chrom)  value: (start, end)

        for chrom in chroms:
            infile = "{}/PGP1_ALL/fragmat/cov1_strict/{}".format(data_dir,chrom)
            with open(infile,'r') as inf:
                for line in inf:
                    if len(line) < 3:
                        continue
                    el = line.strip().split()
                    ID = el[1]
                    el2 = ID.split(':')
                    boundstr = el2[-1]
                    (start, end) = [int(x) for x in boundstr.split('-')]
                    cell = el2[2]
                    chamber = el2[3]
                    chambernum = int(chamber[2:])
                    cell_boundaries[(cell,chambernum,chrom)].append((start, end))

        for cell in cells:
            for chamber in chambers:
                outfile = '{}/{}/fixed_beds/ch{}.bed'.format(data_dir,cell,chamber)
                with open(outfile,'w') as of:
                    for chrom in chroms:
                        cell_boundaries[(cell,chamber,chrom)].sort()
                        for start, end in cell_boundaries[(cell,chamber,chrom)]:
                            line = '{}\t{}\t{}\t{}\tchamber{}'.format(chrom, start, end, end-start, chamber)
                            print(line,file=of)
'''
# COMBINE FRAGMENT MATRICES
rule fix_fragmat:
    params:
        job_name  = "fix_fragmat.{x}",
    input:
        var_vcfs = expand("{v}/{CHR}.vcf",v=variant_vcf_dir,CHR=chroms),
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
    input:  ccf      = 'sissor_project/data/PGP1_ALL.het.whole_genome.ccf',
            vcf      = 'sissor_project/data/PGP1_ALL.het.whole_genome.vcf',
            bounds   = expand('base_calling/fragment_boundary_beds/{ce}/ch{ch}.bed',ch=chambers,ce=cells),
    output: fragmat  = expand('sissor_project/data/PGP1_ALL/augmented_fragmat/{c}',c=chroms),
            vcfs     = expand('sissor_project/data/PGP1_VCFs/{c}.vcf',c=chroms)
    run:
        odir = 'sissor_project/data/PGP1_ALL/augmented_fragmat'

        # generate fragment matrix from sissorhands base calls
        create_hapcut_fragment_matrix_CCF(chamber_call_file=input.ccf, fragment_boundary_files=input.bounds, output_dir=odir)

        #split vcf file into separate chromosome VCF files
        for chrom in chroms:
            shell('grep -P "^{chrom}\t" {input.vcf} > sissor_project/data/PGP1_VCFs/{chrom}.vcf')

rule combine_filtered:
    params: job_name = 'combine_filtered'
    input:  ccf = expand('sissor_project/data/ccf_split/{r}.ccf',r=regions),
            vcf = expand('sissor_project/data/vcf_split/{r}.vcf',r=regions)
    output: ccf = 'sissor_project/data/PGP1_ALL.het.whole_genome.ccf',
            vcf = 'sissor_project/data/PGP1_ALL.het.whole_genome.vcf',
    shell:
        '''cat {input.ccf} > {output.ccf}
           cat {input.vcf} > {output.vcf}'''

rule filter_CCF:
    params: job_name = 'filter_CCF_{chrom}.{start}.{end}'
    input:  ccf      = 'base_calling/output/split/{chrom}.{start}.{end}.out',
            vcf      = 'sissor_project/data/PGP1_VCFs_old/all.vcf'
    output: ccf      = 'sissor_project/data/ccf_split/{chrom}.{start}.{end}.ccf',
            vcf      = 'sissor_project/data/vcf_split/{chrom}.{start}.{end}.vcf',
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

# COMBINE FRAGMENT MATRICES
'''
rule merge_fragmat:
    params:
        job_name  = "merge_fragmat",
    input:
        P21      = expand("{dat}/PGP1_21/augmented_fragmat/{c}",dat=data_dir,c=chroms),
        P22      = expand("{dat}/PGP1_22/augmented_fragmat/{c}",dat=data_dir,c=chroms),
        PA1      = expand("{dat}/PGP1_A1/augmented_fragmat/{c}",dat=data_dir,c=chroms),
        var_vcfs = expand("{v}/{CHR}.vcf",v=variant_vcf_dir,CHR=chroms)
    output:
        P_ALL = expand("{dat}/PGP1_ALL/augmented_fragmat/{c}",dat=data_dir,c=chroms)
    run:
        for i1, i2, i3, o in zip(input.P21, input.P22, input.PA1, output.P_ALL):
            shell('cat {i1} {i2} {i3} > {o}')

# CREATE AUGMENTED FRAGMENT MATRIX FILES WITH HETEROZYGOUS CALL LOCATIONS
# SAME FORMAT AS BEFORE, BUT NOW '2' represents a heterozygous call
rule create_augmented_fragmat:
    params:
        job_name  = "{s}.create_augmented_fragmat",
    input:
        ploidy1_vcfs = expand("{{data_dir}}/{{s}}/ploidy1/ch{ch}.vcf",ch=chambers),
        ploidy2_vcfs = expand("{{data_dir}}/{{s}}/ploidy2/ch{ch}.vcf",ch=chambers),
        beds         = expand("{{data_dir}}/{{s}}/beds/ch{ch}.bed",ch=chambers),
        var_vcfs     = expand("{v}/{CHR}.vcf",v=variant_vcf_dir,CHR=chroms)
    output:
        expand("{{data_dir}}/{{s}}/augmented_fragmat/{c}",c=chroms)
    run:
        output_dir = os.path.join(data_dir,wildcards.s,'augmented_fragmat')
        create_hapcut_fragment_matrices_freebayes(input.ploidy1_vcfs, input.ploidy2_vcfs, input.beds, wildcards.s, chambers_pad, input.var_vcfs, output_dir, hets_in_seq=True)
'''
# simlink data to make path naming scheme consistent between PGP1_21 and PGP1_22

rule simlinks:
    input:
        #expand("{DIR}/PGP1_21_ch{chpad}.freebayes.remap.depth.new.vcf",DIR=config['PGP1_21_ploidy1'],chpad=chambers_pad),
        #expand("{DIR}/PGP1_21_ch{chpad}.ploidy2.freebayes.depth.vcf",DIR=config['PGP1_21_ploidy2'],chpad=chambers_pad),
        expand("{DIR}/PGP1_21_FragmentBoundaryCh{ch}.bed",DIR=config['PGP1_21_beds'],ch=chambers),
        #expand("{DIR}/PGP1_22_ch{chpad}.freebayes.remap.depth.vcf",DIR=config['PGP1_22_ploidy1'],chpad=chambers_pad),
        #expand("{DIR}/PGP1_22_ch{chpad}.ploidy2.freebayes.vcf",DIR=config['PGP1_22_ploidy2'],chpad=chambers_pad),
        expand("{DIR}/PGP1_22_FragmentBoundaryCh{ch}.bed",DIR=config['PGP1_22_beds'],ch=chambers),
        #expand("{DIR}/PGP1_A1_ch{chpad}.freebayes.allpos.vcf",DIR=config['PGP1_A1_ploidy1'],chpad=chambers_pad),
        #expand("{DIR}/PGP1_A1_ch{chpad}.freebayes.ploidy2.allpos.vcf",DIR=config['PGP1_A1_ploidy2'],chpad=chambers_pad),
        expand("{DIR}/PGP1_A1_FragmentBoundaryCh{ch}.bed",DIR=config['PGP1_A1_beds'],ch=chambers),
    run:
        for ch,chpad in zip(chambers,chambers_pad):
            shell('''
            mkdir -p {data_dir}/PGP1_21/ploidy1
            mkdir -p {data_dir}/PGP1_21/ploidy2
            mkdir -p {data_dir}/PGP1_21/beds
            #ln -s {config[PGP1_21_ploidy1]}/PGP1_21_ch{chpad}.freebayes.remap.depth.new.vcf {data_dir}/PGP1_21/ploidy1/ch{ch}.vcf
            #ln -s {config[PGP1_21_ploidy2]}/PGP1_21_ch{chpad}.ploidy2.freebayes.depth.vcf {data_dir}/PGP1_21/ploidy2/ch{ch}.vcf
            ln -s {config[PGP1_21_beds]}/PGP1_21_FragmentBoundaryCh{ch}.bed {data_dir}/PGP1_21/beds/ch{ch}.bed
            mkdir -p {data_dir}/PGP1_22/ploidy1
            mkdir -p {data_dir}/PGP1_22/ploidy2
            mkdir -p {data_dir}/PGP1_22/beds
            #ln -s {config[PGP1_22_ploidy1]}/PGP1_22_ch{chpad}.freebayes.remap.depth.new.vcf {data_dir}/PGP1_22/ploidy1/ch{ch}.vcf
            #ln -s {config[PGP1_22_ploidy2]}/PGP1_22_ch{chpad}.ploidy2.freebayes.vcf {data_dir}/PGP1_22/ploidy2/ch{ch}.vcf
            ln -s {config[PGP1_22_beds]}/PGP1_22_FragmentBoundaryCh{ch}.bed {data_dir}/PGP1_22/beds/ch{ch}.bed
            mkdir -p {data_dir}/PGP1_A1/ploidy1
            mkdir -p {data_dir}/PGP1_A1/ploidy2
            mkdir -p {data_dir}/PGP1_A1/beds
            #ln -s {config[PGP1_A1_ploidy1]}/PGP1_A1_ch{chpad}.freebayes.allpos.vcf {data_dir}/PGP1_A1/ploidy1/ch{ch}.vcf
            #ln -s {config[PGP1_A1_ploidy2]}/PGP1_A1_ch{chpad}.freebayes.ploidy2.allpos.vcf {data_dir}/PGP1_A1/ploidy2/ch{ch}.vcf
            ln -s {config[PGP1_A1_beds]}/PGP1_A1_FragmentBoundaryCh{ch}.bed {data_dir}/PGP1_A1/beds/ch{ch}.bed
            ''')

rule simlinks_new:
    input:
        #expand("{DIR}/PGP1_21_ch{chpad}.freebayes.remap.depth.new.vcf",DIR=config['PGP1_21_ploidy1'],chpad=chambers_pad),
        #expand("{DIR}/PGP1_21_ch{chpad}.ploidy2.freebayes.depth.vcf",DIR=config['PGP1_21_ploidy2'],chpad=chambers_pad),
        expand("{DIR}/ch{ch}.bed",DIR=config['PGP1_21_beds'],ch=chambers_pad),
        #expand("{DIR}/PGP1_22_ch{chpad}.freebayes.remap.depth.vcf",DIR=config['PGP1_22_ploidy1'],chpad=chambers_pad),
        #expand("{DIR}/PGP1_22_ch{chpad}.ploidy2.freebayes.vcf",DIR=config['PGP1_22_ploidy2'],chpad=chambers_pad),
        expand("{DIR}/ch{ch}.bed",DIR=config['PGP1_22_beds'],ch=chambers_pad),
        #expand("{DIR}/PGP1_A1_ch{chpad}.freebayes.allpos.vcf",DIR=config['PGP1_A1_ploidy1'],chpad=chambers_pad),
        #expand("{DIR}/PGP1_A1_ch{chpad}.freebayes.ploidy2.allpos.vcf",DIR=config['PGP1_A1_ploidy2'],chpad=chambers_pad),
        expand("{DIR}/ch{ch}.bed",DIR=config['PGP1_A1_beds'],ch=chambers_pad),
    run:
        for ch,chpad in zip(chambers,chambers_pad):
            shell('''
            mkdir -p {data_dir}/PGP1_21/ploidy1
            mkdir -p {data_dir}/PGP1_21/ploidy2
            mkdir -p {data_dir}/PGP1_21/beds
            #ln -s {config[PGP1_21_ploidy1]}/PGP1_21_ch{chpad}.freebayes.remap.depth.new.vcf {data_dir}/PGP1_21/ploidy1/ch{ch}.vcf
            #ln -s {config[PGP1_21_ploidy2]}/PGP1_21_ch{chpad}.ploidy2.freebayes.depth.vcf {data_dir}/PGP1_21/ploidy2/ch{ch}.vcf
            ln -s {config[PGP1_21_beds]}/ch{ch}.bed {data_dir}/PGP1_21/beds/ch{ch}.bed
            mkdir -p {data_dir}/PGP1_22/ploidy1
            mkdir -p {data_dir}/PGP1_22/ploidy2
            mkdir -p {data_dir}/PGP1_22/beds
            #ln -s {config[PGP1_22_ploidy1]}/PGP1_22_ch{chpad}.freebayes.remap.depth.new.vcf {data_dir}/PGP1_22/ploidy1/ch{ch}.vcf
            #ln -s {config[PGP1_22_ploidy2]}/PGP1_22_ch{chpad}.ploidy2.freebayes.vcf {data_dir}/PGP1_22/ploidy2/ch{ch}.vcf
            ln -s {config[PGP1_22_beds]}/ch{ch}.bed {data_dir}/PGP1_22/beds/ch{ch}.bed
            mkdir -p {data_dir}/PGP1_A1/ploidy1
            mkdir -p {data_dir}/PGP1_A1/ploidy2
            mkdir -p {data_dir}/PGP1_A1/beds
            #ln -s {config[PGP1_A1_ploidy1]}/PGP1_A1_ch{chpad}.freebayes.allpos.vcf {data_dir}/PGP1_A1/ploidy1/ch{ch}.vcf
            #ln -s {config[PGP1_A1_ploidy2]}/PGP1_A1_ch{chpad}.freebayes.ploidy2.allpos.vcf {data_dir}/PGP1_A1/ploidy2/ch{ch}.vcf
            ln -s {config[PGP1_A1_beds]}/ch{ch}.bed {data_dir}/PGP1_A1/beds/ch{ch}.bed
            ''')

rule clean:
    shell:
        '''
        rm -rf sissor_project/plots/bak sissor_project/error_rates/bak sissor_project/experiments/bak sissor_project/data/bak || true 2>/dev/null
        mkdir -p sissor_project/plots/bak sissor_project/error_rates/bak sissor_project/experiments/bak sissor_project/data/bak || true 2>/dev/null
        mv sissor_project/plots/*.png sissor_project/plots/bak || true 2>/dev/null
        mv sissor_project/error_rates/*.p sissor_project/error_rates/bak || true 2>/dev/null
        mv sissor_project/experiments/hapcut2* sissor_project/experiments/bak || true 2>/dev/null
        mv sissor_project/data/PGP1_21 sissor_project/data/bak || true 2>/dev/null
        mv sissor_project/data/PGP1_22 sissor_project/data/bak || true 2>/dev/null
        mv sissor_project/data/PGP1_A1 sissor_project/data/bak || true 2>/dev/null
        mv sissor_project/data/PGP1_ALL sissor_project/data/bak || true 2>/dev/null
        mv sissor_project/data/PGP1_VCFs sissor_project/data/bak || true 2>/dev/null
        '''

rule dummy:
    run:
        pass
