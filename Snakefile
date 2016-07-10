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
#sys.path.append(config['haptools_dir'])
import run_tools
import error_rates
import fileIO
import plot_data
import os
from os.path import join
from create_hapcut_fragment_matrix import create_hapcut_fragment_matrices_freebayes
from fix_chamber_contamination import fix_chamber_contamination
import pickle
# qsub stderr and stdout directories are not automatically created by snakemake!
if not os.path.exists(config['qsub_stdout_dir']):
    os.makedirs(config['qsub_stdout_dir'])


# specify chambers
data_dir = config["data_dir"]
experiments_dir = config["experiments_dir"]
plots_dir = config["plots_dir"]
variant_vcf_dir = join(data_dir,'PGP1_VCFs')
chambers = list(range(1,25))
chambers_pad = ['{0:02d}'.format(c) for c in chambers]
chroms  = ['chr{}'.format(i) for i in range(1,23)]
samples = ['PGP1_ALL','PGP1_21','PGP1_22','PGP1_A1']

rule all:
    input:
        #expand("{P}/sissor_hapcut2.png",P=config['plots_dir']),
        expand("{d}/{s}/augmented_fragmat/cov1/normal/{c}",d=data_dir, s=samples[1:], c=chroms)
        # "new" snipped out runs
        #hapblocks = expand("{E}/hapcut2_{s}/{c}.output",E=experiments_dir,s=runs,c=chroms),
        #runtime = expand("{E}/hapcut2_{s}/{c}.runtime",E=experiments_dir,s=runs,c=chroms)

# PLOT RESULTS
rule plot_hapcut2_results:
    params:
        job_name = "plot_hapcut2_sissor"
    input:
        stats_file  = "{}/hapcut2.stats.p".format(config['error_rates_dir']),
        labels_file = "{}/hapcut2.labels.p".format(config['error_rates_dir'])
    output:
        plot = "%s/sissor_hapcut2.png" % config['plots_dir']
    run:
        data = pickle.load(open(input.stats_file,"rb"))
        labels = pickle.load(open(input.labels_file,"rb"))
        plot_data.plot_experiment_sissor(data,labels,[],output.plot)

exp =['cov1','cov2','cov3']
exp_labels=['All fragments','Cov >= 2','Cov >= 3',]
rule calculate_error_rates:
    params:
        job_name = "hapcut2_error_rates"
    input:
        hapblocks = expand("{E}/hapcut2_PGP1_ALL/{exp}/{c}.output",E=experiments_dir,exp=exp,c=chroms),
        runtimes  = expand("{E}/hapcut2_PGP1_ALL/{exp}/{c}.runtime",E=experiments_dir,exp=exp,c=chroms),
        fragmats  = expand("{d}/PGP1_ALL/fragmat/{exp}/fixed/{c}",d=data_dir,exp=exp,c=chroms),
        var_vcfs  = expand("{v}/{c}.vcf", v=variant_vcf_dir,c=chroms),
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
                frag_file     = "{}/PGP1_ALL/fragmat/{}/fixed/{}".format(data_dir,x,c)
                vcf_file      = "{}/{}.vcf".format(variant_vcf_dir,c)
                truth_file    = "{}/{}.filtered".format(config['BAC_hapblocks'], c)
                err = error_rates.hapblock_hapblock_error_rate(truth_file, assembly_file, frag_file, vcf_file, runtime_file)
                datalist.append(err)

            print("{} results over all chromosomes:".format(x))
            print(sum(datalist,error_rates.error_result(None, 0, 0, 0, 0, 0, 0, 0, None,[],[],[])))

            data.append(datalist)

        pickle.dump(data,open(output.stats_file,"wb"))
        pickle.dump(labels,open(output.labels_file,"wb"))

# RUN HAPCUT2
rule run_hapcut2:
    params:
        job_name = "{s}.{c}.{x}.hapcut",
    input:
        frag_file = "%s/{s}/fragmat/{x}/fixed/{c}" % data_dir,
        vcf_file  = lambda wildcards: expand("{v}/{c}.vcf", v=variant_vcf_dir,c=wildcards.c)
    output:
        hapblocks = "{E}/hapcut2_{s}/{x}/{c}.output",
        runtime = "{E}/hapcut2_{s}/{x}/{c}.runtime"
    run:
        # run hapcut
        runtime = run_tools.run_hapcut2(config['hapcut2'], input.frag_file, input.vcf_file, output.hapblocks, 1000, 25, 1, 0, 5, 0.95,'')
        with open(output.runtime,'w') as rf:
            print(runtime, file=rf)

# COMBINE FRAGMENT MATRICES
rule merge_and_fix_fragmat:
    params:
        job_name  = "merge_and_fix_fragmat.{x}",
    input:
        P21      = expand("{{data_dir}}/PGP1_21/fragmat/cov1/normal/{c}",c=chroms),
        P22      = expand("{{data_dir}}/PGP1_22/fragmat/cov1/normal/{c}",c=chroms),
        PA1      = expand("{{data_dir}}/PGP1_A1/fragmat/cov1/normal/{c}",c=chroms),
        var_vcfs = expand("{{data_dir}}/PGP1_VCFs/{CHR}.vcf",CHR=chroms)
    output:
        fixed = expand("{{data_dir}}/PGP1_ALL/fragmat/{{x}}/fixed/{c}",c=chroms)
    run:
        P_ALL = expand("{dat}/PGP1_ALL/fragmat/{exp}/normal/{c}",dat=data_dir,exp=wildcards.x,c=chroms)
        merged_dir = "{}/PGP1_ALL/fragmat/{}/normal".format(data_dir,wildcards.x)
        shell('mkdir -p {merged_dir}')
        for i1, i2, i3, o in zip(input.P21, input.P22, input.PA1, P_ALL):
            shell('cat {i1} {i2} {i3} > {o}')

        for i, o, v in zip(P_ALL, output.fixed, input.var_vcfs):
            if wildcards.x == 'cov3':
                fix_chamber_contamination(i,o,v,2,3)
            elif wildcards.x == 'cov2':
                fix_chamber_contamination(i,o,v,2,2)
            else:
                fix_chamber_contamination(i,o,v,2,0)

# CREATE FRAGMENT MATRIX FILES FOR EXPERIMENT
rule create_fragmat:
    params:
        job_name  = "{s}.create_fragmat",
    input:
        ploidy1_vcfs         = expand("{{data_dir}}/{{s}}/ploidy1/ch{ch}.vcf",ch=chambers),
        ploidy2_vcfs = expand("{{data_dir}}/{{s}}/ploidy2/ch{ch}.vcf",ch=chambers),
        beds         = expand("{{data_dir}}/{{s}}/beds/ch{ch}.bed",ch=chambers),
        var_vcfs     = expand("{{data_dir}}/PGP1_VCFs/{CHR}.vcf",CHR=chroms)
    output:
        expand("{{data_dir}}/{{s}}/fragmat/cov1/normal/{c}",c=chroms)
    run:
        output_dir = os.path.join(data_dir,wildcards.s,'fragmat','cov1','normal')
        #if wildcards.x == 'old':
        create_hapcut_fragment_matrices_freebayes(input.ploidy1_vcfs, input.ploidy2_vcfs, input.beds, chambers_pad, input.var_vcfs, output_dir)
        #else:
        #    create_hapcut_fragment_matrices_freebayes(input.ploidy1_vcfs, input.beds, chambers_pad, input.var_vcfs, output_dir)


# CREATE AUGMENTED FRAGMENT MATRIX FILES WITH HETEROZYGOUS CALL LOCATIONS
# SAME FORMAT AS BEFORE, BUT NOW '2' represents a heterozygous call
rule create_augmented_fragmat:
    params:
        job_name  = "{s}.create_augmented_fragmat",
    input:
        ploidy1_vcfs = expand("{{data_dir}}/{{s}}/ploidy1/ch{ch}.vcf",ch=chambers),
        ploidy2_vcfs = expand("{{data_dir}}/{{s}}/ploidy2/ch{ch}.vcf",ch=chambers),
        beds         = expand("{{data_dir}}/{{s}}/beds/ch{ch}.bed",ch=chambers),
        var_vcfs     = expand("{{data_dir}}/PGP1_VCFs/{CHR}.vcf",CHR=chroms)
    output:
        expand("{{data_dir}}/{{s}}/augmented_fragmat/cov1/normal/{c}",c=chroms)
    run:
        output_dir = os.path.join(data_dir,wildcards.s,'augmented_fragmat','cov1','normal')
        create_hapcut_fragment_matrices_freebayes(input.ploidy1_vcfs, input.ploidy2_vcfs, input.beds, chambers_pad, input.var_vcfs, output_dir, hets_in_seq=True)

# simlink data to make path naming scheme consistent between PGP1_21 and PGP1_22
rule simlinks:
    input:
        expand("{DIR}/PGP1_21_ch{chpad}.freebayes.remap.depth.new.vcf",DIR=config['PGP1_21_ploidy1'],chpad=chambers_pad),
        expand("{DIR}/PGP1_21_ch{chpad}.ploidy2.freebayes.depth.vcf",DIR=config['PGP1_21_ploidy2'],chpad=chambers_pad),
        expand("{DIR}/PGP1_21_FragmentBoundaryCh{ch}.bed",DIR=config['PGP1_21_beds'],ch=chambers),
        expand("{DIR}/PGP1_22_ch{chpad}.freebayes.remap.depth.vcf",DIR=config['PGP1_22_ploidy1'],chpad=chambers_pad),
        expand("{DIR}/PGP1_22_ch{chpad}.ploidy2.freebayes.vcf",DIR=config['PGP1_22_ploidy2'],chpad=chambers_pad),
        expand("{DIR}/PGP1_22_FragmentBoundaryCh{ch}.bed",DIR=config['PGP1_22_beds'],ch=chambers),
        expand("{DIR}/PGP1_A1_ch{chpad}.freebayes.allpos.vcf",DIR=config['PGP1_A1_ploidy1'],chpad=chambers_pad),
        expand("{DIR}/PGP1_A1_ch{chpad}.freebayes.ploidy2.allpos.vcf",DIR=config['PGP1_A1_ploidy2'],chpad=chambers_pad),
        expand("{DIR}/PGP1_A1_FragmentBoundaryCh{ch}.bed",DIR=config['PGP1_A1_beds'],ch=chambers),
    run:
        for ch,chpad in zip(chambers,chambers_pad):
            shell('''
            mkdir -p {data_dir}/PGP1_21/ploidy1
            mkdir -p {data_dir}/PGP1_21/ploidy2
            mkdir -p {data_dir}/PGP1_21/beds
            ln -s {config[PGP1_21_ploidy1]}/PGP1_21_ch{chpad}.freebayes.remap.depth.new.vcf {data_dir}/PGP1_21/ploidy1/ch{ch}.vcf
            ln -s {config[PGP1_21_ploidy2]}/PGP1_21_ch{chpad}.ploidy2.freebayes.depth.vcf {data_dir}/PGP1_21/ploidy2/ch{ch}.vcf
            ln -s {config[PGP1_21_beds]}/PGP1_21_FragmentBoundaryCh{ch}.bed {data_dir}/PGP1_21/beds/ch{ch}.bed

            mkdir -p {data_dir}/PGP1_22/ploidy1
            mkdir -p {data_dir}/PGP1_22/ploidy2
            mkdir -p {data_dir}/PGP1_22/beds
            ln -s {config[PGP1_22_ploidy1]}/PGP1_22_ch{chpad}.freebayes.remap.depth.new.vcf {data_dir}/PGP1_22/ploidy1/ch{ch}.vcf
            ln -s {config[PGP1_22_ploidy2]}/PGP1_22_ch{chpad}.ploidy2.freebayes.vcf {data_dir}/PGP1_22/ploidy2/ch{ch}.vcf
            ln -s {config[PGP1_22_beds]}/PGP1_22_FragmentBoundaryCh{ch}.bed {data_dir}/PGP1_22/beds/ch{ch}.bed

            mkdir -p {data_dir}/PGP1_A1/ploidy1
            mkdir -p {data_dir}/PGP1_A1/ploidy2
            mkdir -p {data_dir}/PGP1_A1/beds
            ln -s {config[PGP1_A1_ploidy1]}/PGP1_A1_ch{chpad}.freebayes.allpos.vcf {data_dir}/PGP1_A1/ploidy1/ch{ch}.vcf
            ln -s {config[PGP1_A1_ploidy2]}/PGP1_A1_ch{chpad}.freebayes.ploidy2.allpos.vcf {data_dir}/PGP1_A1/ploidy2/ch{ch}.vcf
            ln -s {config[PGP1_A1_beds]}/PGP1_A1_FragmentBoundaryCh{ch}.bed {data_dir}/PGP1_A1/beds/ch{ch}.bed
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
        mv sissor_project/data/PGP1_21_22 sissor_project/data/bak || true 2>/dev/null
        mv sissor_project/data/PGP1_ALL sissor_project/data/bak || true 2>/dev/null
        '''

rule dummy:
    run:
        pass
