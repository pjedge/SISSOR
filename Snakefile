# author: Peter Edge
# 3/29/2016
# email: pedge@eng.ucsd.edu

# this is a snakemake Snakefile, written using snakemake 3.5.5
# configure paths and values in cluster.yaml and config.yaml to your system
# example execution for a TORQUE cluster:
# snakemake -j 200 --cluster-config cluster.yaml --cluster "qsub -A {cluster.group} -V -q {cluster.queue} -o {cluster.qsub_stdout_dir}/{params.job_name}.o -e {cluster.qsub_stdout_dir}/{params.job_name}.e -N {params.job_name} -l nodes=1:ppn={cluster.ppn} -l walltime={cluster.walltime} -M {cluster.email} -m e -d {cluster.working_dir}" --local-cores 1

configfile: "config.yaml"
localrules: all, plot_hapcut2_results, simlink_PGP1_22, simlink_PGP1_21, simlink_PGP1_21_freebayes, simlink_PGP1_22_freebayes

import sys
sys.path.append(config['haptools_dir'])
import run_tools
import error_rates
import fileIO
import plot_data
import os
from os.path import join
from create_hapcut_fragment_matrix import create_hapcut_fragment_matrices_freebayes
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
samples = ['PGP1_21_22','PGP1_21','PGP1_22']

runs = ['PGP1_21_new','PGP1_22_new']
rule all:
    input:
        expand("{P}/sissor_hapcut2.png",P=config['plots_dir']),
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

fix_status=['normal','fixed']
fix_labels=['Normal','Fragment Split Pre-processing']
rule calculate_error_rates:
    params:
        job_name = "hapcut2_error_rates"
    input:
        hapblocks = expand("{E}/hapcut2_PGP1_21_22_{fx}/{p}/{c}.output",E=experiments_dir,fx=fix_status,p=prune_labels,c=chroms),
        runtimes  = expand("{E}/hapcut2_PGP1_21_22_{fx}/{p}/{c}.runtime",E=experiments_dir,fx=fix_status,p=prune_labels,c=chroms),
        fragmats  = expand("{d}/PGP1_21_22/fragmat_{fx}/{c}",d=data_dir,fx=fix_status,s=samples,c=chroms),
        var_vcfs  = expand("{v}/{c}.vcf", v=variant_vcf_dir,c=chroms),
        bac_haps  = expand("{bac}/{c}",bac=config['BAC_hapblocks'], c=chroms)
    output:
        stats_file  = "{}/hapcut2.stats.p".format(config['error_rates_dir']),
        labels_file = "{}/hapcut2.labels.p".format(config['error_rates_dir'])
    run:
        # list of lists of error results,
        data   = [] # index of the list is a 'condition' each inner list has 23 error results (each chrom).
        labels = fix_labels

        for fx in fix_status: # sample
            datalist = []
            for c in chroms: # chromosome
                assembly_file = "{}/hapcut2_PGP1_21_22_{}/{}/{}.output".format(config['experiments_dir'],fx,p,c)
                runtime_file  = "{}/hapcut2_PGP1_21_22_{}/{}/{}.runtime".format(config['experiments_dir'],fx,p,c)
                frag_file     = "{}/PGP1_21_22/fragmat_{}/{}".format(data_dir,fx,c)
                vcf_file      = "{}/{}.vcf".format(variant_vcf_dir,c)
                truth_file    = "{}/{}".format(config['BAC_hapblocks'], c)
                err = error_rates.hapblock_hapblock_error_rate(truth_file, assembly_file, frag_file, vcf_file, runtime_file)
                datalist.append(err)

            print("{} results over all chromosomes:".format(p))
            print(sum(datalist,error_rates.error_result(None, 0, 0, 0, 0, 0, 0, 0, None,[],[],[])))

        data.append(datalist)

        pickle.dump(data,open(output.stats_file,"wb"))
        pickle.dump(labels,open(output.labels_file,"wb"))
'''
# PRUNE HAPCUT2 RESULTS
rule prune_hapcut2:
    params:
        job_name = "PGP1_21_22_prune_hapcut",
    input:
        hapblocks = expand("{{E}}/hapcut2_PGP1_21_22/{c}.output",c=chroms),
        runtime = expand("{{E}}/hapcut2_PGP1_21_22/{c}.runtime",c=chroms)
    output:
        hapblocks = expand("{{E}}/pruned_hapcut2_PGP1_21_22/{p}/{c}.output",p=prune_labels,c=chroms),
        runtime = expand("{{E}}/pruned_hapcut2_PGP1_21_22/{p}/{c}.runtime",p=prune_labels,c=chroms)
    run:
        # run hapcut
        for pl, (sq,rh,sb) in zip(prune_labels, prune_tuples):

            for c in chroms:
                old_runtime_file = "{}/hapcut2_PGP1_21_22/{}.runtime".format(config['experiments_dir'],c)
                runtime_file = "{}/pruned_hapcut2_PGP1_21_22/{}/{}.runtime".format(config['experiments_dir'],pl,c)
                shell("cp {old_runtime_file} {runtime_file}")
                input_hap = "{}/hapcut2_PGP1_21_22/{}.output".format(config['experiments_dir'],c)
                output_hap = "{}/pruned_hapcut2_PGP1_21_22/{}/{}.output".format(config['experiments_dir'],pl,c)
                fileIO.prune_hapblock_file(input_hap, output_hap, sq, sb, rh)
'''
# RUN HAPCUT2
rule run_hapcut2:
    params:
        job_name = "{s}.{c}_hapcut",
    input:
        frag_file = "%s/{s}/fragmat/{fx}/{c}" % data_dir,
        vcf_file  = lambda wildcards: expand("{v}/{c}.vcf", v=variant_vcf_dir,c=wildcards.c)
    output:
        hapblocks = "{E}/hapcut2_{s}_{fx}/{c}.output",
        runtime = "{E}/hapcut2_{s}_{fx}/{c}.runtime"
    run:
        # run hapcut
        runtime = run_tools.run_hapcut2(config['hapcut2'], input.frag_file, input.vcf_file, output.hapblocks, 1000, 100, 1, 0, 5, 0.8)
        with open(output.runtime,'w') as rf:
            print(runtime, file=rf)

# CREATE FRAGMENT MATRIX FILES FOR EXPERIMENT
rule fix_fragmat:
    params:
        job_name  = "{s}_fix_fragmat",
    input:
        samples = expand("{{data_dir}}/{{s}}/fragmat/normal/{c}",c=chroms)
    output:
        fixed = expand("{{data_dir}}/{{s}}/fragmat/fixed/{c}",c=chroms)
    run:
        for i, o in zip(input.samples, output.fixed):

# COMBINE FRAGMENT MATRICES
rule merge_fragmat:
    params:
        job_name  = "merge_fragmat",
    input:
        P21 = expand("{{data_dir}}/PGP1_21/fragmat/normal/{c}",c=chroms),
        P22 = expand("{{data_dir}}/PGP1_22/fragmat/normal/{c}",c=chroms)
    output:
        P21_22 = expand("{{data_dir}}/PGP1_21_22/fragmat/normal/{c}",c=chroms)
    run:
        for i1, i2, o in zip(input.P21, input.P22, output.P21_22):
            shell('cat {i1} {i2} > {o}')

# CREATE FRAGMENT MATRIX FILES FOR EXPERIMENT
rule create_fragmat:
    params:
        job_name  = "{s}_create_fragmat",
    input:
        vcfs     = expand("{{data_dir}}/{{s}}/freebayes/ch{ch}.vcf",ch=chambers),
        beds     = expand("{{data_dir}}/{{s}}/beds/ch{ch}.bed",ch=chambers),
        var_vcfs = expand("{{data_dir}}/PGP1_VCFs/{CHR}.vcf",CHR=chroms)
    output:
        expand("{{data_dir}}/{{s}}/fragmat/normal/{c}",c=chroms)
    run:
        output_dir = os.path.join(data_dir,wildcards.s,'fragmat')
        create_hapcut_fragment_matrices_freebayes(input.vcfs, input.beds, input.var_vcfs, output_dir)

# simlink data to make path naming scheme consistent between PGP1_21 and PGP1_22
rule simlink_PGP1_21_freebayes:
    input:
        vcf = expand("{DIR}/PGP1_21_ch{chpad}.freebayes.remap.depth.vcf",DIR=config['PGP1_21_freebayes'],chpad=chambers_pad),
    run:
        for ch,chpad in zip(chambers,chambers_pad):
            shell('''
            mkdir -p {data_dir}/PGP1_21/freebayes
            ln -s {config[PGP1_21_freebayes]}/PGP1_21_ch{chpad}.freebayes.remap.depth.vcf {data_dir}/PGP1_21/freebayes/ch{ch}.vcf''')

# simlink data to make path naming scheme consistent between PGP1_21 and PGP1_22
rule simlink_PGP1_22_freebayes:
    input:
        vcf = expand("{DIR}/PGP1_22_ch{chpad}.freebayes.remap.depth.vcf",DIR=config['PGP1_22_freebayes'],chpad=chambers_pad),
    run:
        for ch,chpad in zip(chambers,chambers_pad):
            shell('''
            mkdir -p {data_dir}/PGP1_22/freebayes
            ln -s {config[PGP1_22_freebayes]}/PGP1_22_ch{chpad}.freebayes.remap.depth.vcf {data_dir}/PGP1_22/freebayes/ch{ch}.vcf''')

# simlink data to make path naming scheme consistent between PGP1_21 and PGP1_22
rule simlink_PGP1_21:
    input:
        vcf    = expand("{DIR}/PGP1_21_ch{chpad}.sorted.fragment.depth5.vcf",DIR=config['PGP1_21_vcfs_old'],chpad=chambers_pad),
        pileup = expand("{DIR}/PGP1_21_ch{chpad}.sorted.fragment.depth.pileup",DIR=config['PGP1_21_pileups_old'],chpad=chambers_pad),
        bed    = expand("{DIR}/PGP1_21_FragmentBoundaryCh{ch}.bed",DIR=config['PGP1_21_boundaries_old'],ch=chambers),
    run:
        for ch,chpad in zip(chambers,chambers_pad):
            shell('''
            mkdir -p {data_dir}/PGP1_21/VCFs
            mkdir -p {data_dir}/PGP1_21/pileups
            mkdir -p {data_dir}/PGP1_21/beds
            ln -s {config[PGP1_21_vcfs_old]}/PGP1_21_ch{chpad}.sorted.fragment.depth5.vcf {data_dir}/PGP1_21/VCFs/ch{ch}.vcf
            ln -s {config[PGP1_21_pileups_old]}/PGP1_21_ch{chpad}.sorted.fragment.depth.pileup {data_dir}/PGP1_21/pileups/ch{ch}.pileup
            ln -s {config[PGP1_21_boundaries_old]}/PGP1_21_FragmentBoundaryCh{ch}.bed {data_dir}/PGP1_21/beds/ch{ch}.bed''')

# simlink data to make path naming scheme consistent between PGP1_21 and PGP1_22
rule simlink_PGP1_22:
    input:
        vcf    = expand("{DIR}/PGP1_22_ch{chpad}.sorted.fragment.vcf",DIR=config['PGP1_22_old'],chpad=chambers_pad),
        pileup = expand("{DIR}/PGP1_22_ch{chpad}.sorted.fragment.depth.pileup",DIR=config['PGP1_22_old'],chpad=chambers_pad),
        bed    = expand("{DIR}/PGP1_22_FragmentBoundaryCh{ch}.bed",DIR=config['PGP1_22_old'],ch=chambers),
    run:
        for ch,chpad in zip(chambers,chambers_pad):
            shell('''
            mkdir -p {data_dir}/PGP1_22/VCFs
            mkdir -p {data_dir}/PGP1_22/pileups
            mkdir -p {data_dir}/PGP1_22/beds
            ln -s {config[PGP1_22_old]}/PGP1_22_ch{chpad}.sorted.fragment.vcf {data_dir}/PGP1_22/VCFs/ch{ch}.vcf
            ln -s {config[PGP1_22_old]}/PGP1_22_ch{chpad}.sorted.fragment.depth.pileup {data_dir}/PGP1_22/pileups/ch{ch}.pileup
            ln -s {config[PGP1_22_old]}/PGP1_22_FragmentBoundaryCh{ch}.bed {data_dir}/PGP1_22/beds/ch{ch}.bed''')

# NEW MINI-RUNS
# simlink data to make path naming scheme consistent between PGP1_21 and PGP1_22
rule simlink_PGP1_21_new:
    input:
        vcf    = expand("{DIR}/PGP1_21_ch{chpad}.freebayes.remap.depth.new.vcf",DIR=config['PGP1_21_vcfs_new'],chpad=chambers_pad),
        pileup = expand("{DIR}/PGP1_21_ch{chpad}.sorted.fragment.depth.pileup",DIR=config['PGP1_21_pileups_old'],chpad=chambers_pad),
        bed    = expand("{DIR}/PGP1_21_FragmentBoundaryCh{ch}.bed",DIR=config['PGP1_21_boundaries_old'],ch=chambers),
    run:
        for ch,chpad in zip(chambers,chambers_pad):
            shell('''
            mkdir -p {data_dir}/PGP1_21_new/freebayes
            mkdir -p {data_dir}/PGP1_21_new/pileups
            mkdir -p {data_dir}/PGP1_21_new/beds
            ln -s {config[PGP1_21_vcfs_new]}/PGP1_21_ch{chpad}.freebayes.remap.depth.new.vcf {data_dir}/PGP1_21_new/freebayes/ch{ch}.vcf
            ln -s {config[PGP1_21_pileups_old]}/PGP1_21_ch{chpad}.sorted.fragment.depth.pileup {data_dir}/PGP1_21_new/pileups/ch{ch}.pileup
            ln -s {config[PGP1_21_boundaries_old]}/PGP1_21_FragmentBoundaryCh{ch}.bed {data_dir}/PGP1_21_new/beds/ch{ch}.bed''')

# simlink data to make path naming scheme consistent between PGP1_21 and PGP1_22
rule simlink_PGP1_22_new:
    input:
        vcf    = expand("{DIR}/PGP1_22_ch{chpad}.freebayes.remap.depth.new.vcf",DIR=config['PGP1_22_vcfs_new'],chpad=chambers_pad),
        pileup = expand("{DIR}/PGP1_22_ch{chpad}.sorted.fragment.depth.pileup",DIR=config['PGP1_22_old'],chpad=chambers_pad),
        bed    = expand("{DIR}/PGP1_22_FragmentBoundaryCh{ch}.bed",DIR=config['PGP1_22_old'],ch=chambers),
    run:
        for ch,chpad in zip(chambers,chambers_pad):
            shell('''
            mkdir -p {data_dir}/PGP1_22_new/freebayes
            mkdir -p {data_dir}/PGP1_22_new/pileups
            mkdir -p {data_dir}/PGP1_22_new/beds
            ln -s {config[PGP1_22_vcfs_new]}/PGP1_22_ch{chpad}.freebayes.remap.depth.new.vcf {data_dir}/PGP1_22_new/freebayes/ch{ch}.vcf
            ln -s {config[PGP1_22_old]}/PGP1_22_ch{chpad}.sorted.fragment.depth.pileup {data_dir}/PGP1_22_new/pileups/ch{ch}.pileup
            ln -s {config[PGP1_22_old]}/PGP1_22_FragmentBoundaryCh{ch}.bed {data_dir}/PGP1_22_new/beds/ch{ch}.bed''')

rule dummy:
    run:
        pass
