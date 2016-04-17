# author: Peter Edge
# 3/29/2016
# email: pedge@eng.ucsd.edu

# this is a snakemake Snakefile, written using snakemake 3.5.5
# configure paths and values in cluster.yaml and config.yaml to your system
# example execution for a TORQUE cluster:
# snakemake -j 200 --cluster-config cluster.yaml --cluster "qsub -A {cluster.group} -V -q {cluster.queue} -o {cluster.qsub_stdout_dir}/{params.job_name}.o -e {cluster.qsub_stdout_dir}/{params.job_name}.e -N {params.job_name} -l nodes=1:ppn={cluster.ppn} -l walltime={cluster.walltime} -M {cluster.email} -m e -d {cluster.working_dir}" --local-cores 1

configfile: "config.yaml"
localrules: all, plot_hapcut2_results, simlink_PGP1_22, simlink_PGP1_21

import sys
import run_tools
import plot_data
sys.path.append(config['haptools_dir'])

from create_hapcut_fragment_matrix import create_hapcut_fragment_matrix


# qsub stderr and stdout directories are not automatically created by snakemake!
if not os.path.exists(config['qsub_stdout_dir']):
    os.makedirs(config['qsub_stdout_dir'])

# specify chambers
data_dir = config["data_dir"]
experiments_dir = config["experiments_dir"]
plots_dir = config["plots_dir"]
variant_vcf_dir = config["variant_vcf_dir"]
chambers = list(range(1,25))
chambers_pad = ['{0:02d}'.format(c) for c in chambers]
chroms   = list(range(1,23))
samples = ['PGP1_21_22','PGP1_21','PGP1_22']

rule all:
    input:
        expand("{P}/sissor_hapcut2.png",P=config['plots_dir'])

# PLOT RESULTS
rule plot_hapcut2_results:
    params:
        job_name = "plot_hapcut2_sissor"
    input:
        stats_file  = "{}/hapcut2.stats.p".format(config['error_rates_dir']),
        labels_file = "{}/hapcut2.labels.p".format(config['error_rates_dir'])
    output:
        plot = expand("{P}/sissor_hapcut2.png",P=config['plots_dir'])
    run:
        data = pickle.load(open(input.stats_file,"rb"))
        labels = pickle.load(open(input.labels_file,"rb"))
        plot_data.plot_experiment_no_emphasis(data,labels,[],output.plot)

# ERROR RATES
rule calculate_error_rates:
    params:
        job_name = "hapcut2_error_rates"
    input:
        hapblocks = expand("{E}/hapcut_{s}/{c}.output",E=experiments_dir,s=samples,c=chroms),
        runtimes = expand("{E}/hapcut_{s}/{c}.runtime",E=experiments_dir,s=samples,c=chroms),
        expand("{{data_dir}}/{s}/fragmat/chr{c}",c=chroms),
        expand("{v}/{c}.vcf", v=config['variant_vcf_dir'],c=chroms),
        expand("{bac}/haplotype.chr{c}.bac.txt".format(bac=config['BAC_hapblocks'], c=chroms)
    output:
        stats_file  = "{}/hapcut2.stats.p".format(config['error_rates_dir']),
        labels_file = "{}/hapcut2.labels.p".format(config['error_rates_dir'])
    run:
        # list of lists of error results,
        data   = [] # index of the list is a 'condition' each inner list has 23 error results (each chrom).
        labels = ['PGP1_21 + PGP1_22','PGP1_21','PGP1_22']
        for s in samples: # sample
            datalist = []
            for c in chroms: # chromosome
                assembly_file = "{}/hapcut2_{}/{}.output".format(config['experiments_dir'],s,c)
                runtime_file  = "{}/hapcut2_{}/{}.runtime".format(config['experiments_dir'],s,c)
                frag_file     = "{}/{}/fragmat/chr{}".format(d,s,c)
                vcf_file      = "{}/chr{}.vcf".format(config['variant_vcf_dir'],c)
                truth_file    = "{}/haplotype.chr{}.bac.txt".format(config['BAC_hapblocks'], c)
                err = error_rates.hapblock_hapblock_error_rate(truth_file, assembly_file, frag_file, vcf_file, runtime_file)
                datalist.append(err)

            print("{} results over all chromosomes:".format(s))
            print(sum(datalist,error_rates.error_result(None, 0, 0, 0, 0, 0, 0, 0, None,[],[],[])))

            data.append(datalist)

        pickle.dump(data,open(output.stats_file,"wb"))
        pickle.dump(labels,open(output.labels_file,"wb"))

# RUN HAPCUT2
rule run_hapcut2:
    params:
        job_name = "{s}.{c}_hapcut",
    input:
        frag_file = "{data_dir}/{s}/fragmat/chr{c}",
        vcf_file  = lambda wildcards: expand("{v}/{c}.vcf", v=config['variant_vcf_dir'],c=wildcards.c)
    output:
        hapblocks = "{E}/hapcut2_{s}/{c}.output",
        runtime = "{E}/hapcut2_{s}/{c}.runtime"
    run:
        # run hapcut
        runtime = run_tools.run_hapcut2(config['hapcut2'], input.frag_file, input.vcf_file, output.hapblocks, 100, 100, 1, 0, 15, 0.8)
        with open(output.runtime,'w') as rf:
            print(runtime, file=rf)

# COMBINE FRAGMENT MATRICES
rule merge_fragmat:
    params:
        job_name  = "merge_fragmat",
        ppn       = "1",
        walltime  = "1:00:00",
    input:
        P21 = expand("{{data_dir}}/PGP1_21/fragmat/chr{CHR}",CHR=chroms),
        P22 = expand("{{data_dir}}/PGP1_22/fragmat/chr{CHR}",CHR=chroms)
    output:
        P21_22 = expand("{{data_dir}}/PGP1_21_22/fragmat/chr{CHR}",CHR=chroms)
    run:
        for i1, i2, o in zip(P21, P22, P21_22):
            shell('cat {i1} {i2} > {o}')

# CREATE FRAGMENT MATRIX FILES FOR EXPERIMENT
rule create_fragmat:
    params:
        job_name  = "{s}_create_fragmat",
        ppn       = "4",
        walltime  = "20:00:00",
    input:
        vcfs     = expand("{{data_dir}}/{{s}}/ch{C}.vcf",C=chambers),
        pileups  = expand("{{data_dir}}/{{s}}/ch{C}.pileup",C=chambers),
        beds     = expand("{{data_dir}}/{{s}}/ch{C}.bed",C=chambers),
        var_vcfs = expand("{{data_dir}}/PGP1_VCFs/chr{CHR}.vcf",CHR=chroms)
    output:
        expand("{{data_dir}}/{{s}}/fragmat/chr{CHR}",CHR=chroms)
    run:
        output_dir = os.path.join(data_dir,wildcards.s,'fragmat')
        create_hapcut_fragment_matrix(input.vcfs, input.beds, input.pileups, input.var_vcfs, output_dir)

# simlink data to make path naming scheme consistent between PGP1_21 and PGP_22
rule simlink_PGP1_21:
    input:
        expand("{DIR}/PGP1_21_ch{C}.sorted.fragment.depth5.vcf",DIR=config['PGP_21_vcfs_old'],C=chambers_pad),
        expand("{DIR}/PGP1_21_ch{C}.sorted.fragment.depth.pileup",DIR=config['PGP_21_pileups_old'],C=chambers_pad),
        expand("{DIR}/PGP1_21_FragmentBoundaryCh{C}.bed",DIR=config['PGP_21_boundaries_old'],C=chambers),
#    output:
#        expand("{D}/PGP1_21/ch{C}.{F}",D=config["data_dir"],C=chambers,F=config["input_filetypes"])
    run:
        for c,C in zip(chambers,chambers_pad):
            shell('''
            mkdir -p {data_dir}/PGP1_21
            ln -s {config[PGP_21_vcfs_old]}/PGP1_21_ch{C}.sorted.fragment.depth5.vcf {data_dir}/PGP1_21/ch{c}.vcf
            ln -s {config[PGP_21_pileups_old]}/PGP1_21_ch{C}.sorted.fragment.depth.pileup {data_dir}/PGP1_21/ch{c}.pileup
            ln -s {config[PGP_21_boundaries_old]}/PGP1_21_FragmentBoundaryCh{c}.bed {data_dir}/PGP1_21/ch{c}.bed''')

# simlink data to make path naming scheme consistent between PGP1_21 and PGP_22
rule simlink_PGP1_22:
    input:
        expand("{DIR}/PGP1_22_ch{C}.sorted.fragment.vcf",DIR=config['PGP_22_old'],C=chambers_pad),
        expand("{DIR}/PGP1_22_ch{C}.sorted.fragment.depth.pileup",DIR=config['PGP_22_old'],C=chambers_pad),
        expand("{DIR}/PGP1_22_FragmentBoundaryCh{C}.bed",DIR=config['PGP_22_old'],C=chambers),
#    output:
#        expand("{D}/PGP1_22/ch{C}.{F}",D=config["data_dir"],C=chambers,F=config["input_filetypes"])
    run:
        for c,C in zip(chambers,chambers_pad):
            shell('''
            mkdir -p {data_dir}/PGP1_22
            ln -s {config[PGP_22_old]}/PGP1_22_ch{C}.sorted.fragment.vcf {data_dir}/PGP1_22/ch{c}.vcf
            ln -s {config[PGP_22_old]}/PGP1_22_ch{C}.sorted.fragment.depth.pileup {data_dir}/PGP1_22/ch{c}.pileup
            ln -s {config[PGP_22_old]}/PGP1_22_FragmentBoundaryCh{c}.bed {data_dir}/PGP1_22/ch{c}.bed''')
