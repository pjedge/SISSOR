# author: Peter Edge
# 3/29/2016
# email: pedge@eng.ucsd.edu

# this is a snakemake Snakefile, written using snakemake 3.5.5
# configure paths and values in cluster.yaml and config.yaml to your system
# example execution for a TORQUE cluster:
# snakemake -j 200 --cluster-config cluster.yaml --cluster "qsub -A {cluster.group} -V -q {cluster.queue} -o {cluster.qsub_stdout_dir}/{params.job_name}.o -e {cluster.qsub_stdout_dir}/{params.job_name}.e -N {params.job_name} -l nodes=1:ppn={cluster.ppn} -l walltime={cluster.walltime} -M {cluster.email} -m e -d {cluster.working_dir}" --local-cores 1

configfile: "config.yaml"
localrules: simlink_PGP1_22, simlink_PGP1_21

from create_hapcut_fragment_matrix import create_hapcut_fragment_matrix

# qsub stderr and stdout directories are not automatically created by snakemake!
if not os.path.exists(config['qsub_stdout_dir']):
    os.makedirs(config['qsub_stdout_dir'])

# specify chambers
data_dir = config["data_dir"]
variant_vcf_dir = config["variant_vcf_dir"]
chambers = list(range(1,25))
chambers_pad = ['{0:02d}'.format(c) for c in chambers]
chroms   = list(range(1,23))

rule make_PGP1_21_fragmat:
    input:
        expand("{D}/PGP1_21/fragmat_all/chr{CHR}",D=data_dir,CHR=chroms)

rule merge_fragmat:
    params:
        job_name    = "{I}.{CHR}_merge",
        ppn         = "1",
        walltime    = "1:00:00",
    input:
        expand("{{D}}/{{I}}/fragmat_sep/chamber{C}/chr{{CHR}}",C=chambers)
    output:
        "{D}/{I}/fragmat_all/chr{CHR}"
    run:
        for i in input:
            shell("cat {i} >> {output}")

rule create_fragmat:
    params:
        job_name  = "{I}_create_fragmat",
        ppn       = "1",
        walltime  = "20:00:00",
        #V         = config["variant_vcf_dir"],
    input:
        vcf         = expand("{{data_dir}}/{{I}}/ch{C}.vcf",C=chambers),
        pileup      = expand("{{data_dir}}/{{I}}/ch{C}.pileup",C=chambers),
        bed         = expand("{{data_dir}}/{{I}}/ch{C}.bed",C=chambers),
        variant_vcf = "/oasis/tscc/scratch/pedge/sissor/PGP1_VCFs/chr{CHR}.vcf", # hardcoding this in temporarily. config dictionary was acting janky for no apparent reason
    output:
        expand("{{data_dir}}/{{I}}/fragmat_sep/chamber{C}/chr{{CHR}}",C=chambers)
    run:
        for vcf, pileup, bed, variant_vcf, o in zip(input.vcf, input.pileup, input.bed, input.variant_vcf, output):
            create_hapcut_fragment_matrix(vcf, bed, variant_vcf, pileup, o)

'''
rule create_fragmat:
    params:
        job_name  = "{I}.{C}_create_fragmat",
        ppn       = "1",
        walltime  = "20:00:00",
        #V         = config["variant_vcf_dir"],
    input:
        vcf         = "{data_dir}/{I}/ch{C}.vcf",
        pileup      = "{data_dir}/{I}/ch{C}.pileup",
        bed         = "{data_dir}/{I}/ch{C}.bed",
        variant_vcf = "/oasis/tscc/scratch/pedge/sissor/PGP1_VCFs/chr{CHR}.vcf", # hardcoding this in temporarily. config dictionary was acting janky for no apparent reason
    output:
        "{data_dir}/{I}/fragmat_sep/chamber{C}/chr{CHR}"
    run:
        create_hapcut_fragment_matrix(input.vcf, input.bed, input.variant_vcf, input.pileup, output)
'''
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
            ln -s {config[PGP_21_boundaries_old]}/PGP1_21_FragmentBoundaryCh{c}.bed {data_dir}/PGP1_21/ch{c}.bed
            ''')

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
            ln -s {config[PGP_22_old]}/PGP1_22_FragmentBoundaryCh{c}.bed {data_dir}/PGP1_22/ch{c}.bed
            ''')
