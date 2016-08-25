# SISSOR
Scripts for analysis of SISSOR data

Steps to run haplotyping pipeline for someone using UCSD TSCC's Torque/PBS system:
```
git clone https://github.com/pjedge/SISSOR```
cd SISSOR
export SNAKEMAKE=/home/pedge/installed/opt/python/bin/snakemake
```
change email and working_dir at top of cluster.yaml
```
$SNAKEMAKE simlinks # make simlinks to Eric's data
ln -s /home/pedge/git/SISSOR/sissor_project/data/PGP1_VCFs sissor_project/data/PGP1_VCFs  # simlink to CGI heterozygous variant data
```
test that pipeline dependencies are met, should list steps needed and not throw error
```
./test_pipeline.sh
```
run whole haplotyping pipeline, it will manage Qsub jobs for you, wise to do this in a screen session
```
./run_pipeline.sh
```
