# SISSOR (Single-Stranded Sequencing using micrOfluidic Reactors)

## About
Read the pre-print paper on bioarxiv here:
http://biorxiv.org/content/early/2017/05/08/135384

### Abstract (from bioarxiv):
Accurate detection of variants and long-range haplotypes in genomes of single human cells remains very challenging. Common approaches require extensive in vitro amplification of genomes of individual cells using DNA polymerases and high-throughput short-read DNA sequencing. These approaches have two notable drawbacks. First, polymerase replication errors could generate tens of thousands of false positive calls per genome. Second, relatively short sequence reads contain little to no haplotype information. Here we report a method, which is dubbed SISSOR (Single-Stranded Sequencing using micrOfluidic Reactors), for accurate single-cell genome sequencing and haplotyping. A microfluidic processor is used to separate the Watson and Crick strands of the double-stranded chromosomal DNA in a single cell and to randomly partition megabase-size DNA strands into multiple nano-liter compartments for amplification and construction of barcoded libraries for sequencing. The separation and partitioning of large single-stranded DNA fragments of the homologous chromosome pairs allows for the independent sequencing of each of the complementary and homologous strands. This enables the assembly of long haplotypes and reduction of sequence errors by using the redundant sequence information and haplotype-based error removal. We demonstrated the ability to sequence single-cell genomes with error rates as low as 1e-8 and average 500kb long DNA fragments that can be assembled into haplotype contigs with N50 greater than 7Mb. The performance could be further improved with more uniform amplification and more accurate sequence alignment. The ability to obtain accurate genome sequences and haplotype information from single cells will enable applications of genome sequencing for diverse clinical needs.

## Code
This repository contains the majority of the codebase for analysis of SISSOR libraries (with the exception of the read mapping and the HMM-based fragment boundary segmentation steps).

It is important to note that this code is here for documentation purposes and has relatively little applicability except with data from the SISSOR device. Various components may be adapted in the future for general haplotype assembly or MDA-based protocols.

The general workflow, executed by the rules in the Snakefile, is the following:
1. Generate a high-quality BAC reference haplotype from raw BAC reads for PGP1.
2. Generate a high-quality reference dataset of SNV and reference calls for PGP1 by intersecting high-quality SNV/ref datasets from two cell lines (lymphocyte complete genomics, and fibroblast 60x Illumina WGS).
3. Perform a samtools mpileup of 72 SISSOR libraries from PGP1f cells (3 cells x 24 reaction chambers).
4. Estimate various model parameters from the mpileup data.
5. Use custom variant caller (codename SISSORhands) to make SNV/ref calls at every genomic position in every SISSOR library, as well as consensus calls over all libraries.
6. Generate HapCUT2 format haplotype fragments from the SISSORhands variant calls, based on heterozygous SNV positions from 60x WGS of PGP1f cells.
7. Process haplotype fragments to split and/or filter them at (a) clusters of mixed alleles (haploid strands from two haplotypes overlapping in the same chamber to make a diploid mixture) and (b) detectable fragment-fragment switch errors (2 or more switches with respect to other fragments)
8. Assemble the processed haplotype fragments into haplotype blocks with HapCUT2.
9. Calculate accuracy of haplotype blocks with respect to BAC haplotypes from step [1] and also statistics such as N50.
10. Assign haplotype fragments back to haplotype blocks.
11. Calculate the accuracy of various base calling approaches with respect to reference set from step [2]:
    1. consensus calls from step [5]
    2. chamber-specific calls from step [5] that match between strands on the same haplotype:
        1. in the same cell (to get validated cell specific mutations)
        2. in all cells (to get maximum coverage of haplotype-validated calls)
        3. between different cells (cross-cell) (to estimate maximum bound on error rate of the haplotype-strand-matching method against a reference without inflation from cell-specific mutations)

For more details, see the SISSOR manuscript (above) and supplementary materials.

## Requirements
* Snakemake (execution on a cluster, with appropriate cluster configuration, is basically required)
* BAM files for 3 SISSOR libraries, placed in bams/{cells}/ch{ch}.bam
     * {cells} should be named PGP1_21, PGP1_22, PGP1_A1
     * {ch} should be in the range 1-24
* hg38ToHg19 chainfile (hg38ToHg19.over.chain)
* samtools
* BCFtools
* Picard
* CrossMap
* HapCUT2
* Pypy
* Complete Genomics Variant Calls for PGP1 (GFF format)
* BAC library SNV calls for PGP1 (VCF format)
* Fragment boundaries for SISSOR libraries (called with HMM, see methods).
