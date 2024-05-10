# Differential-NGS-Analysis

This repository contains all the necessary code for reproducing the JUND ChIP-seq results from the paper titled: "JUND plays a genome-wide role in the quiescent to contractile switch in the pregnant human myometrium" by Nawrah Khader, Anna Dorogin, Oksana Shynlova, and Jennifer A. Mitchell. 

# Requirements

This script requires the following R packages:

- tibble
- tidyverse
- DiffBind
- DESeq2
- RColorBrewer
- ComplexHeatmap
- circlize
- digest
- cluster

# Usage

To run this analysis, execute the R script (originally written for ChIP-seq analysis but can also be done for ATAC-seq)

# • Input files:
- aligned ChIP/ATAC-seq files (including bam and index [.bai] files)
- peak files in the bed or narrowpeak file format
- .csv file containing all sample information

# • Output files:
- file: tab-delimited file of all consenus peaks
- plot: Venn diagram of the overlap between different treatment/conditions
- plot: Heatmap and PCA clustering of samples
- file: .csv file reporting all differential peaks between different treatment/conditions
- plot: Volcano plot of differential peaks
- file: .csv file of read counts for all ChIP/ATAC replicates at all assessed peaks
- plot: Heatmap (with clustering) of scaled counts for all replicates
- file: .csv file containing scaled coounts with associated cluster number that corresponds the heatmap
