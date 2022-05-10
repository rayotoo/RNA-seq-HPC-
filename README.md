# Single Cell RNA-seq Tutorial
[![GitHub issues](https://img.shields.io/github/issues/rayotoo/RNA-seq-HPC-?style=flat-square)](https://github.com/rayotoo/RNA-seq-HPC-/issues)
[![GitHub stars](https://img.shields.io/github/stars/rayotoo/RNA-seq-HPC-?style=flat-square&color=important)](https://github.com/rayotoo/RNA-seq-HPC-/stargazers)
[![GitHub forks](https://img.shields.io/github/forks/rayotoo/RNA-seq-HPC-?style=flat-square&color=blueviolet)](https://github.com/rayotoo/RNA-seq-HPC-/network/members)
[![LICENSE](https://img.shields.io/github/license/rayotoo/RNA-seq-HPC-?style=flat-square&color=green)](https://github.com/rayotoo/RNA-seq-HPC-/blob/main/LICENSE)

# Adapted from Training material used in Spring 2018 by Computational Biology & Bioinformatics, Spring 2018, Emory University

# RNA-Seq Exercise: raw reads to differential expression

In this exercise we will go through all of the steps of a simple RNA-Seq differential expression analysis using High Performance Computing (HPC), which will include:

## 1. Connect HPC Server via SSH.
SSH allows you to connect to Human Genetics Compute Cluster (HGCC) server securely and perform linux command-line operations.

## 2. Create project directory and sub-directories.
From your home directory i.e /home/<user_name>

## 3. Copy raw data into ~/rnaseq/data and extract .tar.gz file.
Test dataset is at **/scratch** directory. 

## 4. Quality control of raw reads, FastQC.
Create a shell script called **my_fastqc.sh** in your `~/rnaseq/script` sub-directory and submit this job from `~/rnaseq/logs` sub-directory. This script is available for your convenience at [my_fastqc.sh](https://bitbucket.org/adinasarapu/ibs_class/src).

## 5. Mapping raw reads to a reference genome, STAR.
**A).** Download and uncompress Human Reference Genome Sequence data from
[https://support.illumina.com/sequencing/sequencing_software/igenome.html](https://support.illumina.com/sequencing/sequencing_software/igenome.html) using `wget`.

## 6. Quantification of reads against genes, HTSeq.
Create a shell script called **my_htseq_count.sh** in your `~/rnaseq/script` sub-directory and submit this job from `~/rnaseq/logs` sub-directory. This script is available for your convenience [my_htseq_count.sh](https://bitbucket.org/adinasarapu/ibs_class/src).

## 7. Differential expression analysis.
Now using the above merged table perform 
  
 &emsp; a) **Gene filtering**: Keep genes with counts (> 5) in at least 3 (50%) samples. 
  
 &emsp; b) Perform **trimmed mean of M-values (TMM) normalization** (using edgeR, R package) and get the normalization factors. 

The main aim in TMM normalization is to account for library size variation between samples of interest. After filtering, it is a good idea to reset the library sizes. The library size and normalization factors are multiplied together to act as the effective library size.


