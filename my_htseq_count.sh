#!/bin/sh

#$ -N my_htseq
#$ -q b.q
#$ -pe smp 9
#$ -cwd
#$ -j y
#$ -m abe
#$ -M <your_email>@emory.edu

PROJ_DIR=$HOME/rnaseq
OUT_DIR=$PROJ_DIR/out/HTSeq

# Create a sub-directory for output
if [ ! -d $OUT_DIR ]; then
 /bin/mkdir -p $OUT_DIR
fi

# Reference Genome Annotations directory
ANNOTATION_DIR=/sw/hgcc/Data/Illumina/Homo_sapiens/UCSC/hg38/Annotation/Genes

# STAR aligned BAM files  
BAM_FILE_45=$PROJ_DIR/out/MAPPING/SL100145/Aligned.sortedByCoord.out.bam
BAM_FILE_46=$PROJ_DIR/out/MAPPING/SL100146/Aligned.sortedByCoord.out.bam
BAM_FILE_47=$PROJ_DIR/out/MAPPING/SL100147/Aligned.sortedByCoord.out.bam
BAM_FILE_49=$PROJ_DIR/out/MAPPING/SL100149/Aligned.sortedByCoord.out.bam
BAM_FILE_50=$PROJ_DIR/out/MAPPING/SL100150/Aligned.sortedByCoord.out.bam
BAM_FILE_51=$PROJ_DIR/out/MAPPING/SL100151/Aligned.sortedByCoord.out.bam

# Load a bioinformatic package  
module load Anaconda2/4.2.0

# htseq-count from STAR aligned BAM files
htseq-count -f bam -m union -r pos -i gene_id -a 10 -s no $BAM_FILE_45 $ANNOTATION_DIR/genes.gtf > $OUT_DIR/SL100145.counts
htseq-count -f bam -m union -r pos -i gene_id -a 10 -s no $BAM_FILE_46 $ANNOTATION_DIR/genes.gtf > $OUT_DIR/SL100146.counts
htseq-count -f bam -m union -r pos -i gene_id -a 10 -s no $BAM_FILE_47 $ANNOTATION_DIR/genes.gtf > $OUT_DIR/SL100147.counts
htseq-count -f bam -m union -r pos -i gene_id -a 10 -s no $BAM_FILE_49 $ANNOTATION_DIR/genes.gtf > $OUT_DIR/SL100149.counts
htseq-count -f bam -m union -r pos -i gene_id -a 10 -s no $BAM_FILE_50 $ANNOTATION_DIR/genes.gtf > $OUT_DIR/SL100150.counts
htseq-count -f bam -m union -r pos -i gene_id -a 10 -s no $BAM_FILE_51 $ANNOTATION_DIR/genes.gtf > $OUT_DIR/SL100151.counts

# Unload the package
module unload Anaconda2/4.2.0
