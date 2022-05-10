#!/bin/sh

#$ -N my_star
#$ -q b.q
#$ -pe smp 9
#$ -cwd
#$ -j y
#$ -m abe
#$ -M <your_email>@emory.edu

PROJ_DIR=$HOME/rnaseq
DATA_DIR=$PROJ_DIR/data/ibs_class
OUT_DIR=$PROJ_DIR/out

SEQUENCE_DIR=/sw/hgcc/Data/Illumina/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta
ANNOTATION_DIR=/sw/hgcc/Data/Illumina/Homo_sapiens/UCSC/hg38/Annotation/Genes

module load STAR/2.5.3a
module load samtools/1.5

# make sub-directories for output 
if [ ! -d $OUT_DIR/INDEX ]; then
 mkdir -p $OUT_DIR/INDEX
fi

if [ ! -d $OUT_DIR/MAPPING ]; then
 mkdir -p $OUT_DIR/MAPPING/{SL100145,SL100146,SL100147,SL100149,SL100150,SL100151}
fi

# build STAR indexing of reference genome, hg38

STAR --runMode genomeGenerate \
 --genomeDir $OUT_DIR/INDEX \
 --genomeFastaFiles $SEQUENCE_DIR/genome.fa \
 --runThreadN 9

# For sample : SL100145
# Read alignment to refrence genome (or transcriptome)
STAR --genomeDir $OUT_DIR/INDEX \
 --readFilesIn $DATA_DIR/Human_R1_SL100145.fastq.gz $DATA_DIR/Human_R2_SL100145.fastq.gz \
 --runThreadN 9 \
 --readFilesCommand zcat \
 --sjdbGTFfile $ANNOTATION_DIR/genes.gtf \
 --outSAMtype BAM SortedByCoordinate \
 --outFileNamePrefix $OUT_DIR/MAPPING/SL100145/
# create bam index file
samtools index $OUT_DIR/MAPPING/SL100145/Aligned.sortedByCoord.out.bam

# For sample : SL100146
# Read alignment to refrence genome (or transcriptome)
STAR --genomeDir $OUT_DIR/INDEX \
 --readFilesIn $DATA_DIR/Human_R1_SL100146.fastq.gz $DATA_DIR/Human_R2_SL100146.fastq.gz \
 --runThreadN 9 \
 --readFilesCommand zcat \
 --sjdbGTFfile $ANNOTATION_DIR/genes.gtf \
 --outSAMtype BAM SortedByCoordinate \
 --outFileNamePrefix $OUT_DIR/MAPPING/SL100146/
# create bam index file
samtools index $OUT_DIR/MAPPING/SL100146/Aligned.sortedByCoord.out.bam

# For sample : SL100147
# Read alignment to refrence genome (or transcriptome)
STAR --genomeDir $OUT_DIR/INDEX \
 --readFilesIn $DATA_DIR/Human_R1_SL100147.fastq.gz $DATA_DIR/Human_R2_SL100147.fastq.gz \
 --runThreadN 9 \
 --readFilesCommand zcat \
 --sjdbGTFfile $ANNOTATION_DIR/genes.gtf \
 --outSAMtype BAM SortedByCoordinate \
 --outFileNamePrefix $OUT_DIR/MAPPING/SL100147/
# create bam index file
samtools index $OUT_DIR/MAPPING/SL100147/Aligned.sortedByCoord.out.bam

# For sample : SL100149
# Read alignment to refrence genome (or transcriptome)
STAR --genomeDir $OUT_DIR/INDEX \
 --readFilesIn $DATA_DIR/Human_R1_SL100149.fastq.gz $DATA_DIR/Human_R2_SL100149.fastq.gz \
 --runThreadN 9 \
 --readFilesCommand zcat \
 --sjdbGTFfile $ANNOTATION_DIR/genes.gtf \
 --outSAMtype BAM SortedByCoordinate \
 --outFileNamePrefix $OUT_DIR/MAPPING/SL100149/
# create bam index file
samtools index $OUT_DIR/MAPPING/SL100149/Aligned.sortedByCoord.out.bam

# For sample : SL100150
# Read alignment to refrence genome (or transcriptome)
STAR --genomeDir $OUT_DIR/INDEX \
 --readFilesIn $DATA_DIR/Human_R1_SL100150.fastq.gz $DATA_DIR/Human_R2_SL100150.fastq.gz \
 --runThreadN 9 \
 --readFilesCommand zcat \
 --sjdbGTFfile $ANNOTATION_DIR/genes.gtf \
 --outSAMtype BAM SortedByCoordinate \
 --outFileNamePrefix $OUT_DIR/MAPPING/SL100150/
# create bam index file
samtools index $OUT_DIR/MAPPING/SL100150/Aligned.sortedByCoord.out.bam

# For sample : SL100151
# Read alignment to refrence genome (or transcriptome)
STAR --genomeDir $OUT_DIR/INDEX \
 --readFilesIn $DATA_DIR/Human_R1_SL100151.fastq.gz $DATA_DIR/Human_R2_SL100151.fastq.gz \
 --runThreadN 9 \
 --readFilesCommand zcat \
 --sjdbGTFfile $ANNOTATION_DIR/genes.gtf \
 --outSAMtype BAM SortedByCoordinate \
 --outFileNamePrefix $OUT_DIR/MAPPING/SL100151/
# create bam index file
samtools index $OUT_DIR/MAPPING/SL100151/Aligned.sortedByCoord.out.bam

module unload samtools/1.5
module unload STAR/2.5.3a
