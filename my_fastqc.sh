#!/bin/sh

# Note that all qsub directives start with '#$' while comment line starts with '#'

# name for the job.
#$ -N my_fastqc

# The cluster resources (the nodes) are grouped into Queues.
#$ -q b.q

# parallel environment (-pe) requesting 4 slots/cores/threads.
#$ -pe smp 4

# current working directory.
#$ -cwd

# merge standard error with standard output.
#$ -j y

# email options about the status of your job.
#$ -m abe
#$ -M <your_email>@emory.edu

PROJ_DIR=$HOME/rnaseq
SEQ_DATA_DIR=$PROJ_DIR/data/ibs_class
OUT_DIR=$PROJ_DIR/out

module load FastQC/0.11.4

# QC sub-directory will be created 
if [ ! -d $OUT_DIR/QC ]; then
 mkdir -p $OUT_DIR/QC
fi

# QC metrics for all of FASTQs
fastqc -t 9 $SEQ_DATA_DIR/*.fastq.gz -o $OUT_DIR/QC

module unload FastQC/0.11.4
