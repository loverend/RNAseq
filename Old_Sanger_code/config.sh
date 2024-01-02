#!/bin/bash

# This is your configuration file
# It is the only file you need to edit
# Put it in your project's base directory
# When you run the pipeline use the path to this file as an argument
# i.e. sh /well/jknight/RNASeqSTARPipeline/RunPipeline.sh /well/jknight/AbuDhabiRNA/Katie/STAR/config.sh
# It will then use the other scripts in the RNASeqSTARPipeline folder to create scripts for each sample and run all the mapping, filtering, and count steps
# Output files will be put in subdirectories of your base directory

###################################
#
# Paths to your files
#
###################################

# Your project directory
export BASEDIR=/well/immune-rep/shared/MISEQ/VIRAL_SEQ/HUMAN_FIRST
echo "Base directory: "$BASEDIR


# Pathway to your sample information file
# This is made by copying the information in the emails from Core
# e.g.
#FileName	            Sample Name
#WTCHG_384296_201191	    S97
#WTCHG_384296_202179	    S328
#WTCHG_384296_203167	    S81

export KEY=/gpfs2/well/immune-rep/shared/CODE/RNASeqSTARPipeline/SAMPLE_FILES/EBV_samples_KATIEPIPELINE.txt
echo "Sample information file: "$KEY


# Pathway to directory containing your fastq files (or softlinks to them)
# If you have downloaded batches of data and have them stored in subdirectories make a new folder with softlinks to all:
# mkdir /well/jknight/MyProject/fastq
# ln -s /well/jknight/MyProject/*/*fastq.gz /well/jknight/MyProject/fastq

export /well/immune-rep/shared/PUBLISHED_DATASETS/BULK_RNA/EBV_LCL/Merged_Data
echo "Fastq directory: "$FASTQ


###################################
#
# Parameters
#
###################################

# What genome build do you want to use?
export GENOME=GRCh38
#export GENOME=hg19
echo "Genome build: "$GENOME

export READLENGTH=100
echo "Read length: "$READLENGTH

# Do you want to remove duplicates?
#export DEDUP=Y
export DEDUP=N
echo "Deduplication: "$DEDUP
