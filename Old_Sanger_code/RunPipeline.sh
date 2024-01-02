#!/bin/bash

# Run pipeline

module load R/3.2.2

######### Read in your configuration file

CONFIG=$1

echo "Using configuration file" $CONFIG

echo "Setting up pipeline using your configuration file" $CONFIG

source $CONFIG



######## Genome references

if [ "$GENOME" == "GRCh38" ]
then
	echo "The pipeline will use build 38"
	export GTF=/well/immune-rep/shared/CODE/RNASeqSTARPipeline/data/GRCh38/annotation/Homo_sapiens.GRCh38.99.gtf
	export GENOME=/well/immune-rep/shared/CODE/RNASeqSTARPipeline/data/GRCh38/star_indices_overhang74
fi

if [ "$READLENGTH" -ne "75" ]
then
  echo "The pipeline will be run using the generic value of 100 for the genome index. According to the STAR manual, this will work well but is not optimised for your read length."
  export GENOME=/well/immune-rep/shared/CODE/RNASeqSTARPipeline/data/GRCh38/star_indices_overhang100
fi

if [ "$GENOME" == "hg19" ]
then
	echo "The pipeline will use build 19"
	export GTF=/well/jknight/reference/mapping/GRCh37/gencode.v19.chr_patch_hapl_scaff.annotation.gtf
	export GENOME=/well/immune-rep/shared/CODE/RNASeqSTARPipeline/data/hg19/star_indices_overhang74
fi

if [ "$DEDUP" == "N" ]
then
	echo "Duplicates will not be removed"
	export FORCOUNTS=dups.marked
fi

if [ "$DEDUP" == "Y" ]
then
	echo "Duplicates will be removed"
	export FORCOUNTS=nodup
fi

##### move to your base directory containing the pipeline scripts and your 
# inital sample-fastq information key

cd $BASEDIR

##### make a key to link sample ID to multiple fastq files
echo "Making sample key"

Rscript --vanilla /well/immune-rep/shared/CODE/RNASeqSTARPipeline/1.MakeSampleKey.R $KEY

#### make scripts
echo "Making scripts for each sample"

# Make a script for each sample to align reads to the genome using STAR
sh /well/immune-rep/shared/CODE/RNASeqSTARPipeline/2.MakeSTARMappingScripts.sh

# Make a script for each sample to filter out duplicates and output some mapping metrics using Picard
sh /well/immune-rep/shared/CODE/RNASeqSTARPipeline/3.MakeFilteringScripts.sh

# Make a script for each sample to count reads for each feature
sh /well/immune-rep/shared/CODE/RNASeqSTARPipeline/4.MakeCountScripts.sh




#### submit jobs
echo "Submitting jobs for alignment"

sh $BASEDIR/Mapping/submit_jobs_for_mapping.sh
# rm $BASEDIR/Mapping/submit_jobs_for_mapping.sh

# echo "Mapping QC"  # hold on mapping, need to edit to input $BASEDIR?

# qsub -hold_jid "map.*" /well/jknight/RNASeqSTARPipeline/MappingQC.sh

echo "Submitting jobs for filtering"

sh $BASEDIR/Filtering/submit_jobs_for_filtering.sh
#rm $BASEDIR/Filtering/submit_jobs_for_filtering.sh

echo "Submitting jobs for counting"

sh $BASEDIR/Counts/submit_jobs_for_counting.sh
#rm $BASEDIR/Counts/submit_jobs_for_counting.sh

echo "Merge counts"

qsub -hold_jid "count.*" /well/immune-rep/shared/CODE/RNASeqSTARPipeline/5.MergeCounts.sh
