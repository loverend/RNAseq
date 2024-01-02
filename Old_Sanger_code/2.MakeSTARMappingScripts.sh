#!/bin/bash

### The sample key made by the 1.MakeSampleKey.R script
FILE=$BASEDIR"/mapping.info.txt"

# Make a directory for the mapping step and move into it
if [[ ! -e "Mapping" ]]; then

mkdir ./Mapping
fi

cd Mapping

OUTPUT_DIR=$BASEDIR"/Mapping/"

# Make a directory for the mapping scripts
if [[ ! -e "Mappingscripts" ]]; then

mkdir ./Mappingscripts
fi

# Start a script that will submit all the per-sample mapping scripts
rst="./submit_jobs_for_mapping.sh"
echo "" > $rst

# read in the sample key
# for each sample input the fastq files to the STAR mapping script
# this makes a script per sample to run the alignment
while read -r line
do
	sample="$(echo "$line" | cut -f 1)"
	echo $sample
	
	DIR_SAMPLE_NAME=$OUTPUT_DIR$sample
	if [[ ! -e $DIR_SAMPLE_NAME ]]; then
                mkdir -p $DIR_SAMPLE_NAME
  fi

	fastq="$(echo "$line" | cut -f 2)"
	fastq1="$(echo "$fastq" | sed -e "s|,|_R1.fastq.gz,$FASTQ\/|g")"
	fastq2="$(echo "$fastq" | sed -e "s|,|_R2.fastq.gz,$FASTQ\/|g")"
	fastq1=$fastq1"_R1.fastq"
	fastq2=$fastq2"_R2.fastq"
	fastq1=$FASTQ/$fastq1
	fastq2=$FASTQ/$fastq2

	echo $fastq1
	echo $fastq2

	exc_name=$OUTPUT_DIR"Mappingscripts/mapping."$sample.sh
	cmd="sh /well/immune-rep/shared/CODE/RNASeqSTARPipeline/2.STARmapping.sh $sample $fastq1 $fastq2"
	$cmd > $exc_name

	echo "qsub "$exc_name >> $rst
done < $FILE

cd ..
