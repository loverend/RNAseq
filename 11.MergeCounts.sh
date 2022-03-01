#!/bin/bash
#$ -cwd
#$ -pe shmem 1 -N mergecounts
#$ -q short.qc

echo "------------------------------------------------"
echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}
echo SGE_TASK_FIRST=${SGE_TASK_FIRST}, SGE_TASK_LAST=${SGE_TASK_LAST}, SGE_TASK_STEPSIZE=${SGE_TASK_STEPSIZE}
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "------------------------------------------------"

COUNT_DIR=$1
CODE_DIRECTORY=$2
MAPPING_DIR=$3
KEY=$4
cd ${COUNT_DIR}

##Load Module 
module load R/3.6.2-foss-2019b

echo "********************************************************"
echo "* Job Parameters"
echo "********************************************************"
echo "KEY    		  : ${KEY}"
echo "FASTQDIR     	  : ${FASTQDIR}"
echo "MAPPING_DIR     : ${MAPPING_DIR}"
echo "QC_DIR          : ${QC_DIR}"
echo "GTF     		  : ${GTF}"
echo "FASTQ1          : ${FASTQ1}"
echo "FASTQ2          : ${FASTQ2}"
echo "GENOME          : ${GENOME}"
echo "********************************************************"
echo "Running stage 12. Merging Read Counts"


## Run R script 
Rscript --vanilla ${CODE_DIRECTORY}11.MergeCounts.R

## Now we want to delete BAMS to save memory!
SAMPLE_NAME=$(cat ${KEY} | tail -n+${SGE_TASK_ID} | head -1 | cut -f1 )
echo ${SAMPLE_NAME}

# move to sample directory
DIR_SAMPLE_NAME=${MAPPING_DIR}${SAMPLE_NAME}

## Remove large mapped bam and star files
## Leave the unmapped reads for further analysis 
cd ${DIR_SAMPLE_NAME}
rm *.bam*
rm -r *STARgenome 
rm *out 
rm *tab
## Compress the unmapped reads! 
gzip * 

## Remove the trimmed fastq file 
cd ${MAPPING_DIR}AdapterTrimmed
rm ${SAMPLE_NAME}*_trimmed.fq.gz


# Done 
echo
echo "********************************************************"
echo "["`date`"] Done"
echo "********************************************************"
exit 0


