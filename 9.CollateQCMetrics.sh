#!/bin/bash
#$ -cwd
#$ -pe shmem 1 -N mergeqc
#$ -q short.qc

echo "------------------------------------------------"
echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}
echo SGE_TASK_FIRST=${SGE_TASK_FIRST}, SGE_TASK_LAST=${SGE_TASK_LAST}, SGE_TASK_STEPSIZE=${SGE_TASK_STEPSIZE}
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "------------------------------------------------"

QC_DIR=$1
CODE_DIRECTORY=$2
cd ${QC_DIR}

## Load Module 
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
echo "Running stage 9. Collate QC Metrics"

## Run Rscript within the QC directory
Rscript --vanilla ${CODE_DIRECTORY}9.CollateQCMetrics.R

# Done 
echo
echo "********************************************************"
echo "["`date`"] Done"
echo "********************************************************"
exit 0
