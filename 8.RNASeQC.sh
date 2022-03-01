#!/bin/bash
#$ -cwd
#$ -pe shmem 1 -N rnaseqc
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
MAPPING_DIR=$2
KEY=$3
GTF=$4
CODE_DIRECTORY=$5
BASEDIR=$6
RNASEQDIR=$7


# Get sample ID for this task
SAMPLE_NAME=$(cat ${KEY} | tail -n+${SGE_TASK_ID} | head -1 | cut -f1 )
echo ${SAMPLE_NAME}

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
echo "Running stage 8. RNASEQ QC"


# Make sample directory
DIR_SAMPLE_NAME=${QC_DIR}${SAMPLE_NAME}
if [[ ! -d ${DIR_SAMPLE_NAME} ]]; then
              echo "Creating Directory"
			  mkdir -p ${DIR_SAMPLE_NAME}
fi

## Note will need to provide location of rnaseqc
## RUN COMMAND
${RNASEQDIR}rnaseqc.v2.3.5.linux ${GTF} ${MAPPING_DIR}${SAMPLE_NAME}/${SAMPLE_NAME}.Aligned.sortedByCoord.out.bam ${DIR_SAMPLE_NAME}/


NWCMD="echo ${SAMPLE_NAME} >> ${BASEDIR}/job_rnaseqc.txt"
eval "${NWCMD}"


# Done 
echo
echo "********************************************************"
echo "["`date`"] Done"
echo "********************************************************"
exit 0
