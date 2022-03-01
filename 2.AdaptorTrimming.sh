#!/bin/bash
#$ -cwd
#$ -pe shmem 5 -N trimming
#$ -q short.qc
# Katie Burnham
# Adapted Lauren Overend

echo "------------------------------------------------"
echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}
echo SGE_TASK_FIRST=${SGE_TASK_FIRST}, SGE_TASK_LAST=${SGE_TASK_LAST}, SGE_TASK_STEPSIZE=${SGE_TASK_STEPSIZE}
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "------------------------------------------------"

## Arguments 
KEY=$1
FASTQDIR=$2
MAPPING_DIR=$3
BASEDIR=$4
PAIRED=$5

# Set parameters
module load Trim_Galore/0.6.2-GCCcore-8.2.0-Java-11

FASTQ=$(cat $KEY | tail -n+${SGE_TASK_ID} | head -1 | cut -f1 )
echo ${FASTQ}
FASTQ1=${FASTQDIR}/${FASTQ}_R1_001.fastq.gz
FASTQ2=${FASTQDIR}/${FASTQ}_R2_001.fastq.gz
READS="${FASTQ1} ${FASTQ2}"

if [[ ${PAIRED} == "FALSE" ]]; then
	READS="${FASTQ1}"
fi 


echo "********************************************************"
echo "* Job Parameters"
echo "********************************************************"
echo "KEY    		  : ${KEY}"
echo "FASTQDIR     	  : ${FASTQDIR}"
echo "MAPPING_DIR     : ${MAPPING_DIR}"
echo "FASTQ1          : ${FASTQ1}"
echo "FASTQ2          : ${FASTQ2}"
echo "PAIRED          : ${PAIRED}"
echo "********************************************************"


## Check if outdirectory Needs to be created!
TRIMDIR=${MAPPING_DIR}/AdapterTrimmed
if [[ ! -d ${TRIMDIR} ]]; then
echo "Creating Mapping Directory"
mkdir ${TRIMDIR}
else 
echo "Directory Present"
fi

## Run Command
if [[ ${PAIRED} == "TRUE" || -z "$PAIRED" ]]; then
trim_galore --paired --gzip --output_dir ${TRIMDIR} ${READS}
fi

if [[ ${PAIRED} == "FALSE" ]]; then
trim_galore  --gzip --output_dir ${TRIMDIR} ${READS}
fi 


NWCMD="echo ${FASTQ} >> ${BASEDIR}/job_trimming.txt"
eval "${NWCMD}"


# Done 
echo
echo "********************************************************"
echo "["`date`"] Done"
echo "********************************************************"
exit 0