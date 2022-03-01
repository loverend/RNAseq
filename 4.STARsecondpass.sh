#!/bin/bash
#$ -cwd
#$ -pe shmem 6 -N star2ndpass
#$ -q short.qc

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
GTF=$4
GENOME=$5
BASEDIR=$6
PAIRED=$7

# Set parameters
module load STAR/2.7.3a-GCC-8.3.0

# Get sample ID for this task
SAMPLE_NAME=$(cat ${KEY} | tail -n+${SGE_TASK_ID} | head -1 | cut -f1 )
echo ${SAMPLE_NAME}

# Make sample directory
DIR_SAMPLE_NAME=${MAPPING_DIR}${SAMPLE_NAME}
if [[ ! -d ${DIR_SAMPLE_NAME} ]]; then
    echo "Creating Sample Directory"
	mkdir -p ${DIR_SAMPLE_NAME}
fi

# Get fastq file names
FASTQ=$(cat $KEY | tail -n+${SGE_TASK_ID} | head -1 | cut -f1 )
FASTQ1=${MAPPING_DIR}AdapterTrimmed/${FASTQ}_R1_001_val_1.fq.gz
FASTQ2=${MAPPING_DIR}AdapterTrimmed/${FASTQ}_R2_001_val_2.fq.gz

if [[ ${PAIRED} == "FALSE" ]]; then
	FASTQ1=${MAPPING_DIR}AdapterTrimmed/${FASTQ}_R1_001_trimmed.fq.gz
	FASTQ2=""
fi 


echo "********************************************************"
echo "* Job Parameters"
echo "********************************************************"
echo "KEY    		  : ${KEY}"
echo "FASTQDIR     	  : ${FASTQDIR}"
echo "MAPPING_DIR     : ${MAPPING_DIR}"
echo "GTF     		  : ${GTF}"
echo "FASTQ1          : ${FASTQ1}"
echo "FASTQ2          : ${FASTQ2}"
echo "GENOME          : ${GENOME}"
echo "PAIRED          : ${PAIRED}"
echo "********************************************************"

STAR --genomeDir ${GENOME} \
      --runThreadN 6 \
      --readFilesIn ${FASTQ1} ${FASTQ2} \
      --readFilesCommand gunzip -c \
      --outFileNamePrefix ${DIR_SAMPLE_NAME}/${SAMPLE_NAME}. \
      --outSAMtype BAM SortedByCoordinate \
      --limitSjdbInsertNsj 10000000 \
      --sjdbGTFfile ${GTF} \
      --outFilterMultimapNmax 20 \
      --alignSJoverhangMin 8 \
      --alignSJDBoverhangMin 1 \
      --outFilterMismatchNmax 999 \
      --outFilterMismatchNoverReadLmax 0.04 \
      --alignIntronMin 20 \
      --alignIntronMax 1000000 \
      --alignMatesGapMax 1000000 \
      --sjdbFileChrStartEnd ${MAPPING_DIR}*.SJ.out.tab \
      --outFilterType BySJout \
      --outReadsUnmapped Fastx

rm ${DIR_SAMPLE_NAME}/${SAMPLE_NAME}.Aligned.out.bam


NWCMD="echo ${SAMPLE_NAME} >> ${BASEDIR}/job_star2ndpass.txt"
eval "${NWCMD}"


# Done 
echo
echo "********************************************************"
echo "["`date`"] Done"
echo "********************************************************"
exit 0
