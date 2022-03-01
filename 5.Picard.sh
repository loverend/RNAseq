#!/bin/bash
#$ -cwd
#$ -pe shmem 8 -N markdups
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
QC_DIR=$4
GTF=$5
GENOME=$6
BASEDIR=$7

# Set parameters
module load SAMtools/1.9-foss-2018b
MARKDUP=/apps/htseq/gatk-4.0.3.0/gatk

# Get sample ID for this task
SAMPLE_NAME=$(cat ${KEY} | tail -n+${SGE_TASK_ID} | head -1 | cut -f1 )
echo ${SAMPLE_NAME}

# Make sample directory
DIR_SAMPLE_NAME=${QC_DIR}${SAMPLE_NAME}
if [[ ! -d ${DIR_SAMPLE_NAME} ]]; then
              echo "Creating QC Sample directory"
			  mkdir -p ${DIR_SAMPLE_NAME}
fi

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
echo "Running stage 5. PICARD / SAMTOOLS INDEXING"

# Index bam file
samtools index ${MAPPING_DIR}${SAMPLE_NAME}/${SAMPLE_NAME}.Aligned.sortedByCoord.out.bam

# Quantify duplicates
# $MARKDUP --java-options -Xmx8g MarkDuplicates -I $MAPPING_DIR$SAMPLE_NAME/${SAMPLE_NAME}.Aligned.sortedByCoord.out.bam -O $DIR_SAMPLE_NAME/$SAMPLE_NAME.dups.marked.bam -M $DIR_SAMPLE_NAME/${SAMPLE_NAME}.dup_metrics.txt

# grep Unknown $DIR_SAMPLE_NAME/${SAMPLE_NAME}.dup_metrics.txt > $DIR_SAMPLE_NAME/$SAMPLE_NAME.dup.counts.txt
# rm $DIR_SAMPLE_NAME/$SAMPLE_NAME.dups.marked.bam
# rm $DIR_SAMPLE_NAME/$SAMPLE_NAME.dup_metrics.txt

NWCMD="echo ${SAMPLE_NAME} >> ${BASEDIR}/job_markdups.txt"
eval "${NWCMD}"


# Done 
echo
echo "********************************************************"
echo "["`date`"] Done"
echo "********************************************************"
exit 0

