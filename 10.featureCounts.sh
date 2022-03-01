#!/bin/bash
#$ -cwd
#$ -pe shmem 6 -N counts
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
MAPPING_DIR=$2
KEY=$3
GTF=$4
PAIRED=$5
BASEDIR=$6

# Set parameters
module load Subread/1.6.4-foss-2018b

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
echo "PAIRED          : ${PAIRED}"
echo "********************************************************"
echo "Running stage 10. FEATURE COUNTS"

# Make sample directory
DIR_SAMPLE_NAME=${COUNT_DIR}${SAMPLE_NAME}
if [[ ! -d ${DIR_SAMPLE_NAME} ]]; then
    echo "Creating Directory"
	mkdir -p ${DIR_SAMPLE_NAME}
fi

# Add path to config file
if [[ ${PAIRED} == "TRUE" || -z "$PAIRED" ]]; then
featureCounts -T 4 -a ${GTF} -g gene_id \
      -o ${DIR_SAMPLE_NAME}/${SAMPLE_NAME}.counts.txt -p \
      -s 2 ${MAPPING_DIR}${SAMPLE_NAME}/${SAMPLE_NAME}.Aligned.sortedByCoord.out.bam
cut -f 1,7 ${DIR_SAMPLE_NAME}/${SAMPLE_NAME}.counts.txt > ${DIR_SAMPLE_NAME}/${SAMPLE_NAME}.reduced.counts.txt   #  This
mv ${DIR_SAMPLE_NAME}/${SAMPLE_NAME}.reduced.counts.txt ${DIR_SAMPLE_NAME}/${SAMPLE_NAME}.counts.txt             #  reduces the file size from ~ 30M to ~1M
fi 


#-p   If specified, fragments (or templates) will  be  counted  instead  of  reads.  This option is only applicable for paired-end reads.
			  			  
if [[ ${PAIRED} == "FALSE" ]]; then
featureCounts -T 9 -a ${GTF} -g gene_id \
      -o ${DIR_SAMPLE_NAME}/${SAMPLE_NAME}.counts.txt \
      -s 2 ${MAPPING_DIR}${SAMPLE_NAME}/${SAMPLE_NAME}.Aligned.sortedByCoord.out.bam
cut -f 1,7 ${DIR_SAMPLE_NAME}/${SAMPLE_NAME}.counts.txt > ${DIR_SAMPLE_NAME}/${SAMPLE_NAME}.reduced.counts.txt   #  This
mv ${DIR_SAMPLE_NAME}/${SAMPLE_NAME}.reduced.counts.txt ${DIR_SAMPLE_NAME}/${SAMPLE_NAME}.counts.txt             #  reduces the file size from ~ 30M to ~1M
fi 

NWCMD="echo ${SAMPLE_NAME} >> ${BASEDIR}/job_counting.txt"
eval "${NWCMD}"

# Done 
echo
echo "********************************************************"
echo "["`date`"] Done"
echo "********************************************************"
exit 0
