#!/bin/bash
## RNAseq Pipeline 
## Katie Burnham Sanger Institute - Cambridge
## Adapted Lauren Overend
## lauren.overend@live.co.uk

##################
# Job script arguments
BASEDIR=$1
STAGE=$2
KEY=$3
PAIRED=$4


######### Set parameters and paths
echo "Base directory: ${BASEDIR}"

## Locations of Files to USE
CODE_DIRECTORY=/well/immune-rep/shared/CODE/RNAseqPipeline/
RNASEQ_QC_DIR=/gpfs2/well/immune-rep/shared/CODE/rnaseqc/

##### make Required subdirectories
MAPPING_DIR=${BASEDIR}Mapping/
if [[ ! -d ${MAPPING_DIR} ]]; then
echo "Creating Mapping Directory"
mkdir ${MAPPING_DIR}
else 
echo "Directory Present"
fi

QC_DIR=${BASEDIR}QC/
if [[ ! -d ${QC_DIR} ]]; then
echo "Creating QC Directory"
mkdir ${QC_DIR}
else 
echo "Directory Present"
fi

COUNT_DIR=${BASEDIR}Counts/
if [[ ! -d ${COUNT_DIR} ]]; then
echo "Creating Counts Directory"
mkdir ${COUNT_DIR}
else 
echo "Directory Present"
fi


LOG_DIR=${BASEDIR}LOGFILES/
if [[ ! -d ${LOG_DIR} ]]; then
echo "Creating LOG file Directory"
mkdir ${LOG_DIR}
else 
echo "Directory Present"
fi


# Sample information file
echo "Sample information file: ${KEY}"

# Fastq files (or softlinks to them)
FASTQDIR=${BASEDIR}fastqs
echo "Fastq directory: ${FASTQDIR}"

######## Genome references
echo "The pipeline will use build 38"
GTF=/well/immune-rep/shared/REFERENCES/LAUREN_REFERENCES/DATA/Homo_sapiens.GRCh38.100.gtf
GTF2=/well/immune-rep/shared/REFERENCES/LAUREN_REFERENCES/DATA/Homo_sapiens.GRCh38.100.collapsed.gtf
GENOME=/well/immune-rep/shared/REFERENCES/LAUREN_REFERENCES/GRCh38/star_indices
echo "using GTF: ${GTF}"
echo "using STAR GENOME BUILD: ${GENOME}"

## How many samples?
TASKS=$(cat ${KEY} | wc -l)
TASKS=$((TASKS+1))
echo "No. Sampels: ${TASKS}"

#########################################################################
## Pipeline Stages! 
########################################################################

## IF STAGE =SK we will make a sample KEY!
if [[ "$STAGE" == "SK" ]]; then
echo "Making sample key"
module load R/3.6.2-foss-2019b
Rscript --vanilla ${CODE_DIRECTORY}0.MakeSampleKey.R ${KEY}
module purge
fi 

## IF STAGE = G we will make a star genome references
if [[ "$STAGE" == "G" ]]; then
echo "Making STAR index"
qsub -P ${CODE_DIRECTORY}1.MakeGenomeIndex.sh
fi 

## STAGE = R we will run the pipeline trim adaptor sequences
if [[ "$STAGE" == "R" ]]; then
	
	## Adaptor Trimming
	JOB_ID=$(qsub -t 1-${TASKS} -terse -e ${LOG_DIR} -o ${LOG_DIR} ${CODE_DIRECTORY}2.AdaptorTrimming.sh ${KEY} ${FASTQDIR} ${MAPPING_DIR} ${BASEDIR} ${PAIRED})
	echo "Submitting jobs for trimming"
	JOB_ID=${JOB_ID%.*}
	
	## Run STAR PASS 1 
	JOB_ID=$(qsub -hold_jid ${JOB_ID} -terse -e ${LOG_DIR} -o ${LOG_DIR} -t 1-${TASKS} ${CODE_DIRECTORY}3.STARfirstpass.sh ${KEY} ${FASTQDIR} ${MAPPING_DIR} ${GTF} ${GENOME} ${BASEDIR} ${PAIRED})
	echo "Submitting jobs for Star pass 1 "
	JOB_ID=${JOB_ID%.*}
	
	## Run STAR PASS 2
	JOB_ID=$(qsub -hold_jid ${JOB_ID} -terse -e ${LOG_DIR} -o ${LOG_DIR} -t 1-${TASKS}  ${CODE_DIRECTORY}4.STARsecondpass.sh ${KEY} ${FASTQDIR} ${MAPPING_DIR} ${GTF} ${GENOME} ${BASEDIR} ${PAIRED})
	echo "Submitting jobs for Star pass 2 "
	JOB_ID=${JOB_ID%.*}
	JOB_ID_2=${JOB_ID} 
	
	## Run SAMTOOLS INDEXING 
	echo "Submitting jobs for QC"
	JOB_ID=$(qsub -hold_jid ${JOB_ID} -terse -e ${LOG_DIR} -o ${LOG_DIR} -t 1-${TASKS}  ${CODE_DIRECTORY}5.Picard.sh ${KEY} ${FASTQDIR} ${MAPPING_DIR} ${QC_DIR} ${GTF} ${GENOME} ${BASEDIR})
	echo "Submitting jobs for INDEXING "
	JOB_ID=${JOB_ID%.*}
	
	## Collating STAR LOGS
	JOB_ID=$(qsub -hold_jid ${JOB_ID_2} -e ${LOG_DIR} -o ${LOG_DIR} -terse ${CODE_DIRECTORY}6.CollateSTARLogs.sh ${MAPPING_DIR} ${CODE_DIRECTORY})
	echo "Submitting jobs for COLLATING STAR LOGS "
	JOB_ID=${JOB_ID%.*}
	
	### Collating Duplicates - Not Needed 
	#qsub -hold_jid "markdups"  /well/immune-rep/shared/CODE/RNASeqSTARPipeline/SangerCombatUpdate/7.CollateDups.sh
	
	## QC RNAseq 
	JOB_ID=$(qsub -hold_jid ${JOB_ID_2} -terse -t 1-${TASKS} -e ${LOG_DIR} -o ${LOG_DIR} ${CODE_DIRECTORY}8.RNASeQC.sh ${QC_DIR} ${MAPPING_DIR} ${KEY} ${GTF2} ${CODE_DIRECTORY} ${BASEDIR} ${RNASEQ_QC_DIR})
	echo "Submitting jobs for QC"
	JOB_ID=${JOB_ID%.*}
	
	## COLLATE METRICS 
	JOB_ID=$(qsub -hold_jid ${JOB_ID} -terse -e ${LOG_DIR} -o ${LOG_DIR} ${CODE_DIRECTORY}9.CollateQCMetrics.sh ${QC_DIR} ${CODE_DIRECTORY})
	echo "Summarising QC_Metrics"
	
	## RNA seq counting
	JOB_ID=$(qsub -hold_jid ${JOB_ID_2} -terse -t 1-${TASKS} -e ${LOG_DIR} -o ${LOG_DIR} ${CODE_DIRECTORY}10.featureCounts.sh ${COUNT_DIR} ${MAPPING_DIR} ${KEY} ${GTF} ${PAIRED} ${BASEDIR})
	echo "Submitting jobs for counts"
	JOB_ID3=${JOB_ID%.*}
		
	## Merging Counts
	JOB_ID=$(qsub -hold_jid ${JOB_ID3} -terse -e ${LOG_DIR} -o ${LOG_DIR} ${CODE_DIRECTORY}11.MergeCounts.sh ${COUNT_DIR} ${CODE_DIRECTORY} ${MAPPING_DIR} ${KEY})
	echo "Summarising read counts"
	JOB_ID=${JOB_ID%.*}
fi 

