#\!/bin/bash
#$ -cwd -pe shmem 10 -V -N fastQC
#$ -P jknight.prjc -q short.qc
#$ -o /well/jknight/michael/rnaseq_pandit/PANDIT_PERFORM_STAR/QC/fastQC.stdout
#$ -e /well/jknight/michael/rnaseq_pandit/PANDIT_PERFORM_STAR/QC/fastQC.sterr

# Your base directory
cd /well/jknight/michael/rnaseq_pandit/PANDIT_PERFORM_STAR/

module load java/1.8.0_latest

# Change output path, location of fastq files as appropriate
/apps/htseq/FastQC/fastqc -t 10 -o /well/jknight/michael/rnaseq_pandit/PANDIT_PERFORM_STAR/QC/ /well/jknight/michael/rnaseq_pandit/PANDIT_PERFORM_STAR/fastq_data/*.fastq.gz > /well/jknight/michael/rnaseq_pandit/PANDIT_PERFORM_STAR/QC/fastQC.log

