#\!/bin/bash
#$ -cwd -pe shmem 1 -V 
#$ -N mappingQC
#$ -P jknight.prjc -q short.qc
#$ -o RNASeQC.stdout
#$ -e RNASeQC.sterr
#$ -t 1:223:1

IND=$(sed "$SGE_TASK_ID"'q;d' sample_list.txt)
echo $IND

/well/jknight/RNASeqSTARPipeline/rnaseqc/rnaseqc.v2.3.5.linux /well/jknight/RNASeqSTARPipeline/data/GRCh38/annotation/gencode.v7.annotation_goodContig.gtf Mapping/$IND/$IND.Aligned.sortedByCoord.out.bam RNASeQC/


