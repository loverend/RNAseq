#\!/bin/bash
#$ -cwd -pe shmem 1 -V -N index
#$ -P jknight.prjc -q short.qc

echo Start analysis at: `date`

cd /well/jknight/RNASeqSTARPipeline/data/GRCh38/sequence

java -jar /apps/well/picard-tools/1.111/CreateSequenceDictionary.jar R=/well/jknight/RNASeqSTARPipeline/data/GRCh38/sequence/GRCh38_r92.all.fa O=/well/jknight/RNASeqSTARPipeline/data/GRCh38/sequence/GRCh38_r92.all.dict

echo Finished analysis at: `date`
