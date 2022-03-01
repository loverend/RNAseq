#!/bin/bash
#$ -cwd
#$ -pe shmem 4 -N index
#$ -q short.qc 

echo "------------------------------------------------"
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "------------------------------------------------"

mkdir -p /well/immune-rep/shared/REFERENCES/LAUREN_REFERENCES/GRCh38/star_indices
module load STAR/2.7.3a-GCC-8.3.0

STAR --runThreadN 4 --runMode genomeGenerate --genomeDir /well/immune-rep/shared/REFERENCES/LAUREN_REFERENCES/GRCh38/star_indices --genomeFastaFiles /well/immune-rep/shared/CODE/RNASeqSTARPipeline/data/GRCh38/sequence/Homo_sapiens.GRCh38.dna.primary_assembly.fa --sjdbGTFfile /well/immune-rep/shared/CODE/RNASeqSTARPipeline/data/GRCh38/annotation/Homo_sapiens.GRCh38.100.gtf --sjdbOverhang 99

echo "Finished at: "`date`
