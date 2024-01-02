#!/bin/bash
#$ -cwd
#$ -pe shmem 4
#$  -N index149
#$ -P jknight.prjc -q short.qc

echo "------------------------------------------------"
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "------------------------------------------------"

mkdir -p /well/jknight/RNASeqSTARPipeline/data/GRCh38/star_indices_overhang149

/well/jknight/RNASeqSTARPipeline/STAR-2.7.3a/bin/Linux_x86_64/STAR --runThreadN 4 --runMode genomeGenerate --genomeDir /well/jknight/RNASeqSTARPipeline/data/GRCh38/star_indices_overhang149 --genomeFastaFiles /well/jknight/RNASeqSTARPipeline/data/GRCh38/sequence/Homo_sapiens.GRCh38.dna.primary_assembly.fa --sjdbGTFfile /well/jknight/RNASeqSTARPipeline/data/GRCh38/annotation/Homo_sapiens.GRCh38.100.gtf --sjdbOverhang 149

echo "Finished at: "`date`

