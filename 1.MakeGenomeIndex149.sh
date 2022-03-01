#!/bin/bash
#$ -cwd
#$ -pe shmem 4 -N index100
#$ -q short.qc -P combat.prjc

echo "------------------------------------------------"
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "------------------------------------------------"

mkdir -p /well/combat/projects/rnaseq/P200133/mapping/GRCh38/star_indices
module load STAR/2.7.3a-GCC-8.3.0

STAR --runThreadN 4 --runMode genomeGenerate --genomeDir /well/combat/projects/rnaseq/P200133/mapping/GRCh38/star_indices_overhang149 --genomeFastaFiles /well/combat/projects/rnaseq/P200133/mapping/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa --sjdbGTFfile /well/combat/projects/rnaseq/P200133/mapping/GRCh38/Homo_sapiens.GRCh38.100.gtf --sjdbOverhang 99

echo "Finished at: "`date`
