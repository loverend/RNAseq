#\!/bin/bash
#$ -cwd -pe shmem 4 -V -N index38
#$ -P jknight.prjc -q short.qc

echo Start analysis at: `date`

# mkdir -p data/hg19/sequence
#cd data/hg19/sequence/
#wget ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.{1..22}.fa.gz
# wget ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.{MT,X,Y}.fa.gz
# gunzip -c Homo_sapiens.GRCh38.dna.chromosome.* > GRCh38_r92.all.fa
# cd ../../../
# 
# mkdir -p data/GRCh38/annotation
# cd data/GRCh38/annotation/
# wget ftp://ftp.ensembl.org/pub/release-92/gtf/homo_sapiens/Homo_sapiens.GRCh38.92.gtf.gz
# gunzip Homo_sapiens.GRCh38.92.gtf.gz
# cd ../../../

mkdir -p /well/jknight/RNASeqSTARPipeline/data/GRCh38/star_indices_overhang74
/well/jknight/RNASeqSTARPipeline/STAR-2.7.3a/bin/Linux_x86_64/STAR --runThreadN 4 --runMode genomeGenerate --genomeDir /well/jknight/RNASeqSTARPipeline/data/GRCh38/star_indices_overhang74 --genomeFastaFiles /well/jknight/RNASeqSTARPipeline/data/GRCh38/sequence/GRCh38_r92.all.fa

echo Finished analysis at: `date`
