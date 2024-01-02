#!/bin/bash

###################################################
#
#  SETUP PATHS TO TOOLS AND REFERENCE
#
###################################################

export STAR=/well/immune-rep/shared/CODE/RNASeqSTARPipeline/STAR-2.7.3a/bin/Linux_x86_64/STAR

###################################################
#
# PROJECT DIRECTORY AND INPUT DATA
#
###################################################

OUTPUT_DIR=$BASEDIR"/Mapping/"

SAMPLE_NAME=$1

INPUT_READ1=$2
INPUT_READ2=$3

###################################################
#
# STAR MAPPING
#
###################################################

DIR_SAMPLE_NAME=$OUTPUT_DIR$SAMPLE_NAME

STAR_MAPPING(){
        echo
	      echo "echo Start mapping at: \`date\`"
        echo

	tmp="$STAR --genomeDir $GENOME \
	          --runThreadN 12 \
	          --readFilesIn $INPUT_READ1 $INPUT_READ2 \
		  --readFilesCommand gunzip -c \
		  --outFileNamePrefix $DIR_SAMPLE_NAME"/"$SAMPLE_NAME \
	          --outSAMtype BAM Unsorted \
	          --outFilterMismatchNoverLmax 0.04 \
	        --limitOutSJcollapsed 2000000  
		--outSAMprimaryFlag OneBestScore \
	          2> $DIR_SAMPLE_NAME/stderr.$SAMPLE_NAME.STAR.txt"
	echo $tmp

	echo "echo MAPPING DONE."
}



###################################################
#
# WRITE THE SCRIPT
#
###################################################

WRITE_A_SCRIPT(){

        echo "#\!/bin/bash"
        echo "#$ -cwd -pe shmem 16 -V -N map."$1
        echo "#$ -q short.qc"
        echo "#$ -o $DIR_SAMPLE_NAME/$1.stdout"
        echo "#$ -e $DIR_SAMPLE_NAME/$1.sterr"
	echo ""

	STAR_MAPPING $1 $2 $3
	
        echo "echo Finished analysis at: \`date\`"

}

WRITE_A_SCRIPT $1 $2 $3
