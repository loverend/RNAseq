#!/bin/bash

###################################################
#
#  SETUP PATHS TO TOOLS AND REFERENCE
#
###################################################

export PATH=$PATH:/apps/well/samtools/1.2/bin
SAMTOOLS="/apps/well/samtools/1.2/bin/samtools"
export PATH=$PATH:/apps/well/bedtools/2.24.0-18-gb0bc5b7/bin/
BAMTOOLS="/apps/well/bamtools/2.3.0/bin/bamtools"
export PATH=$PATH:/apps/well/bamtools/2.3.0/bin
export PATH=$PATH:/apps/well/bamtools/2.3.0/lib/bamtools
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/apps/well/bamtools/2.3.0/lib/bamtools/
MARKDUP="/apps/htseq/gatk-4.0.3.0/gatk"


###################################################
#
# PROJECT DIRECTORY AND INPUT DATA
#
###################################################

INPUT_DIR=$BASEDIR"/Mapping/"

OUTPUT_DIR=$BASEDIR"/Filtering/"

SAMPLE_NAME=$1


###################################################
#
# FILTERING
#
###################################################

DIR_SAMPLE_NAME=$OUTPUT_DIR$SAMPLE_NAME

FILTERING(){
  
  echo
	echo "echo Start filtering at: \`date\`"
  echo

	tmp="$SAMTOOLS sort $INPUT_DIR$SAMPLE_NAME/${SAMPLE_NAME}Aligned.out.bam $INPUT_DIR$SAMPLE_NAME/${SAMPLE_NAME}Aligned.sortedByCoord.out"
	echo $tmp

	tmp="rm $INPUT_DIR$SAMPLE_NAME/${SAMPLE_NAME}Aligned.out.bam"
	echo $tmp

	tmp="$SAMTOOLS index $INPUT_DIR$SAMPLE_NAME/${SAMPLE_NAME}Aligned.sortedByCoord.out.bam"
	echo $tmp

	tmp="$MARKDUP --java-options -Xmx8g MarkDuplicates -I $INPUT_DIR$SAMPLE_NAME/${SAMPLE_NAME}Aligned.sortedByCoord.out.bam -O $DIR_SAMPLE_NAME/$SAMPLE_NAME.dups.marked.bam -M $DIR_SAMPLE_NAME/$SAMPLE_NAME.dup_metrix.txt"
	echo $tmp
	
	tmp="$MARKDUP --java-options -Xmx8g MarkDuplicates -I $INPUT_DIR$SAMPLE_NAME/${SAMPLE_NAME}Aligned.sortedByCoord.out.bam -O $DIR_SAMPLE_NAME/$SAMPLE_NAME.nodup.bam -REMOVE_SEQUENCING_DUPLICATES true -M $DIR_SAMPLE_NAME/$SAMPLE_NAME.dup_metrix.txt"
	echo $tmp

	tmp="$SAMTOOLS sort -n $DIR_SAMPLE_NAME/$SAMPLE_NAME.$FORCOUNTS.bam $DIR_SAMPLE_NAME/${SAMPLE_NAME}.$FORCOUNTS.namesorted"
        echo $tmp

	echo "echo FILTERING DONE."
}

###################################################
#
# ENVIRONMENT SETUP
#
###################################################

ENVIRONMENT_SETUP(){

        echo "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/apps/well/bamtools/2.3.0/lib/bamtools/"
        echo "echo "

}


###################################################
#
# WRITE A SCRIPT
#
###################################################

WRITE_A_SCRIPT(){

        echo "#\!/bin/bash"
        echo "#$ -cwd -pe shmem 1 -V -N filter."$1
        echo "#$ -hold_jid map."$1
        echo "#$ -q short.qc"
        echo "#$ -o $DIR_SAMPLE_NAME/$1.stdout"
        echo "#$ -e $DIR_SAMPLE_NAME/$1.sterr"
	
        
	ENVIRONMENT_SETUP
	FILTERING $1
	
        echo "echo Finished analysis at: \`date\`"

}

WRITE_A_SCRIPT $1
