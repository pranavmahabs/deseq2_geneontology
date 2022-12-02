#!/bin/bash

#######################
# pipeline.sh
# Pranav Mahableshwarkar
# Last Mod. 9/10/22
# Purpose: Take fastq files and produce DESeq2 and Gene Ontology results from produced counts tables. 
#######################

if [ $# -ne 4 ]; then
	echo $0: "Usage: './1_preprocess.sh 1 2 3 4'
            1) Data Input Directory (created in createfolderstruct.sh)
            2) Results Output Directory (created in createfolderstruct.sh)
            3) Reference Directory
            4) Experimental Condition (Will be used for Output Naming)
            
            FASTQ folder should ONLY contain directories for each sample (which contain the fastq files).
            (i.e.: ~/fastq/... where '...' is the UNIQUE name for each replicate."
        exit 1
fi

# run: "conda info | grep -i 'base environment'" to get the path for everything before /etc below.
source /gpfs/runtime/opt/miniconda/4.10/etc/profile.d/conda.sh
conda activate deseq2pipeline

# Variable Assignment from Input Parameters. 
DATA_DIR=$1
RES_DIR=$2
REF_DIR=$3
OUTNAME=$4
FASTQF="${DATA_DIR}/fastq"

echo "Make sure that you have edited preprocess.sh in runcomponents to note the number of samples you are running."

# Job One: Run the Preprocessing Scripts
# TODO: Go to preprocess.sh and change the array size to the number of samples being analyzed. 
jobID_1=$(sbatch runcomponents/preprocess.sh $DATA_DIR $RES_DIR $REF_DIR $OUTNAME | cut -f 4 -d' ')

# Job Two: Run the fastcounts Script to Generate Count Tables
sbatch --dependency=afterok:$jobID_1 runcomponents/gen_counts.sh $DATA_DIR $RES_DIR $REF_DIR $OUTNAME

conda deactivate

