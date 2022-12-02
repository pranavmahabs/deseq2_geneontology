#!/bin/bash

if [ $# -ne 2 ]; then
	echo $0: "Usage: 'bash createfolderstruct.sh'
            1) Directory to store all data and results. (i.e. ~/CLAMP_RNAi_Experiment)
            2) Experiment Name"
        exit 1
fi

INDIR=$1
EXPNAME=$2

echo "Making ${EXPNAME} folder within the directory."
mkdir "${INDIR}/${EXPNAME}"

echo "Making data and results folders."
echo "Please put all input fastq files into inputdir/${EXPNAME}/data/fastq."
echo "These fastq files should be in folders labeled with replicate name (no _R* in naming)."
mkdir "${INDIR}/${EXPNAME}/data"
mkdir "${INDIR}/${EXPNAME}/results"
mkdir "${INDIR}/${EXPNAME}/data/fastq"
mkdir "${INDIR}/${EXPNAME}/data/count"
mkdir "${INDIR}/${EXPNAME}/results/preprocess"
mkdir "${INDIR}/${EXPNAME}/results/analysis"

mkdir "${INDIR}/${EXPNAME}/results/preprocess/fastqc"
mkdir "${INDIR}/${EXPNAME}/results/preprocess/trim_galore_fastqc"
mkdir "${INDIR}/${EXPNAME}/results/preprocess/alignment"

echo "Prior to running the scripts, make sure to format all data correctly. Please place all 
      reference files in a reference folder. The working file path for that directory is 
      needed for the scripts."

echo "createfolderstruct.sh Done."
