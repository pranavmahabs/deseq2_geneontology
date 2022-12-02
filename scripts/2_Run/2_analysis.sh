#!/bin/bash

#SBATCH -J DeSeqGO

#SBATCH -n 4
#SBATCH --mem=8G
#SBATCH -t 1:00:00

#SBATCH -o output/analysis.out

#SBATCH --mail-type=END
#SBATCH --mail-user=pranav_mahableshwarkar@brown.edu

#######################
# rundeseq_go.sh
# Pranav Mahableshwarkar
# Last Mod. 9/11/22
# Purpose: Take featurecounts Count Matrices and produce DESeq2 and Gene Ontology results.
#######################

module load R/4.2.0
module load gcc/10.2 pcre2/10.35 intel/2020.2 texlive/2018



# run: "conda info | grep -i 'base environment'" to get the path for everything before /etc below.
# source /gpfs/runtime/opt/miniconda/4.10/etc/profile.d/conda.sh
# conda activate deseq_gopipe

if [ $# -ne 5 ]; then
	echo $0: "Usage: 'sbatch preproc_countTab.sh' or 'bash preproc_countTab.sh'
            2) Results Output Directory (created in createfolderstruct.sh)
            3) Experimental Condition
            4) Number of Control Samples
            5) Number of Mutant Samples
            6) Counts Table"
        exit 1
fi

RES_DIR=$1
OUTNAME=$2
NUM_CON=$3
NUM_MUT=$4
MATRIX=$5

# Run the RScript
# Job Three: Run DESeq2 and Gene Ontology (GProfiler2) 
Rscript runcomponents/rundeseq_go.r -r ${MATRIX} -c ${NUM_CON} -e ${NUM_MUT} -g null -p ${OUTNAME} -o ${RES_DIR}/analysis