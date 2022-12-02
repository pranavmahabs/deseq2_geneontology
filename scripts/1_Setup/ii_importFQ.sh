#!/bin/bash

#SBATCH -J LoadSRAsDump

#SBATCH -n 4
#SBATCH --mem=8G
#SBATCH -t 3:00:00

if [ $# -ne 1 ]; then
	echo $0: "Usage: 'bash createfolderstruct.sh'
            1) Directory to store all data and results. (i.e. ~/CLAMP_RNAi_Experiment)"
        exit 1
fi

#SBATCH -o output/output.out
# SBATCH  --array=1-2
module load sratoolkit/2.11.0

# wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR8895530/SRR8895530
# wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos3/sra-pub-run-20/SRR8895531/SRR8895531.1

# wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR8895523/SRR8895523
# wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos3/sra-pub-run-21/SRR8895524/SRR8895524.1

# wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR8895523/SRR8895523
# wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos3/sra-pub-run-21/SRR8895524/SRR8895524.1


fastq-dump -I --split-files --gzip *SRR*