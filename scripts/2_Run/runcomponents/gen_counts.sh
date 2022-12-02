#!/bin/bash
#SBATCH -J GenerateCounts

#SBATCH -n 3
#SBATCH --mem=8G
#SBATCH -t 1:00:00

#SBATCH -o output/counts.out

#SBATCH --mail-type=END
#SBATCH --mail-user=pranav_mahableshwarkar@brown.edu

module load samtools/1.13 
module load gcc/8.3
module load subread/1.6.2

if [ $# -ne 4 ]; then
	echo $0: "Usage: 'sbatch gen_counts.sh' or 'bash gen_counts.sh'
            1) Data Input Directory (created in createfolderstruct.sh)
            2) Results Output Directory (created in createfolderstruct.sh)
            3) Reference Directory
            4) Output Name for Counts Table."
        exit 1
fi

echo "The Time2Splice folder is assumed to be in the same folder as this script."

# Variable Assignment from Input Parameters. 
DATA_DIR=$1
RES_DIR=$2
REF_DIR=$3
OUTNAME=$4

# Dmel references the genome of D Melanogaster that is needed for BowTie2 Alignment. 
DMELGTF="${REF_DIR}/dmel-all-r6.46.gtf"
if [ ! -e $DMELGTF ] 
then
    echo $0:"${DMELGTF} not found. Please make sure file path is correct."
    exit 1
fi
echo "Reference file found. Continuing..."

# 5. Run featureCounts to generate WITHOUT THE SBATCH ARRAY JOB!
featureCounts -t gene -C -T 8 \
  -a ${DMELGTF} \
  -o ${DATA_DIR}/count/${OUTNAME}_featurecounts.txt \
  ${RES_DIR}/preprocess/alignment/*/out.sorted.bam