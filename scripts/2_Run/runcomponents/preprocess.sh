#!/bin/bash
#SBATCH -J Preprocess_ArrayJobs

#SBATCH -n 2
#SBATCH --mem=16G
#SBATCH -t 24:00:00

# Create an array SBATCH Job for all files
### ONLY CHANGE THE FOLLOWING LINE IN THIS FILE FOR THE NUMBER OF SAMPLES (1-x) ####
#SBATCH  --array=1-16

#SBATCH -o output/output-%a.out

# Modules needed for preprocessing the data. 
module load fastqc/0.11.5 
module load trimgalore/0.5.0
module load cutadapt/1.14 
module load bowtie2/2.4.2 
module load samtools/1.13 
module load gcc/8.3
module load subread/1.6.2

if [ $# -ne 4 ]; then
	echo $0: "Usage: 'sbatch preproc_countTab.sh' or 'bash preproc_countTab.sh'
            1) Data Input Directory (created in createfolderstruct.sh)
            2) Results Output Directory (created in createfolderstruct.sh)
            3) Reference Directory
            4) Output Name for Counts Table.
            
            FASTQ folder should ONLY contain directories for each sample (which contain the fastq files).
            (i.e.: ~/fastq/... where '...' is the UNIQUE name for each replicate.

            Please open this file to change the end of line 10 to be 1-x where x is the number of replicates."
        exit 1
fi

echo "The Time2Splice folder is assumed to be in the same folder as this script."

# Variable Assignment from Input Parameters. 
DATA_DIR=$1
RES_DIR=$2
REF_DIR=$3
T2SPL="time2splice/scripts"
OUTNAME=$4
FASTQF="${DATA_DIR}/fastq"

# BOWTIE references the genome of D Melanogaster that is needed for BowTie2 Alignment. 
BOWTIE="${REF_DIR}/bowtie2_index/dmel_all_chromosome"
if [ -d "${REF_DIR}/bowtie2_index/" ] 
then
    echo "dmel_all_chromsome found in reference folder. Continuing..."
else 
    echo $0:"${BOWTIE} NOT found. There should be 6 files in the bowtie2 index in the *.bt2 format.
            Is your dmel_all_chromsome in '${REF_DIR}/reference/bowtie2_index/'?"
    exit 1
fi

# Array of all directories being run at a given moment. (MUST BE IN SAME FOLDER!)
readarray -t ALLDATA < <(find ${FASTQF}/. -maxdepth 1 -type d -printf '%P\n')
declare -p ALLDATA
echo "These are the samples being checked. The first one is meant to be null/empty. Please check."

# IMPORTANT: SAMPLE NAMES SHOULD NOT CONTAIN '_R1' or '_R2'. This will cause ERRORS in trim galore.
for SAMPLE in "${ALLDATA[@]}"
do
    if [[ "$SAMPLE" == *"_R"* ]]; then
        echo "Sample title including '_R*' was detected. Please change the naming schema."
        exit 1
    fi
done

# Setting up the individaul sample for each Array Job in SLURM. 
DIREC=${ALLDATA[$SLURM_ARRAY_TASK_ID]}
echo "Preparations All Made. Starting FastQC..."

# 1. Run the FASTQC on all files in all folders and make new folders in results.
./${T2SPL}/preprocess/2_run_fastQC.sh 4 $DATA_DIR/fastq/${DIREC} $RES_DIR/preprocess/fastqc/${DIREC}
echo "Completed FastQC. Starting Trim_Galore..."

# 2. Run TrimGalore 
./${T2SPL}/preprocess/3_run_trim_galore.sh 4 $DATA_DIR/fastq/${DIREC} $RES_DIR/preprocess/trim_galore_fastqc/${DIREC} 30 cutadapt 0
echo "Completed Trim_Galore. Starting BowTie2..."

# 3. Run BowTie2 to make bam files for featureCounts.   
./${T2SPL}/preprocess/4_run_Bowtie2.sh ${RES_DIR}/preprocess/trim_galore_fastqc/${DIREC} ${RES_DIR}/preprocess/alignment/${DIREC} 1 ${BOWTIE} dme 0
echo "completed BowTie2. Starting featureCounts..."

