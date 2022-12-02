### Process Analyze RNA-Seq Data and Produce Gene Ontology Analysis

###  Author: Pranav Mahableshwarkar 
###  Justin Currie's DESeq2 Script and Ashley Conard's Time2Splice are implemented in this pipeline!

1. Package Installation:
# Oscar Modules
This package installation guide is largely meant for Oscar usage as the installation of R packages is a little
more complex! The following modules are necessary for the R package installation shell script:

module load miniconda/4.10
module load R/4.2.0
module load gcc/10.2 pcre2/10.35 intel/2020.2 texlive/2018

You can add these to the .modules file in your home directory of Oscar so that they are loaded automatically
when you open a new Oscar terminal/access the login node. 

# R Packages
Once all of the above modules are installed, you have to enter RStudio (in the login node) and run the package
installations that are documented in the rinstall.txt file. 

This will take a while, but it will install all of the necessary R packages. For some reason, this takes much longer on Oscar but once this is done, you will not have to do it again!

When you run the installation, you will be prompted three times. The first two prompts will ask you about local package installation - enter "y" for both (you will never have to do this again). Then it will ask you to select a CRAN mirror. Choose the mirror for US - [OH] - which was number 76 for me.

# R Packages
There is a conda environment for the R packages but for the timebeing, the above R package installation is more appropriate as there are no python packages needed!

2. Folder Setup
Go to scripts/1_setup and run the folder structure create script. The specific instructions are provided in the shell script. 

If you need to upload FastQ files to Oscar, you can use the ii_importFQ.sh file to handle that. 

3. Run the Analysis
# Part One: Preprocess
Following the rules outlined in the preprocess.sh file you can run the preprocess script. The preprocess.sh script is designed to automatically make CountsTables for 16 total samples (4 replicates for 4 sample types maximum). However, the individual scripts for each step of the preprocessing can be found in 
2_Run/runcomponents. 

# Part Two: Analysis
Once the CountsTable is made, you will have to edit the CountsTable produced before it can be run in analysis. 
For some reason, featurecounts keeps the entire filepath for each .bam file. Remove the extraneous information so the headings look something like this:

Geneid	Chr	Start	End	Strand	Length  KC_CLAMP_1  KC_CLAMP_2	KC_CLAMP_3	KC_GFP_1    KC_GFP_2	KC_GFP_3	S2_CLAMP_1	S2_CLAMP_2	S2_CLAMP_3	S2_GFP_1	S2_GFP_2	S2_GFP_3

Now that we have it cleaner, there is one more change that must be made for the script to recognize the mutant vs control groups. The script will look at the first x columns where x is the number of control samples and designate them the control. Our script alphabetically orders the sample columns - so the control and the mutant groups must have the same prefix (regardless of sex as this script does not do sex-specific analysis).
This change would look something like this:

Geneid  Chr Start   End Strand  Length  RNAiCLAMP_1  RNAiCLAMP_2   RNAiCLAMP_3  Control_1    Control_2   Control_3     RNAiCLAMP_4   RNAiCLAMP_5    RNAiCLAMP_6  Control_4   Control_5   Control_6

Notice how alphabetically the Control Samples would come first once the sorting is performed! Also note that
the sex-differentiation in the naming is gone as both sexes are merged for each sample type in DESeq2.

## EXAMPLE:
sbatch 2_analysis.sh /users/pmahable/data/pmahable/rna_analysis/CLAMPRNAi_CellLine/data
\ /users/pmahable/data/pmahable/rna_analysis/CLAMPRNAi_CellLine/results RNAiCLAMP 6 6 
\ /users/pmahable/data/pmahable/rna_analysis/CLAMPRNAi_CellLine/data/count/CellLineCLAMPRNAi_featurecounts.txt

Then, you are good to go! Run the analysis.sh script to complete the analysis and everything will end up in its appropriate folders. 


