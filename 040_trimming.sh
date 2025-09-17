#!/bin/bash
#SBATCH --time=48:00:00 # cpu partition max time
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --job-name=trim_adapters
#SBATCH --array=0-30 # !!! modify to number of unique samples being processing !!!
#SBATCH --output=/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/04_trimmed/logs/trim_adapters_%A_%a.out
#SBATCH --error=/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/04_trimmed/logs/trim_adapters_%A_%a.err

#---------# 
# purpose #
#---------#

# this file uses trimgalore to quality trim and remove illumina adaptors from raw sequencing reads (fastq files)
# trim galore produces read QC reports after trimming 

#--------------------# 
# set up environment #
#--------------------#

# change dir to where conda envs are 
cd /users/k2587336

# load shell
source ~/.bashrc

# activate conda env
source activate trim_galore

# root dir where raw fastq files are stored
ROOT_DIR="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/03_merged_fastq"

# directory to save files
OUT_DIR="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/04_trimmed"

#---------------------------------------------------# 
# create list of files corresponding to each sample #
#---------------------------------------------------#

# get sorted list of unique raw samples (R1 files)
SAMPLES=($(find "${ROOT_DIR}" -maxdepth 1 -type f -name "*_R1_001.fastq.gz" | sort))

# pick sample based on job array index
R1=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
R2="${R1/_R1_001.fastq.gz/_R2_001.fastq.gz}"

#---------------------------------# 
# trim adapters using trim galore #
#---------------------------------#

echo "Processing sample: ${R1} and ${R2}"

# run trim galore
# default trimming below phred 20
trim_galore --gzip \
    --paired \
    --fastqc \
    --output_dir ${OUT_DIR} \
    ${R1} \
    ${R2}