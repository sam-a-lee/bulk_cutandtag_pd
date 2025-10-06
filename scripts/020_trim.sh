#!/bin/bash
#SBATCH --time=48:00:00 # cpu partition max time
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --job-name=cutadapt
#SBATCH --array=0-29 # !!! modify to number of unique samples being processing !!!
#SBATCH --output=/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/020_trimmed/logs/020_cutadapt_%A_%a.out
#SBATCH --error=/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/020_trimmed/logs/020_cutadapt_%A_%a.err

#---------# 
# purpose #
#---------#

# this file uses cutadapt to quality trim and remove illumina adaptors from raw sequencing reads (fastq files)

#--------------------# 
# set up environment #
#--------------------#

# activate cutadapt
module load py-cutadapt/4.4-gcc-13.2.0-python-3.11.6

# root dir where raw fastq files are stored
ROOT_DIR="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/010_merged"

# directory to save files
OUT_DIR="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/020_trimmed"

mkdir -p $OUT_DIR

#---------------------------------------------------# 
# create list of files corresponding to each sample #
#---------------------------------------------------#

# get sorted list of unique raw samples (R1 files)
SAMPLES=($(find "${ROOT_DIR}" -maxdepth 1 -type f -name "*_R1_001.fastq.gz" | sort))

# pick sample based on job array index
R1=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
R2="${R1/_R1_001.fastq.gz/_R2_001_fastq.gz}"

BASE=$(basename "${R1}" _R1_001.fastq.gz)

#---------------------------------# 
# trim adapters using trim galore #
#---------------------------------#

echo "Processing sample: ${R1} and ${R2}"

# for consistency with cell ranger
# no quality trimming
# min lenght of 30 bp
# trim tn5 (CTGTCTCTTATACACATCT) and illumina univesal adapters (AGATCGGAAGAGC)
# max length of 50 since 4 samples seq with 75 bp instead of 50 bp and read length can affect peak identificaiton 
cutadapt -m 30 -l 50 -n 2 \
  -a AGATCGGAAGAGC \
  -A AGATCGGAAGAGC \
  -a CTGTCTCTTATACACATCT \
  -A CTGTCTCTTATACACATCT \
  -o "${OUT_DIR}/${BASE}_R1_001_cutadapt.fastq.gz" \
  -p "${OUT_DIR}/${BASE}_R2_001_cutadapt.fastq.gz" \
  ${R1} ${R2}
  