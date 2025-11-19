#!/bin/bash
#SBATCH --time=48:00:00 # cpu partition max time per kcl policy
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 # could specify more here if wanted 
#SBATCH --mem=1G # shouldnt require more than a gb per file
#SBATCH --hint=nomultithread # prefer physical cores for better throughput
#SBATCH --job-name=fastqc
#SBATCH --array=0-29 # modify based on number of files
#SBATCH --output=/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/050_deduplicated/logs/051_fastqc_%A_%a.out
#SBATCH --error=/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/050_deduplicated/logs/051_fastqc_%A_%a.err

#---------#
# purpose #
#---------#

# this files examines quality of raw fastq files prior to merging of technical duplicates

#--------------------#
# set up environment #
#--------------------#

# activate conda env
module load fastqc/0.12.1-gcc-13.2.0

# input directory of fastq files
IN_DIR="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/050_deduplicated"

# output directory for fastqc files
OUT_DIR="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/050_deduplicated/051_fastqc"

mkdir -p ${OUT_DIR}

#----------------------------------------#
# get list of files for fastqc and array #
#----------------------------------------#

SAMPLES=($(find "${IN_DIR}" -maxdepth 1 -type f -name "*_picard_dedup.sam"))

# specify sample based on array number 
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

#---------------------------#
# fastqc on raw fastq files #
#---------------------------#

# (Optional) force headless Java; harmless if inputs are valid
export JAVA_TOOL_OPTIONS="-Djava.awt.headless=true"

fastqc \
    -o "${OUT_DIR}" \
    "${SAMPLE}"
