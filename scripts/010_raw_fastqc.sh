#!/bin/bash
#SBATCH --time=48:00:00 # cpu partition max time per kcl policy
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 # could specify more here if wanted 
#SBATCH --mem=1G # shouldnt require more than a gb per file
#SBATCH --hint=nomultithread # prefer physical cores for better throughput
#SBATCH --job-name=fastqc
#SBATCH --array=0-47 # modify based on number of files
#SBATCH --output=/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/010_fastqc/logs/010_fastqc_%A_%a.out
#SBATCH --error=/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/010_fastqc/logs/010_fastqc_%A_%a.err

#---------#
# purpose #
#---------#

# this files examines quality of raw fastq files prior to merging of technical duplicates

#--------------------#
# set up environment #
#--------------------#

# load fastq module
module load fastqc/0.12.1-gcc-13.2.0

# output directory for fastqc files
OUT_DIR="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/02_fastqc"

# sample sheet with sample IDs corresponding to each sample and cell type
SAMPLE_SHEET="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/resources/metadata/sample_sheet.txt"

#----------------------------------------#
# get list of files for fastqc and array #
#----------------------------------------#

# put r1 and r2 files into array all_files
mapfile -t ALL_FILES < <(awk -F'\t' 'NF && $1!~/^#/{ if($4!="")print $4; if($5!="")print $5 }' "$SAMPLE_SHEET")

# specify sample based on array number 
SAMPLE=${ALL_FILES[$SLURM_ARRAY_TASK_ID]}

#---------------------------#
# fastqc on raw fastq files #
#---------------------------#

# (Optional) force headless Java; harmless if inputs are valid
export JAVA_TOOL_OPTIONS="-Djava.awt.headless=true"

fastqc \
    -o "${OUT_DIR}" \
    "${SAMPLE}"
