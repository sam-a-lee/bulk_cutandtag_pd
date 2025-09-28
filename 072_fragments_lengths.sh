#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --hint=nomultithread
#SBATCH --job-name=fragment_lengths
#SBATCH --output=/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/07_filtered_reads/logs/fragment_lengths_%A_%a.out
#SBATCH --error=/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/07_filtered_reads/logs/fragment_lengths_%A_%a.err
#SBATCH --array=1-30   # update or see dynamic sizing note below

set -euo pipefail

#--------------------#
# env setup
#--------------------#
cd /users/k2587336

source ~/.bashrc

module load samtools

IN_DIR="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/07_filtered_reads"

OUT_DIR="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/07_filtered_reads/fragment_lengths"

#-----------------#
# build file list #
#-----------------#

pattern="*_filtered_namesorted.sam"

# Make a robust, sorted array (handles spaces via -print0)
mapfile -d '' -t SAMPLES_IN < <(find "${IN_DIR}" -maxdepth 1 -type f -name "${pattern}" -print0 | sort -z)

SAMPLE=${SAMPLES_IN[$SLURM_ARRAY_TASK_ID]}

# Derive sample name by stripping the exact SUFFIX
BASE="$(basename "${SAMPLE}")"

SUFFIX="_filtered_namesorted.sam"

SAMPLE_NAME="${BASE%${SUFFIX}}"

echo "[INFO] SAMPLE_NAME=${SAMPLE_NAME}"

#----------------------------------------#
# extract fragment lengths from sam file #
#----------------------------------------#

# get fragment lengths from bedpe
samtools view \
    -f 0x42 \
    "${SAMPLE}" \
  | awk 'BEGIN{OFS="\t"} { L = $9; if (L<0) L = -L; if (L>0) print L }' \
  | sort -n \
  | uniq -c \
  | awk -v OFS="\t" '{print $2,$1}' \
  > "${OUT_DIR}/${SAMPLE_NAME}_fragment_lengths.txt"