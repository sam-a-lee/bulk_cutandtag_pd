#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --hint=nomultithread
#SBATCH --job-name=fragments
#SBATCH --output=/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/08_fragments/logs/fragments_%A_%a.out
#SBATCH --error=/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/08_fragments/logs/fragments_%A_%a.err
#SBATCH --array=0-29   # update or see dynamic sizing note below

set -euo pipefail

#--------------------#
# env setup
#--------------------#
cd /users/k2587336

source ~/.bashrc

module load samtools

source activate bedtools

IN_DIR="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/07_filtered_reads"

OUT_DIR="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/08_fragments"

#-----------------#
# build file list #
#-----------------#

pattern="*_filtered_namesorted.bam"

# Make a robust, sorted array (handles spaces via -print0)
mapfile -d '' -t SAMPLES_IN < <(find "${IN_DIR}" -maxdepth 1 -type f -name "${pattern}" -print0 | sort -z)

SAMPLE=${SAMPLES_IN[$SLURM_ARRAY_TASK_ID]}

# Derive sample name by stripping the exact SUFFIX
BASE="$(basename "${SAMPLE}")"

SUFFIX="_filtered_namesorted.bam"

SAMPLE_NAME="${BASE%${SUFFIX}}"

echo "[INFO] SAMPLE_NAME=${SAMPLE_NAME}"

#------------------#
# convert to BEDPE #
#------------------#

# keep properly paired, both mapped, drop secondary/supp/dup/QC-fail, MAPQ ≥30
samtools view -b -f 0x2 -F 3852 -q 30 -o "${OUT_DIR}/${SAMPLE_NAME}_hq.bam" "${SAMPLE}"
samtools sort -n -o "${OUT_DIR}/${SAMPLE_NAME}_hq_namesorted.bam" "${OUT_DIR}/${SAMPLE_NAME}_hq.bam"

# to BEDPE (name-sorted input)
bedtools bamtobed -bedpe -i "${OUT_DIR}/${SAMPLE_NAME}_hq_namesorted.bam" > "${OUT_DIR}/${SAMPLE_NAME}_hq_namesorted.bedpe"

#-----------------#
# filter BED file #
#-----------------#
awk '($1==$4) && ($1!="chrM" && $1!="MT") {
  start = ($2 < $5 ? $2 : $5);
  end   = ($3 > $6 ? $3 : $6);
  if ((end - start) < 1000) print
}' "${OUT_DIR}/${SAMPLE_NAME}_hq_namesorted.bedpe" > "${OUT_DIR}/${SAMPLE_NAME}_hq_namesorted_clean.bedpe"

#----------------------------#
# extract and sort fragments #
#----------------------------#

cut -f1,2,6 "${OUT_DIR}/${SAMPLE_NAME}_hq_namesorted_clean.bedpe" > "${OUT_DIR}/${SAMPLE_NAME}_fragments.bed"
sort -k1,1 -k2,2n -k3,3n "${OUT_DIR}/${SAMPLE_NAME}_fragments.bed" > "${OUT_DIR}/${SAMPLE_NAME}_fragments_sorted.bed"

#---------------------------------------------#
# extract fragment lengths to sam TLEN format #
#---------------------------------------------#

# get fragment lengths from bedpe
# to sam file TLEN (col 9) format
awk 'BEGIN{OFS="\t"}
     /^#/ {next}
     { s = ($2 < $5 ? $2 : $5);   
       e = ($3 > $6 ? $3 : $6);   
       L = e - s + 1;             
       if (L > 0) print L
     }' "${OUT_DIR}/${SAMPLE_NAME}_hq_namesorted_clean.bedpe" \
  | sort -n \
  | uniq -c \
  | awk -v OFS="\t" '{print $2,$1}' \
  > "${OUT_DIR}/${SAMPLE_NAME}_fragment_lengths.txt"

