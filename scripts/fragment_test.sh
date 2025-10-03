#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
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
samtools view -h -f 0x2 -F 3852 "${SAMPLE}" \
| awk 'BEGIN{OFS="\t"}
  /^@/ { print; next }                           # pass headers
  function is_mito(c,   lc){ lc=tolower(c); return (lc=="chrm" || lc=="mt" || lc=="m") }
  {
    # buffer pairs by qname (name-sorted, primary only → 2 lines per qname)
    if ($1 != qn) {
      qn = $1; l1 = $0; split($0,f,"\t")
      r1chr=f[3]; r1mapq=f[5]+0; r1tlen=f[9]+0
      have1=1; next
    } else {
      l2 = $0; split($0,g,"\t")
      r2chr=g[3]; r2mapq=g[5]+0

      abs_tlen = (r1tlen<0)? -r1tlen : r1tlen
      pass = (r1mapq>=30 && r2mapq>=30) \
             && (!is_mito(r1chr) && !is_mito(r2chr)) \
             && (r1chr == r2chr) \
             && (abs_tlen <= 1000)

      if (pass) { print l1; print l2 }
      have1=0
    }
  }' \
| samtools view -b -o "${OUT_DIR}/${SAMPLE_NAME}_hq.bam" -

# re sort by name before bedpe converstion
samtools sort -n -o "${OUT_DIR}/${SAMPLE_NAME}_hq_namesorted.bam" "${OUT_DIR}/${SAMPLE_NAME}_hq.bam"

# to BEDPE (name-sorted input)
bedtools bamtobed -bedpe -i "${OUT_DIR}/${SAMPLE_NAME}_hq_namesorted.bam" > "${OUT_DIR}/${SAMPLE_NAME}_hq_namesorted.bed"

#----------------------------#
# extract and sort fragments #
#----------------------------#

cut -f1,2,6 "${OUT_DIR}/${SAMPLE_NAME}_hq_namesorted.bed" > "${OUT_DIR}/${SAMPLE_NAME}_fragments.bed"
sort -k1,1 -k2,2n -k3,3n "${OUT_DIR}/${SAMPLE_NAME}_fragments.bed" > "${OUT_DIR}/${SAMPLE_NAME}_fragments_sorted.bed"

#----------------------------------------#
# extract fragment lengths from bam file #
#----------------------------------------#

# get fragment lengths from bedpe
samtools view \ 
    -f 0x42 \ # require proper pair AND read1 to about double counting
    "${OUT_DIR}/${SAMPLE_NAME}_hq_namesorted.bam" \ 
    | awk 'BEGIN{OFS="\t"} { L = $9; if (L<0) L = -L; if (L>0) print L }' \ 
    | sort -n \ 
    | uniq -c \ 
    | awk -v OFS="\t" '{print $2,$1}' \ 
    > "${OUT_DIR}/${SAMPLE_NAME}_fragment_lengths.txt"