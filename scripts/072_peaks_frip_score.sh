#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --job-name=frip_score
#SBATCH --output=/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/070_peaks/logs/072_frip_score_%A_%a.out
#SBATCH --error=/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/070_peaks/logs/072_frip_score_%A_%a.err
#SBATCH --array=0-29

#--------------------#
# set up environment #
#--------------------#

# initiate conda
CONDA_ROOT="/software/spackages_v0_21_prod/apps/linux-ubuntu22.04-zen2/gcc-13.2.0/anaconda3-2022.10-5wy43yh5crcsmws4afls5thwoskzarhe"
if [ -f "${CONDA_ROOT}/etc/profile.d/conda.sh" ]; then
  . "${CONDA_ROOT}/etc/profile.d/conda.sh"
else
  export PATH="${CONDA_ROOT}/bin:$PATH"
  eval "$(${CONDA_ROOT}/bin/conda shell.bash hook)"
fi

# activate deeptools
conda activate deeptools

# load samtools module
module load samtools/1.17-gcc-13.2.0-python-3.11.6

# in directory for bam files
IN_DIR="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/060_filtered"

# out directory for frip scores
FRIP_DIR="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/070_peaks/072_frip_scores"
mkdir -p "${FRIP_DIR}"

# peak dir (after removing blacklisted)
PEAK_DIR="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/070_peaks/071_peaks_blacklist_removed"

# assign to array
mapfile -t BAM_FILES < <(find "${IN_DIR}" -maxdepth 1 -type f -name "*_filtered_namesorted.bam" | sort)

# get file path
FILE="${BAM_FILES[$SLURM_ARRAY_TASK_ID]}"

# extract sample name
SAMPLE_NAME="$(basename "$FILE" | cut -d'_' -f1)"
echo "Processing sample: ${SAMPLE_NAME}"

#------------------------------------#
# coordinate sort and index bam file #
#------------------------------------#

# paths for coord-sorted BAM and possible index filenames
BAM_COORD="${IN_DIR}/${SAMPLE_NAME}_filtered_coordsorted.bam"
BAM_INDEX1="${BAM_COORD}.bai"             # e.g., file.bam.bai
BAM_INDEX2="${BAM_COORD%.bam}.bai"        # e.g., file.bai

# coordinate sort
samtools sort -@ 4 -o "${BAM_COORD}" "${FILE}"

# index only if not already indexed
if [[ -f "${BAM_INDEX1}" || -f "${BAM_INDEX2}" ]]; then
  echo "Index for ${BAM_COORD} already exists. Skipping indexing."
else
  echo "Index not found for ${BAM_COORD}. Creating index..."
  samtools index -@ 4 "${BAM_COORD}"
fi

#-----------------#
# get frip scores #
#-----------------#

plotEnrichment \
  -b "${BAM_COORD}" \
  --BED "${PEAK_DIR}/${SAMPLE_NAME}_macs3_q_0_00001_h3k27ac_peaks_clean.narrowPeak" \
  -o "${FRIP_DIR}/${SAMPLE_NAME}_sample_enrichment.pdf" \
  --labels "${SAMPLE_NAME}" \
  --blackListFileName "/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/resources/nordin_hg38_problematic_clean.bed" \
  --outRawCounts "${FRIP_DIR}/${SAMPLE_NAME}_frip_score.txt"
