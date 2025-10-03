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

# in directory for bam files
IN_DIR="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/060_filtered"

# out directory for frip scores
FRIP_DIR="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/070_peaks/072_frip_scores"

mkdir -p ${FRIP_DIR}

# peak dir (after removing blacklisted)
PEAK_DIR="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/070_peaks/071_peaks_blacklist_removed"

# assaign to array
mapfile -t BAM_FILES < <(find "${IN_DIR}" -maxdepth 1 -type f -name "*_filtered_namesorted.bam" | sort)

# get file path
FILE="${BAM_FILES[$SLURM_ARRAY_TASK_ID]}"

# extract sample name
SAMPLE_NAME="$(basename "$FILE" | cut -d'_' -f1)"

echo "Processing sample: ${SAMPLE_NAME}"

#------------------------------------#
# coordinate sort and index bam file #
#------------------------------------#

# load samtools module
 
# coordinate sor
samtools sort -@ 4 -o "${IN_DIR}/${SAMPLE_NAME}_filtered_coordsorted.bam" "${FILE}"

# index
samtools index -@ 4 "${IN_DIR}/${SAMPLE_NAME}_filtered_coordsorted.bam"

#-----------------#
# get frip scores #
#-----------------#

# activate deeptools
conda activate deeptools    

plotEnrichment \
  -b "${IN_DIR}/${SAMPLE_NAME}_filtered_coordsorted.bam" \
  --BED "${PEAK_DIR}/${SAMPLE_NAME}_macs3_q_0_00001_h3k27ac_peaks_clean.narrowPeak" \
  -o "${FRIP_DIR}/${SAMPLE_NAME}_sample_enrichment.pdf" \
  --labels "${SAMPLE_NAME}" \
  --blackListFileName "/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/resources/nordin_hg38_problematic_clean.bed" \
  --outRawCounts "${FRIP_DIR}/${SAMPLE_NAME}_frip_score.txt"
