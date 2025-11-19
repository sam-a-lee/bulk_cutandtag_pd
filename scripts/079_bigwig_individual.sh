#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1                 
#SBATCH --mem=8G
#SBATCH --job-name=bigwig
#SBATCH --output=/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/070_peaks/logs/big_wig_%A_%a.out
#SBATCH --error=/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/070_peaks/logs/big_wig_%A_%a.err
#SBATCH --array=0-29                     

#---------#
# purpose #
#---------#

# Create bigWig files for CUT&Tag data using deepTools bamCoverage with literature-consistent defaults:
#   - CPM normalization
#   - MAPQ ≥ 30
#   - ignore duplicates
#   - blacklist masking
#   - fragment-centric coverage: --extendReads (typical for PE CUT&Tag) + --centerReads
#   - fine binning + smoothing (defaults below for histone marks; see MARK_CLASS)

#--------------------#
# set up environment #
#--------------------#

CONDA_ROOT="/software/spackages_v0_21_prod/apps/linux-ubuntu22.04-zen2/gcc-13.2.0/anaconda3-2022.10-5wy43yh5crcsmws4afls5thwoskzarhe"

if [ -f "${CONDA_ROOT}/etc/profile.d/conda.sh" ]; then
  . "${CONDA_ROOT}/etc/profile.d/conda.sh"
else
  export PATH="${CONDA_ROOT}/bin:$PATH"
  eval "$(${CONDA_ROOT}/bin/conda shell.bash hook)"
fi

# load custom deeptools
conda activate deeptools

# load samtools module
module load samtools/1.17-gcc-13.2.0-python-3.11.6

#--------------------#
# user configuration #
#--------------------#

# Input directory containing aligned, filtered BAMs
IN_DIR="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/060_filtered"

# Output directory for bigWigs
OUT_DIR="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/070_peaks/079_bigwigs"

mkdir -p "${OUT_DIR}"

# ENCODE-/project-specific blacklist (leave empty to disable)
BLACKLIST="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/resources/nordin_hg38_problematic_clean.bed"

#---------------------#
# prep & file listing #
#---------------------#

# Find BAMs
mapfile -t BAM_FILES < <(find "${IN_DIR}" -maxdepth 1 -type f -name "*_filtered_coordsorted.bam" | sort)
if [[ ${#BAM_FILES[@]} -eq 0 ]]; then
  echo "No BAM files found in ${IN_DIR}" >&2
  exit 1
fi

# Get BAM for this array index
FILE="${BAM_FILES[$SLURM_ARRAY_TASK_ID]}"
BASENAME="$(basename "$FILE")"

# Derive sample name (adjust this parsing if your naming differs)
SAMPLE_NAME="$(basename "$FILE" | cut -d'_' -f1)"

echo "[$(date)] Processing sample: ${SAMPLE_NAME}"
echo "BAM: ${FILE}"

#------------------#
# ensure bam index #
#------------------#

if command -v samtools >/dev/null 2>&1; then
  if [[ ! -f "${FILE}.bai" && ! -f "${FILE%.*}.bai" ]]; then
    echo "Index not found for ${FILE}. Creating index..."
    samtools index -@ "${THREADS}" "${FILE}"
  else
    echo "Index present for ${FILE}."
  fi
else
  echo "WARNING: samtools not found; cannot index BAMs."
fi

#------------------#
# run bamCoverage  #
#------------------#

if ! command -v bamCoverage >/dev/null 2>&1; then
  echo "bamCoverage not found in PATH." >&2
  exit 1
fi

# output sample name 
OUT_BW="${OUT_DIR}/${SAMPLE_NAME}_CPM_bs25.bw"

echo "Running bamCoverage → ${OUT_BW}"

# Bin/smoothing by class
# H3K4me3 (narrow):  BINSIZE=10 SMOOTH=30 
# H3K27ac/H3K4me1 (enhancer): BINSIZE=25 SMOOTH=75 
# H3K27me3/H3K9me3 (broad): BINSIZE=50 SMOOTH=150 

bamCoverage \
  --bam "${FILE}" \
  --outFileName "${OUT_BW}" \
  --outFileFormat bigwig \
  --binSize 25 \
  --smoothLength 75 \
  --minMappingQuality 30 \
  --ignoreDuplicates \
  --centerReads \
  --extendReads \
  --skipNonCoveredRegions \
  --normalizeUsing CPM \
  --blackListFileName "${BLACKLIST}"

echo "[$(date)] Done: ${OUT_BW}"
