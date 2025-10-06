#!/bin/bash
#SBATCH --time=48:00:00 # cpu partition max time
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 # bedtools not multithreaded
#SBATCH --mem=4G
#SBATCH --job-name=bigwig
#SBATCH --output=/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/060_peaks/logs/big_wig_%A_%a.out
#SBATCH --error=/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/060_peaks/logs/big_wig_%A_%a.err
#SBATCH --array=0-29 # !!! ADJUST ME AS NEEDED BASED ON FILE NUMBER !!!

set -euo pipefail

#---------#
# purpose #
#---------#
# Create bigWig files for acetylation peaks using deepTools bamCoverage.
# Scale by DESeq2 size factors: bamCoverage --scaleFactor = 1 / (DESeq2 size factor).
# Only samples present in the SCALE_FACTORS file are processed.

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

# load personal conda env
conda activate deeptools

#--------------------#
# user configuration #
#--------------------#
# dir containing input filtered bam files
IN_DIR="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/050_filtered"

# output dir for bigWigs
OUT_DIR="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/060_peaks/069_bigwigs"

# blacklist region
BLACKLIST="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/resources/nordin_hg38_problematic_clean.bed"

# scale factors from DESeq2 (TSV: <sample>\t<size_factor>)
SCALE_FACTORS="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_vanderkant/data_out/070_differential_analysis/size_factors.tsv"

# parameters
BINSIZE=25
MIN_MAPQ=30
THREADS="${SLURM_CPUS_PER_TASK:-1}"

#---------------------#
# prep & file listing #
#---------------------#
mkdir -p "${OUT_DIR}"

# get list of bam files of fragments without duplicates
mapfile -t BAM_FILES < <(find "${IN_DIR}" -maxdepth 1 -type f -name "*_filtered_coordsorted.bam" | sort)

if [[ ${#BAM_FILES[@]} -eq 0 ]]; then
  echo "No BAM files found in ${IN_DIR}" >&2
  exit 1
fi

# get sample based on array index
FILE="${BAM_FILES[$SLURM_ARRAY_TASK_ID]}"
SAMPLE_NAME="$(basename "$FILE" | cut -d'_' -f1)"

echo "[$(date)] Processing sample: ${SAMPLE_NAME}"
echo "BAM: ${FILE}"

#------------------------------#
# get inverse DESeq2 scale fac #
#------------------------------#
if [[ ! -s "${SCALE_FACTORS}" ]]; then
  echo "Scale factors file missing or empty: ${SCALE_FACTORS}" >&2
  exit 1
fi

# Compute 1 / size_factor for this sample (supports scientific notation)
set +e
SCALE_FACTOR="$(
  awk -v s="${SAMPLE_NAME}" 'BEGIN{FS="\t"; sf=""}
    $1==s { sf=$2 }
    END {
      if (sf=="" || sf==0) { exit 1 }
      # awk handles floating and scientific notation natively
      printf("%.12f", 1/sf)
    }' "${SCALE_FACTORS}"
)"
rc=$?
set -e

if [[ $rc -ne 0 || -z "${SCALE_FACTOR}" ]]; then
  echo "No valid (non-zero) DESeq2 size factor for '${SAMPLE_NAME}' in ${SCALE_FACTORS}. Skipping."
  exit 0
fi

echo "Using bamCoverage --scaleFactor = ${SCALE_FACTOR}  (i.e., 1 / DESeq2 size factor)"

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
  echo "WARNING: samtools not found; cannot index BAMs. bamCoverage may be slower or fail."
fi

#------------------#
# run bamCoverage  #
#------------------#
if ! command -v bamCoverage >/dev/null 2>&1; then
  echo "bamCoverage not found in PATH." >&2
  exit 1
fi

OUT_BW="${OUT_DIR}/${SAMPLE_NAME}_deeptools_deseq2_sf_norm.bw"

echo "Running bamCoverage → ${OUT_BW}"
bamCoverage \
  --bam "${FILE}" \
  --outFileName "${OUT_BW}" \
  --outFileFormat bigwig \
  --binSize "${BINSIZE}" \
  --scaleFactor "${SCALE_FACTOR}" \
  --normalizeUsing None \
  --numberOfProcessors "${THREADS}" \
  --minMappingQuality "${MIN_MAPQ}" \
  --ignoreDuplicates \
  --blackListFileName "${BLACKLIST}"

echo "[$(date)] Done: ${OUT_BW}"
