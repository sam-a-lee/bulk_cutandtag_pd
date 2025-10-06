#!/bin/bash
#SBATCH --time=48:00:00 # cpu partition max time
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --job-name=merge_groups_bw
#SBATCH --output=/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/070_peaks/logs/merge_bw_%A_%a.out
#SBATCH --error=/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/070_peaks/logs/merge_bw_%A_%a.err
# Tip: set the array after we know how many cell types exist, e.g. --array=0-<N-1>

set -euo pipefail

#---------#
# purpose #
#---------#
# For each cell type (from metadata col3), merge BAMs within each group (col4: control/disease),
# then create bigWig tracks from the merged BAMs using deepTools bamCoverage.
# No DESeq2 size factors are used here.

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

conda activate deeptools

#--------------------#
# user configuration #
#--------------------#
# Input filtered BAMs directory (files like <sample>_filtered_coordsorted.bam)
IN_DIR="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/060_filtered"

# Metadata (TSV): col1=sample, col3=cell_type, col4=group(control|disease). Header allowed.
SAMPLE_META="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_vanderkant/data_out/080_differential_analysis/size_factors.tsv"

# Output directory
OUT_DIR="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/070_peaks/079_bigwigs_groups"

# Blacklist BED
BLACKLIST="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/resources/nordin_hg38_problematic_clean.bed"

# bamCoverage parameters
BINSIZE=25
MIN_MAPQ=30
THREADS="${SLURM_CPUS_PER_TASK:-1}"

# Normalization for visualization of merged groups:
# RPGC requires an effective genome size. For hg38:
NORM="RPGC"
EFFECTIVE_GENOME_SIZE=2913022398   # hg38 effective genome size

mkdir -p "${OUT_DIR}"

#------------------------------#
# collect BAMs and sample map  #
#------------------------------#
# Build a map from sample -> BAM path by matching prefix before first underscore.
declare -A SAMPLE2BAM
while IFS= read -r -d '' bam; do
  base="$(basename "$bam")"
  # sample name is the prefix before the first underscore
  sample="${base%%_*}"
  SAMPLE2BAM["$sample"]="$bam"
done < <(find "${IN_DIR}" -maxdepth 1 -type f -name "*_filtered_coordsorted.bam" -print0 | sort -z)

if [[ ${#SAMPLE2BAM[@]} -eq 0 ]]; then
  echo "No BAM files found in ${IN_DIR}" >&2
  exit 1
fi

#-------------------------#
# parse metadata, get CTs #
#-------------------------#
if [[ ! -s "${SAMPLE_META}" ]]; then
  echo "Metadata TSV missing or empty: ${SAMPLE_META}" >&2
  exit 1
fi

# Get unique cell types (col3), ignore header and empty
mapfile -t CELL_TYPES < <(
  awk -F'\t' 'NR==1 && ($3 ~ /cell|Cell|CELL/) {next} {if($3!="") print $3}' "${SAMPLE_META}" | sort -u
)

if [[ ${#CELL_TYPES[@]} -eq 0 ]]; then
  echo "No cell types found in column 3 of ${SAMPLE_META}" >&2
  exit 1
fi

# Select the cell type for this array task
IDX="${SLURM_ARRAY_TASK_ID:-0}"
if (( IDX < 0 || IDX >= ${#CELL_TYPES[@]} )); then
  echo "SLURM_ARRAY_TASK_ID ${IDX} out of range (0..$(( ${#CELL_TYPES[@]} - 1 )))" >&2
  echo "Found ${#CELL_TYPES[@]} cell types. Set --array=0-$(( ${#CELL_TYPES[@]} - 1 ))" >&2
  exit 1
fi

CELL_TYPE="${CELL_TYPES[$IDX]}"
echo "[$(date)] Processing cell type: ${CELL_TYPE}"

#------------------------------#
# gather sample lists per group#
#------------------------------#
# Get lists of samples for this cell type by group (col4)
mapfile -t CONTROL_SAMPLES < <(
  awk -F'\t' -v ct="${CELL_TYPE}" '
    NR==1 && ($3 ~ /cell|Cell|CELL/ || $4 ~ /group|Group|GROUP/) {next}
    tolower($3)==tolower(ct) && tolower($4) ~ /^control$/ {print $1}
  ' "${SAMPLE_META}" | sort -u
)

mapfile -t DISEASE_SAMPLES < <(
  awk -F'\t' -v ct="${CELL_TYPE}" '
    NR==1 && ($3 ~ /cell|Cell|CELL/ || $4 ~ /group|Group|GROUP/) {next}
    tolower($3)==tolower(ct) && tolower($4) ~ /^disease$/ {print $1}
  ' "${SAMPLE_META}" | sort -u
)

echo "Found ${#CONTROL_SAMPLES[@]} control and ${#DISEASE_SAMPLES[@]} disease samples for ${CELL_TYPE}"

# Turn sample lists into BAM lists (skip missing)
gather_bams() {
  local -n samples_ref=$1
  local out=()
  for s in "${samples_ref[@]}"; do
    if [[ -n "${SAMPLE2BAM[$s]:-}" ]]; then
      out+=("${SAMPLE2BAM[$s]}")
    else
      echo "WARNING: No BAM found for sample '${s}' (expected prefix match in ${IN_DIR})"
    fi
  done
  printf '%s\0' "${out[@]}"  # NUL-separated for safety
}

CONTROL_BAMS=()
DISEASE_BAMS=()
while IFS= read -r -d '' p; do CONTROL_BAMS+=("$p"); done < <(gather_bams CONTROL_SAMPLES)
while IFS= read -r -d '' p; do DISEASE_BAMS+=("$p"); done < <(gather_bams DISEASE_SAMPLES)

# Output subdir per cell type
CT_OUT="${OUT_DIR}/${CELL_TYPE}"
mkdir -p "${CT_OUT}"

#------------------------------#
# function: merge and bigwig   #
#------------------------------#
merge_and_make_bw() {
  local label="$1"         # control or disease
  shift
  local -a bams=( "$@" )

  if [[ ${#bams[@]} -eq 0 ]]; then
    echo "No ${label} BAMs for ${CELL_TYPE}; skipping ${label}."
    return 0
  fi

  local merged_bam="${CT_OUT}/${CELL_TYPE}_${label}_merged.bam"
  local merged_bw="${CT_OUT}/${CELL_TYPE}_${label}_merged_${NORM,,}.bw"

  echo "Merging ${#bams[@]} ${label} BAM(s) → ${merged_bam}"
  if [[ -f "${merged_bam}" ]]; then
    echo "Merged BAM exists, will overwrite."
    rm -f "${merged_bam}" "${merged_bam}.bai"
  fi

  # samtools merge
  samtools merge -@ "${THREADS}" -f "${merged_bam}" "${bams[@]}"
  echo "Indexing ${merged_bam}"
  samtools index -@ "${THREADS}" "${merged_bam}"

  # bamCoverage
  echo "Creating bigWig → ${merged_bw}"
  bamCoverage \
    --bam "${merged_bam}" \
    --outFileName "${merged_bw}" \
    --outFileFormat bigwig \
    --binSize "${BINSIZE}" \
    --normalizeUsing "${NORM}" \
    --effectiveGenomeSize "${EFFECTIVE_GENOME_SIZE}" \
    --numberOfProcessors "${THREADS}" \
    --minMappingQuality "${MIN_MAPQ}" \
    --ignoreDuplicates \
    --blackListFileName "${BLACKLIST}"

  echo "Done: ${merged_bw}"
}

#------------------------------#
# run per group                #
#------------------------------#
merge_and_make_bw "control" "${CONTROL_BAMS[@]}"
merge_and_make_bw "disease" "${DISEASE_BAMS[@]}"

echo "[$(date)] Completed cell type: ${CELL_TYPE}"
