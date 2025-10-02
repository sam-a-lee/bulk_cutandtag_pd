#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --job-name=consensus_peaks_by_celltype
#SBATCH --output=/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/070_peaks/logs/074_consensus_peaks_%j.out
#SBATCH --error=/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/070_peaks/logs/074_consensus_peaks_%j.err

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

# load macs conda env
conda activate macs3

# load module 
module load samtools/1.17-gcc-13.2.0-python-3.11.6

# BAMs live here; filenames like IGF136815_filtered_namesorted.bam
BAM_DIR="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/060_filtered"

# sample sheet (tab-delimited). Col1: PDCxxx_celltype; Col4: FASTQ path containing IGF######
SAMPLE_SHEET="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/resources/metadata/sample_sheet.txt"

# output directories
PEAK_DIR="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/070_peaks/074_peaks_consensus_celltype"

# tmp dir for organizing bams by cell type
TMP_DIR="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/070_peaks/074_peaks_consensus_celltype/tmp_$$"

mkdir -p "${PEAK_DIR}" "${TMP_DIR}"

# Exclude any sample IDs starting with these prefixes
EXCLUDE_PREFIXES=("IGF136815" "IGF136816" "IGF136817")

is_excluded() {
  local id="$1"
  for p in "${EXCLUDE_PREFIXES[@]}"; do
    [[ "$id" == ${p}* ]] && return 0
  done
  return 1
}

#----------------------------#
# parse sheet -> ct,sampleid #
#----------------------------#

echo "[INFO] Parsing cell types and sample IDs from ${SAMPLE_SHEET}"

MAP_TSV="${TMP_DIR}/celltype_sample.tsv"

# Produces lines: "<celltype>\t<IGF#######>"
awk -F'\t' '
  NR>1 {
    if ($1 ~ /^PD(C)?[0-9]+_/) {
      split($1, a, /_/);
      ct = (length(a)>1 ? a[2] : "");
      id = "";
      if (match($4, /(IGF[0-9]+)/, m)) { id = m[1]; }
      if (ct != "" && id != "") {
        print ct "\t" id;
      }
    }
  }' "${SAMPLE_SHEET}" | sort -u > "${MAP_TSV}"

if [ ! -s "${MAP_TSV}" ]; then
  echo "[ERROR] No (celltype, sampleID) pairs parsed. Check the sample sheet format." >&2
  exit 1
fi

# unique cell types
mapfile -t CELLTYPES < <(cut -f1 "${MAP_TSV}" | sort -u)

echo "[INFO] Found cell types: ${CELLTYPES[*]}"

#----------------------------------#
# per-celltype: gather BAMs, merge #
#----------------------------------#

for ct in "${CELLTYPES[@]}"; do
  echo "[INFO] ==== Cell type: ${ct} ===="

  # sample IDs for this cell type
  mapfile -t SIDS < <(awk -v ct="${ct}" '$1==ct{print $2}' "${MAP_TSV}" | sort -u)

  if [ "${#SIDS[@]}" -eq 0 ]; then
    echo "[WARN] No samples for ${ct}; skipping."
    continue
  fi

  # build list of BAMs, skipping excluded prefixes
  BAM_LIST="${TMP_DIR}/bams_${ct}.list"
  : > "${BAM_LIST}"

  missing=0
  excluded=0
  included=0
  for sid in "${SIDS[@]}"; do
    if is_excluded "$sid"; then
      excluded=$((excluded+1))
      continue
    fi
    bam="${BAM_DIR}/${sid}_filtered_namesorted.bam"
    if [ -s "${bam}" ]; then
      echo "${bam}" >> "${BAM_LIST}"
      included=$((included+1))
    else
      echo "[WARN] Missing BAM for ${sid}: ${bam}"
      missing=$((missing+1))
    fi
  done

  if [ "${included}" -eq 0 ]; then
    echo "[WARN] No BAMs present for ${ct} after exclusions (excluded=${excluded}, missing=${missing}); skipping."
    continue
  fi
  echo "[INFO] Using ${included} BAM(s) for ${ct} (excluded=${excluded}, missing=${missing})."

  # paths for merged outputs
  MERGED_NAME="${PEAK_DIR}/${ct}_merged_namesorted.bam"

  # merge name-sorted BAMs (inputs are *_namesorted.bam) -> name-sorted merged BAM
  echo "[INFO] Merging with samtools merge -n ..."
  samtools merge -n -@ 4 -b "${BAM_LIST}" -o "${MERGED_NAME}"

  #--------------------#
  # macs3 peak calling #
  #--------------------#

  echo "[INFO] Calling MACS3 peaks for ${ct}"

  macs3 callpeak \
    -t "${MERGED_NAME}" \
    -f BAMPE \
    -g "hs" \
    -q "0.00001" \
    --keep-dup all \
    --verbose 3 \
    --nolambda \
    --nomodel \
    --outdir "${PEAK_DIR}" \
    -n "consensus_${ct}_macs3_q_0_00001"

  echo "[INFO] Finished ${ct}"

done

echo "[INFO] All cell types processed."
echo "[INFO] Peak outputs in: ${PEAK_DIR}"
