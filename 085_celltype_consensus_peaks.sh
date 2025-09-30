#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --job-name=consensus_peaks_by_celltype
#SBATCH --output=/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/08_peaks/logs/consensus_peaks_%j.out
#SBATCH --error=/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/08_peaks/logs/consensus_peaks_%j.err

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

# load modules
module load samtools

# BAMs live here; filenames like IGF136815_filtered_namesorted.bam
BAM_DIR="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/07_filtered_reads"

# sample sheet (tab-delimited). Col1: PDCxxx_celltype; Col4: FASTQ path containing IGF######
SAMPLE_SHEET="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/resources/metadata/sample_sheet.txt"

# output directories
PEAK_DIR="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/08_peaks/consensus_peaks"

# tmp dir for organizing bams by cell type
TMP_DIR="$/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/08_peaks/tmp_$$"

mkdir -p "${TMP_DIR}"

#----------------------------#
# parse sheet -> ct,sampleid #
#----------------------------#

echo "[INFO] Parsing cell types and sample IDs from ${SAMPLE_SHEET}"

MAP_TSV="${TMP_DIR}/celltype_sample.tsv"

# Produces lines: "<celltype>\t<IGF#######>"
awk -F'\t' '
  NR>1 {
    # Enforce: first column must start with PDC or PD + digits + underscore
    if ($1 ~ /^PD(C)?[0-9]+_/) {
      # cell type = everything after first underscore
      split($1, a, /_/);
      ct = (length(a)>1 ? a[2] : "");
      # sample id from col4: IGF followed by digits
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

  # build list of BAMs
  BAM_LIST="${TMP_DIR}/bams_${ct}.list"
  : > "${BAM_LIST}"

  missing=0
  for sid in "${SIDS[@]}"; do
    bam="${BAM_DIR}/${sid}_filtered_namesorted.bam"
    if [ -s "${bam}" ]; then
      echo "${bam}" >> "${BAM_LIST}"
    else
      echo "[WARN] Missing BAM for ${sid}: ${bam}"
      missing=$((missing+1))
    fi
  done

  total=$(wc -l < "${BAM_LIST}" || echo 0)
  if [ "${total}" -eq 0 ]; then
    echo "[WARN] No BAMs present for ${ct}; skipping."
    continue
  fi
  echo "[INFO] Using ${total} BAM(s) for ${ct} (${missing} missing)."

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
    --bw "401" \
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
