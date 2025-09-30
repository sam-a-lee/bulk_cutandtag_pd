#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --job-name=feature_counts
#SBATCH --output=/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/08_peaks/logs/counts_%j.out
#SBATCH --error=/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/08_peaks/logs/counts_%j.err

#-------------#
# Environment #
#-------------#

CONDA_ROOT="/software/spackages_v0_21_prod/apps/linux-ubuntu22.04-zen2/gcc-13.2.0/anaconda3-2022.10-5wy43yh5crcsmws4afls5thwoskzarhe"
if [ -f "${CONDA_ROOT}/etc/profile.d/conda.sh" ]; then
  . "${CONDA_ROOT}/etc/profile.d/conda.sh"
else
  export PATH="${CONDA_ROOT}/bin:$PATH"
  eval "$(${CONDA_ROOT}/bin/conda shell.bash hook)"
fi

# load samtools module
module load samtools

# activate custom bedtools env
conda activate bedtools

# sample sheet
SAMPLE_SHEET="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/resources/metadata/sample_sheet.txt"

# bam files to read in 
BAM_DIR="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/07_filtered_reads"

# dir with cell type specific consensus peaks
CONSENSUS_DIR="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/08_peaks/consensus_peaks"

# out dir
OUT_DIR="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/08_peaks/cell_type_feature_counts"

# tmp dir for processing
TMP_DIR="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/08_peaks/tmp_counts_$$"

mkdir -p "${OUT_DIR}" "${TMP_DIR}"

# num cpus
CPUS="${SLURM_CPUS_PER_TASK:-8}"

ts() { date "+%Y-%m-%d %H:%M:%S"; }

#--------------------------------------------#
# parse sample sheet to get sample/cell type #
#--------------------------------------------#

echo "[INFO $(ts)] Parsing sample sheet: ${SAMPLE_SHEET}"

# Map of celltype -> sampleID (IGF...); parse col1 (PD/PDC..._) and col4 (IGF####### path)
MAP_TSV="${TMP_DIR}/celltype_sample.tsv"
awk -F'\t' '
  NR>1 {
    if ($1 ~ /^PD(C)?[0-9]+_/) {
      split($1,a,/_/); ct=(length(a)>1?a[2]:"");
      if (match($4, /(IGF[0-9]+)/, m)) id=m[1];
      if (ct!="" && id!="") print ct "\t" id;
    }
  }' "${SAMPLE_SHEET}" | sort -u > "${MAP_TSV}"

if [ ! -s "${MAP_TSV}" ]; then
  echo "[ERROR $(ts)] No (celltype, sampleID) pairs parsed. Check sample sheet format." >&2
  exit 1
fi

mapfile -t CELLTYPES < <(cut -f1 "${MAP_TSV}" | sort -u)

echo "[INFO $(ts)] Cell types: ${CELLTYPES[*]}"

#------------------------------#
# For each cell type / sample  #
#------------------------------#

for ct in "${CELLTYPES[@]}"; do
  echo "[INFO $(ts)] === Cell type: ${ct} ==="

  # consensus peaks (assumes MACS3 narrowPeak filename pattern)
  PEAK_BED="$(ls "${CONSENSUS_DIR}/consensus_${ct}_macs2_q_0_00001_peaks_clean.narrowPeak" 2>/dev/null | head -n1 || true)"
  if [ -z "${PEAK_BED}" ]; then
    echo "[WARN $(ts)] No consensus peaks for ${ct}; skipping."
    continue
  fi
  echo "[INFO $(ts)] Using peaks: ${PEAK_BED}"

  # detect peak columns for header (BED6 vs narrowPeak10; fallback generic)
  PEAK_NCOLS=$(awk -F'\t' 'NF>0{print NF; exit}' "${PEAK_BED}")
  case "${PEAK_NCOLS}" in
    6)  PEAK_HEADER="chrom\tstart\tend\tname\tscore\tstrand" ;;
    10) PEAK_HEADER="chrom\tstart\tend\tname\tscore\tstrand\tsignalValue\tpValue\tqValue\tsummit" ;;
    *)
      PEAK_HEADER=$(awk -v n="${PEAK_NCOLS}" 'BEGIN{for(i=1;i<=n;i++){printf (i==1?"peak_col1":"\tpeak_col"i)} }')
      ;;
  esac

  # sample IDs for this cell type
  mapfile -t SIDS < <(awk -v ct="${ct}" '$1==ct{print $2}' "${MAP_TSV}" | sort -u)
  if [ "${#SIDS[@]}" -eq 0 ]; then
    echo "[WARN $(ts)] No samples mapped for ${ct}; skipping."
    continue
  fi

  CT_TMP="${TMP_DIR}/${ct}"
  mkdir -p "${CT_TMP}"

  for sid in "${SIDS[@]}"; do
    inbam="${BAM_DIR}/${sid}_filtered_namesorted.bam"
    if [ ! -s "${inbam}" ]; then
      echo "[WARN $(ts)] Missing BAM for ${sid}: ${inbam} — skipping."
      continue
    fi

    csbam="${CT_TMP}/${sid}.cs.bam"
    if [ ! -s "${csbam}" ]; then
      echo "[INFO $(ts)] Sorting to coordinate order: ${sid}"
      samtools sort -@ "${CPUS}" -o "${csbam}" "${inbam}"
      samtools index "${csbam}" >/dev/null 2>&1 || true
    fi

    OUT_SAMPLE="${OUT_DIR}/${sid}_${ct}_multicov.tsv"
    echo "[INFO $(ts)] Counting → ${OUT_SAMPLE}"
    bedtools multicov -bams "${csbam}" -bed "${PEAK_BED}" > "${OUT_SAMPLE}"

  done

  echo "[INFO $(ts)] Finished cell type: ${ct}"

done

echo "[INFO $(ts)] Done. Per-sample TSVs in: ${OUT_DIR}"

