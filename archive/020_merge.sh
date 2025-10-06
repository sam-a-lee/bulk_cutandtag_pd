#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --job-name=merge_tech_dups
#SBATCH --array=0-59
#SBATCH --output=/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/020_merged/logs/020_merge_%A_%a.out
#SBATCH --error=/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/020_merged/logs/020_merge_%A_%a.err

set -euo pipefail

# set default once (20 reads) if not provided via env
: "${DEBUG_PREVIEW_READS:=20}"


#---------#
# purpose #
#---------#
# Merge FASTQ files of the same library sequenced across lanes/runs to increase depth.
# Forward (R1) and reverse (R2) reads are kept separate.
# Merge early (before filtering/dedup) so reads from the same sample are considered together.

#--------------------#
# set up environment #
#--------------------#

OUT_DIR="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/020_merged"
SAMPLE_SHEET="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/resources/metadata/sample_sheet.txt"

echo "[INFO] Host: $(hostname)"
echo "[INFO] JobID: ${SLURM_JOB_ID:-na}  ArrayTaskID: ${SLURM_ARRAY_TASK_ID:-na}"
echo "[INFO] OUT_DIR: $OUT_DIR"
echo "[INFO] SAMPLE_SHEET: $SAMPLE_SHEET"

mkdir -p "$OUT_DIR"

#-----------------------------------------------#
# identify files to combine and assign to array #
#-----------------------------------------------#

# Build PAIRS as "<out_base>\t<full_path>" after trimming CR/LF and whitespace.
mapfile -t PAIRS < <(awk -F'\t' '
  function trim(s){ gsub(/^[ \t\r\n]+|[ \t\r\n]+$/, "", s); return s }
  BEGIN { OFS="\t" }
  NR==1 { next }                  # skip header
  NF {
    for (i=4; i<=5; i++) {        # columns: fastq_1, fastq_2
      f = trim($i)
      if (f != "") {
        base = f
        sub(/^.*\//, "", base)           # basename only
        gsub(/_L00[0-9]+_/, "_", base)   # collapse lane markers L001/L002/...
        base = trim(base)
        if (base != "") print base, f
      }
    }
  }
' "$SAMPLE_SHEET")

# Unique output basenames (expected merged files)
mapfile -t KEYS < <(printf '%s\n' "${PAIRS[@]}" | awk -F'\t' 'NF && $1!=""{print $1}' | sort -u)

N=${#KEYS[@]}
echo "[INFO] Merged outputs to produce: $N"
if (( N == 0 )); then
  echo "[ERROR] No merge targets parsed from sample sheet (check tabs and columns 4–5)."
  exit 1
fi
echo "[INFO] Recommended --array=0-$((N-1))  (current script has --array=0-59)"

# Debug: list keys with indices so mapping to array is obvious
echo "[DEBUG] KEYS (index : key):"
nl -ba < <(printf '%s\n' "${KEYS[@]}") | sed 's/^/  /'

# Guard array bounds
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" || ${SLURM_ARRAY_TASK_ID} -lt 0 || ${SLURM_ARRAY_TASK_ID} -ge ${N} ]]; then
  echo "[ERROR] SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID:-unset} is out of range 0..$((N-1))."
  exit 2
fi

KEY="${KEYS[$SLURM_ARRAY_TASK_ID]}"
if [[ -z "$KEY" ]]; then
  echo "[ERROR] Computed empty KEY for task ${SLURM_ARRAY_TASK_ID}."
  exit 3
fi

OUT="${OUT_DIR}/${KEY}"

# Collect files corresponding to this KEY
mapfile -t FILES < <(printf '%s\n' "${PAIRS[@]}" \
  | awk -F'\t' -v k="$KEY" 'NF && $1==k && $2!=""{print $2}')

echo "[INFO] KEY: $KEY"
printf "[INFO] Inputs (%d):\n" "${#FILES[@]}"
for f in "${FILES[@]}"; do printf "  [%s]\n" "$(printf "%q" "$f")"; done

#----------------------#
# Merge fastq.gz files #
#----------------------#

if (( ${#FILES[@]} == 0 )); then
  echo "[ERROR] No input files found for KEY='$KEY'."
  exit 4
fi

# Safety: forbid OUT to equal any input
for f in "${FILES[@]}"; do
  if [[ "$f" == "$OUT" ]]; then
    echo "[ERROR] Output path equals an input: $f"
    exit 5
  fi
done

# Input sanity: existence, non-zero, gzip integrity
bad=0
for f in "${FILES[@]}"; do
  if [[ ! -e "$f" ]]; then
    echo "[ERROR] Missing input: $(printf "%q" "$f")"; bad=1; continue
  fi
  if [[ ! -s "$f" ]]; then
    echo "[ERROR] ZERO-BYTE: $(printf "%q" "$f")"; bad=1; continue
  fi
  if ! gzip -t "$f" >/dev/null 2>&1; then
    echo "[ERROR] CORRUPT GZIP: $(printf "%q" "$f")"; bad=1; continue
  fi
done
if (( bad )); then
  echo "[ERROR] Aborting due to bad inputs."
  exit 6
fi

tmp="${OUT}.tmp.$$"
echo "[INFO] Merging -> $tmp"
cat "${FILES[@]}" > "$tmp"

# Validate merged gzip
if ! gzip -t "$tmp" >/dev/null 2>&1; then
  echo "[ERROR] Merged file failed gzip test: $tmp"
  rm -f "$tmp"
  exit 7
fi

_preview_reads=$DEBUG_PREVIEW_READS
if (( _preview_reads > 0 )); then
  echo "[DEBUG] Showing first ${_preview_reads} reads (${_preview_reads}*4 lines):"
  echo "-------- BEGIN PREVIEW (${KEY}) --------"
  set +o pipefail
  { gzip -cd "$tmp" | { head -n $((_preview_reads * 4)); cat >/dev/null; }; } \
    | sed -e $'s/\r$//' \
    | nl -ba -w6 -s': ' \
    || true
  set -o pipefail
  echo "--------- END PREVIEW (${KEY}) ---------"
fi



# Optional: quick count of full reads (decompresses whole file)
reads_in_merged=$(gzip -cd "$tmp" | awk 'END{print int(NR/4)}')
echo "[INFO] Merged file contains at least ${reads_in_merged} full FASTQ records."

mv -f "$tmp" "$OUT"
echo "[INFO] Wrote: $OUT"
