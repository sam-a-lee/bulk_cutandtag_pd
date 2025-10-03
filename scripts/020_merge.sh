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

#---------#
# purpose #
#---------#

# this script merges fastq files of the same library that was sequences across lanes/runs to increase depth
# forward (R1) and reverse (R2) reads are kept separate 
# nb: merging should be done as early as possible and before filtering and deduplication as reads 
# come from the same sample and should be considered together

#--------------------#
# set up environment #
#--------------------#

# Output directory for merged files (fixed missing '/')
OUT_DIR="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/020_merged"

# Sample sheet (tab-delimited: sample, group, replicate, fastq_1, fastq_2, control)
SAMPLE_SHEET="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/resources/metadata/sample_sheet.txt"

echo "[INFO] Host: $(hostname)"
echo "[INFO] JobID: ${SLURM_JOB_ID:-na}  ArrayTaskID: ${SLURM_ARRAY_TASK_ID:-na}"
echo "[INFO] OUT_DIR: $OUT_DIR"
echo "[INFO] SAMPLE_SHEET: $SAMPLE_SHEET"


#-----------------------------------------------#
# identify files to combine and assign to array #
#-----------------------------------------------#

# key = output basename with lane removed (L001/L002 collapsed)
# note: array size should be equal to the number of expected files after merging

# Pairs look like: "<out_base>\t<full_path>"
mapfile -t PAIRS < <(awk -F'\t' '
  NR==1 { next }      # skip header
  NF {
    for (i=4;i<=5;i++) if ($i!="") {
      f=$i
      gsub(/\r/,"",f)
      base=f
      sub(/^.*\//,"",base)          # basename only
      gsub(/_L00[12]_/, "_", base)  # collapse lane
      print base "\t" f
    }
  }' "$SAMPLE_SHEET")

# Unique output basenames (each corresponds to one merged file)
mapfile -t KEYS < <(printf '%s\n' "${PAIRS[@]}" | cut -f1 | sort -u)

N=${#KEYS[@]}
echo "[INFO] Merged outputs to produce: $N"
if (( N == 0 )); then
  echo "[WARN] Nothing to do. Check the sample sheet."
  exit 0
fi

KEY=${KEYS[$SLURM_ARRAY_TASK_ID]}
OUT="${OUT_DIR}/${KEY}"

# Collect all files for this KEY (one or two lanes; occasionally more)
mapfile -t FILES < <(printf '%s\n' "${PAIRS[@]}" | awk -v k="$KEY" -F'\t' '$1==k{print $2}')

echo "[INFO] KEY: $KEY"
printf "[INFO] Inputs (%d):\n" "${#FILES[@]}"
printf "  %s\n" "${FILES[@]}"

# Sanity checks: files exist
for f in "${FILES[@]}"; do
  [[ -f "$f" ]] || { echo "[ERROR] Missing input: $f"; exit 1; }
done

#----------------------#
# Merge fastq.gz files #
#----------------------#

cat "${FILES[@]}" > "$OUT"
echo "[INFO] Wrote: $OUT"
