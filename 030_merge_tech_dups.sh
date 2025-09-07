#!/bin/bash
#SBATCH --time=48:00:00 
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 
#SBATCH --mem=2G
#SBATCH --job-name=merge_tech_dups
#SBATCH --array=0-91 # could mcodify based on number of files - %20 throttles to 20 at a time to be polite
#SBATCH --output=/scratch/prj/bcn_pd_pesticides/Files-From-Imperial/analysis_cutandtag_pd_bulk/data_out/03_merged_fastq/logs/merge_%A_%a.out
#SBATCH --error=/scratch/prj/bcn_pd_pesticides/Files-From-Imperial/analysis_cutandtag_pd_bulk/data_out/03_merged_fastq/logs/merge_%A_%a.err

#---------#
# purpose #
#---------#

# this files merges library technical duplicates into single files of forward and reverse reads for a given sample/cell type

#--------------------#
# set up environment #
#--------------------#

# change dir to where conda envs are 
cd /users/k2587336

# load shell 
source ~/.bashrc

# input directory of fastq files
IN_DIR="/scratch/prj/bcn_pd_pesticides/Files-From-Imperial/bulkCutandTag_PD/AACNGMFHV"

# output directory for merged files
OUT_DIR="/scratch/prj/bcn_pd_pesticides/Files-From-Imperial/analysis_cutandtag_pd_bulk/data_out/03_merged_fastq"

# sample sheet with sample IDs corresponding to each sample and cell type
SAMPLE_SHEET="/scratch/prj/bcn_pd_pesticides/Files-From-Imperial/analysis_cutandtag_pd_bulk/resources/metadata/sample_sheet.txt"

echo "[INFO] Host: $(hostname)"
echo "[INFO] JobID: ${SLURM_JOB_ID:-na}  ArrayTaskID: ${SLURM_ARRAY_TASK_ID:-na}"
echo "[INFO] IN_DIR:  $IN_DIR"
echo "[INFO] OUT_DIR: $OUT_DIR"
echo "[INFO] SAMPLE_SHEET: $SAMPLE_SHEET"

#-----------------------------------#
# get list of libraries in lane one #
#-----------------------------------#

# read first column (tab-delimited) as SAMPLE_IDs 
mapfile -t SAMPLE_ID < <(awk -F'\t' 'NF && $1!~/^#/{gsub(/\r/,"",$1); print $1}' "$SAMPLE_SHEET")
echo "[INFO] Loaded ${#SAMPLE_ID[@]} sample IDs from sheet."

# build robust -name args for find
name_args=()
for id in "${SAMPLE_ID[@]}"; do
  name_args+=(-o -name "*${id}_S[0-9]*_L001_R[12]_001.fastq.gz")
done

mapfile -t ALL_L001_SAMPLES < <(find "$IN_DIR" -type f \( -false "${name_args[@]}" \) -print | sort)

N=${#ALL_L001_SAMPLES[@]}
echo "[INFO] Found ${N} L001 files matching patterns."
if (( N == 0 )); then
  echo "[WARN] No L001 matches. Check that filenames look like *_L001_R1_001.fastq.gz and *_L001_R2_001.fastq.gz."
  echo "[WARN] Also check for Windows CR in the sample sheet; we strip it, but verify IDs."
  exit 0
fi

if (( SLURM_ARRAY_TASK_ID >= N )); then
  echo "[INFO] No work for task ${SLURM_ARRAY_TASK_ID} (N=${N}); exiting."
  exit 0
fi

L001_SAMPLE=${ALL_L001_SAMPLES[$SLURM_ARRAY_TASK_ID]}

L002_SAMPLE="${L001_SAMPLE/\/1\/16_NA\//\/2\/16_NA\/}"
L002_SAMPLE="${L002_SAMPLE/_L001_/_L002_}"

base=$(basename -- "$L001_SAMPLE")
NO_LANES_ID=$(sed -E 's/_L00[12]_/_/' <<<"$base")

echo "[INFO] L001: $L001_SAMPLE"
echo "[INFO] L002: $L002_SAMPLE"
echo "[INFO] Out:  $OUT_DIR/$NO_LANES_ID"

# sanity check
[[ -n "$NO_LANES_ID" && "$NO_LANES_ID" != "." ]]
[[ -f "$L001_SAMPLE" ]] || { echo "[ERROR] Missing L001: $L001_SAMPLE"; exit 1; }
[[ -f "$L002_SAMPLE" ]] || { echo "[ERROR] Missing L002: $L002_SAMPLE"; exit 1; }


#-----------------------------#
# merge lane one and lane two #
#-----------------------------#

cat ${L001_SAMPLE} ${L002_SAMPLE} > "${OUT_DIR}/${NO_LANES_ID}"
