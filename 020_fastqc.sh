#!/bin/bash
#SBATCH --time=48:00:00 # cpu partition max time per kcl policy
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 # could specify more here if wanted 
#SBATCH --mem=1G # shouldnt require more than a gb per file
#SBATCH --hint=nomultithread # prefer physical cores for better throughput
#SBATCH --job-name=fastqc
#SBATCH --array=0-111 # modify based on number of files
#SBATCH --output=/scratch/prj/bcn_pd_pesticides/Files-From-Imperial/analysis_cutandtag_pd_bulk/data_out/02_fastqc/logs/fastqc_%A_%a.out
#SBATCH --error=/scratch/prj/bcn_pd_pesticides/Files-From-Imperial/analysis_cutandtag_pd_bulk/data_out/02_fastqc/logs/fastqc_%A_%a.err

#---------#
# purpose #
#---------#

# this files examines quality of raw fastq files prior to merging of technical duplicates

#--------------------#
# set up environment #
#--------------------#

# change dir to where conda envs are 
cd /users/k2587336

# load shell
source ~/.bashrc

# activate conda env
module load fastqc/0.12.1-gcc-13.2.0

# input directory of fastq files
IN_DIR="/scratch/prj/bcn_pd_pesticides/Files-From-Imperial/bulkCutandTag_PD/AACNGMFHV"

# output directory for fastqc files
OUT_DIR="/scratch/prj/bcn_pd_pesticides/Files-From-Imperial/analysis_cutandtag_pd_bulk/data_out/02_fastqc"

# sample sheet with sample IDs corresponding to each sample and cell type
SAMPLE_SHEET="/scratch/prj/bcn_pd_pesticides/Files-From-Imperial/analysis_cutandtag_pd_bulk/resources/metadata/sample_sheet.txt"

#----------------------------------------#
# get list of files for fastqc and array #
#----------------------------------------#

# read first column (tab-delimited) as SAMPLE_IDs 
SAMPLE_ID=($(awk -F'\t' 'NF && $1!~/^#/{print $1}' "$SAMPLE_SHEET"))

# build -name patterns safely (no literal quotes), then find & sort
args=()
for id in "${SAMPLE_ID[@]}"; do
  args+=(-o -name "*${id}_S[0-9]*_L00[12]_R[12]_001.fastq.gz")
done

# bollect matches (null-delimited for safety)
readarray -d '' ALL_SAMPLES < <(find "${IN_DIR}" -type f \( -false "${args[@]}" \) -print0 | sort -z)

# guard against empty list or bad array index (prevents GUI fallback)
if (( ${#ALL_SAMPLES[@]} == 0 )); then
  echo "ERROR: No FASTQ files matched under ${IN_DIR} using IDs from ${SAMPLE_SHEET}" >&2
  exit 2
fi
if (( SLURM_ARRAY_TASK_ID < 0 || SLURM_ARRAY_TASK_ID >= ${#ALL_SAMPLES[@]} )); then
  echo "ERROR: SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} out of range (0..$(( ${#ALL_SAMPLES[@]} - 1 )))" >&2
  exit 3
fi

# specify sample based on array number 
SAMPLE=${ALL_SAMPLES[$SLURM_ARRAY_TASK_ID]}

#---------------------------#
# fastqc on raw fastq files #
#---------------------------#

# (Optional) force headless Java; harmless if inputs are valid
export JAVA_TOOL_OPTIONS="-Djava.awt.headless=true"

fastqc \
    -o "${OUT_DIR}" \
    "${SAMPLE}"
