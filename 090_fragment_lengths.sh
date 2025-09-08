#!/bin/bash 
#SBATCH --time=48:00:00 # cpu partition max time
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --job-name=get_frag_len
#SBATCH --output=/scratch/prj/bcn_pd_pesticides/Files-From-Imperial/analysis_cutandtag_pd_bulk/data_out/09_fragment_lengths/logs/get_frag_len_%A_%a.out
#SBATCH --error=/scratch/prj/bcn_pd_pesticides/Files-From-Imperial/analysis_cutandtag_pd_bulk/data_out/09_fragment_lengths/logs/get_frag_len_%A_%a.err
#SBATCH --array=0-27 # !!! ADJUST ME AS NEEDED BASED ON FILE NUMBER !!! 

#--------------------#
# environment set up #
#--------------------#

# change dir to where conda envs are 
cd /users/k2587336

# load shell
source ~/.bashrc

# load samtools module (v1.17-gcc-13.2.0-python-3.11.6_ 
# libncurse 5.0 req by conda samtools by 6.5 installed
module load samtools 

# set root dir for working  
ROOT_DIR="/scratch/prj/bcn_pd_pesticides/Files-From-Imperial/analysis_cutandtag_pd_bulk/data_out"

# dir containing input files
IN_DIR="${ROOT_DIR}/06_duplicates_removed"

# output dir for fragment length files
OUT_DIR="${ROOT_DIR}/09_fragment_lengths"

#------------------------#
# sample and file naming #
#------------------------#

# get list of bam files 
mapfile -t SAM_FILES_IN < <(find "$IN_DIR" -maxdepth 1 -type f -name "*_picard_dup_removed.sam" | sort)

# sanity
# get num of files found and echo
N=${#SAM_FILES_IN[@]}
echo "[INFO] Found $N BAMs in $BAM_DIR"

# fail gracefully if no files found or more files than specified in array
if (( N == 0 )); then echo "[WARN] No inputs; exiting."; exit 0; fi
if (( SLURM_ARRAY_TASK_ID >= N )); then echo "[INFO] No work for task $SLURM_ARRAY_TASK_ID (N=$N)"; exit 0; fi

# get sample based on master sample list and array index
FILE="${SAM_FILES_IN[$SLURM_ARRAY_TASK_ID]}"

# get sample name by stripping basename and file ext
BN=$(basename -- "$FILE")
SAMPLE_NAME="${BN%_picard_dup_removed.sam}"

# out file name
OUT_FILE="${OUT_DIR}/${SAMPLE_NAME}_samtools_frag_len.txt"

#--------------------------#
# extract fragment lengths #
#--------------------------#

# 9th column in sam file contains fragment length - extract this! 
samtools view -@ 4 -F 0x04 ${FILE} | \
    awk -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($9)}' | \
    sort | \
    uniq -c | \
    awk -v OFS="\t" '{print $2, $1/2}' > \
    ${OUT_FILE}