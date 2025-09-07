#!/bin/bash 
#SBATCH --time=48:00:00 # cpu partition max time
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --job-name=get_frag_len
#SBATCH --output=/scratch/users/k2587336/test/data_out/06_fragment_lengths/get_frag_len_%A_%a.out
#SBATCH --error=/scratch/users/k2587336/test/data_out/06_fragment_lengths/get_frag_len_%A_%a.err
#SBATCH --array=0-0 # !!! ADJUST ME AS NEEDED BASED ON FILE NUMBER !!! 

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
ROOT_DIR="/scratch/users/k2587336/test/data_out"


#------------------------#
# sample and file naming #
#------------------------#

# get list of files without duplicates
SAM_FILES_IN=($(find "$ROOT_DIR/03_duplicates_removed" -maxdepth 1 -type f -name "*_picard_dup_removed.sam" | sort))

# get sample based on master sample list and array index
FILE="${SAM_FILES_IN[$SLURM_ARRAY_TASK_ID]}"

# get sample name by stripping basename and file ext
SAMPLE_NAME=$(basename "$FILE" | cut -d'_' -f1)

# out file name
OUT_FILE="$ROOT_DIR/06_fragment_lengths/${SAMPLE_NAME}_samtools_frag_len.txt"


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