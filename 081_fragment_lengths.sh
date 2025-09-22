#!/bin/bash 
#SBATCH --time=48:00:00 # cpu partition max time
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --job-name=get_frag_len
#SBATCH --output=/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/08_fragments/logs/get_frag_len_%A_%a.out
#SBATCH --error=/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/08_fragments/logs/get_frag_len_%A_%a.err
#SBATCH --array=0-29 # !!! ADJUST ME AS NEEDED BASED ON FILE NUMBER !!! 

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

# dir containing input files
IN_DIR="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/08_fragments"

# output dir for fragment length files
OUT_DIR="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/08_fragments/lengths"

#------------------------#
# sample and file naming #
#------------------------#

# get list of bam files 
mapfile -t FILES_IN < <(find "$IN_DIR" -maxdepth 1 -type f -name "*_filtered_noMT_namesorted_hq.bed" | sort)

# get sample based on master sample list and array index
FILE=${FILES_IN[$SLURM_ARRAY_TASK_ID]}

# get sample name by stripping basename and file ext
BN=$(basename -- "$FILE")
SAMPLE_NAME="${BN%_filtered_noMT_namesorted_hq.bed}"

# out file name
OUT_FILE="${OUT_DIR}/${SAMPLE_NAME}_fragment_lengths.txt"

#--------------------------#
# extract fragment lengths #
#--------------------------#

# get fragment lengths from bed (not - bedpe)
# after final filtering step
awk '{L=$3-$2; print L}' "$FILE" \
  | sort -n \
  | uniq -c \
  | awk -v OFS="\t" '{print $2,$1}' \
  > "$OUT_FILE"
