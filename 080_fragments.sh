#!/bin/bash
#SBATCH --time=48:00:00 # cpu partition max time
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --hint=nomultithread 
#SBATCH --job-name=fragments
#SBATCH --output=/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/08_fragments/logs/fragments_%A_%a.out
#SBATCH --error=/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/08_fragments/logs/fragments_%A_%a.err
#SBATCH --array=0-29 # !!! change this as needed

#---------# 
# purpose #
#---------#

# this script converts .sam files to .bam files
# only mapped reads are retained in the output .bam files

#--------------------# 
# set up environment #
#--------------------#

# change dir to where conda envs are 
cd /users/k2587336

# load shell
source ~/.bashrc

# load samtools module (v1.17-gcc-13.2.0-python-3.11.6_ 
# libncurse 5.0 req by conda samtools by 6.5 installed
module load samtools 

source activate bedtools 

# directory of input files
IN_DIR="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/07_filtered_reads"

#out dir
OUT_DIR="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/08_fragments"

#----------------------------#
# array list and file naming #
#----------------------------#

# get list of samples without duplicates
SAMPLES_IN=($(find "${IN_DIR}" -maxdepth 1 -type f -name "*_filtered_noMT_namesorted.bam" | sort))

# get sample from samples list based on array index
SAMPLE="${SAMPLES_IN[$SLURM_ARRAY_TASK_ID]}"

# get sample name by stripping basename and file ext
SAMPLE_NAME=$(basename "$SAMPLE" | cut -d'_' -f1)

#------------------#
# convert to BEDPE #
#------------------#

bedtools bamtobed -bedpe \
  -i "${SAMPLE}" \
  > "${OUT_DIR}/${SAMPLE_NAME}_filtered_namesorted.bed"

#-----------------#
# filter BED file #
#-----------------#

# keep only pairs where BOTH mates have MAPQ ≥ 30 and NEITHER mate maps to chrM/MT AND
# keep only pairs on the same chromosome and where fragment length < 1000 bp
awk '($8>=30) && ($1==$4) && ($6-$2<1000) && ($1!="chrM" && $1!="MT")' \
  "${OUT_DIR}/${SAMPLE_NAME}_filtered_namesorted.bed" \
  > "${OUT_DIR}/${SAMPLE_NAME}_filtered_namesorted_clean.bed"

#----------------------------#
# extract and sort fragments #
#----------------------------#

# Extract fragment intervals
cut -f1,2,6 "${OUT_DIR}/${SAMPLE_NAME}_filtered_namesorted_clean.bed" \
  > "${OUT_DIR}/${SAMPLE_NAME}_fragments.bed"

# sort fragment intervals
sort -k1,1 -k2,2n -k3,3n \
  "${OUT_DIR}/${SAMPLE_NAME}_fragments.bed" \
  > "${OUT_DIR}/${SAMPLE_NAME}_fragments_sorted.bed"

