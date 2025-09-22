#!/bin/bash
#SBATCH --time=48:00:00 # cpu partition max time
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --hint=nomultithread 
#SBATCH --job-name=filter_reads
#SBATCH --output=/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/07_filtered_reads/logs/filter_reads_%A_%a.out
#SBATCH --error=/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/07_filtered_reads/logs/filter_reads_%A_%a.err
#SBATCH --array=0-30 # !!! change this as needed

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

# root working directory
ROOT_DIR="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out"

# directory of input files
IN_DIR="${ROOT_DIR}/06_duplicates_removed"

#out dir
OUT_DIR="${ROOT_DIR}/07_filtered_reads"

#----------------------------#
# array list and file naming #
#----------------------------#

# get list of samples without duplicates
SAMPLES_IN=($(find "${IN_DIR}" -maxdepth 1 -type f -name "*_picard_dup_rm.sam" | sort))

# get sample from samples list based on array index
SAMPLE="${SAMPLES_IN[$SLURM_ARRAY_TASK_ID]}"

# get sample name by stripping basename and file ext
SAMPLE_NAME=$(basename "$SAMPLE" | cut -d'_' -f1)

#---------------------------#
# filter and convert to bam #
#---------------------------#

#  name sort reads
samtools sort -n -@ 4 -O SAM "${SAMPLE}" -o "${OUT_DIR}/${SAMPLE_NAME}_namesorted.sam"

# filter sorted
# -f 2 require properly paired reads (both mates mapped in a proper pair) - includes dovetail for bwamem2
# -F 256 drop secondary alignments
# -F 2048 drop supplementary alignments
# -q 30 keep only reads with mapq > 30 (done on bedpe)
# -F 12 drop unmapped read (0x4) and unmapped mate (0x8) - both mates mapped
samtools view -@ 4 -h -f 2 -F 12 -F 256 -F 2048 \
  "${OUT_DIR}/${SAMPLE_NAME}_namesorted.sam" \
  -o "${OUT_DIR}/${SAMPLE_NAME}_filtered.sam"

# sort again
samtools sort -n -@ 4 -O SAM "${OUT_DIR}/${SAMPLE_NAME}_filtered.sam" \
  -o "${OUT_DIR}/${SAMPLE_NAME}_filtered_namesorted.sam"

# convert to BAM
samtools view -@ 4 -b \
  -o "${OUT_DIR}/${SAMPLE_NAME}_filtered_namesorted.bam" \
  "${OUT_DIR}/${SAMPLE_NAME}_filtered_namesorted.sam"