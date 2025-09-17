#!/bin/bash
#SBATCH --time=48:00:00 # cpu partition max time
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --hint=nomultithread 
#SBATCH --job-name=sam_to_bam
#SBATCH --output=/scratch/prj/bcn_pd_pesticides/Files-From-Imperial/analysis_cutandtag_pd_bulk/data_out/07_bam/logs/sam_to_bam_%A_%a.out
#SBATCH --error=/scratch/prj/bcn_pd_pesticides/Files-From-Imperial/analysis_cutandtag_pd_bulk/data_out/07_bam/logs/sam_to_bam_%A_%a.err
#SBATCH --array=0-27 # !!! change this as needed

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

# root working directory
ROOT_DIR="/scratch/prj/bcn_pd_pesticides/Files-From-Imperial/analysis_cutandtag_pd_bulk/data_out"

# directory of input files
IN_DIR="${ROOT_DIR}/06_duplicates_removed"


#----------------------------#
# array list and file naming #
#----------------------------#

# get list of samples without duplicates
SAMPLES_IN=($(find "${IN_DIR}" -maxdepth 1 -type f -name "*_picard_dup_removed.sam" | sort))

# get sample from samples list based on array index
SAMPLE="${SAMPLES_IN[$SLURM_ARRAY_TASK_ID]}"

# get sample name by stripping basename and file ext
SAMPLE_NAME=$(basename "$SAMPLE" | cut -d'_' -f1)

# naming of output file based on sample and suffix
OUT_BAM="$ROOT_DIR/07_bam/${SAMPLE_NAME}_samtools_dup_removed_mapped_reads.bam"


#----------------#
# convert to bam #
#----------------#

# output BAM directly from samtools sort and keep only mapped reads, using 3 threads (+ 1 for I/O)
samtools sort -n -@ 4 -O BAM "$SAMPLE" \
| samtools view -b -F 0x04 -@ 4 - \
> "$OUT_BAM"
