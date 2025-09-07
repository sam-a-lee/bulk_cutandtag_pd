#!/bin/bash 
#SBATCH --time=48:00:00 # cpu partition max time
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 # bedtools not multithreaded
#SBATCH --mem=16G
#SBATCH --job-name=bedgraph
#SBATCH --output=/scratch/users/k2587336/test/data_out/08_bedgraph/bedgraph_%A_%a.out
#SBATCH --error=/scratch/users/k2587336/test/data_out/07_bedgraph/bedgraph_%A_%a.err
#SBATCH --array=0-1 # !!! ADJUST ME AS NEEDED BASED ON FILE NUMBER !!! 


##### PURPOSE #####

#--------------------#
# set up environment #
#--------------------#

# change dir to where conda envs are 
cd /users/k2587336

# load shell
source ~/.bashrc

# load bedtools conda env 
source activate bedtools

# specify root dir (working dir)
ROOT_DIR="/scratch/users/k2587336/test/data_out"

# specify location of hg19 chromosome sizes 
CHROM_SIZE="/scratch/users/k2587336/test/hg19/hg19.chrom.sizes"

#------------------------------------#
# create array for list of bed files #
#------------------------------------#

# get list of bam files of fragments without duplicates
FRAG_BED_FILES=($(find "${ROOT_DIR}/05_bed" -maxdepth 1 -type f -name "*_dup_removed_mapped_reads_clean_fragments.bed" | sort))

# get sample based on sample list and array index
FRAG_BED_FILE="${FRAG_BED_FILES[$SLURM_ARRAY_TASK_ID]}"

# get sample name by stripping basename and file ext
SAMPLE_NAME=$(basename "$FRAG_BED_FILE" | cut -d'_' -f1)

echo "Processing sample: ${SAMPLE_NAME}"

#------------------#
# create bed graph #
#------------------# 

bedtools genomecov -bg \
    -i ${FRAG_BED_FILE} \
    -g $chromSize > \
    ${ROOT_DIR}/08_begraph/${SAMPLE_NAME}_dup_removed_clean_fragments.bedgraph