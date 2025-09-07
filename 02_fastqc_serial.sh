#!/bin/bash
#SBATCH --time=48:00:00 # cpu partition max time per kcl policy
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 # could specify more here if wanted 
#SBATCH --mem=2G # shouldnt require more than a gb per file
#SBATCH --hint=nomultithread # prefer physical cores for better throughput
#SBATCH --job-name=fastqc
#SBATCH --output=/scratch/prj/bcn_pd_pesticides/Files-From-Imperial/analysis_cutandtag_pd_bulk/data_out/02_fastqc/logs/fastqc_%j.out
#SBATCH --error=/scratch/prj/bcn_pd_pesticides/Files-From-Imperial/analysis_cutandtag_pd_bulk/data_out/02_fastqc/logs/fastqc_%j.err

#---------#
# purpose #
#---------#

# this files examines quality of raw fastq files prior to merging of technical duplicates

#--------------------#
# set up environment #
#--------------------#

# change dir to where conda envs are 
cd /users/k2587336

#source ~/conda/etc/profile.d/conda.sh
source ~/.bashrc 

# activate fastqc module
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

# find all forward (R1) and reverse (R2) read files for those SAMPLE_IDs from lanes one (L001) and two (L002)
# sort to ensure consistent order
ALL_SAMPLES=($(find "${IN_DIR}" -type f \( -false $(for id in "${SAMPLE_ID[@]}"; do printf " -o -name '*%s_S[0-9]*_L00[12]_R[12]_001.fastq.gz'" "$id"; done) \) | sort))

#---------------------------#
# fastqc on raw fastq files #
#---------------------------#

for SAMPLE in "${ALL_SAMPLES[@]}"; do
  fastqc \
      -o ${OUT_DIR} \
      ${SAMPLE}
done
