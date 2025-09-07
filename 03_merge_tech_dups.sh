#!/bin/bash
#SBATCH --time=48:00:00 
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 
#SBATCH --mem=2G
#SBATCH --job-name=merge_tech_dups
#SBATCH --array=0-91%20 # modify based on number of files - %20 throttles to 20 at a time to be polite
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
IN_DIR="/rds/prj/bcn_marzi_lab/Files-From-Imperial/data_cutandtag_pd_bulk/AACNGMFHV"

# output directory for merged files
OUT_DIR="/scratch/prj/bcn_pd_pesticides/Files-From-Imperial/analysis_cutandtag_pd_bulk/data_out/03_merged_fastq"

# sample sheet with sample IDs corresponding to each sample and cell type
SAMPLE_SHEET="/scratch/prj/bcn_pd_pesticides/Files-From-Imperial/analysis_cutandtag_pd_bulk/resources/metadata/sample_sheet.txt"

#-----------------------------------#
# get list of libraries in lane one #
#-----------------------------------#

# read first column (tab-delimited) as SAMPLE_IDs 
SAMPLE_ID=($(awk -F'\t' 'NF && $1!~/^#/{print $1}' "$SAMPLE_SHEET"))

# get list of all samples in lane 1
ALL_L001_SAMPLES=($(find "${IN_DIR}" -type f \( -false $(for id in "${SAMPLE_ID[@]}"; do printf " -o -name '*%s_S[0-9]*_L001_R[12]_001.fastq.gz'" "$id"; done) \) | sort))

# get lane one sample per array
L001_SAMPLE=${ALL_L001_SAMPLES[$SLURM_ARRAY_TASK_ID]}

# get lane two sample based on lane one
L002_SAMPLE="${L001_SAMPLE/_L001_/_L002_}"

# get new sample id from lane one sample by removing lane info 
NO_LANES_ID=$(basename -- "${L001_SAMPLE}" | sed -E 's/_L00[12]_/_/')

#-----------------------------#
# merge lane one and lane two #
#-----------------------------#

cat ${L001_SAMPLE} ${L002_SAMPLE} > ${OUT_DIR}/${NO_LANES_ID}
