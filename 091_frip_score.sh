#!/bin/bash
#SBATCH --time=48:00:00 # cpu partition max time
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 
#SBATCH --mem=4G
#SBATCH --job-name=frip_score
#SBATCH --output=/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/09_peaks/logs/frip_score_%A_%a.out
#SBATCH --error=/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/09_peaks/logs/frip_score_%A_%a.err
#SBATCH --array=0-29 # !!! ADJUST ME AS NEEDED BASED ON FILE NUMBER !!! 

#--------------------#
# set up environment #
#--------------------#

# change dir to where conda envs are 
cd /users/k2587336

# load shell
load module miniforge3

# load macs conda env 
mamba activate deeptools

# dir with input bed files
IN_DIR="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/09_peaks"

# out dir for peak calling files
OUT_DIR="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/09_peaks/frip_scores"

#--------------------------------#
# create array of BAMPE files #
#---------------------------------#

# get list of bam files of fragments without duplicates
BAM_FILES=($(find "${IN_DIR}" -maxdepth 1 -type f -name "*_filtered_namesorted.bam" | sort))

# get sample based on sample list and array indexs
FILE="${BAMPE_FILES[$SLURM_ARRAY_TASK_ID]}"

# get sample name by stripping basename and file ext
SAMPLE_NAME=$(basename "$FILE" | cut -d'_' -f1)

echo "Processing sample: ${SAMPLE_NAME}"

#---------------------------------#
# h3k27ac peak calling with seacr #
#---------------------------------#

# macs threshold based on doi.org/10.1038/s41467-025-58137-2
# bandwith = 401 for consistency with cellranger-atac seq
macs3 callpeak \
    -t ${FILE} \
    -f BAMPE \
    -g hs \
    -q 0.00001 \
    --bw 401 \
    --keep-dup all \
    --verbose 3 \
    --nolambda \
    --nomodel \
    --outdir ${OUT_DIR} \
    -n ${SAMPLE_NAME}_macs2_q_0_00001_h3k27ac