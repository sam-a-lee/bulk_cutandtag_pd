#!/bin/bash
#SBATCH --time=48:00:00 # cpu partition max time
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 
#SBATCH --mem=4G
#SBATCH --job-name=rm_blacklisted
#SBATCH --output=/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/070_peaks/logs/075_rm_blacklisted_%A_%a.out
#SBATCH --error=/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/070_peaks/logs/075_rm_blacklisted_%A_%a.err
#SBATCH --array=0-29 # !!! ADJUST ME AS NEEDED BASED ON FILE NUMBER !!! 

#--------------------#
# set up environment #
#--------------------#

# load bedtools module 
module load bedtools2/2.31.0-gcc-12.3.0-python-3.11.6

# directory with macs3 peaks
PEAKS_IN="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/070_peaks/073_peaks_consensus_all"

# output dir
OUT_DIR="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/070_peaks/073_peaks_consensus_all/075_peaks_consensus_all_blacklist_removed"

mkdir -p ${OUT_DIR}

# blacklist file
# clean file has "chr" stripped from chromosome number to match macs3 peak bed files
BLACKLIST="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/resources/nordin_hg38_problematic_clean.bed"

#------------------------#
# index files into array #
#------------------------#

# assaign to array
mapfile -t PEAK_FILES < <(find "${PEAKS_IN}" -maxdepth 1 -type f -name "*_macs3_q_0_00001_peaks.narrowPeak" | sort)

# get file path
PEAKS="${PEAK_FILES[$SLURM_ARRAY_TASK_ID]}"

# extract sample name
SAMPLE_NAME="$(basename "$PEAKS" | cut -d'_' -f1)"

echo "Removing blacklisted peaksf from sample: ${SAMPLE_NAME}"

#--------------------------#
# remove blacklisted peaks #
#--------------------------#

# both bed files are already sorted by chr and start pos

bedtools intersect -v \
  -a "${PEAKS}" \
  -b "${BLACKLIST}" \
  > "${OUT_DIR}/${SAMPLE_NAME}_macs3_q_0_00001_peaks_clean.narrowPeak"
