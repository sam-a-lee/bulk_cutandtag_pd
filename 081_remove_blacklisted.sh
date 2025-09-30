#!/bin/bash
#SBATCH --time=48:00:00 # cpu partition max time
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 
#SBATCH --mem=4G
#SBATCH --job-name=rm_blacklisted
#SBATCH --output=/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/08_peaks/logs/rm_blacklisted_%A_%a.out
#SBATCH --error=/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/08_peaks/logs/rm_blacklisted_%A_%a.err
#SBATCH --array=0-29 # !!! ADJUST ME AS NEEDED BASED ON FILE NUMBER !!! 

#--------------------#
# set up environment #
#--------------------#

# initiate conda
CONDA_ROOT="/software/spackages_v0_21_prod/apps/linux-ubuntu22.04-zen2/gcc-13.2.0/anaconda3-2022.10-5wy43yh5crcsmws4afls5thwoskzarhe"
if [ -f "${CONDA_ROOT}/etc/profile.d/conda.sh" ]; then
  . "${CONDA_ROOT}/etc/profile.d/conda.sh"
else
  export PATH="${CONDA_ROOT}/bin:$PATH"
  eval "$(${CONDA_ROOT}/bin/conda shell.bash hook)"
fi

# load custom env
conda activate bedtools

# directory with macs3 peaks
PEAKS_IN="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/08_peaks/raw_peaks"

# output dir
OUT_DIR="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/08_peaks/clean_peaks"

# blacklist file
# clean file has "chr" stripped from chromosome number to match macs3 peak bed files
BLACKLIST="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/resources/nordin_hg38_problematic_clean.bed"

#------------------------#
# index files into array #
#------------------------#

# assaign to array
mapfile -t PEAK_FILES < <(find "${PEAKS_IN}" -maxdepth 1 -type f -name "*_macs3_q_0_00001_h3k27ac_peaks.narrowPeak" | sort)

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
  > "${OUT_DIR}/${SAMPLE_NAME}_macs3_q_0_00001_h3k27ac_peaks_clean.narrowPeak"
