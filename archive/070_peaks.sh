#!/bin/bash
#SBATCH --time=48:00:00 # cpu partition max time
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 
#SBATCH --mem=4G
#SBATCH --job-name=macs3
#SBATCH --output=/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/070_peaks/logs/07_macs3_%A_%a.out
#SBATCH --error=/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/070_peaks/logs/07_macs3_%A_%a.err
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

# load macs conda env 
conda activate macs3

# dir with input bed files
IN_DIR="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/060_filtered"

# out dir for peak calling files
OUT_DIR="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/070_peaks"

mkdir -p ${OUT_DIR}
mkdir -p "${OUT_DIR}/logs"

#-----------------------------#
# create array of BAMPE files #
#-----------------------------#

# get list of bam files of fragments without duplicates
BAMPE_FILES=($(find "${IN_DIR}" -maxdepth 1 -type f -name "*_filtered_namesorted.bam" | sort))

# get sample based on sample list and array indexs
FILE="${BAMPE_FILES[$SLURM_ARRAY_TASK_ID]}"

# get sample name by stripping basename and file ext
SAMPLE_NAME=$(basename "$FILE" | cut -d'_' -f1)

echo "Processing sample: ${SAMPLE_NAME}"

#---------------------------------#
# h3k27ac peak calling with macs3 #
#---------------------------------#

# macs threshold based on doi.org/10.1038/s41467-025-58137-2
# dont adjust bandwidth - will call peaks from cellranger counts data using macs3
macs3 callpeak \
    -t ${FILE} \
    -f BAMPE \
    -g hs \
    -q 0.00001 \
    --keep-dup all \
    --verbose 3 \
    --nolambda \
    --nomodel \
    --outdir ${OUT_DIR} \
    -n ${SAMPLE_NAME}_macs3_q_0_00001_h3k27ac