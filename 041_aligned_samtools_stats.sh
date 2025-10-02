#!/bin/bash -l
#SBATCH --time=48:00:00
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --hint=nomultithread
#SBATCH --job-name=samtools_stats
#SBATCH --array=0-29
#SBATCH --output=/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/040_aligned/logs/041_samtools_stats_%A_%a.out
#SBATCH --error=/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/040_aligned/logs/041_samtools_stats_%A_%a.err

#--------------------#
# set up environment #
#--------------------#

# load module
module load samtools/1.17-gcc-13.2.0-python-3.11.6

# in dir with bams
IN_DIR="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/040_aligned"

# out dir
OUT_DIR="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/040_aligned/041_samtools_stats"

mkdir -p ${OUT_DIR}

#------------------------------#
# get BAM list and pick sample #
#------------------------------#

# assign bam, files to array 
mapfile -t SAMPLES < <(find "$IN_DIR" -maxdepth 1 -type f -name "*_bwamem2.sam" | sort)

# get sample from array
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

# get file name
SAMPLE_BASENAME=$(basename "$SAMPLE")

# get sample name
SAMPLE_NAME="${SAMPLE_BASENAME%.sam}"

echo "Processing sample: ${SAMPLE}"

#----------#
# samtools #
#----------#

# stats
samtools stats "${SAMPLE}" \
  > "${OUT_DIR}/${SAMPLE_NAME}_samtools_stats.txt"

echo "Done: ${SAMPLE_NAME}"