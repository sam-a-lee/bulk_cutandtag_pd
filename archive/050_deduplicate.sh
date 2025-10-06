#!/bin/bash 
#SBATCH --time=48:00:00 # cpu partition max time
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --job-name=picard_dedup
#SBATCH --output=/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/050_deduplicated/logs/050_picard_dedup_%j.out
#SBATCH --error=/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/050_deduplicated/logs/050_picard_dedup_%j.err
#SBATCH --array=0-29 # !!! change to correct number

# NB: picard is single threaded
# garbage collection can parallelize - up to 4 threads for max returns
# total threads requested = 5 (1 sort + 4 gc)

#---------# 
# purpose #
#---------#

# this script takes trimmed sequences in .sam format, sorts them and then either
# marks duplicate reads without moving them
# or removes duplicate reads
# two output files are generated in different file locations depending on whether 
# duplicates are marked or discarded
# from here, reads without duplicates should be used for processing (https://doi.org/10.1038/s41467-025-58137-2)

#--------------------# 
# set up environment #
#--------------------#

# load picard module
module load picard/2.26.2-gcc-13.2.0

# directory containing input files (trimmed + sorted sam files )
IN_DIR="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/040_aligned"

# directory to store output files in
OUT_DIR="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/050_deduplicated"

# temporary directory for file processing 
TEMP_DIR="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/050_deduplicated/temp"

mkdir -p ${OUT_DIR} 
mkdir -p "${OUT_DIR}/logs"
mkdir -p ${TEMP_DIR}

#---------------------------------------------#
# create list of samples for array and naming #
#---------------------------------------------#

# get list of unique trimmed + aligned samples
# NB: files names must be explicit or this will become infinite loop
SAMPLES=($(find "${IN_DIR}" -maxdepth 1 -type f -name "*_bwamem2.sam"))

# select the sample for this array task
SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"

# get sample name (remove basename and extension)
SAMPLE_NAME=$(basename "$SAMPLE" | cut -d'_' -f1)

# set java environment options (memory, gc threads, gc time limit)
# only need to do once per shell
export _JAVA_OPTIONS="-Xmx12G -XX:ParallelGCThreads=1 -XX:GCTimeLimit=50"

#-----------------------------------------#
# mark and remove duplicates using picard #
#-----------------------------------------#

# create tmp file for coordinate sorting (required for dedup)
# and remove after dedup
TMP_SORTED="$(mktemp -p "$TMP_DIR" "${SAMPLE_NAME}_picard_sorted.XXXXXX.sam")"

cleanup() {
  rm -f "$TMP_SORTED"
}

# clean temp files on exit (success or failure)
trap cleanup EXIT INT TERM

# coordinate sort sam fil
picard SortSam \
    I=$SAMPLE \
    O=${TMP_SORTED} \
    SORT_ORDER=coordinate \
    TMP_DIR=$TEMP_DIR 

# remove duplicates
picard MarkDuplicates \
    I=${TMP_SORTED} \
    O=${OUT_DIR}/${SAMPLE_NAME}_picard_dedup.sam \
    TMP_DIR=$TEMP_DIR \
    METRICS_FILE=${OUT_DIR}/${SAMPLE_NAME}_picard_dedup_summary.txt \
    REMOVE_DUPLICATES=true
