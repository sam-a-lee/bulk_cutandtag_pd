#!/bin/bash 
#SBATCH --time=48:00:00 # cpu partition max time
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --job-name=rm_dups
#SBATCH --output=/scratch/prj/bcn_pd_pesticides/Files-From-Imperial/analysis_cutandtag_pd_bulk/data_out/06_duplicates_removed/logs/remove_dups_%j.out
#SBATCH --error=/scratch/prj/bcn_pd_pesticides/Files-From-Imperial/analysis_cutandtag_pd_bulk/data_out/06_duplicates_removed/logs/remove_dups_%j.err
#SBATCH --array=0-27 # !!! change to correct number

# NB: picard is single threaded
# garbage collection can parallelize - up to 4 threads for max returns
# total threads requested = 5 (1 sort + 4 gc)

##### PURPOSE #####

# this script takes trimmed sequences in .sam format, sorts them and then either
# marks duplicate reads without moving them
# or removes duplicate reads
# two output files are generated in different file locations depending on whether 
# duplicates are marked or discarded
# from here, reads without duplicates should be used for processing (https://doi.org/10.1038/s41467-025-58137-2)


#--------------------# 
# set up environment #
#--------------------#

# change dir to where conda envs are 
cd /users/k2587336

# load shell
source ~/.bashrc

# activate conda env
source activate picard

# directory containing input files (trimmed + sorted sam files )
IN_DIR="/scratch/prj/bcn_pd_pesticides/Files-From-Imperial/analysis_cutandtag_pd_bulk/data_out/05_aligned"

# directory to store output files in
OUT_DIR="/scratch/prj/bcn_pd_pesticides/Files-From-Imperial/analysis_cutandtag_pd_bulk/data_out/06_duplicates_removed"

# temporary directory
TEMP_DIR="/scratch/prj/bcn_pd_pesticides/Files-From-Imperial/analysis_cutandtag_pd_bulk/data_out/06_duplicates_removed/temp"


#---------------------------------------------#
# create list of samples for array and naming #
#---------------------------------------------#

# get list of unique trimmed + aligned samples
# NB: files names must be explicit or this will become infinite loop
SAMPLES=($(find "${IN_DIR}" -maxdepth 1 -type f -name "*_bowtie2_local_vsensitive_nomixed_nodiscord.sam"))

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

# sort by coordinate
picard SortSam \
    I=$SAMPLE \
    O=${OUT_DIR}/${SAMPLE_NAME}_picard_sorted.sam \
    SORT_ORDER=coordinate \
    TMP_DIR=$TEMP_DIR

# mark duplicates (don't remove)
picard MarkDuplicates \
    I=${OUT_DIR}/${SAMPLE_NAME}_picard_sorted.sam \
    O=${OUT_DIR}/${SAMPLE_NAME}_picard_dup_marked.sam \
    TMP_DIR=$TEMP_DIR \
    METRICS_FILE=${OUT_DIR}/${SAMPLE_NAME}_picard_dup_marked_summary.txt 

# remove duplicates
picard MarkDuplicates \
    I=${OUT_DIR}/${SAMPLE_NAME}_picard_sorted.sam \
    O=${OUT_DIR}/${SAMPLE_NAME}_picard_dup_removed.sam \
    TMP_DIR=$TEMP_DIR \
    METRICS_FILE=${OUT_DIR}/${SAMPLE_NAME}_picard_dup_removed_summary.txt \
    REMOVE_DUPLICATES=true
