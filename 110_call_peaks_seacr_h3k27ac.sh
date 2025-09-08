#!/bin/bash
#SBATCH --time=48:00:00 # cpu partition max time
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 
#SBATCH --mem=4G
#SBATCH --job-name=seacr_h3k27ac
#SBATCH --output=/scratch/prj/bcn_pd_pesticides/Files-From-Imperial/analysis_cutandtag_pd_bulk/data_out/11_seacr/h3k27ac/logs/seacr_h3k27ac_%A_%a.out
#SBATCH --error=/scratch/prj/bcn_pd_pesticides/Files-From-Imperial/analysis_cutandtag_pd_bulk/data_out/11_seacr/h3k27ac/logs/seacr_h3k27ac_%A_%a.err
#SBATCH --array=0-27 # !!! ADJUST ME AS NEEDED BASED ON FILE NUMBER !!! 

#--------------------#
# set up environment #
#--------------------#

# change dir to where conda envs are 
cd /users/k2587336

# load shell
source ~/.bashrc

# load seacr conda env 
source activate seacr

# root working directory
ROOT_DIR="/scratch/prj/bcn_pd_pesticides/Files-From-Imperial/analysis_cutandtag_pd_bulk/data_out"

# dir with input bed files
IN_DIR="${ROOT_DIR}/10_bedgraph"

# out dir for peak calling files
OUT_DIR="${ROOT_DIR}/11_seacr/h3k27ac"

#--------------------------------#
# create array of bedgraph files #
#---------------------------------#

# get list of bam files of fragments without duplicates
BEDGRAPH_FILES=($(find "${IN_DIR}" -maxdepth 1 -type f -name "*_dup_removed_clean_fragments.bedgraph" | sort))

# sanity
# get num of files found and echo
N=${#BEDGRAPH_FILES[@]}
echo "[INFO] Found $N bedgraph files in $IN_DIR"

# get sample based on sample list and array index
FILE="${BEDGRAPH_FILES[$SLURM_ARRAY_TASK_ID]}"

# get sample name by stripping basename and file ext
SAMPLE_NAME=$(basename "$FILE" | cut -d'_' -f1)

echo "Processing sample: ${SAMPLE_NAME}"

#---------------------------------#
# h3k27ac peak calling with seacr #
#---------------------------------#

# set threshold for seacr
# based on doi.org/10.1038/s41467-025-58137-2
THRESHOLD="0.01"

## call peaks
SEACR_1.3.sh \
    ${FILE} \
    $THRESHOLD \
    non \
    stringent \
    ${OUT_DIR}/${SAMPLE_NAME}_seacr_stringent_${THRESHOLD}_peaks