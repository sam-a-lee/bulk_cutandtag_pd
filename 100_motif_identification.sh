#!/bin/bash
#SBATCH --time=48:00:00 # cpu partition max time
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 
#SBATCH --mem=4G
#SBATCH --job-name=homer_motifs
#SBATCH --output=/scratch/prj/bcn_pd_pesticides/Files-From-Imperial/analysis_cutandtag_pd_bulk/data_out/12_motifs/h3k27ac/logs/homer_motifs_h3k27ac_%A_%a.out
#SBATCH --error=/scratch/prj/bcn_pd_pesticides/Files-From-Imperial/analysis_cutandtag_pd_bulk/data_out/12_motifs/h3k27ac/logs/homer_motifs_h3k27ac_%A_%a.err
#SBATCH --array=0-27 # !!! ADJUST ME AS NEEDED BASED ON FILE NUMBER !!! 


#--------------------#
# set up environment #
#--------------------#

# change dir to where conda envs are 
cd /users/k2587336

# load shell
source ~/.bashrc

# load seacr conda env 
source activate homer

# root working directory
ROOT_DIR="/scratch/prj/bcn_pd_pesticides/Files-From-Imperial/analysis_cutandtag_pd_bulk/data_out"

# dir with input bed files
IN_DIR="${ROOT_DIR}/11_seacr/h3k27ac"

# out dir for peak calling files
OUT_DIR="${ROOT_DIR}/12_motifs/h3k27ac"

# ref genome file


#-----------------------------------#
# create array of h3k27ac bed files #
#-----------------------------------#

# get list of bam files of fragments without duplicates
PEAK_FILES=($(find "${IN_DIR}" -maxdepth 1 -type f -name "*_seacr_stringent_0.01_peaks.stringent.bed" | sort))

# sanity
# get num of files found and echo
N=${#PEAK_FILES[@]}
echo "[INFO] Found $N peak ped files in $IN_DIR"

# get sample based on sample list and array index
FILE="${PEAK_FILES[$SLURM_ARRAY_TASK_ID]}"

# get sample name by stripping basename and file ext
SAMPLE_NAME=$(basename "$FILE" | cut -d'_' -f1)

echo "Processing sample: ${SAMPLE_NAME}"

#----------------------------------#
# motify identification with homer #
#----------------------------------#

findMotifsGenome.pl \
    ${FILE} \
    hg38 \
    ${OUT_DIR}/${SAMPLE_NAME} \
    -size 1000
