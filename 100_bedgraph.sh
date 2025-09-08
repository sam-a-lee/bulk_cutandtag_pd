#!/bin/bash 
#SBATCH --time=48:00:00 # cpu partition max time
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 # bedtools not multithreaded
#SBATCH --mem=4G
#SBATCH --job-name=bedgraph
#SBATCH --output=/scratch/prj/bcn_pd_pesticides/Files-From-Imperial/analysis_cutandtag_pd_bulk/data_out/10_bedgraph/logs/bedgraph_%A_%a.out
#SBATCH --error=/scratch/prj/bcn_pd_pesticides/Files-From-Imperial/analysis_cutandtag_pd_bulk/data_out/10_bedgraph/logs/bedgraph_%A_%a.err
#SBATCH --array=0-27 # !!! ADJUST ME AS NEEDED BASED ON FILE NUMBER !!! 

#--------------------#
# set up environment #
#--------------------#

# change dir to where conda envs are 
cd /users/k2587336

# load shell
source ~/.bashrc

# load bedtools conda env 
source activate bedtools

# set root dir for working  
ROOT_DIR="/scratch/prj/bcn_pd_pesticides/Files-From-Imperial/analysis_cutandtag_pd_bulk/data_out"

# dir containing input files
IN_DIR="${ROOT_DIR}/08_bed"

# output dir for fragment length files
OUT_DIR="${ROOT_DIR}/10_bedgraph"

# specify location of hg19 chromosome sizes 
CHROM_SIZE="/scratch/prj/bcn_pd_pesticides/Files-From-Imperial/analysis_cutandtag_pd_bulk/resources/ref_genome/hg38.chrom.sizes"

#------------------------------------#
# create array for list of bed files #
#------------------------------------#

# get list of bam files of fragments without duplicates
FRAG_BED_FILES=($(find "${IN_DIR}" -maxdepth 1 -type f -name "*_dup_removed_mapped_reads_clean_fragments.bed" | sort))

# sanity
# get num of files found and echo
N=${#FRAG_BED_FILES[@]}
echo "[INFO] Found $N BAMs in $IN_DIR"

# get sample based on sample list and array index
FRAG_BED_FILE="${FRAG_BED_FILES[$SLURM_ARRAY_TASK_ID]}"

# get sample name by stripping basename and file ext
SAMPLE_NAME=$(basename "$FRAG_BED_FILE" | cut -d'_' -f1)

echo "Processing sample: ${SAMPLE_NAME}"

#------------------#
# create bed graph #
#------------------# 

bedtools genomecov -bg \
    -i ${FRAG_BED_FILE} \
    -g ${CHROM_SIZE} > \
    ${OUT_DIR}/${SAMPLE_NAME}_dup_removed_clean_fragments.bedgraph