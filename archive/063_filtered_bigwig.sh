#!/bin/bash 
#SBATCH --time=48:00:00 # cpu partition max time
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 # bedtools not multithreaded
#SBATCH --mem=4G
#SBATCH --job-name=bedgraph
#SBATCH --output=/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/08_fragments/logs/bedgraph_%A_%a.out
#SBATCH --error=/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/08_fragments/logs/bedgraph_%A_%a.err
#SBATCH --array=0-29 # !!! ADJUST ME AS NEEDED BASED ON FILE NUMBER !!! 

#--------------------#
# set up environment #
#--------------------#

# load bedtools module 
module load bedtools2/2.31.0-gcc-12.3.0-python-3.11.6

# dir containing input files
IN_DIR="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/08_fragments"

# output dir for fragment length files
OUT_DIR="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/08_fragments/bedgraphs"

# specify location of hg19 chromosome sizes 
CHROM_SIZE="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/resources/ref_genome/Homo_sapiens.GRCh38.dna.primary_assembly.chrom.sizes"

#------------------------------------#
# create array for list of bed files #
#------------------------------------#
# get list of bam files of fragments without duplicates
FRAG_FILES=($(find "${IN_DIR}" -maxdepth 1 -type f -name "*_fragments_sorted.bed" | sort))

# get sample based on sample list and array index
FILE="${FRAG_FILES[$SLURM_ARRAY_TASK_ID]}"

# get sample name by stripping basename and file ext
SAMPLE_NAME=$(basename "$FILE" | cut -d'_' -f1)

echo "Processing sample: ${SAMPLE_NAME}"

#------------------#
# create bed graph #
#------------------# 

bedtools genomecov -bg \
    -i ${FILE} \
    -g ${CHROM_SIZE} > \
    ${OUT_DIR}/${SAMPLE_NAME}_fragments_sorted.bedgraph