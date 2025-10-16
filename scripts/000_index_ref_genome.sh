#!/bin/bash 
#SBATCH --time=48:00:00 # cpu partition max time
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --job-name=make_index
#SBATCH --output=/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/resources/ref_genome/logs/make_index_job_%j.out
#SBATCH --error=/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/resources/ref_genome/logs/make_index_job_%j.err

# before running this script download the appropriate reference files using wget (or something similar)
# in the directory containing the reference, i have noted the version and date of download in a readme

#---------# 
# purpose #
#---------#

# this script indexes the desired reference genome

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

# activate bowtie2 conda env
# source activate bowtie2
conda activate bwa-mem2

#-------------# 
# build index #
#-------------#

#bowtie2-build --threads 8 \
#    /scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/resources/ref_genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
#   /scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/resources/ref_genome/grch38_primary_assembly_index \

bwa-mem2 index /scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/resources/ref_genome/Homo_sapiens.GRCh38.dna.chromosomes_only.fa