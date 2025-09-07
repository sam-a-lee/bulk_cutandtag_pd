#!/bin/bash 
#SBATCH --time=48:00:00 # cpu partition max time
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --job-name=make_index
#SBATCH --output=/scratch/prj/bcn_pd_pesticides/Files-From-Imperial/analysis_cutandtag_pd_bulk/resources/ref_genome/logs/make_index_job_%j.out
#SBATCH --error=/scratch/prj/bcn_pd_pesticides/Files-From-Imperial/analysis_cutandtag_pd_bulk/resources/ref_genome/logs/make_index_job_%j.err

# before running this script download the appropriate reference files using wget (or something similar)
# in the directory containing the reference, i have noted the version and date of download in a readme


#--------------------# 
# set up environment #
#--------------------#

# change dir to where conda envs are 
cd /users/k2587336

# load shell
source ~/.bashrc

# activate bowtie2 conda env
source activate bowtie2


#-------------# 
# build index #
#-------------#

bowtie2-build --threads 8 \
    /scratch/prj/bcn_pd_pesticides/Files-From-Imperial/analysis_cutandtag_pd_bulk/resources/ref_genome/hg38.fa \
    /scratch/prj/bcn_pd_pesticides/Files-From-Imperial/analysis_cutandtag_pd_bulk/resources/ref_genome/hg38_indexed
