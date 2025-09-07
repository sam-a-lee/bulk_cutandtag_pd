#!/bin/bash
#SBATCH --time=48:00:00 # cpu partition max time per kcl policy
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 # could specify more here if wanted 
#SBATCH --mem=4G # shouldnt require more than 2 gb for all files ss
#SBATCH --hint=nomultithread # prefer physical cores for better throughput
#SBATCH --job-name=multiqc
#SBATCH --output=/scratch/prj/bcn_pd_pesticides/Files-From-Imperial/analysis_cutandtag_pd_bulk/data_out/02_fastqc/logs/multiqc_%j.out
#SBATCH --error=/scratch/prj/bcn_pd_pesticides/Files-From-Imperial/analysis_cutandtag_pd_bulk/data_out/02_fastqc/logs/multiqc_%j.err

#---------#
# purpose #
#---------#

# this files creates a summary of all fastqc files to make qc easier and faster

#--------------------#
# set up environment #
#--------------------#

# change dir to where conda envs are 
cd /users/k2587336

# load shell
source ~/.bashrc

# activate conda env
module load py-multiqc/1.15-gcc-13.2.0-python-3.11.6

OUT_DIR="/scratch/prj/bcn_pd_pesticides/Files-From-Imperial/analysis_cutandtag_pd_bulk/data_out/02_fastqc"
IN_DIR="/scratch/prj/bcn_pd_pesticides/Files-From-Imperial/analysis_cutandtag_pd_bulk/data_out/02_fastqc"

#---------------------------------#
# run multiqc in fastqc directory #
#---------------------------------#

multiqc \
    --outdir ${OUT_DIR} \
    --verbose \
    ${IN_DIR}

