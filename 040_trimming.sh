#!/bin/bash
#SBATCH --time=48:00:00 # cpu partition max time
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --job-name=trim_adapters
#SBATCH --array=0-27 # !!! modify to number of unique samples being processing !!!
#SBATCH --output=/scratch/prj/bcn_pd_pesticides/Files-From-Imperial/analysis_cutandtag_pd_bulk/data_out/04_trimmed/logs/trim_adapters_%A_%a.out
#SBATCH --error=/scratch/prj/bcn_pd_pesticides/Files-From-Imperial/analysis_cutandtag_pd_bulk/data_out/04_trimmed/logs/trim_adapters_%A_%a.err

##### PURPOSE #####

# this file uses fastp to quality trim and remove illumina adaptors from raw sequencing reads (fastq files)
#--------------------# 
# set up environment #
#--------------------#

# change dir to where conda envs are 
cd /users/k2587336

# load shell
source ~/.bashrc

# activate conda env
source activate fastp

# root dir where raw fastq files are stored
ROOT_DIR="/scratch/prj/bcn_pd_pesticides/Files-From-Imperial/analysis_cutandtag_pd_bulk/data_out/03_merged_fastq"

# directory to save files
OUT_DIR="/scratch/prj/bcn_pd_pesticides/Files-From-Imperial/analysis_cutandtag_pd_bulk/data_out/04_trimmed"

#---------------------------------------------------# 
# create list of files corresponding to each sample #
#---------------------------------------------------#

# get sorted list of unique raw samples (R1 files)
SAMPLES=($(find "${ROOT_DIR}" -maxdepth 1 -type f -name "*_R1_001.fastq.gz" | sort))

# pick sample based on job array index
R1=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
R2="${R1/_R1_001.fastq.gz/_R2_001.fastq.gz}"

# get sample names for r1 and r2
SAMPLE_R1=$(basename "${R1}" "_001.fastq.gz")
SAMPLE_R2=$(basename "${R2}" "_001.fastq.gz")

# create output names
OUT_R1="${OUT_DIR}/${SAMPLE_R1}_fastp_polyg_20bp.fq.gz"
OUT_R2="${OUT_DIR}/${SAMPLE_R2}_fastp_polyg_20bp.fq.gz"
#---------------------------------# 
# trim adapters using trim galore #
#---------------------------------#

echo "Processing sample: ${R1} and ${R2}"

# run trim galore
# default trimming below phred 20
#trim_galore --gzip \
#    --paired \
#    --fastqc \
#    --output_dir ${OUT_DIR} \
#    ${R1} \
#    ${R2}

fastp \
    --length_required 20 \
    --trim_poly_g \
    -i ${R1} \
    -I ${R2} \
    -o ${OUT_R1} \
    -O ${OUT_R2}