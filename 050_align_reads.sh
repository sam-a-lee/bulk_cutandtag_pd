#!/bin/bash 
#SBATCH --time=48:00:00 # cpu partition max time
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --hint=nomultithread 
#SBATCH --job-name=seq_align
#SBATCH --array=0-27 # adjust based on number of samples
#SBATCH --output=/scratch/users/k2587336/test/data_out/02_aligned/seq_align_%A_%a.out
#SBATCH --error=/scratch/users/k2587336/test/data_out/02_aligned/seq_align_%A_%a.err

##### PURPOSE #####

# this script aligns trimmed reads to a reference genome (in this case human hg19)


#--------------------# 
# set up environment #
#--------------------#

# change dir to where conda envs are 
cd /users/k2587336

# load shell
source ~/.bashrc

# activate conda env
source activate bowtie2

# directory where trimmed files are located
IN_DIR="/scratch/prj/bcn_pd_pesticides/Files-From-Imperial/analysis_cutandtag_pd_bulk/data_out/04_trimmed"

# directory to save output files to
OUT_DIR="/scratch/prj/bcn_pd_pesticides/Files-From-Imperial/analysis_cutandtag_pd_bulk/data_out/05_aligned"

# directory where reference genome (hg19) is
REF="/scratch/prj/bcn_pd_pesticides/Files-From-Imperial/analysis_cutandtag_pd_bulk/resources/ref_genome/hg38_indexed"

# memory samtools to use for sorting 
SAMTOOLS_SORT_MEM="1G" 

#---------------------------------------------------# 
# create list of samples to get corresponding files #
#---------------------------------------------------#

# get sorted list of unique trimmed R1 fastqs (sorted for consistency)
SAMPLES=($(find "${IN_DIR}" -maxdepth 1 -type f -name "*_R1_001_val_1.fq.gz" | sort))

# pick sample for this array index
R1=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
R2="${R1/_R1_001_val_1.fq.gz/_R2_001_val_2.fq.gz}"

# extract clean sample name (strip suffixes)
SAMPLE=$(basename "$R1" | sed -E 's/_R[12]_001_val_[12]\.fq\.gz$//')

#-----------------------------------#
# align samples to ref with bowtie2 #
#-----------------------------------#

echo "Processing sample: ${sample}"

# run bowtie2
bowtie2 --local \
    --very-sensitive \
    --no-mixed \
    --no-discordant \
    --phred33 \
    -I 10 \
    -X 700 \
    -p 8 \
    -x ${REF} \
    -1 ${R1} \
    -2 ${R2} \
    -S ${OUT_DIR}/${SAMPLE}_bowtie2_local_vsensitive_nomixed_nodiscord.sam \
    &> ${OUT_DIR}/logs/${SAMPLE}_bowtie2_local_vsensitive_nomixed_nodiscord_summary.txt
