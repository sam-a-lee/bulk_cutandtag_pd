#!/bin/bash 
#SBATCH --time=48:00:00 # cpu partition max time
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --hint=nomultithread 
#SBATCH --job-name=seq_align
#SBATCH --array=0-30 # adjust based on number of samples
#SBATCH --output=/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/02_aligned/seq_align_%A_%a.out
#SBATCH --error=/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/02_aligned/seq_align_%A_%a.err

#---------# 
# purpose #
#---------#

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
IN_DIR="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/04_trimmed"

# directory to save output files to
OUT_DIR="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/05_aligned"

# directory where reference genome (hg19) is
REF="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/resources/ref_genome/grch38_primary_assembly_index"

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
    --dovetail \
    --phred33 \
    -I 10 \
    -X 700 \
    -p 8 \
    -x ${REF} \
    -1 ${R1} \
    -2 ${R2} \
    -S ${OUT_DIR}/${SAMPLE}_bowtie2_local_vsensitive_nomixed_nodiscord.sam \
    &> ${OUT_DIR}/logs/${SAMPLE}_bowtie2_local_vsensitive_nomixed_nodiscord_summary.txt

# local = local alignment
# very sensitivte =  -D 20 -R 3 -N 0 -L 20 -i S,1,0.50 
# where -D = max number of mismatches in a string; 
# -R = for reads with repetitive seeds, try X seeds;
# -L= length of seed substring; 
# -i =  interval between seed substrings w/r/t read length
# no mixed = suppress unpaired alignments for paired reads
# no-discordant = suppress discordant alignments for paired reads
# dovetail = concordant when mates extend past each other (per nf-core cutandrun pipeline)
# -I = minimum read length
# -X = maximum fragment length (10-20ish above max expected length)
# -p = threads
