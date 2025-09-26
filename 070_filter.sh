#!/bin/bash
#SBATCH --time=48:00:00 # cpu partition max time
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --hint=nomultithread 
#SBATCH --job-name=filter_reads
#SBATCH --output=/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/07_filtered_reads/logs/filter_reads_%A_%a.out
#SBATCH --error=/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/07_filtered_reads/logs/filter_reads_%A_%a.err
#SBATCH --array=0-30 # !!! change this as needed

#---------# 
# purpose #
#---------#

# this script converts .sam files to .bam files
# only mapped reads are retained in the output .bam files

#--------------------# 
# set up environment #
#--------------------#

# change dir to where conda envs are 
cd /users/k2587336

# load shell
source ~/.bashrc

# load samtools module (v1.17-gcc-13.2.0-python-3.11.6_ 
# libncurse 5.0 req by conda samtools by 6.5 installed
module load samtools 

source activate bedtools 

# root working directory
ROOT_DIR="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out"

# directory of input files
IN_DIR="${ROOT_DIR}/06_duplicates_removed"

#out dir
OUT_DIR="${ROOT_DIR}/07_filtered_reads"

#----------------------------#
# array list and file naming #
#----------------------------#

# get list of samples without duplicates
mapfile -t SAMPLES_IN < <(find "$IN_DIR" -maxdepth 1 -type f -name "*_picard_dup_rm.sam" | sort)

# get sample from samples list based on array index
SAMPLE="${SAMPLES_IN[$SLURM_ARRAY_TASK_ID]}"

# get sample name by stripping basename and file ext
SAMPLE_NAME=$(basename "$SAMPLE" | cut -d'_' -f1)

#---------------------------#
# filter and convert to bam #
#---------------------------#

#  name sort reads
samtools sort -n -@ 4 -O SAM "${SAMPLE}" -o "${OUT_DIR}/${SAMPLE_NAME}_namesorted.sam"

#samtools view -@ 4 -h -f 2 -F 12 -F 256 -F 2048 \
# "${OUT_DIR}/${SAMPLE_NAME}_namesorted.sam" \
# -o "${OUT_DIR}/${SAMPLE_NAME}_filtered.sam"

# filter name-sorted reads in a pair-aware way
# -f 2        require properly paired reads (both mates mapped; includes dovetail for bwamem2)
# -F 4        drop read unmapped
# -F 8        drop mate unmapped
# -F 256      drop secondary alignments
# -F 2048     drop supplementary alignments
# drop reads pairs where one or both are mapped to MT
# drop read pairs where one or both have MAPQ <30
# drop read pairs where reads are on different chromosomes
# drop read pairs where one or more has multiple alignments (only want primary)
samtools view -h -@ 4 \
  -f 2 \
  -F 4 \
  -F 8 \
  -F 256 \
  -F 2048 \
  "${OUT_DIR}/${SAMPLE_NAME}_namesorted.sam" \
| awk 'BEGIN{
         OFS="\t"
       }
       function is_mito(chr,   lc){
         lc=tolower(chr)
         return (lc=="chrm" || lc=="mt" || lc=="m" || lc=="chrmt")
       }
  /^@/ { print; next }
  {
    if ($1 != qn) {
      qn = $1; a = $0
      split($0, A, "\t")
      a_chr  = A[3]
      a_mapq = A[5]+0
      a_xa   = (a ~ /\tXA:Z:/)
      a_sa   = (a ~ /\tSA:Z:/)
      a_nh   = match(a, /\tNH:i:([0-9]+)/, m) ? m[1]+0 : -1
      next
    } else {
      b = $0
      split($0, B, "\t")
      b_chr  = B[3]
      b_mapq = B[5]+0
      b_xa   = (b ~ /\tXA:Z:/)
      b_sa   = (b ~ /\tSA:Z:/)
      b_nh   = match(b, /\tNH:i:([0-9]+)/, n) ? n[1]+0 : -1

      passA = (a_mapq>=30) && (!a_xa) && (!a_sa) && (a_nh==-1 || a_nh==1) && (!is_mito(a_chr))
      passB = (b_mapq>=30) && (!b_xa) && (!b_sa) && (b_nh==-1 || b_nh==1) && (!is_mito(b_chr))
      sameC = (a_chr == b_chr)

      if (passA && passB && sameC) { print a; print b }
    }
  }' \
| samtools view -@ 4 -b -o "${OUT_DIR}/${SAMPLE_NAME}_filtered.sam" -

# sort again
samtools sort -n -@ 4 -O SAM "${OUT_DIR}/${SAMPLE_NAME}_filtered.sam" \
  -o "${OUT_DIR}/${SAMPLE_NAME}_filtered_namesorted.sam"

# convert to BAM
samtools view -@ 4 -b \
  -o "${OUT_DIR}/${SAMPLE_NAME}_filtered_namesorted.bam" \
  "${OUT_DIR}/${SAMPLE_NAME}_filtered_namesorted.sam"