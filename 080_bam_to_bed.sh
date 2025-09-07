#!/bin/bash
#SBATCH --time=48:00:00 # cpu partition max time
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 # bedtools not multithreaded 
#SBATCH --mem=2G
#SBATCH --job-name=bam_to_bed
#SBATCH --output=/scratch/prj/bcn_pd_pesticides/Files-From-Imperial/analysis_cutandtag_pd_bulk/data_out/08_bed/logs/bam_to_bed_%A_%a.out
#SBATCH --error=/scratch/prj/bcn_pd_pesticides/Files-From-Imperial/analysis_cutandtag_pd_bulk/data_out/08_bed/logs/bam_to_bed_%A_%a.err
#SBATCH --array=0-27 # change as needed

##### PURPOSE #####

# this script converts .bam files to .bed files
# after conversion, paired reads are filtered for those on the same chromosome and those that are less than 1000 bp
# after filtering, fragment-related columns are extracted into a new .bed file


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

# activate conda bedtools env
source activate bedtools

# root working directory
ROOT_DIR="/scratch/prj/bcn_pd_pesticides/Files-From-Imperial/analysis_cutandtag_pd_bulk/data_out"

# dir with bam files
BAM_DIR="${ROOT_DIR}/07_bam"

# out dir for bed files
OUT_DIR="${ROOT_DIR}/08_bed"

#------------------------------------------------# 
# create master list of files without duplicates #
#------------------------------------------------#

# get list of bam files 
mapfile -t BAM_IN < <(find "$BAM_DIR" -maxdepth 1 -type f -name "*_samtools_dup_removed_mapped_reads.bam" | sort)

# sanity
# get num of files found and echo
N=${#BAM_IN[@]}
echo "[INFO] Found $N BAMs in $BAM_DIR"

# fail gracefully if no files found or more files than specified in array
if (( N == 0 )); then echo "[WARN] No inputs; exiting."; exit 0; fi
if (( SLURM_ARRAY_TASK_ID >= N )); then echo "[INFO] No work for task $SLURM_ARRAY_TASK_ID (N=$N)"; exit 0; fi

# get file
SAMPLE="${BAM_IN[$SLURM_ARRAY_TASK_ID]}"
# get basename
BN=$(basename -- "$SAMPLE")
# keep full sample prefix by stripping the exact suffix
SAMPLE_NAME="${BN%_samtools_dup_removed_mapped_reads.bam}"

echo "[INFO] SAMPLE: $SAMPLE"
echo "[INFO] SAMPLE_NAME: $SAMPLE_NAME"

#--------------------------------------#
# convert to get and extract fragments #
#--------------------------------------#

# Convert bam -> bed 
samtools sort -n -@ 1 -o - "$SAMPLE" \
  | bedtools bamtobed -bedpe -i - \
  > "$OUT_DIR/${SAMPLE_NAME}_dup_removed_mapped_reads.bed"

# same chromosome + fragment length < 1000 bp
awk '$1==$4 && $6-$2 < 1000' \
  "$OUT_DIR/${SAMPLE_NAME}_dup_removed_mapped_reads.bed" \
  > "$OUT_DIR/${SAMPLE_NAME}_dup_removed_mapped_reads_clean.bed"

# extract fragment intervals (chrom, start, end), sorted
cut -f1,2,6 "$OUT_DIR/${SAMPLE_NAME}_dup_removed_mapped_reads_clean.bed" \
  | sort -k1,1 -k2,2n -k3,3n \
  > "$OUT_DIR/${SAMPLE_NAME}_dup_removed_mapped_reads_clean_fragments.bed"

