#!/bin/bash -l
#SBATCH --time=48:00:00 # cpu partition max time
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --hint=nomultithread 
#SBATCH --job-name=bwamem2
#SBATCH --array=0-29 # adjust based on number of samples
#SBATCH --output=/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/040_aligned/logs/040_bwamem2_%A_%a.out
#SBATCH --error=/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/040_aligned/logs/040_bwamem2_%A_%a.err

#---------# 
# purpose #
#---------#

# this script aligns trimmed reads to a reference genome (in this case human hg19)

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

# activate conda env
conda activate bwa-mem2

# load samtools
module load samtools/1.17-gcc-13.2.0-python-3.11.6

# directory where trimmed files are located
IN_DIR="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/030_trimmed"

# directory to save output files to
OUT_DIR="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/040_aligned"

# directory where reference genome (hg19) is
REF="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/resources/ref_genome/Homo_sapiens.GRCh38.dna.chromosomes_only.fa"

# memory samtools to use for sorting 
SAMTOOLS_SORT_MEM="1G" 

#---------------------------------------------------# 
# create list of samples to get corresponding files #
#---------------------------------------------------#

# get sorted list of unique trimmed R1 fastqs (sorted for consistency)
SAMPLES=($(find "${IN_DIR}" -maxdepth 1 -type f -name "*_R1_001_cutadapt.fastq.gz" | sort))

# pick sample for this array index
R1=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
R2="${R1/_R1_001_cutadapt.fastq.gz/_R2_001_cutadapt.fastq.gz}"

# extract clean sample name (strip suffixes)
SAMPLE=$(basename "$R1" | sed -E 's/_R[12]_001_cutadapt\.fastq\.gz$//')

#-----------------------------------#
# align samples to ref with bowtie2 #
#-----------------------------------#

echo "Processing sample: ${SAMPLE}"

bwa-mem2 mem -t 8 \
  "${REF}" \
  "${R1}" \
  "${R2}" \
  > "${OUT_DIR}/${SAMPLE}_bwamem2.sam"

echo "Done: ${SAMPLE}"


