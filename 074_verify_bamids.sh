#!/bin/bash -l
#SBATCH --time=48:00:00
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --hint=nomultithread
#SBATCH --job-name=verify_bamid
#SBATCH --array=0-29
#SBATCH --output=/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/07_filtered_reads/logs/verify_bamid_%A_%a.out
#SBATCH --error=/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/07_filtered_reads/logs/verify_bamid_%A_%a.err

#---------#
# purpose #
#---------#

# examine dna contamination in individual samples based on snp proporitons with verifybamid2

#--------------------#
# environment set up #
#--------------------#

# set up conda (not through bashrc)

CONDA_ROOT="/software/spackages_v0_21_prod/apps/linux-ubuntu22.04-zen2/gcc-13.2.0/anaconda3-2022.10-5wy43yh5crcsmws4afls5thwoskzarhe"
if [ -f "${CONDA_ROOT}/etc/profile.d/conda.sh" ]; then
  . "${CONDA_ROOT}/etc/profile.d/conda.sh"
else
  export PATH="${CONDA_ROOT}/bin:$PATH"
  eval "$(${CONDA_ROOT}/bin/conda shell.bash hook)"
fi

# activate custom env
conda activate verifybamid2

# dir with bam files
BAM_DIR="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/07_filtered_reads"

# out dir for results
OUT_DIR="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/07_filtered_reads/verify_bamid"

#-----------------#
# index bam files #
#-----------------#

# Collect BAM paths 
mapfile -t BAM_PATHS < <(find "${BAM_DIR}" -maxdepth 1 -type f -name "*_filtered_namesorted.bam" -print | sort)

# select the sample for this array task
BAM="${BAM_PATHS[$SLURM_ARRAY_TASK_ID]}"

# get sample name (remove basename and extension)
SAMPLE_NAME=$(basename "$BAM" | cut -d'_' -f1)

#---------------#
# verify bamids #
#---------------#

cd $OUT_DIR 

verifybamid2 \
  --SVDPrefix $(VERIFY_BAM_ID_HOME)/resource/1000g.100k.b38.vcf.gz.dat \
  --Reference "/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/resources/verifybamid_ref//GRCh38_full_analysis_set_plus_decoy_hla.fa" \
  --BamFile ${BAM}