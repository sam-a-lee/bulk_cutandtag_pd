#!/bin/bash
#SBATCH --time=48:00:00 # cpu partition max time
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --job-name=feature_counts
#SBATCH --output=/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/08_peaks/logs/feature_counts_%j.out
#SBATCH --error=/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/08_peaks/logs/feature_counts_%j.err


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

conda activate bedtools

BAM_IN="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/07_filtered_reads"

BED_IN="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/08_peaks/consensus_peaks/all_samples_macs2_q_0_00001_h3k27ac_peaks.narrowPeak"





bedtools multicov -bams s1.bam s2.bam ... -bed union_peaks.bed > counts.tsv