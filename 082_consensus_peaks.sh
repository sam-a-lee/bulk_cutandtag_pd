#!/bin/bash
#SBATCH --time=48:00:00 # cpu partition max time
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --job-name=consensus_peaks
#SBATCH --output=/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/08_peaks/logs/consensus_peaks_%j.out
#SBATCH --error=/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/08_peaks/logs/consensus_peaks_%j.err

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

# load macs conda env 
conda activate macs3

# load module
module load samtools

# dir with named sorted bams to merge
BAM_DIR="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/07_filtered_reads"

# outdir for consensus peaks
PEAK_DIR="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/08_peaks/consensus_peaks"

#---------------------#
# merge all bam files #
#---------------------#

# make list of bam files
ls "${BAM_DIR}"/*.bam > "${BAM_DIR}/bam.list"

# Merge; writes a coordinate-sorted BAM if inputs are also coord-sorted.
samtools merge -@ 8 -b "${BAM_DIR}/bam.list" -o "${BAM_DIR}/all_samples.bam" 

samtools sort -@ 8 -o "${BAM_DIR}/all_samples_namesorted.bam" "${BAM_DIR}/all_samples.bam" 

#--------------------#
# macs3 peak calling #
#--------------------#

macs3 callpeak \
    -t "${BAM_DIR}/all_samples_namesorted.bam" \
    -f BAMPE \
    -g hs \
    -q 0.00001 \
    --bw 401 \
    --keep-dup all \
    --verbose 3 \
    --nolambda \
    --nomodel \
    --outdir "${PEAK_DIR}" \
    -n all_samples_macs2_q_0_00001_h3k27ac