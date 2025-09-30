#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --job-name=consensus_peaks_all
#SBATCH --output=/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/08_peaks/logs/consensus_peaks_all_%j.out
#SBATCH --error=/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/08_peaks/logs/consensus_peaks_all_%j.err

#-------------#
# Environment #
#-------------#

# initiate conda 
CONDA_ROOT="/software/spackages_v0_21_prod/apps/linux-ubuntu22.04-zen2/gcc-13.2.0/anaconda3-2022.10-5wy43yh5crcsmws4afls5thwoskzarhe"
if [ -f "${CONDA_ROOT}/etc/profile.d/conda.sh" ]; then
  . "${CONDA_ROOT}/etc/profile.d/conda.sh"
else
  export PATH="${CONDA_ROOT}/bin:$PATH"
  eval "$(${CONDA_ROOT}/bin/conda shell.bash hook)"
fi

# activate custom env
conda activate macs3

# load module 
module load samtools

# dir with all bams
BAM_DIR="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/07_filtered_reads"

# dir for writing consensus peaks
PEAK_DIR="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/08_peaks/consensus_peaks"

# tmp dir for processing
TMP_DIR="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/08_peaks/tmp_$$"

# mkdir if they dont exist
mkdir -p "${PEAK_DIR}" "${TMP_DIR}"

# requested cpus
CPUS="${SLURM_CPUS_PER_TASK:-8}"

# date for output logs
ts(){ date "+%Y-%m-%d %H:%M:%S"; }

echo "[INFO $(ts)] Building ALL-sample consensus peaks from ${BAM_DIR}"

# Prevent too-many-open-files issues
ulimit -n 4096 2>/dev/null || true

#--------------#
# Collect BAMs #
#--------------#

# get list of bam files
ALL_LIST="${TMP_DIR}/bams_all.list"

# create array
find "${BAM_DIR}" -maxdepth 1 -type f -name "*_filtered_namesorted.bam" | sort > "${ALL_LIST}" || true

# sanity checks
TOTAL=$(wc -l < "${ALL_LIST}" || echo 0)
if [ "${TOTAL}" -eq 0 ]; then
  echo "[ERROR $(ts)] No BAMs matching *_filtered_namesorted.bam found in ${BAM_DIR}."
  exit 1
fi
echo "[INFO $(ts)] Found ${TOTAL} BAM(s)."

#---------------------#
# Merge (name-sorted) #
#---------------------#

MERGED_BAM="${PEAK_DIR}/all_samples_merged_namesorted.bam"

echo "[INFO $(ts)] Merging with samtools merge -n → ${MERGED_BAM}"

samtools merge -n -@ "${CPUS}" -b "${ALL_LIST}" -o "${MERGED_BAM}"

#--------------------#
# MACS3 callpeak     #
#--------------------#

echo "[INFO $(ts)] Calling MACS3 (BAMPE) on merged BAM"

macs3 callpeak \
  -t "${MERGED_BAM}" \
  -f BAMPE \
  -g hs \
  -q 0.00001 \
  --bw 401 \
  --keep-dup all \
  --nolambda \
  --nomodel \
  --verbose 3 \
  --outdir "${PEAK_DIR}" \
  -n "consensus_all_samples_macs3_q_0_00001"

echo "[INFO $(ts)] Done. Outputs in: ${PEAK_DIR}"
