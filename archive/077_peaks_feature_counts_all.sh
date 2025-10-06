#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --job-name=feature_counts_all
#SBATCH --output=/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/070_peaks/logs/077_feature_counts_all_%A_%a.out
#SBATCH --error=/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/070_peaks/logs/077_feature_counts_all_%A_%a.err
#SBATCH --array=0-29
#-------------#
# Environment #
#-------------#

CONDA_ROOT="/software/spackages_v0_21_prod/apps/linux-ubuntu22.04-zen2/gcc-13.2.0/anaconda3-2022.10-5wy43yh5crcsmws4afls5thwoskzarhe"
if [ -f "${CONDA_ROOT}/etc/profile.d/conda.sh" ]; then
  . "${CONDA_ROOT}/etc/profile.d/conda.sh"
else
  export PATH="${CONDA_ROOT}/bin:$PATH"
  eval "$(${CONDA_ROOT}/bin/conda shell.bash hook)"
fi

# load module 
module load samtools/1.17-gcc-13.2.0-python-3.11.6

# load custom bedtools env
conda activate bedtools

# dir with filtered bam files
BAM_DIR="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/060_filtered"

# dir with consensus bam file and peaks
CONSENSUS_PEAKS="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/070_peaks/073_peaks_consensus_all/075_peaks_consensus_all_blacklist_removed/consensus_macs3_q_0_00001_peaks_clean.narrowPeak"

# out dir for samples vs consensus
OUT_DIR="/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_bulk/data_out/070_peaks/077_peaks_feature_count_all"

mkdir -p "${OUT_DIR}"

# cpus requested
CPUS="${SLURM_CPUS_PER_TASK:-4}"

# foramtting for log files
ts(){ date "+%Y-%m-%d %H:%M:%S"; }

echo "[INFO $(ts)] Per-sample counts vs ALL-consensus peaks (no sample sheet)"

#----------------------------#
# Discover BAMs & sample IDs #
#----------------------------#

# Collect BAM paths and aligned sample IDs
mapfile -t BAM_PATHS < <(find "${BAM_DIR}" -maxdepth 1 -type f -name "*_filtered_namesorted.bam" -print | sort)

# Pick the one file for this array task (0-based)
IDX=${SLURM_ARRAY_TASK_ID}

BAM_IN="${BAM_PATHS[$IDX]}"

# Derive sample ID from filename
BASENAME="${BAM_IN##*/}"
SAMPLE_NAME="${BASENAME%_filtered_namesorted.bam}"

# Make coord-sorted temp (only for this sample)
BAM_CS="${BAM_DIR}/${SAMPLE_NAME}_filtered_coordsorted.bam"

[[ -s "$BAM_CS" ]] || samtools sort -@ "$CPUS" -o "$BAM_CS" "$BAM_IN"

# Count vs ALL-consensus peaks
bedtools multicov -bams "${BAM_CS}" -bed "${CONSENSUS_PEAKS}" > "${OUT_DIR}/${SAMPLE_NAME}_feature_counts_all_multicov.tsv"
