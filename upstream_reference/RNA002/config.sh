#!/usr/bin/env bash
# Central configuration for the HEK293_WT/HMEC_WT RNA002 NanoRemix workflow.
# Edit only this file after copying the whole pipeline directory to the HPC.

set -euo pipefail

PIPELINE_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

# Inputs are never modified. Every generated file is written under RUN_ROOT.
FOUR_DATASET_ROOT="${FOUR_DATASET_ROOT:-/gpfs/work/bio/zhixingwu16/Nanopore/four_dataset}"
RUN_ROOT="${RUN_ROOT:-/gpfs/work/bio/zhixingwu16/Nanopore/nanoremix5mer_m6A_HEK_HMEC}"
REFERENCE_FASTA="${REFERENCE_FASTA:-/gpfs/work/bio/yuxinzhang17/DirectMultiRM/Homo_sapiens.GRCh38.cdna.all.fa}"

SAMPLES=(HEK293_WT HMEC_WT)

# RNA002 basecalling and Tombo. These paths were taken from the existing,
# read-only HEK/HMEC and canonical NanoRemix scripts on the cluster.
GUPPY_BASECALLER="${GUPPY_BASECALLER:-/gpfs/work/bio/yuxinzhang17/software/ont-guppy/bin/guppy_basecaller}"
FLOWCELL="${FLOWCELL:-FLO-MIN106}"
KIT="${KIT:-SQK-RNA002}"
FAST5_API_ENV="${FAST5_API_ENV:-/gpfs/work/bio/yuxinzhang17/environment/.conda/envs/ont-env2}"
MULTI_TO_SINGLE="${MULTI_TO_SINGLE:-${FAST5_API_ENV}/bin/multi_to_single_fast5}"
TOMBO_ENV="${TOMBO_ENV:-/gpfs/work/bio/yuxinzhang17/environment/.conda/envs/tombo1.0-env}"
TOMBO="${TOMBO:-${TOMBO_ENV}/bin/tombo}"
GRAPHMAP="${GRAPHMAP:-/gpfs/work/bio/yuxinzhang17/software/graphmap/bin/Linux-x64/graphmap}"
TOMBO_BASECALL_GROUP="${TOMBO_BASECALL_GROUP:-Basecall_1D_000}"

# The original 5-mer extractor used this ONT environment. The replacement
# extractor only requires Python, NumPy and h5py, and preserves the original
# signal reversal plus median/MAD normalization.
FEATURE_ENV="${FEATURE_ENV:-/gpfs/work/bio/jiayili1905/environment/.conda/envs/ont-env}"
FEATURE_PYTHON="${FEATURE_PYTHON:-${FEATURE_ENV}/bin/python}"

# NanoRemix inference and downstream Python environment.
NANOREMIX_ENV="${NANOREMIX_ENV:-/gpfs/work/bio/zhixingwu16/.conda/envs/nanoremix_model}"
NANOREMIX5="${NANOREMIX5:-${NANOREMIX_ENV}/bin/nanoremix5mer}"
POSTPROCESS_PYTHON="${POSTPROCESS_PYTHON:-${NANOREMIX_ENV}/bin/python}"
NANOREMIX_MODELS="${NANOREMIX_MODELS:-/gpfs/work/bio/zhixingwu16/nanoRemix_benchmark/5mer_tool/5mer_models}"

# Primary inference specification. Do not choose these after looking at the
# final HEK/HMEC results.
NANOREMIX_STRATEGY="${NANOREMIX_STRATEGY:-strategy3}"
NANOREMIX_CONDITION="${NANOREMIX_CONDITION:-full}"
READ_M6A_THRESHOLD="${READ_M6A_THRESHOLD:-0.5}"
SITE_MIN_READS="${SITE_MIN_READS:-10}"
SITE_MOD_RATIO_THRESHOLD="${SITE_MOD_RATIO_THRESHOLD:-0.5}"
KDE_BANDWIDTH="${KDE_BANDWIDTH:-0.04}"
PLOT_BINS="${PLOT_BINS:-50}"

# For a transcriptome FASTA, valid DRS reads are expected on the transcript
# sense (+). Keeping this at 1 prevents antisense alignment artefacts from
# entering the biological 5'-to-3' density plot.
PLUS_STRAND_ONLY="${PLUS_STRAND_ONLY:-1}"

# Existing m6Anet outputs are optional comparison inputs, never NanoRemix
# inputs. HMEC is included only if the file exists when the plotting job runs.
HEK_M6ANET_SITE="${HEK_M6ANET_SITE:-/gpfs/work/bio/zhixingwu16/Nanopore/HEK293WT_m6anet/results/m6anet/data.site_proba.csv}"
HMEC_M6ANET_SITE="${HMEC_M6ANET_SITE:-/gpfs/work/bio/zhixingwu16/Nanopore/HMEC_m6anet/results/m6anet/data.site_proba.csv}"
M6ANET_PROB_THRESHOLD="${M6ANET_PROB_THRESHOLD:-0.9}"
M6ANET_MIN_READS="${M6ANET_MIN_READS:-10}"

# Scheduler settings copied from the partitions/QoS used by the existing jobs.
GPU_QOS="${GPU_QOS:-4gpus}"
GPU_PARTITION="${GPU_PARTITION:-gpu3090}"
CPU_QOS="${CPU_QOS:-32cores}"
CPU_PARTITION="${CPU_PARTITION:-cpu6348}"
HUGE_QOS="${HUGE_QOS:-46cores}"
HUGE_PARTITION="${HUGE_PARTITION:-cpuhuge1t}"

sample_from_index() {
    local index="${1:?array index is required}"
    if (( index < 0 || index >= ${#SAMPLES[@]} )); then
        echo "Invalid sample index: ${index}" >&2
        return 2
    fi
    printf '%s\n' "${SAMPLES[$index]}"
}

sample_input_dir() {
    local sample="${1:?sample is required}"
    printf '%s/%s/multi_fast5\n' "$FOUR_DATASET_ROOT" "$sample"
}

sample_run_dir() {
    local sample="${1:?sample is required}"
    printf '%s/%s\n' "$RUN_ROOT" "$sample"
}
