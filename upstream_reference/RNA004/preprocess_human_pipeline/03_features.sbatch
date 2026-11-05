#!/usr/bin/env bash
#SBATCH --account=jiameng
#SBATCH --partition=cpuhuge1t
#SBATCH --qos=46cores
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=7-00:00:00
#SBATCH --job-name=nr004-hs-feature
#SBATCH --output=/gpfs/work/bio/zhixingwu16/NanoRemix_fused_tool/validation_RNA004/preprocess_human_pipeline/logs/%x-%A_%a.out
#SBATCH --error=/gpfs/work/bio/zhixingwu16/NanoRemix_fused_tool/validation_RNA004/preprocess_human_pipeline/logs/%x-%A_%a.err
#SBATCH --array=0-5%4

set -euo pipefail
source "/gpfs/work/bio/zhixingwu16/NanoRemix_fused_tool/validation_RNA004/preprocess_human_pipeline/config.sh"

split="split${SLURM_ARRAY_TASK_ID}"
marker="$SAMPLE_ROOT/.features.$split.complete"

test -f "$SAMPLE_ROOT/.align.$split.complete"
test -s "$REGIONS"
test -s "$KMER_LEVELS"
if [[ ! -f "$marker" ]]; then
    "$REMORA_PYTHON" "$FEATURE_SCRIPT" \
        --pod5_dir "$SAMPLE_ROOT/pod5" \
        --bam "$SAMPLE_ROOT/aligned" \
        --reg "$REGIONS" \
        --level "$KMER_LEVELS" \
        --output "$SAMPLE_ROOT/features" \
        --threads "${SLURM_CPUS_PER_TASK}" \
        --run-id "${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}" \
        --splits "${SLURM_ARRAY_TASK_ID}" \
        --kmer 9 \
        --step 1 \
        --resquiggle True
    test -s "$SAMPLE_ROOT/features/$split.npz"
    test -s "$SAMPLE_ROOT/features/$split.csv"
    rows="$(awk 'END {print NR - 1}' "$SAMPLE_ROOT/features/$split.csv")"
    test "$rows" -gt 0
    printf 'feature_rows\t%s\n' "$rows" > "$SAMPLE_ROOT/features/$split.read_counts.tsv.tmp"
    mv "$SAMPLE_ROOT/features/$split.read_counts.tsv.tmp" \
        "$SAMPLE_ROOT/features/$split.read_counts.tsv"
    touch "$marker"
fi
