#!/usr/bin/env bash
#SBATCH --account=jiameng
#SBATCH --partition=gpu3090,gpu4090
#SBATCH --qos=4gpus
#SBATCH --gpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2-00:00:00
#SBATCH --job-name=nr004-hs-basecall
#SBATCH --output=/gpfs/work/bio/zhixingwu16/NanoRemix_fused_tool/validation_RNA004/preprocess_human_pipeline/logs/%x-%A_%a.out
#SBATCH --error=/gpfs/work/bio/zhixingwu16/NanoRemix_fused_tool/validation_RNA004/preprocess_human_pipeline/logs/%x-%A_%a.err
#SBATCH --array=0-5%4

set -euo pipefail
source "/gpfs/work/bio/zhixingwu16/NanoRemix_fused_tool/validation_RNA004/preprocess_human_pipeline/config.sh"
module load cuda/11.8
module load samtools

split="split${SLURM_ARRAY_TASK_ID}"
input="$SAMPLE_ROOT/pod5/$split"
output="$SAMPLE_ROOT/basecall/$split.bam"
legacy_output="$output.tmp"
marker="$SAMPLE_ROOT/.basecall.$split.complete"

test -f "$SAMPLE_ROOT/.extract.complete"
test -d "$input"
test -d "$DORADO_MODEL"
test -x "$DORADO"
if [[ ! -f "$marker" ]]; then
    if [[ -s "$legacy_output" ]] && samtools quickcheck -q -u "$legacy_output"; then
        echo "Recovering complete unaligned BAM from $legacy_output"
        mv "$legacy_output" "$output"
    elif [[ -e "$legacy_output" ]]; then
        echo "ERROR: incomplete legacy basecall output exists: $legacy_output" >&2
        exit 1
    elif [[ -s "$output" ]] && samtools quickcheck -q -u "$output"; then
        echo "Recovering complete unaligned BAM from $output"
    else
        tmp="$output.tmp.${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
        test ! -e "$tmp"
        "$DORADO" basecaller \
            --device cuda:0 \
            --batchsize 512 \
            --chunksize 2000 \
            --emit-moves \
            --recursive \
            "$DORADO_MODEL" "$input" > "$tmp"
        samtools quickcheck -v -u "$tmp"
        mv "$tmp" "$output"
    fi
    samtools quickcheck -v -u "$output"
    samtools view -c "$output" > "$SAMPLE_ROOT/basecall/$split.read_count.txt"
    touch "$marker"
fi
