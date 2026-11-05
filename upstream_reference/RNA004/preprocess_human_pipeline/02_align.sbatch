#!/usr/bin/env bash
#SBATCH --account=jiameng
#SBATCH --partition=cpuhuge1t
#SBATCH --qos=46cores
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=2-00:00:00
#SBATCH --job-name=nr004-hs-align
#SBATCH --output=/gpfs/work/bio/zhixingwu16/NanoRemix_fused_tool/validation_RNA004/preprocess_human_pipeline/logs/%x-%A_%a.out
#SBATCH --error=/gpfs/work/bio/zhixingwu16/NanoRemix_fused_tool/validation_RNA004/preprocess_human_pipeline/logs/%x-%A_%a.err
#SBATCH --array=0-5%4

set -euo pipefail
source "/gpfs/work/bio/zhixingwu16/NanoRemix_fused_tool/validation_RNA004/preprocess_human_pipeline/config.sh"
module load samtools

split="split${SLURM_ARRAY_TASK_ID}"
input="$SAMPLE_ROOT/basecall/$split.bam"
output="$SAMPLE_ROOT/aligned/${split}_sorted.bam"
marker="$SAMPLE_ROOT/.align.$split.complete"

test -f "$SAMPLE_ROOT/.basecall.$split.complete"
test -s "$REFERENCE_FASTA"
if [[ ! -f "$marker" ]]; then
    "$DORADO" aligner "$REFERENCE_FASTA" "$input" \
        | samtools view -@ "${SLURM_CPUS_PER_TASK}" -bh -F 2308 -q 10 - \
        | samtools sort -@ "${SLURM_CPUS_PER_TASK}" -o "$output.tmp" -
    samtools quickcheck -v "$output.tmp"
    mv "$output.tmp" "$output"
    samtools index -@ "${SLURM_CPUS_PER_TASK}" "$output"
    {
        printf 'all_basecalled\t'; samtools view -c "$input"
        printf 'primary_mapped_mapq10\t'; samtools view -c "$output"
    } > "$SAMPLE_ROOT/aligned/$split.read_counts.tsv"
    touch "$marker"
fi
