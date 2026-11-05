#!/usr/bin/env bash

set -euo pipefail

PIPELINE_DIR="${NANOREMIX_PIPELINE_DIR:-${SLURM_SUBMIT_DIR:-}}"
[[ -f "$PIPELINE_DIR/config.sh" ]] || {
    echo "cannot locate pipeline config.sh; PIPELINE_DIR=$PIPELINE_DIR" >&2
    exit 1
}
# shellcheck source=../config.sh
source "$PIPELINE_DIR/config.sh"

sample="$(sample_from_index "${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}")"
input_dir="$(sample_input_dir "$sample")"
run_dir="$(sample_run_dir "$sample")"
out_dir="$run_dir/01_guppy_rna002"
success="$out_dir/_SUCCESS"

if [[ -f "$success" ]]; then
    echo "[$sample] basecalling already complete: $success"
    exit 0
fi
if [[ -d "$out_dir" ]] && [[ -n "$(find "$out_dir" -mindepth 1 -print -quit)" ]]; then
    echo "[$sample] refusing to mix with a partial basecalling directory: $out_dir" >&2
    exit 1
fi

mkdir -p "$out_dir"

"$GUPPY_BASECALLER" \
    -i "$input_dir" \
    -s "$out_dir" \
    --flowcell "$FLOWCELL" \
    --kit "$KIT" \
    --recursive \
    --fast5_out \
    --device cuda:0 \
    --num_callers "${SLURM_CPUS_PER_TASK:-1}" \
    --cpu_threads_per_caller 1

if [[ -z "$(find "$out_dir/workspace" -type f -name '*.fast5' -print -quit)" ]]; then
    echo "[$sample] Guppy completed but no FAST5 was found in $out_dir/workspace" >&2
    exit 1
fi
touch "$success"
echo "[$sample] RNA002 Guppy basecalling complete"
