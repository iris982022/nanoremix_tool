#!/usr/bin/env bash

set -euo pipefail

PIPELINE_DIR="${NANOREMIX_PIPELINE_DIR:-${SLURM_SUBMIT_DIR:-}}"
[[ -f "$PIPELINE_DIR/config.sh" ]] || {
    echo "cannot locate pipeline config.sh; PIPELINE_DIR=$PIPELINE_DIR" >&2
    exit 1
}
source "$PIPELINE_DIR/config.sh"

sample="$(sample_from_index "${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}")"
run_dir="$(sample_run_dir "$sample")"
input_dir="$run_dir/01_guppy_rna002/workspace"
out_dir="$run_dir/02_single_fast5"
success="$out_dir/_SUCCESS"

[[ -f "$run_dir/01_guppy_rna002/_SUCCESS" ]] || {
    echo "[$sample] missing basecalling success marker" >&2
    exit 1
}
if [[ -f "$success" ]]; then
    echo "[$sample] single-FAST5 conversion already complete"
    exit 0
fi
if [[ -d "$out_dir" ]] && [[ -n "$(find "$out_dir" -mindepth 1 -print -quit)" ]]; then
    echo "[$sample] refusing to mix with a partial single-FAST5 directory: $out_dir" >&2
    exit 1
fi

mkdir -p "$out_dir"
"$MULTI_TO_SINGLE" \
    -i "$input_dir" \
    -s "$out_dir" \
    -t "${SLURM_CPUS_PER_TASK:-1}" \
    --recursive

if [[ -z "$(find "$out_dir" -type f -name '*.fast5' -print -quit)" ]]; then
    echo "[$sample] no single-read FAST5 was produced" >&2
    exit 1
fi
touch "$success"
echo "[$sample] single-FAST5 conversion complete"
