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
fast5_dir="$run_dir/02_single_fast5"
out_dir="$run_dir/04_step4"
fast5_list="$out_dir/$sample.fast5.list"
success="$out_dir/_SUCCESS"

[[ -f "$run_dir/03_tombo_resquiggle.SUCCESS" ]] || {
    echo "[$sample] missing Tombo success marker" >&2
    exit 1
}
if [[ -f "$success" ]]; then
    echo "[$sample] step4 extraction already complete"
    exit 0
fi
if [[ -d "$out_dir" ]] && [[ -n "$(find "$out_dir" -mindepth 1 -print -quit)" ]]; then
    echo "[$sample] refusing to mix with a partial step4 directory: $out_dir" >&2
    exit 1
fi

mkdir -p "$out_dir"
find "$fast5_dir" -type f -name '*.fast5' -print | LC_ALL=C sort > "$fast5_list"

export PYTHONDONTWRITEBYTECODE=1
"$FEATURE_PYTHON" "$PIPELINE_DIR/python/extract_tombo_step4.py" \
    --fast5-list "$fast5_list" \
    --output-prefix "$out_dir/$sample.step4" \
    --workers "${SLURM_CPUS_PER_TASK:-1}" \
    --clip 10 \
    --progress-every 10000 \
    --rna

touch "$success"
echo "[$sample] step4 features and coordinates complete"
