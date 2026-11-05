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
success="$run_dir/03_tombo_resquiggle.SUCCESS"
reference_dir="$run_dir/03_tombo_reference"
local_reference="$reference_dir/transcriptome.fa"
failed_reads="$run_dir/03_tombo_failed_reads.${SLURM_JOB_ID:-manual}_${SLURM_ARRAY_TASK_ID}.tsv"

[[ -f "$fast5_dir/_SUCCESS" ]] || {
    echo "[$sample] missing single-FAST5 success marker" >&2
    exit 1
}
if [[ -f "$success" ]]; then
    echo "[$sample] Tombo resquiggle already complete"
    exit 0
fi

# GraphMap may need to build an index next to the reference. The shared input
# reference is intentionally read-only, so give each sample a private writable
# copy and avoid index races between the two array tasks.
mkdir -p "$reference_dir"
if [[ ! -s "$local_reference" ]]; then
    cp -p "$REFERENCE_FASTA" "$local_reference"
fi

# Tombo writes RawGenomeCorrected_000 into the COPIES under RUN_ROOT. The
# source multi-FAST5 directories are never modified.
export PYTHONDONTWRITEBYTECODE=1
export HDF5_USE_FILE_LOCKING=FALSE
"$TOMBO" resquiggle \
    --basecall-group "$TOMBO_BASECALL_GROUP" \
    --basecall-subgroups BaseCalled_template \
    --failed-reads-filename "$failed_reads" \
    "$fast5_dir" \
    "$local_reference" \
    --graphmap-executable "$GRAPHMAP" \
    --processes "${SLURM_CPUS_PER_TASK:-1}" \
    --align-processes "${SLURM_CPUS_PER_TASK:-1}" \
    --overwrite \
    --include-event-stdev

"$FEATURE_PYTHON" "$PIPELINE_DIR/python/validate_tombo_fast5.py" \
    --fast5-dir "$fast5_dir" \
    --sample-size 2000 \
    --min-ready 10 \
    --min-ready-fraction 0.01

touch "$success"
echo "[$sample] Tombo resquiggle complete"
