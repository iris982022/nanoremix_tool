#!/usr/bin/env bash

set -euo pipefail

PIPELINE_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
VALIDATION_ROOT="$(cd -- "$PIPELINE_DIR/.." && pwd)"
PROJECT_ROOT="$(cd -- "$VALIDATION_ROOT/.." && pwd)"
SAMPLE="HeLa"
SAMPLE_ROOT="$VALIDATION_ROOT/$SAMPLE"
ARCHIVE="/gpfs/work/bio/zhixingwu16/nanopore_9mer_bio/HeLa_directRNA.tar.gz"
N_SPLITS=6

DORADO_ROOT="/gpfs/work/bio/jiayili1905/Benchmark/Tools/dorado-1.3.0-linux-x64"
DORADO="$DORADO_ROOT/bin/dorado"
DORADO_MODEL="$DORADO_ROOT/models/rna004_130bps_hac@v5.2.0"
REFERENCE_FASTA="/gpfs/work/bio/yuxinzhang17/DirectMultiRM/hg38.fa"
KMER_LEVELS="/gpfs/work/bio/yuxinzhang17/DirectMultiRM/ONT/9mer_levels_v1.txt"
FEATURE_SCRIPT="$VALIDATION_ROOT/preprocess_pipeline/feature_extraction_rna004.py"
REMORA_PYTHON="/gpfs/work/bio/zhixingwu16/.conda/envs/remora-env/bin/python"
UTIL_PYTHON="$REMORA_PYTHON"
RESERVED="$VALIDATION_ROOT/../reserved"

MODEL_ENV="/gpfs/work/bio/zhixingwu16/.conda/envs/nanoremix_env"
MODEL_PYTHON="$MODEL_ENV/bin/python"
BINARY_ROOT="$PROJECT_ROOT/Binary/nanoremix_binary_fused_RNA004"
MULTI_ROOT="$PROJECT_ROOT/Multi/nanoremix_multi_fused_RNA004"
SECOND_ROOT="$MULTI_ROOT/models/second-gen"
CONTEXT_SCRIPT="$VALIDATION_ROOT/mouse_model_pipeline/add_splice_aware_context.py"
FUSED_ROOT="$SAMPLE_ROOT/fused_model_results"

REGIONS="$PIPELINE_DIR/human_hg38_candidate_regions.csv"
ANNOTATIONS="$PIPELINE_DIR/human_hg38_candidate_annotations.csv"
GROUND_TRUTH="$PIPELINE_DIR/human_hg38_ground_truth_sites.csv"
GROUND_TRUTH_CONFLICTS="$PIPELINE_DIR/human_hg38_ground_truth_conflicts.csv"
GROUND_TRUTH_SUMMARY="$PIPELINE_DIR/human_hg38_ground_truth_summary.json"

export PIPELINE_DIR VALIDATION_ROOT PROJECT_ROOT SAMPLE SAMPLE_ROOT ARCHIVE N_SPLITS
export DORADO_ROOT DORADO DORADO_MODEL REFERENCE_FASTA KMER_LEVELS FEATURE_SCRIPT
export REMORA_PYTHON UTIL_PYTHON RESERVED REGIONS ANNOTATIONS GROUND_TRUTH
export GROUND_TRUTH_CONFLICTS GROUND_TRUTH_SUMMARY
export MODEL_ENV MODEL_PYTHON BINARY_ROOT MULTI_ROOT SECOND_ROOT CONTEXT_SCRIPT FUSED_ROOT
