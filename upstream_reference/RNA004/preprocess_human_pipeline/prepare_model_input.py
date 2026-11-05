#!/usr/bin/env python3
"""Select truth-centered step6 rows and prepare audited RNA004 model input."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from collections import Counter
from pathlib import Path

import pandas as pd


METRICS = ("mean", "std", "median", "length")
NUMERIC_COLUMNS = [f"{metric}_output{position}" for position in range(1, 10) for metric in METRICS]
STEP6_REQUIRED = {
    "read_id", "ref_name", "ref_start", "ref_end", "kmer_length", "strand",
    "sequence_input", *NUMERIC_COLUMNS,
}
TRUTH_REQUIRED = {
    "site_id", "chr", "position", "position0", "strand", "true_label",
    "reference_sequence",
}
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def reverse_signal_positions(frame: pd.DataFrame, mask: pd.Series) -> None:
    if not mask.any():
        return
    original = frame.loc[mask, NUMERIC_COLUMNS].copy()
    for destination in range(1, 10):
        source = 10 - destination
        for metric in METRICS:
            frame.loc[mask, f"{metric}_output{destination}"] = original[
                f"{metric}_output{source}"
            ].to_numpy()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--step6", required=True, type=Path)
    parser.add_argument("--truth", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--binary-output", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--chunksize", type=int, default=500_000)
    args = parser.parse_args()
    if args.chunksize < 1:
        raise ValueError("chunksize must be positive")
    for output in (args.output, args.binary_output):
        if output.exists() or output.with_suffix(output.suffix + ".tmp").exists():
            raise FileExistsError(f"Refusing to overwrite model input: {output}")

    truth = pd.read_csv(args.truth, dtype={"chr": str, "strand": str})
    missing = TRUTH_REQUIRED - set(truth.columns)
    if missing:
        raise ValueError(f"Ground truth lacks columns: {sorted(missing)}")
    truth["position0"] = pd.to_numeric(truth["position0"], errors="raise").astype("int64")
    if truth.duplicated(["chr", "position0", "strand"]).any() or truth["site_id"].duplicated().any():
        raise ValueError("Ground truth must contain one non-conflicting row per genomic site")
    truth["expected9"] = truth["reference_sequence"].str.upper().str.replace("U", "T").str[6:15]
    truth["context17"] = truth["reference_sequence"].str.upper().str.replace("U", "T").str[2:19]
    if truth["expected9"].str.len().ne(9).any() or truth["context17"].str.len().ne(17).any():
        raise ValueError("Ground-truth reference sequences must be centered 21-mers")
    truth_payload = truth.rename(columns={"chr": "ref_name", "strand": "truth_strand"})

    header = pd.read_csv(args.step6, nrows=0).columns.tolist()
    missing = STEP6_REQUIRED - set(header)
    if missing:
        raise ValueError(f"step6 lacks columns: {sorted(missing)}")

    output_columns = [
        "record_key", "sample_id", "site_id", "true_label", "read_id",
        "reference_name", "center_pos0", "strand", "center_ref_base",
        "9mer_final", "reference_sequence_21bp", "sequence_input",
        "signal_orientation", *NUMERIC_COLUMNS, "supportNum", "cellList",
        "GeneID", "transcriptID", "motifscore", "source_file", "source_line",
    ]
    temp = args.output.with_suffix(args.output.suffix + ".tmp")
    binary_temp = args.binary_output.with_suffix(args.binary_output.suffix + ".tmp")
    counters = Counter()
    seen_keys: set[str] = set()
    wrote_header = False

    truth_columns = [
        "ref_name", "position0", "truth_strand", "site_id", "true_label",
        "reference_sequence", "expected9", "context17", "supportNum", "cellList",
        "GeneID", "transcriptID", "motifscore", "source_file", "source_line",
    ]
    for chunk in pd.read_csv(args.step6, chunksize=args.chunksize, low_memory=False):
        counters["step6_rows"] += len(chunk)
        chunk["ref_start"] = pd.to_numeric(chunk["ref_start"], errors="raise").astype("int64")
        chunk["center_pos0"] = chunk["ref_start"] + 4
        chunk["strand"] = chunk["strand"].astype(str)
        selected = chunk.merge(
            truth_payload[truth_columns],
            how="inner",
            left_on=["ref_name", "center_pos0", "strand"],
            right_on=["ref_name", "position0", "truth_strand"],
            validate="many_to_one",
            sort=False,
        )
        counters["truth_center_rows"] += len(selected)
        if selected.empty:
            continue

        selected["sequence_input"] = selected["sequence_input"].astype(str).str.upper().str.replace("U", "T")
        direct = selected["sequence_input"].eq(selected["expected9"])
        reverse = selected["sequence_input"].map(reverse_complement).eq(selected["expected9"])
        ambiguous = direct & reverse
        reverse_only = reverse & ~direct
        mismatch = ~(direct | reverse)
        counters["context_direct"] += int(direct.sum())
        counters["context_reverse_complemented"] += int(reverse_only.sum())
        counters["context_palindromic"] += int(ambiguous.sum())
        counters["context_mismatch_excluded"] += int(mismatch.sum())
        selected = selected.loc[~mismatch].copy()
        reverse_only = reverse_only.loc[selected.index]
        reverse_signal_positions(selected, reverse_only)
        selected.loc[reverse_only, "sequence_input"] = selected.loc[reverse_only, "sequence_input"].map(
            reverse_complement
        )
        selected["signal_orientation"] = "catalog"
        selected.loc[reverse_only, "signal_orientation"] = "reverse_complemented_to_catalog"
        if not selected["sequence_input"].eq(selected["expected9"]).all():
            raise ValueError("Internal context-orientation failure")

        selected["record_key"] = "HeLa|" + selected["read_id"].astype(str) + "|" + selected["site_id"]
        duplicate_mask = selected["record_key"].isin(seen_keys) | selected["record_key"].duplicated()
        counters["duplicate_read_site_rows_excluded"] += int(duplicate_mask.sum())
        selected = selected.loc[~duplicate_mask].copy()
        if selected.empty:
            continue
        seen_keys.update(selected["record_key"])
        selected["sample_id"] = "HeLa"
        selected["reference_name"] = selected["ref_name"]
        selected["center_ref_base"] = selected["context17"].str[8].str.replace("T", "U")
        selected["9mer_final"] = selected["context17"].map(lambda value: "_".join(value.replace("T", "U")))
        selected["reference_sequence_21bp"] = selected["reference_sequence"]
        selected["strand"] = selected["truth_strand"]
        selected[output_columns].to_csv(
            temp, mode="a", header=not wrote_header, index=False, quoting=csv.QUOTE_MINIMAL
        )
        binary = selected[selected["true_label"].isin({"m6A", "m5C", "pseudoU"})]
        if not binary.empty:
            binary[output_columns].to_csv(
                binary_temp, mode="a", header=not binary_temp.exists(), index=False,
                quoting=csv.QUOTE_MINIMAL,
            )
            counters["binary_prepared_rows"] += len(binary)
        wrote_header = True
        counters["prepared_rows"] += len(selected)

    if not wrote_header or counters["prepared_rows"] == 0:
        raise RuntimeError("No truth-centered, context-valid model records were produced")
    if not binary_temp.exists() or counters["binary_prepared_rows"] == 0:
        raise RuntimeError("No binary-study model records were produced")
    temp.replace(args.output)
    binary_temp.replace(args.binary_output)
    for output in (args.output, args.binary_output):
        with output.open("rb") as source, gzip.open(str(output) + ".gz.tmp", "wb", compresslevel=6) as target:
            for block in iter(lambda: source.read(1024 * 1024), b""):
                target.write(block)
        Path(str(output) + ".gz.tmp").replace(Path(str(output) + ".gz"))

    summary = {
        "sample": "HeLa",
        "chemistry": "RNA004",
        "assembly": "hg38",
        "ground_truth_scope": "human reserved positive-site catalogs",
        "coordinate_rule": "center_pos0 = step6 ref_start + 4",
        "context_rule": "step6 centered 9-mer must match the catalog's strand-oriented 21-mer",
        "input_step6": str(args.step6),
        "input_step6_sha256": sha256(args.step6),
        "ground_truth": str(args.truth),
        "ground_truth_sha256": sha256(args.truth),
        "output": str(args.output),
        "output_sha256": sha256(args.output),
        "binary_output": str(args.binary_output),
        "binary_output_sha256": sha256(args.binary_output),
        "counts": dict(counters),
        "truth_labels_in_prepared": {},
    }
    label_counts = pd.read_csv(args.output, usecols=["true_label"])["true_label"].value_counts().sort_index()
    summary["truth_labels_in_prepared"] = label_counts.to_dict()
    args.summary.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
