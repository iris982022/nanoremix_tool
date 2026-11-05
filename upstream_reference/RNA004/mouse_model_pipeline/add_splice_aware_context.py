#!/usr/bin/env python3
"""Attach read-specific, splice-aware 81-nt reference contexts to mouse rows.

The aligned BAM is genomic.  For each read we follow its CIGAR path, retain
reference-consuming exonic operations (M/D/=/X), skip introns (N), orient the
path in RNA strand direction, and take 40 reference bases on either side of
the truth site.  Every retained sequence is checked against the catalog 21-mer
already validated during input preparation.  Rows lacking complete flanks are
reported and excluded rather than padded or replaced by genomic context.
"""

from __future__ import annotations

import argparse
import bisect
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd
import pysam


COMPLEMENT = str.maketrans("ACGT", "TGCA")
BINARY_TRUTH = {"A", "C", "U", "m6A", "m5C", "pseudoU"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def reference_path(alignment: pysam.AlignedSegment) -> list[int]:
    """Return exonic reference coordinates in ascending genomic order."""
    position = alignment.reference_start
    coordinates: list[int] = []
    for operation, length in alignment.cigartuples or []:
        if operation in (0, 2, 7, 8):  # M, D, =, X: transcript reference bases
            coordinates.extend(range(position, position + length))
            position += length
        elif operation == 3:  # N: spliced intron, absent from mature RNA
            position += length
        elif operation in (1, 4, 5, 6):  # I/S/H/P: no reference base
            continue
        else:
            raise ValueError(f"Unsupported CIGAR operation {operation}")
    if position != alignment.reference_end:
        raise ValueError("CIGAR/reference-end disagreement")
    return coordinates


def fetch_oriented(
    fasta: pysam.FastaFile,
    reference_name: str,
    coordinates: list[int],
    reverse: bool,
) -> str:
    """Fetch possibly discontinuous reference coordinates in RNA orientation."""
    # Coordinates are at most 81 bases.  Per-base fetches are cached by the
    # caller and keep the logic explicit at splice junctions.
    sequence = "".join(
        fasta.fetch(reference_name, position, position + 1).upper()
        for position in coordinates
    )
    return sequence.translate(COMPLEMENT) if reverse else sequence


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument(
        "--bam", required=True, type=Path, action="append",
        help="Aligned BAM; repeat this option when reads are split across BAMs",
    )
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--reference-index", required=True, type=Path)
    parser.add_argument("--multi-output", required=True, type=Path)
    parser.add_argument("--binary-output", required=True, type=Path)
    parser.add_argument("--audit", required=True, type=Path)
    args = parser.parse_args()

    for output in (args.multi_output, args.binary_output, args.audit):
        if output.exists():
            raise FileExistsError(output)
        output.parent.mkdir(parents=True, exist_ok=True)

    frame = pd.read_csv(args.input, low_memory=False)
    required = {
        "record_key", "read_id", "reference_name", "center_pos0", "strand",
        "reference_sequence_21bp", "true_label",
    }
    missing = required - set(frame)
    if missing:
        raise ValueError(f"Input lacks columns: {sorted(missing)}")
    if frame["record_key"].duplicated().any():
        raise ValueError("Duplicate record_key values")
    frame["center_pos0"] = pd.to_numeric(frame["center_pos0"], errors="raise").astype(int)

    targets: dict[str, list[int]] = defaultdict(list)
    for index, read_id in enumerate(frame["read_id"].astype(str)):
        targets[read_id].append(index)

    contexts: list[str | None] = [None] * len(frame)
    reasons: list[str | None] = [None] * len(frame)
    seen_alignments: set[str] = set()
    context_cache: dict[tuple[str, bool, tuple[int, ...]], str] = {}
    counters = Counter()

    fasta = pysam.FastaFile(str(args.reference), filepath_index=str(args.reference_index))
    for bam_path in args.bam:
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            for alignment in bam.fetch(until_eof=True):
                read_id = alignment.query_name
                indices = targets.get(read_id)
                if indices is None:
                    continue
                if read_id in seen_alignments:
                    raise ValueError(f"Multiple retained alignments for read {read_id}")
                seen_alignments.add(read_id)
                coordinates = reference_path(alignment)
                for index in indices:
                    row = frame.iloc[index]
                    if alignment.reference_name != row.reference_name:
                        reasons[index] = "reference_name_mismatch"
                        continue
                    expected_reverse = row.strand == "-"
                    if row.strand not in ("+", "-") or alignment.is_reverse != expected_reverse:
                        reasons[index] = "strand_mismatch"
                        continue
                    center = int(row.center_pos0)
                    offset = bisect.bisect_left(coordinates, center)
                    if offset >= len(coordinates) or coordinates[offset] != center:
                        reasons[index] = "center_not_on_reference_path"
                        continue
                    if offset < 40 or offset + 40 >= len(coordinates):
                        reasons[index] = "insufficient_exonic_flanks"
                        continue
                    window = coordinates[offset - 40:offset + 41]
                    if expected_reverse:
                        window = window[::-1]
                    cache_key = (alignment.reference_name, expected_reverse, tuple(window))
                    sequence = context_cache.get(cache_key)
                    if sequence is None:
                        sequence = fetch_oriented(
                            fasta, alignment.reference_name, window, expected_reverse
                        )
                        context_cache[cache_key] = sequence
                    if len(sequence) != 81 or set(sequence) - set("ACGT"):
                        reasons[index] = "invalid_reference_sequence"
                        continue
                    expected21 = str(row.reference_sequence_21bp).upper().replace("U", "T")
                    if sequence[30:51] != expected21:
                        reasons[index] = "catalog_21mer_mismatch"
                        continue
                    contexts[index] = sequence
                    reasons[index] = "ok"

    fasta.close()
    for index, reason in enumerate(reasons):
        if reason is None:
            reasons[index] = "alignment_not_found"
    counters.update(reasons)

    valid = pd.Series(reasons).eq("ok").to_numpy()
    retained = frame.loc[valid].copy()
    retained["reference_sequence_81bp"] = [
        contexts[index] for index in range(len(frame)) if valid[index]
    ]
    if retained.empty:
        raise ValueError("No rows have valid splice-aware 81-nt context")
    if not retained["reference_sequence_81bp"].str[30:51].eq(
        retained["reference_sequence_21bp"].str.upper().str.replace("U", "T")
    ).all():
        raise AssertionError("Final 21-mer validation failed")

    retained.to_csv(args.multi_output, index=False, compression="gzip")
    binary = retained[retained["true_label"].isin(BINARY_TRUTH)].copy()
    binary.to_csv(args.binary_output, index=False, compression="gzip")
    audit = {
        "method": (
            "Read-specific mature-RNA reference path from genomic BAM CIGAR; "
            "M/D/=/X retained, N skipped, negative strand reverse-complemented"
        ),
        "input": str(args.input),
        "input_sha256": sha256(args.input),
        "bams": [str(path) for path in args.bam],
        "bam_sha256": {str(path): sha256(path) for path in args.bam},
        "reference": str(args.reference),
        "reference_index": str(args.reference_index),
        "input_rows": len(frame),
        "input_unique_reads": int(frame["read_id"].nunique()),
        "seen_alignments": len(seen_alignments),
        "status_counts": dict(sorted(counters.items())),
        "retained_multi_rows": len(retained),
        "retained_multi_sites": int(retained["site_id"].nunique()),
        "retained_binary_rows": len(binary),
        "retained_binary_sites": int(binary["site_id"].nunique()),
        "unique_81nt_contexts": int(retained["reference_sequence_81bp"].nunique()),
        "validation": "seq[30:51] equals catalog strand-oriented 21-mer for every row",
        "excluded_rows": int((~valid).sum()),
    }
    args.audit.write_text(json.dumps(audit, indent=2) + "\n")
    print(json.dumps(audit, indent=2), flush=True)


if __name__ == "__main__":
    main()
