#!/usr/bin/env python3
"""Extract strand-oriented RNA004 9-mer signal features at candidate sites."""

from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pod5
import tqdm
from remora import io, refine_signal_map, util


COORD_COLUMNS = [
    "read_id", "ref_name", "ref_start", "ref_end", "kmer_length", "strand"
]


def extract_base_call_error(io_read):
    """Return reference-anchored basecall-error arrays in read orientation."""
    ref_len = len(io_read.ref_seq)
    quality_text = io_read.full_align.get("qual")
    if quality_text in (None, "*"):
        seq_quality = np.zeros(io_read.seq_len, dtype=np.float32)
    else:
        seq_quality = np.fromiter(
            (ord(char) - 33 for char in quality_text), dtype=np.float32
        )
        if io_read.ref_reg.strand == "-":
            seq_quality = seq_quality[::-1]
        quality_sd = float(seq_quality.std())
        if quality_sd > 0:
            seq_quality = (seq_quality - float(seq_quality.mean())) / quality_sd
        else:
            seq_quality = np.zeros_like(seq_quality)

    insertion = np.zeros(ref_len, dtype=np.float32)
    deletion = np.zeros(ref_len, dtype=np.float32)
    mismatch = np.zeros(ref_len, dtype=np.float32)
    quality = np.zeros(ref_len, dtype=np.float32)
    ref_index = 0
    query_index = 0

    for operation, length in io_read.cigar:
        if operation in (0, 7, 8):  # M, =, X consume query and reference.
            ref_slice = slice(ref_index, ref_index + length)
            query_slice = slice(query_index, query_index + length)
            quality[ref_slice] = seq_quality[query_slice]
            ref_bases = io_read.ref_seq[ref_slice]
            query_bases = io_read.seq[query_slice]
            if len(ref_bases) == len(query_bases):
                mismatch[ref_slice] = np.fromiter(
                    (left != right for left, right in zip(ref_bases, query_bases)),
                    dtype=np.float32,
                    count=length,
                )
            ref_index += length
            query_index += length
        elif operation == 1:  # I consumes query only; anchor it to one ref base.
            if ref_len:
                insertion[min(ref_index, ref_len - 1)] = 1
            query_index += length
        elif operation in (2, 3):  # D and N consume reference only.
            deletion[ref_index:ref_index + length] = 1
            ref_index += length
        elif operation == 4:  # S consumes query only.
            query_index += length
        elif operation in (5, 6):  # H and P consume neither stored sequence.
            continue
        else:
            raise ValueError(f"Unsupported CIGAR operation: {operation}")

    if ref_index != ref_len:
        raise ValueError(
            f"CIGAR/reference length mismatch for {io_read.read_id}: "
            f"cigar={ref_index} reference={ref_len}"
        )
    return np.column_stack((insertion, deletion, mismatch, quality))


def quantiles(signal, sequence_to_signal):
    segments = [
        signal[sequence_to_signal[index]:sequence_to_signal[index + 1]]
        for index in range(len(sequence_to_signal) - 1)
    ]
    segments = [segment if segment.size else np.nan for segment in segments]
    return {
        "min": np.asarray([np.min(segment) for segment in segments], dtype=np.float32),
        "q1": np.asarray([np.quantile(segment, 0.25) for segment in segments], dtype=np.float32),
        "median": np.asarray([np.median(segment) for segment in segments], dtype=np.float32),
        "q3": np.asarray([np.quantile(segment, 0.75) for segment in segments], dtype=np.float32),
        "max": np.asarray([np.max(segment) for segment in segments], dtype=np.float32),
    }


def extract_signal_stats(io_read):
    region = io.RefRegion(
        ctg=io_read.ref_reg.ctg,
        strand=io_read.ref_reg.strand,
        start=io_read.ref_reg.start,
        end=io_read.ref_reg.end,
    )
    stats = io_read.compute_per_base_metric(
        "dwell_mean_sd", ref_anchored=True, region=region
    )
    stats.update(
        io_read.compute_per_base_metric(
            metric_func=quantiles, ref_anchored=True, region=region
        )
    )
    return np.column_stack(
        (
            stats["mean"], stats["sd"], stats["dwell"], stats["min"],
            stats["q1"], stats["median"], stats["q3"], stats["max"],
        )
    )


def one_hot_encoding(sequence):
    encoded = np.fromiter(("ACGT".index(base) for base in sequence), dtype=np.int8)
    one_hot = np.zeros((encoded.size, 4), dtype=np.float32)
    one_hot[np.arange(encoded.size), encoded] = 1
    return one_hot


def extract_kmer_features(io_read):
    key = (io_read.ref_reg.ctg, io_read.ref_reg.strand)
    starts = candidate_starts.get(key)
    if starts is None:
        return [], [], [], [], []

    kmer_size = FLAGS.kmer
    first = np.searchsorted(starts, io_read.ref_reg.start, side="left")
    last = np.searchsorted(
        starts, io_read.ref_reg.end - kmer_size, side="right"
    )
    starts = starts[first:last:FLAGS.step]
    if starts.size == 0:
        return [], [], [], [], []

    basecall_errors = extract_base_call_error(io_read)
    signal_stats = extract_signal_stats(io_read)
    sequence = one_hot_encoding(io_read.ref_seq)
    model_levels = sig_map_refiner.extract_levels(util.seq_to_int(io_read.ref_seq))

    expected_length = io_read.ref_reg.end - io_read.ref_reg.start
    arrays = (sequence, signal_stats, model_levels, basecall_errors)
    if any(len(array) != expected_length for array in arrays):
        lengths = [len(array) for array in arrays]
        raise ValueError(
            f"Feature/reference length mismatch for {io_read.read_id}: "
            f"expected={expected_length} observed={lengths}"
        )

    # Remora stores a negative-strand read in target-strand order (high to low
    # genomic coordinate). Convert all arrays to ascending genomic order for
    # coordinate slicing, then reverse each selected window back to the target
    # strand. Sequence and every signal/BSE row therefore stay paired.
    is_reverse = io_read.ref_reg.strand == "-"
    if is_reverse:
        sequence, signal_stats, model_levels, basecall_errors = (
            array[::-1] for array in arrays
        )

    sequences, stats, levels, errors, coordinates = [], [], [], [], []
    for start in starts:
        offset = int(start) - io_read.ref_reg.start
        stop = offset + kmer_size
        kmer_arrays = [
            sequence[offset:stop],
            signal_stats[offset:stop],
            model_levels[offset:stop],
            basecall_errors[offset:stop],
        ]
        if any(len(array) != kmer_size for array in kmer_arrays):
            raise ValueError(f"Incomplete k-mer window for {io_read.read_id} at {start}")
        if is_reverse:
            kmer_arrays = [array[::-1] for array in kmer_arrays]
        if np.isnan(kmer_arrays[1]).any():
            continue
        sequences.append(kmer_arrays[0])
        stats.append(kmer_arrays[1])
        levels.append(kmer_arrays[2])
        errors.append(kmer_arrays[3])
        coordinates.append(
            [
                io_read.read_id,
                io_read.ref_reg.ctg,
                int(start),
                int(start) + kmer_size,
                kmer_size,
                io_read.ref_reg.strand,
            ]
        )
    return sequences, stats, levels, errors, coordinates


def iterate_pod5(pod5_reader, bam_reader, read_ids):
    aggregate = [[], [], [], [], []]
    counts = Counter()
    error_examples = []
    for read_id in read_ids:
        try:
            pod5_read = pod5_reader.get_read(read_id)
            bam_read = bam_reader.get_first_alignment(read_id)
            io_read = io.Read.from_pod5_and_alignment(pod5_read, bam_read)
            if io_read.seq_len <= 500:
                counts["short_read"] += 1
                continue
            if FLAGS.resquiggle:
                io_read.set_refine_signal_mapping(sig_map_refiner, ref_mapping=True)
            values = extract_kmer_features(io_read)
            for destination, source in zip(aggregate, values):
                destination.extend(source)
            counts["processed_read"] += 1
            counts["feature_rows"] += len(values[-1])
        except Exception as error:  # A bad raw read must not discard its split.
            counts["read_error"] += 1
            if len(error_examples) < 20:
                error_examples.append(f"{read_id}: {type(error).__name__}: {error}")
    for message in error_examples:
        print(f"READ_ERROR {message}", file=sys.stderr)
    print("read_counts " + " ".join(f"{key}={counts[key]}" for key in sorted(counts)))
    return (*aggregate, counts)


def iterate_split(split):
    bam_path = Path(FLAGS.bam) / f"{split}_sorted.bam"
    if not bam_path.is_file() or not Path(str(bam_path) + ".bai").is_file():
        raise FileNotFoundError(f"Missing indexed BAM: {bam_path}")
    bam_reader = io.ReadIndexedBam(str(bam_path))

    bam_read_ids = set()
    for region in regions.itertuples(index=False):
        for alignment in bam_reader.fetch(
            ctg=region.seqnames,
            strand=region.strand,
            start=int(region.start),
            end=int(region.end),
        ):
            bam_read_ids.add(alignment.query_name)
    print(f"split={split} candidate_bam_reads={len(bam_read_ids)}")

    pod5_dir = Path(FLAGS.pod5_dir) / split
    pod5_files = sorted(pod5_dir.rglob("*.pod5"))
    if not pod5_files:
        raise FileNotFoundError(f"No POD5 files found in {pod5_dir}")

    aggregate = [[], [], [], [], []]
    total_counts = Counter()
    for pod5_path in tqdm.tqdm(pod5_files, desc=split):
        with pod5.DatasetReader(pod5_path) as pod5_reader:
            pod5_read_ids = set(pod5_reader.read_ids)
            read_ids = sorted(pod5_read_ids.intersection(bam_read_ids), key=str)
            values = iterate_pod5(pod5_reader, bam_reader, read_ids)
        for destination, source in zip(aggregate, values[:5]):
            destination.extend(source)
        total_counts.update(values[5])

    sequences, stats, levels, errors, coordinates = aggregate
    if not coordinates:
        raise RuntimeError(f"Feature extraction produced zero rows for {split}")
    if not all(len(values) == len(coordinates) for values in aggregate):
        raise RuntimeError(f"Feature array row-count mismatch for {split}")

    output_dir = Path(FLAGS.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    npz_path = output_dir / f"{split}.npz"
    csv_path = output_dir / f"{split}.csv"
    npz_temp = output_dir / f".{split}.npz.{FLAGS.run_id}.tmp"
    csv_temp = output_dir / f".{split}.csv.{FLAGS.run_id}.tmp"
    for temp in (npz_temp, csv_temp):
        if temp.exists():
            raise FileExistsError(f"Refusing to overwrite temporary output: {temp}")
    with npz_temp.open("wb") as handle:
        np.savez(
            handle,
            seq=np.asarray(sequences),
            stat=np.asarray(stats),
            level=np.asarray(levels),
            bse=np.asarray(errors),
        )
    pd.DataFrame(coordinates, columns=COORD_COLUMNS).to_csv(csv_temp, index=False)
    npz_temp.replace(npz_path)
    csv_temp.replace(csv_path)
    print(
        f"split={split} feature_rows={len(coordinates)} "
        f"read_errors={total_counts['read_error']} output={npz_path}"
    )


def parse_args(argv):
    parser = argparse.ArgumentParser("extract RNA004 event-signal features")
    parser.add_argument("--pod5_dir", required=True)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--reg", required=True)
    parser.add_argument("--level", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("-t", "--threads", type=int, default=1)
    parser.add_argument("--splits", nargs="+", required=True)
    parser.add_argument("--kmer", type=int, default=9)
    parser.add_argument("--step", type=int, default=1)
    parser.add_argument(
        "--resquiggle",
        type=lambda value: str(value).lower() == "true",
        default=True,
    )
    parser.add_argument("--run-id", default="manual")
    parsed = parser.parse_args(argv)
    if parsed.kmer < 1 or parsed.step < 1 or parsed.threads < 1:
        parser.error("--kmer, --step, and --threads must be positive")
    return parsed


if __name__ == "__main__":
    FLAGS = parse_args(sys.argv[1:])
    sig_map_refiner = refine_signal_map.SigMapRefiner(
        kmer_model_filename=FLAGS.level,
        do_rough_rescale=True,
        scale_iters=0,
        do_fix_guage=True,
    )
    regions = pd.read_csv(FLAGS.reg)
    required_columns = {"seqnames", "start", "end", "strand"}
    missing = required_columns - set(regions.columns)
    if missing:
        raise ValueError(f"Candidate-region file lacks columns: {sorted(missing)}")
    regions["start"] = pd.to_numeric(regions["start"], errors="raise").astype("int64")
    regions["end"] = pd.to_numeric(regions["end"], errors="raise").astype("int64")
    if (regions["start"] < 0).any() or (regions["end"] <= regions["start"]).any():
        raise ValueError("Candidate regions contain invalid coordinates")
    if not regions["strand"].isin(["+", "-"]).all():
        raise ValueError("Candidate regions contain invalid strands")
    candidate_starts = {
        key: np.unique(
            np.concatenate(
                [np.arange(row.start, row.end, dtype=np.int64) for row in frame.itertuples()]
            )
        )
        for key, frame in regions.groupby(["seqnames", "strand"], sort=False)
    }

    if len(FLAGS.splits) == 2:
        split_numbers = range(int(FLAGS.splits[0]), int(FLAGS.splits[1]) + 1)
    else:
        split_numbers = [int(value) for value in FLAGS.splits]
    for split_number in split_numbers:
        iterate_split(f"split{split_number}")
