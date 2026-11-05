#!/usr/bin/env python3
"""Extract NanoRemix 5-mer step4 features and Tombo coordinates.

This is a streaming, coordinate-aware replacement for the historical
``extract_raw_and_feature_fast.py``. It intentionally preserves the original
NanoRemix 5-mer feature definition:

* RNA raw signal is reversed.
* Signal is normalized by median/MAD over unique ADC values.
* Five consecutive Tombo events form one feature row.
* The first/last ``clip`` events plus two k-mer flanks are excluded.

Tombo Events are read-centric. Tombo's own helper reverses Events on a minus
alignment to obtain genome-centric order, so the reference coordinate for
event ``i`` is ``mapped_start + i`` on plus and
``mapped_start + n_events - 1 - i`` on minus.
"""

from __future__ import print_function

import argparse
import json
import multiprocessing as mp
import os
import sys
import time
import traceback

import h5py
import numpy as np


MAD_SCALE = 0.6744897501960817
VALID_BASES = set("ACGT")
WORKER_CONFIG = None


def decode_text(value):
    if isinstance(value, bytes):
        return value.decode("utf-8", "replace")
    if hasattr(value, "decode"):
        try:
            return value.decode("utf-8", "replace")
        except Exception:
            pass
    return str(value)


def clean_field(value):
    return decode_text(value).replace("\t", " ").replace("\n", " ").strip()


def fmt_float(value):
    return format(float(value), ".17g")


def robust_mad_unique(signal):
    unique = np.unique(np.asarray(signal, dtype=np.float64))
    if unique.size < 2:
        raise ValueError("raw signal has fewer than two unique values")
    center = np.median(unique)
    mad = np.median(np.abs(unique - center)) / MAD_SCALE
    if not np.isfinite(mad) or mad <= 0:
        raise ValueError("raw signal MAD is zero or non-finite")
    return center, mad


def raw_read_group(handle):
    reads = handle["/Raw/Reads"]
    names = list(reads.keys())
    if len(names) != 1:
        raise ValueError("expected one read in single-FAST5; observed %d" % len(names))
    return reads[names[0]]


def event_bases(events):
    output = []
    for value in events["base"]:
        base = decode_text(value).upper().replace("U", "T")
        if len(base) != 1:
            raise ValueError("invalid Tombo event base: %r" % base)
        output.append(base)
    return output


def init_worker(config):
    global WORKER_CONFIG
    WORKER_CONFIG = config


def extract_one_fast5(path):
    config = WORKER_CONFIG
    try:
        with h5py.File(path, "r") as handle:
            raw_group = raw_read_group(handle)
            raw_signal = raw_group["Signal"][()]
            read_id = clean_field(raw_group.attrs.get("read_id", os.path.basename(path)))
            read_id = read_id.replace("|", "_")

            corr_path = "/Analyses/{}/{}".format(
                config["corrected_group"], config["basecall_subgroup"]
            )
            corr = handle[corr_path]
            status = clean_field(corr.attrs.get("status", "success")).lower()
            if status not in ("success", "none", ""):
                raise ValueError("Tombo status is not success: %s" % status)

            events = corr["Events"][()]
            if len(events) < (2 * config["clip"] + 5):
                raise ValueError("too few Tombo events: %d" % len(events))
            required_fields = set(("start", "length", "base"))
            if not required_fields.issubset(set(events.dtype.names or ())):
                raise ValueError("Tombo Events lacks start/length/base fields")

            alignment = corr["Alignment"].attrs
            transcript_id = clean_field(alignment["mapped_chrom"])
            strand = clean_field(alignment["mapped_strand"])
            if strand not in ("+", "-"):
                raise ValueError("invalid mapped strand: %r" % strand)
            mapped_start = int(alignment["mapped_start"])
            mapped_end = int(alignment["mapped_end"])
            expected_end = mapped_start + len(events)
            if mapped_end != expected_end:
                raise ValueError(
                    "mapped span/event mismatch: start=%d end=%d events=%d"
                    % (mapped_start, mapped_end, len(events))
                )

            read_start = int(corr["Events"].attrs["read_start_rel_to_raw"])
            starts = np.asarray(events["start"], dtype=np.int64) + read_start
            lengths = np.asarray(events["length"], dtype=np.int64)
            bases = event_bases(events)

            # The historical extractor always reversed RNA signal. Tombo stores
            # the RNA flag too, but --rna is explicit to avoid metadata drift.
            if config["rna"]:
                raw_signal = raw_signal[::-1]
            center, mad = robust_mad_unique(raw_signal)
            normalized = (np.asarray(raw_signal, dtype=np.float64) - center) / mad

            feature_lines = []
            coordinate_lines = []
            first_index = config["clip"] + 2
            stop_index = len(events) - config["clip"] - 2

            for index in range(first_index, stop_index):
                kmer = "".join(bases[index - 2:index + 3])
                if len(kmer) != 5 or any(base not in VALID_BASES for base in kmer):
                    continue

                segments = []
                for event_index in range(index - 2, index + 3):
                    start = int(starts[event_index])
                    length = int(lengths[event_index])
                    end = start + length
                    if length <= 0 or start < 0 or end > len(normalized):
                        raise ValueError(
                            "event signal bounds invalid at event %d: %d:%d of %d"
                            % (event_index, start, end, len(normalized))
                        )
                    segment = normalized[start:end]
                    if segment.size == 0:
                        raise ValueError("empty event segment at event %d" % event_index)
                    segments.append(segment)

                means = [fmt_float(np.mean(segment)) for segment in segments]
                stds = [fmt_float(np.std(segment)) for segment in segments]
                medians = [fmt_float(np.median(segment)) for segment in segments]
                event_lengths = [str(int(lengths[i])) for i in range(index - 2, index + 3)]

                row_id = "{}|{}|{}".format(read_id, index, strand)
                feature_lines.append(
                    "\t".join(
                        (
                            row_id,
                            bases[index],
                            kmer,
                            "|".join(means),
                            "|".join(stds),
                            "|".join(medians),
                            "|".join(event_lengths),
                            kmer,
                        )
                    )
                    + "\n"
                )

                if strand == "+":
                    transcript_position = mapped_start + index
                else:
                    transcript_position = mapped_start + len(events) - 1 - index
                coordinate_lines.append(
                    "\t".join(
                        (
                            row_id,
                            read_id,
                            str(index),
                            strand,
                            transcript_id,
                            str(transcript_position),
                            str(mapped_start),
                            str(mapped_end),
                            bases[index],
                            kmer,
                            clean_field(path),
                        )
                    )
                    + "\n"
                )

        return {
            "ok": True,
            "path": path,
            "feature_text": "".join(feature_lines),
            "coordinate_text": "".join(coordinate_lines),
            "rows": len(feature_lines),
            "strand": strand,
        }
    except Exception as exc:
        return {
            "ok": False,
            "path": path,
            "error": "{}: {}".format(type(exc).__name__, clean_field(exc)),
            "traceback": traceback.format_exc(),
            "rows": 0,
        }


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fast5-list", required=True)
    parser.add_argument("--output-prefix", required=True)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--clip", type=int, default=10)
    parser.add_argument("--corrected-group", default="RawGenomeCorrected_000")
    parser.add_argument("--basecall-subgroup", default="BaseCalled_template")
    parser.add_argument("--rna", action="store_true", help="Reverse raw signal as required for DRS")
    parser.add_argument(
        "--progress-every",
        type=int,
        default=10000,
        help="Write a machine-readable progress message every N processed FAST5 files",
    )
    return parser.parse_args(argv)


def emit_progress(stats, total, started):
    processed = stats["successful_fast5"] + stats["failed_fast5"]
    elapsed = max(time.time() - started, 0.0)
    rate = processed / elapsed if elapsed > 0 else None
    eta = (total - processed) / rate if rate and processed < total else 0.0
    payload = {
        "processed_fast5": processed,
        "input_fast5": total,
        "progress_percent": 100.0 * processed / total if total else 0.0,
        "successful_fast5": stats["successful_fast5"],
        "failed_fast5": stats["failed_fast5"],
        "feature_rows": stats["feature_rows"],
        "elapsed_seconds": elapsed,
        "rate_fast5_per_hour": rate * 3600.0 if rate is not None else None,
        "eta_seconds": eta,
    }
    print("[EXTRACT PROGRESS] " + json.dumps(payload, sort_keys=True), file=sys.stderr)
    sys.stderr.flush()


def main(argv=None):
    args = parse_args(argv)
    if args.workers < 1:
        raise SystemExit("--workers must be >= 1")
    if args.clip < 0:
        raise SystemExit("--clip must be >= 0")
    if args.progress_every < 1:
        raise SystemExit("--progress-every must be >= 1")

    with open(args.fast5_list, "r") as handle:
        paths = [line.strip() for line in handle if line.strip()]
    if not paths:
        raise SystemExit("FAST5 list is empty: %s" % args.fast5_list)

    output_dir = os.path.dirname(os.path.abspath(args.output_prefix))
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    feature_path = args.output_prefix + ".feature.tsv"
    coordinate_path = args.output_prefix + ".coordinates.tsv"
    error_path = args.output_prefix + ".errors.tsv"
    summary_path = args.output_prefix + ".summary.json"
    for path in (feature_path, coordinate_path, error_path, summary_path):
        if os.path.exists(path):
            raise SystemExit("refusing to overwrite existing output: %s" % path)

    config = {
        "clip": args.clip,
        "corrected_group": args.corrected_group,
        "basecall_subgroup": args.basecall_subgroup,
        "rna": bool(args.rna),
    }
    stats = {
        "input_fast5": len(paths),
        "successful_fast5": 0,
        "failed_fast5": 0,
        "feature_rows": 0,
        "plus_fast5": 0,
        "minus_fast5": 0,
        "clip": args.clip,
        "rna": bool(args.rna),
        "corrected_group": args.corrected_group,
        "basecall_subgroup": args.basecall_subgroup,
    }

    feature_header = "id\tbase\tsequence_input\tmean\tstd\tmedian\tlength\tkmer_raw\n"
    coordinate_header = (
        "id\tread_id\tfast5_order\tmapped_strand\ttranscript_id\t"
        "transcript_position0\tmapped_start0\tmapped_end0\tevent_base\t"
        "event_kmer5\tfast5_path\n"
    )

    started = time.time()
    pool = mp.Pool(args.workers, initializer=init_worker, initargs=(config,))
    try:
        with open(feature_path, "w") as feature_out, open(
            coordinate_path, "w"
        ) as coordinate_out, open(error_path, "w") as error_out:
            feature_out.write(feature_header)
            coordinate_out.write(coordinate_header)
            error_out.write("fast5_path\terror\n")
            for result in pool.imap_unordered(extract_one_fast5, paths, chunksize=1):
                if result["ok"]:
                    feature_out.write(result["feature_text"])
                    coordinate_out.write(result["coordinate_text"])
                    stats["successful_fast5"] += 1
                    stats["feature_rows"] += int(result["rows"])
                    stats[result["strand"] == "+" and "plus_fast5" or "minus_fast5"] += 1
                else:
                    error_out.write(
                        "{}\t{}\n".format(clean_field(result["path"]), result["error"])
                    )
                    stats["failed_fast5"] += 1
                processed = stats["successful_fast5"] + stats["failed_fast5"]
                if processed % args.progress_every == 0:
                    feature_out.flush()
                    coordinate_out.flush()
                    error_out.flush()
                    emit_progress(stats, len(paths), started)
    finally:
        pool.close()
        pool.join()

    processed = stats["successful_fast5"] + stats["failed_fast5"]
    if processed % args.progress_every:
        emit_progress(stats, len(paths), started)

    with open(summary_path, "w") as handle:
        json.dump(stats, handle, indent=2, sort_keys=True)
        handle.write("\n")

    if stats["successful_fast5"] == 0 or stats["feature_rows"] == 0:
        raise SystemExit("no valid NanoRemix features were extracted; inspect %s" % error_path)
    print(json.dumps(stats, sort_keys=True))
    print("feature_tsv=%s" % feature_path)
    print("coordinates_tsv=%s" % coordinate_path)


if __name__ == "__main__":
    main()
