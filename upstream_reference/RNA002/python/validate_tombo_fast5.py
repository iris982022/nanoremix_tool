#!/usr/bin/env python3
"""Fail a Tombo stage when a deterministic FAST5 sample has no valid output.

Read-only validation. This exists because legacy Tombo may exit with status 0
even when every read failed to align.
"""

from __future__ import print_function

import argparse
import json
import os
import sys

import h5py


CORRECTED = "/Analyses/RawGenomeCorrected_000/BaseCalled_template"


def decode_text(value):
    if isinstance(value, bytes):
        return value.decode("utf-8", "replace")
    return str(value)


def iter_fast5(root):
    for directory, dirnames, filenames in os.walk(root):
        dirnames.sort()
        for name in sorted(filenames):
            if name.endswith(".fast5"):
                yield os.path.join(directory, name)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fast5-dir", required=True)
    parser.add_argument("--sample-size", type=int, default=2000)
    parser.add_argument("--min-ready", type=int, default=10)
    parser.add_argument("--min-ready-fraction", type=float, default=0.01)
    args = parser.parse_args()

    if args.sample_size < 1 or args.min_ready < 1:
        raise SystemExit("sample-size and min-ready must be positive")

    checked = 0
    ready = 0
    corrected_incomplete = 0
    failed_status = 0
    absent = 0
    unreadable = 0
    for path in iter_fast5(args.fast5_dir):
        if checked >= args.sample_size:
            break
        checked += 1
        try:
            with h5py.File(path, "r") as handle:
                if CORRECTED not in handle:
                    absent += 1
                    continue
                group = handle[CORRECTED]
                status = decode_text(group.attrs.get("status", "success")).strip().lower()
                if status not in ("success", "none", ""):
                    failed_status += 1
                    continue
                if "Events" in group and "Alignment" in group and len(group["Events"]) > 0:
                    ready += 1
                else:
                    corrected_incomplete += 1
        except Exception:
            unreadable += 1

    fraction = float(ready) / checked if checked else 0.0
    payload = {
        "checked_fast5": checked,
        "ready_fast5": ready,
        "ready_fraction": fraction,
        "corrected_incomplete": corrected_incomplete,
        "failed_tombo_status": failed_status,
        "corrected_absent": absent,
        "unreadable": unreadable,
        "sample_size_limit": args.sample_size,
        "minimum_ready": args.min_ready,
        "minimum_ready_fraction": args.min_ready_fraction,
    }
    print(json.dumps(payload, sort_keys=True))
    if checked == 0:
        raise SystemExit("ERROR: no FAST5 files found")
    if ready < args.min_ready or fraction < args.min_ready_fraction:
        raise SystemExit(
            "ERROR: Tombo validation failed: %d/%d ready (%.4f); refusing success marker"
            % (ready, checked, fraction)
        )
    print("[TOMBO VALIDATION PASS] %d/%d sampled FAST5 ready" % (ready, checked))


if __name__ == "__main__":
    main()
