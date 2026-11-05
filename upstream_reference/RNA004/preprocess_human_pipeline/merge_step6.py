#!/usr/bin/env python3
"""Deterministically convert all feature NPZ/CSV pairs into one step6 CSV."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
import pandas as pd


def one_hot_to_sequence(one_hot: np.ndarray) -> str:
    bases = np.asarray(list("ACGU"))
    return "".join(bases[np.argmax(one_hot, axis=1)])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    npz_files = sorted(args.input_dir.glob("split*.npz"))
    if not npz_files:
        raise RuntimeError(f"No NPZ files found in {args.input_dir}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temp = args.output.with_suffix(args.output.suffix + ".tmp")
    feature_names = [
        f"{metric}_output{position}"
        for position in range(1, 10)
        for metric in ("mean", "std", "median", "length")
    ]
    header = [
        "read_id", "ref_name", "ref_start", "ref_end", "kmer_length", "strand",
        "sequence_input", *feature_names,
    ]
    total = 0
    with temp.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        for npz_path in npz_files:
            sidecar = npz_path.with_suffix(".csv")
            if not sidecar.exists():
                raise FileNotFoundError(f"Missing coordinate sidecar: {sidecar}")
            with np.load(npz_path, allow_pickle=False) as data:
                seq = data["seq"]
                stat = data["stat"]
            coords = pd.read_csv(sidecar)
            if len(coords) != len(seq) or len(seq) != len(stat):
                raise ValueError(
                    f"Row mismatch for {npz_path}: coords={len(coords)} seq={len(seq)} stat={len(stat)}"
                )
            if seq.ndim != 3 or seq.shape[1:] != (9, 4):
                raise ValueError(f"Unexpected sequence shape in {npz_path}: {seq.shape}")
            if stat.ndim != 3 or stat.shape[1] != 9 or stat.shape[2] < 6:
                raise ValueError(f"Unexpected statistic shape in {npz_path}: {stat.shape}")
            for index in range(len(seq)):
                values = []
                for position in range(9):
                    values.extend([
                        stat[index, position, 0],
                        stat[index, position, 1],
                        stat[index, position, 5],
                        stat[index, position, 2],
                    ])
                writer.writerow([
                    *coords.iloc[index].tolist()[:6],
                    one_hot_to_sequence(seq[index]),
                    *values,
                ])
                total += 1
    temp.replace(args.output)
    if total == 0:
        raise RuntimeError("Feature extraction produced zero rows")
    print(f"wrote_rows={total} output={args.output}")


if __name__ == "__main__":
    main()
