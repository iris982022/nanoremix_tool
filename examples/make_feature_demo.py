"""Generate artificial features for model-loading checks, not accuracy evaluation."""

import argparse
import csv
from pathlib import Path
import random


def generate(output):
    output = Path(output)
    if output.exists():
        raise ValueError("Demo output exists; choose a new directory")
    output.mkdir(parents=True)
    rng = random.Random(1041)
    sequences = []
    for base in "ACGUACGU":
        sequence = [rng.choice("ACGT") for _ in range(81)]
        sequence[40] = base.replace("U", "T")
        sequences.append("".join(sequence))
    for chemistry, k in (("RNA002", 5), ("RNA004", 9)):
        rows = []
        for i, sequence in enumerate(sequences):
            row = {
                "record_key": f"artificial-{i:02d}", "sample_id": "artificial-feature-demo",
                "read_id": f"artificial-read-{i:02d}", "chemistry": chemistry,
                "center_ref_base": sequence[40].replace("T", "U"),
                "sequence81": sequence,
                "9mer_final": "_".join(sequence[41-k:40+k].replace("T", "U")),
            }
            for position in range(1, k + 1):
                for stat in ("mean", "median", "std", "length"):
                    row[f"{stat}_output{position}"] = round(rng.uniform(.1, 2), 6)
            rows.append(row)
        with (output / f"{chemistry}_features.csv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True)
    generate(parser.parse_args().output)
