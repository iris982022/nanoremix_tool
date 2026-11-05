"""Generate synthetic probability fixtures, not biological observations."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from nanoremix.labels import BINARY, NATIVE, SECOND


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=False)
    rng = np.random.default_rng(5417)
    for chemistry in NATIVE:
        for task in ("binary", "multi"):
            rows = []
            length = 9 if chemistry == "RNA002" else 17
            for i, base in enumerate("ACGU" * 4):
                context = "A" * (length // 2) + base + "U" * (length // 2)
                row = {
                    "sample_id": "SYNTHETIC_NOT_EXPERIMENTAL",
                    "record_key": f"demo-{i:03d}", "read_id": f"read-{i:03d}",
                    "site_id": f"site-{i:03d}", "chemistry": chemistry,
                    "9mer_final": "_".join(context), "center_ref_base": base,
                    "true_label": base,
                }
                for strategy in ("s1", "s2"):
                    if task == "binary":
                        row.update({f"{strategy}_prob_{c}": rng.uniform(.05, .95) for c in BINARY})
                    else:
                        values = rng.dirichlet(np.ones(len(NATIVE[chemistry])))
                        row.update({f"{strategy}_prob_{c}": p for c, p in zip(NATIVE[chemistry], values)})
                second = rng.dirichlet(np.ones(len(SECOND)))
                row.update({f"second_prob_{c}": p for c, p in zip(SECOND, second)})
                rows.append(row)
            pd.DataFrame(rows).to_csv(output / f"{chemistry}_{task}.csv", index=False)
    print(f"Four synthetic probability tables created in {output}")


if __name__ == "__main__":
    main()
