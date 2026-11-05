import numpy as np
import pandas as pd
import pytest

from nanoremix.labels import BASE, BINARY, NATIVE, SECOND


def fixture_frame(chemistry, task="multi"):
    rng = np.random.default_rng(5417)
    rows = []
    n = 9 if chemistry == "RNA002" else 17
    for i, base in enumerate("ACGU" * 3):
        sequence = "A" * (n // 2) + base + "U" * (n // 2)
        row = {"record_key": f"read-{i:03}", "site_id": f"site-{i:03}", "read_id": f"read-{i:03}",
               "chemistry": chemistry, "center_ref_base": base, "9mer_final": "_".join(sequence),
               "true_label": base}
        for strategy in ("s1", "s2"):
            if task == "binary":
                for mod in BINARY:
                    row[f"{strategy}_prob_{mod}"] = float(rng.uniform(.05, .95))
            else:
                p = rng.uniform(.001, 1, len(NATIVE[chemistry]))
                p /= p.sum()
                for c, value in zip(NATIVE[chemistry], p):
                    row[f"{strategy}_prob_{c}"] = value
        q = rng.uniform(.01, 1, len(SECOND))
        q /= q.sum()
        for c, value in zip(SECOND, q):
            row[f"second_prob_{c}"] = value
        rows.append(row)
    return pd.DataFrame(rows)


@pytest.fixture
def input_factory(tmp_path):
    def make(chemistry="RNA002", task="multi", modify=None):
        frame = fixture_frame(chemistry, task)
        if modify:
            frame = modify(frame)
        path = tmp_path / f"{chemistry}_{task}.csv"
        frame.to_csv(path, index=False)
        return path, frame
    return make
