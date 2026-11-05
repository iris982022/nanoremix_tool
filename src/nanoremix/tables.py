import gzip
import hashlib
from pathlib import Path

import numpy as np
import pandas as pd

from .labels import BASE, canonical_context, label


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def separator(path):
    suffixes = Path(path).suffixes
    return "\t" if ".tsv" in suffixes or ".tab" in suffixes else ","


def read_table(path, **kwargs):
    return pd.read_csv(path, sep=separator(path), dtype=str, keep_default_na=False, **kwargs)


def detect_chemistry(frame, requested):
    declared = set(frame["chemistry"].str.upper()) if "chemistry" in frame else set()
    if declared - {"RNA002", "RNA004"} or len(declared) > 1:
        raise ValueError("Input chemistry must be consistently RNA002 or RNA004")
    lengths = set()
    if "9mer_final" in frame:
        lengths = {len(x.split("_")) if "_" in x else len(x) for x in frame["9mer_final"].astype(str)}
    inferred = {9: "RNA002", 17: "RNA004"}.get(next(iter(lengths))) if len(lengths) == 1 else None
    if requested != "auto":
        if declared and declared != {requested}:
            raise ValueError("Declared input chemistry conflicts with --chemistry")
        return requested
    result = next(iter(declared)) if declared else inferred
    if result is None or (inferred and inferred != result):
        raise ValueError("Cannot unambiguously infer chemistry; supply --chemistry RNA002 or RNA004")
    return result


def prepare_frame(frame, chemistry, offset, source_hash, task, selected, input_kind, condition):
    frame = frame.reset_index(drop=True).copy()
    frame["input_row"] = np.arange(offset, offset + len(frame))
    if "record_key" not in frame:
        frame["record_key"] = [f"{source_hash}:{i}" for i in frame.input_row]
    if frame.record_key.eq("").any() or frame.record_key.duplicated().any():
        raise ValueError("record_key must be nonempty and unique")
    if "chemistry" in frame and not frame.chemistry.str.upper().eq(chemistry).all():
        raise ValueError("Input chemistry mismatch")
    if "condition" in frame and not frame.condition.eq(condition).all():
        raise ValueError("Input signal condition mismatch")
    families, sequences, statuses = [], [], []
    for _, row in frame.iterrows():
        try:
            sequence = canonical_context(row["9mer_final"], chemistry) if "9mer_final" in frame else ""
            if sequence:
                base = sequence[len(sequence) // 2]
                if "center_ref_base" in frame:
                    given = label(row.center_ref_base)
                    if given not in "ACGU" or given != base:
                        raise ValueError("reference/context disagreement")
            else:
                base = label(row.get("center_ref_base", ""))
                if base not in "ACGU":
                    raise ValueError("noncanonical reference base")
            families.append(base)
            sequences.append(sequence)
            statuses.append("ok" if base in {BASE[c] for c in selected} else "outside_selected_family")
        except ValueError:
            families.append("")
            sequences.append("")
            statuses.append("invalid_context_or_reference")
    frame["center_ref_base"] = families
    frame["canonical_context"] = sequences
    frame["_status"] = statuses
    if input_kind == "features":
        if "9mer_final" not in frame:
            raise ValueError("Feature inference requires 9mer_final (9 or 17 context tokens)")
        frame["9mer_final"] = ["_".join(s) for s in sequences]
        k = 5 if chemistry == "RNA002" else 9
        stats = ("mean", "median", "std", "length") if condition == "full" else (condition,)
        numeric = [f"{stat}_output{i}" for i in range(1, k + 1) for stat in stats]
        missing = set(numeric) - set(frame)
        if missing:
            raise ValueError(f"Missing signal columns: {sorted(missing)}")
        for column in numeric:
            values = pd.to_numeric(frame[column], errors="coerce")
            bad = ~np.isfinite(values)
            if column.startswith(("std_", "length_")):
                bad |= values < 0
            frame.loc[bad & frame._status.eq("ok"), "_status"] = "invalid_signal"
            frame[column] = values
    return frame


def read_fasta(path):
    result, name = {}, None
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                tokens = line[1:].split()
                if not tokens:
                    raise ValueError("FASTA header has no identifier")
                name = tokens[0]
                if name in result:
                    raise ValueError(f"Duplicate FASTA identifier: {name}")
                result[name] = ""
            elif name is None:
                raise ValueError("FASTA sequence appears before its identifier")
            else:
                result[name] += line.upper().replace("U", "T")
    if not result:
        raise ValueError("Empty FASTA")
    return result


def sequence81(frame, sequence_column=None, fasta=None, id_column=None, center_column=None,
               strand_column=None, coordinate_base=None):
    sequences, statuses = [], []
    for _, row in frame.iterrows():
        status, sequence = "ok", ""
        try:
            if sequence_column:
                sequence = str(row[sequence_column]).strip().upper().replace("U", "T")
            else:
                name = str(row[id_column])
                if name not in fasta:
                    raise ValueError("missing_fasta_id")
                coordinate = float(row[center_column])
                if not np.isfinite(coordinate) or coordinate != int(coordinate):
                    raise ValueError("invalid_center_coordinate")
                center = int(coordinate) - coordinate_base
                full = fasta[name]
                if center < 40 or center + 40 >= len(full):
                    raise ValueError("insufficient_81nt_flanks")
                sequence = full[center - 40:center + 41]
                strand = str(row[strand_column]) if strand_column else "+"
                if strand not in {"+", "-"}:
                    raise ValueError("invalid_strand")
                if strand == "-":
                    sequence = sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]
            if len(sequence) != 81 or set(sequence) - set("ACGT"):
                raise ValueError("invalid_81nt_sequence")
            if sequence[40].replace("T", "U") != row.center_ref_base:
                raise ValueError("81nt_center_mismatch")
            short = row.canonical_context.replace("U", "T")
            if short and sequence[40 - len(short) // 2:41 + len(short) // 2] != short:
                raise ValueError("81nt_context_mismatch")
        except (ValueError, KeyError, TypeError, OverflowError) as error:
            known = {
                "missing_fasta_id", "invalid_center_coordinate", "insufficient_81nt_flanks",
                "invalid_strand", "invalid_81nt_sequence", "81nt_center_mismatch", "81nt_context_mismatch",
            }
            status = "missing_sequence_metadata" if isinstance(error, KeyError) else (
                str(error) if str(error) in known else "invalid_81nt_input"
            )
            sequence = ""
        sequences.append(sequence)
        statuses.append(status)
    return np.array(sequences), np.array(statuses)


def probability_columns(frame, prefixes, labels, aliases=None):
    aliases = aliases or {}
    columns = []
    for name in labels:
        candidates = list(dict.fromkeys(
            prefix + candidate for prefix in prefixes for candidate in (name, aliases.get(name, name))
        ))
        present = [x for x in candidates if x in frame]
        if len(present) != 1:
            raise ValueError(f"Expected exactly one probability column for {name}; found {present}, candidates {candidates}")
        columns.append(present[0])
    return frame[columns].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
