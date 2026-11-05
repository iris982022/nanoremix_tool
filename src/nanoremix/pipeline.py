import json
import os
import platform
import shutil
import tempfile
from collections import Counter
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

import numpy as np
import pandas as pd

from . import __version__
from .engines import Engines
from .labels import BASE, BINARY, NATIVE, SECOND, code_strategy, label, resolve_labels
from .probability import (
    binary_distribution, binary_product, condition, fuse_multi, row_validity, second_pairs,
)
from .tables import detect_chemistry, prepare_frame, read_table, sha256


def software_versions():
    packages = {}
    for name in ("numpy", "pandas", "scipy", "scikit-learn", "xgboost", "joblib", "torch"):
        try:
            packages[name] = version(name)
        except PackageNotFoundError:
            packages[name] = None
    return {"python": platform.python_version(), **packages}


class Metrics:
    def __init__(self, labels):
        self.labels = labels
        self.counts = np.zeros((len(labels), len(labels)), dtype=np.int64)
        self.totals = Counter()

    def update(self, frame, truth_column):
        self.totals.update(frame.status)
        if not truth_column:
            return
        if truth_column not in frame:
            raise ValueError(f"Missing truth column: {truth_column}")
        for _, row in frame.iterrows():
            try:
                truth = label(row[truth_column])
            except ValueError:
                self.totals["truth_unknown"] += 1
                continue
            if truth not in self.labels:
                self.totals["truth_outside_selected_scope"] += 1
                continue
            self.totals["truth_in_scope"] += 1
            if row.status == "ok":
                self.totals["evaluated_retained"] += 1
                self.counts[self.labels.index(truth), self.labels.index(row.predicted_label)] += 1

    def report(self):
        cm = self.counts
        total = int(cm.sum())
        correct = int(np.trace(cm))
        rows, columns = cm.sum(1), cm.sum(0)
        precision = np.divide(np.diag(cm), columns, out=np.zeros(len(cm)), where=columns > 0)
        recall = np.divide(np.diag(cm), rows, out=np.zeros(len(cm)), where=rows > 0)
        f1 = np.divide(2 * precision * recall, precision + recall, out=np.zeros(len(cm)), where=precision + recall > 0)
        return {
            "status_counts_and_evaluation_counts": dict(self.totals),
            "labels": self.labels, "confusion_counts": cm.tolist(),
            "evaluation_scope": "selected truth labels on retained rows; other truth labels are counted as exclusions",
            "conditional_accuracy": correct / total if total else None,
            "accuracy_with_in_scope_rejections_as_errors": correct / self.totals["truth_in_scope"] if self.totals["truth_in_scope"] else None,
            "per_class": [
                {"label": c, "support_retained": int(rows[i]), "precision": float(precision[i]) if columns[i] else None,
                 "recall": float(recall[i]) if rows[i] else None,
                 "f1": float(f1[i]) if rows[i] or columns[i] else None}
                for i, c in enumerate(self.labels)
            ],
            "note": "Inference performs no fitting. Conditional subset probabilities do not measure all-modification prevalence.",
        }


def output_modes(args):
    if args.data_type and args.mode:
        raise ValueError("Choose --data-type or the advanced --mode override, not both")
    if (args.also_thirdgen or args.also_secondgen) and args.data_type != "biological":
        raise ValueError("--also-thirdgen and --also-secondgen require --data-type biological")
    if args.data_type == "synthetic":
        return ("thirdgen",)
    if args.data_type == "biological":
        modes = ["fused"]
        if args.also_thirdgen:
            modes.append("thirdgen")
        if args.also_secondgen:
            modes.append("secondgen")
        return tuple(modes)
    if args.mode:
        return ("thirdgen", "fused") if args.mode == "both" else (args.mode,)
    raise ValueError("Specify --data-type biological|synthetic, or an explicit advanced --mode")


def selection(args, chemistry, modes):
    return resolve_labels(chemistry, args.task, modes, args.fusion_scope,
                          args.mods, args.targets, args.canonical)


def labels_by_mode(selected, modes):
    result = {mode: selected for mode in modes}
    if "secondgen" in modes and len(modes) > 1:
        supported = []
        for base in "ACGU":
            requested = [c for c in selected if BASE[c] == base]
            shared = [c for c in requested if c in SECOND]
            has_modifications = any(c not in "ACGU" for c in requested)
            # Do not turn an unsupported modification into a certain canonical call.
            if has_modifications and (len(shared) < 2 or not any(c not in "ACGU" for c in shared)):
                continue
            supported.extend(shared)
        result["secondgen"] = tuple(c for c in selected if c in supported)
    return result


def metadata_frame(frame, args):
    if args.keep_features:
        columns = [c for c in frame if not c.startswith("_") and not c.startswith(
            ("s1_prob_", "s2_prob_", "s3_prob_", "second_prob_", "probability_", "prob_")
        )]
    else:
        desired = ["record_key", "input_row", "sample_id", "read_id", "site_id", "reference_name",
                   "chrom", "pos", "strand", "center_ref_base", "canonical_context", args.truth_column]
        columns = list(dict.fromkeys(c for c in desired if c and c in frame))
    return frame[columns].copy()


def produce(frame, args, engines, selected, modes, output_labels):
    indices = np.flatnonzero(frame._status.eq("ok"))
    work = frame.iloc[indices].reset_index(drop=True)
    families = work.center_ref_base.to_numpy()
    strategies = ("s1", "s2") if args.strategy == "both" else (args.strategy,)
    need_third = any(m != "secondgen" for m in modes)
    need_second = any(m != "thirdgen" for m in modes)
    second, second_status = None, np.array([], dtype=object)
    if need_second and len(work):
        second, second_status = engines.second(work)
    raw = {}
    if need_third and len(work):
        for strategy in strategies:
            raw[strategy] = (engines.binary(work, strategy, selected) if args.task == "binary"
                             else engines.multi(work, strategy))
    for mode in modes:
        current_labels = output_labels[mode]
        for strategy in (("none",) if mode == "secondgen" else strategies):
            q = np.full((len(work), len(current_labels)), np.nan)
            selected_mass = np.zeros(len(work))
            second_mass = np.full(len(work), np.nan)
            valid = np.zeros(len(work), dtype=bool)
            status = np.full(len(work), "ok", dtype=object)
            evidence = np.full(len(work), mode, dtype=object)
            if len(work):
                if mode != "thirdgen":
                    status[second_status != "ok"] = second_status[second_status != "ok"]
                if mode == "secondgen":
                    unsupported = ~np.isin(families, [BASE[c] for c in current_labels])
                    status[unsupported] = "unsupported_by_second_generation"
                    evidence[unsupported] = "no_supported_selected_comparison"
                    if current_labels != selected:
                        evidence[~unsupported] = "secondgen_supported_selected_subset"
                if args.task == "binary":
                    pairs = raw.get(strategy, {})
                    if mode != "thirdgen":
                        second_q, masses = second_pairs(second, current_labels)
                        second_mass[:] = 1
                        for mod, mass in masses.items():
                            rows = families == BASE[mod]
                            family_cols = [i for i, c in enumerate(SECOND) if BASE[c] == BASE[mod]]
                            full_mass = second[:, family_cols].sum(1)
                            second_mass[rows] = np.divide(mass, full_mass, out=np.zeros(len(work)), where=full_mass > 0)[rows]
                    if mode == "secondgen":
                        pairs = second_q
                    elif mode == "fused":
                        combined = {}
                        for mod, third_p in pairs.items():
                            values = np.full(len(work), np.nan)
                            relevant = families == BASE[mod]
                            safe = (relevant & (status == "ok") & np.isfinite(third_p)
                                    & (third_p >= 0) & (third_p <= 1) & np.isfinite(second_q[mod]))
                            values[safe] = binary_product(
                                [second_q[mod][safe], third_p[safe]],
                                [args.second_weight, args.third_weight], args.binary_bias,
                            )
                            combined[mod] = values
                        pairs = combined
                    q, valid = binary_distribution(pairs, current_labels, families)
                    selected_mass[:] = 1
                elif mode == "thirdgen":
                    q, selected_mass, valid = condition(raw[strategy], NATIVE[engines.chemistry], current_labels, families)
                elif mode == "secondgen":
                    q, selected_mass, valid = condition(second, SECOND, current_labels, families)
                    second_mass[:] = selected_mass
                else:
                    q, selected_mass, second_mass, valid = fuse_multi(
                        raw[strategy], second, NATIVE[engines.chemistry], selected, families,
                        args.second_weight, args.third_weight,
                    )
                    for base in "ACGU":
                        shared = [c for c in selected if BASE[c] == base and c in SECOND]
                        exclusive = [c for c in selected if BASE[c] == base and c not in SECOND]
                        if len(shared) < 2:
                            evidence[families == base] = "thirdgen_only_in_selected_family"
                        elif exclusive:
                            evidence[families == base] = "fused_shared_classes_with_thirdgen_only_classes"
                status[(status == "ok") & ~valid] = "invalid_or_zero_probability_support"
                valid &= status == "ok"
                q[~valid] = np.nan
                status[(status == "ok") & (selected_mass < args.min_selected_mass)] = "below_selected_mass"
                if mode != "thirdgen":
                    status[(status == "ok") & (second_mass < args.min_second_mass)] = "below_second_mass"
            confidence = np.full(len(work), np.nan)
            candidates = np.full(len(work), "", dtype=object)
            good = row_validity(q)
            if good.any():
                confidence[good] = q[good].max(1)
                candidates[good] = np.asarray(current_labels)[q[good].argmax(1)]
            status[(status == "ok") & (confidence < args.min_confidence)] = "below_confidence"
            output = metadata_frame(frame, args)
            output["chemistry"] = engines.chemistry
            output["task"] = args.task
            output["data_type"] = args.data_type or "unspecified"
            output["mode"] = mode
            output["strategy"] = strategy
            output["status"] = frame._status.to_numpy()
            output["candidate_label"] = ""
            output["predicted_label"] = ""
            output["confidence"] = np.nan
            output["selected_family_mass"] = np.nan
            output["second_selected_family_mass"] = np.nan
            output["evidence_source"] = ""
            output.loc[indices, "status"] = status
            output.loc[indices, "candidate_label"] = candidates
            output.loc[indices, "predicted_label"] = np.where(status == "ok", candidates, "")
            output.loc[indices, "confidence"] = confidence
            output.loc[indices, "selected_family_mass"] = selected_mass
            output.loc[indices, "second_selected_family_mass"] = second_mass
            output.loc[indices, "evidence_source"] = evidence
            for i, c in enumerate(current_labels):
                output[f"prob_{c}"] = np.nan
                output.loc[indices, f"prob_{c}"] = q[:, i]
            if args.on_invalid == "error":
                invalid = output.status.str.startswith(("invalid", "missing", "insufficient", "81nt"))
                if invalid.any():
                    examples = output.loc[invalid, ["input_row", "status"]].head(5).to_dict("records")
                    raise ValueError(f"Invalid selected-family input rows: {examples}; use --on-invalid mark to retain explicit failures")
            yield mode, strategy, output


def validate_options(args):
    modes = output_modes(args)
    reserved = {"predicted_label", "candidate_label", "status", "confidence", "input_row",
                "evidence_source", "selected_family_mass", "second_selected_family_mass"}
    if args.truth_column in reserved or (args.truth_column and args.truth_column.startswith("prob_")):
        raise ValueError("Truth must use a separate input annotation column, not a generated prediction field")
    if args.task == "multi" and args.condition != "full":
        raise ValueError("Multi uses the trained full-feature schemas; signal ablations are Binary-only")
    if args.task == "multi" and args.binary_bias != 0:
        raise ValueError("--binary-bias is Binary-only")
    if args.binary_bias != 0 and "fused" not in modes:
        raise ValueError("--binary-bias applies only to fused output")
    for name in ("min_confidence", "min_selected_mass", "min_second_mass"):
        if not np.isfinite(getattr(args, name)) or not 0 <= getattr(args, name) <= 1:
            raise ValueError(f"--{name.replace('_', '-')} must be within [0,1]")
    for name in ("second_weight", "third_weight"):
        if not np.isfinite(getattr(args, name)) or getattr(args, name) <= 0:
            raise ValueError("Fusion weights must be finite and positive")
    if not np.isfinite(args.binary_bias):
        raise ValueError("Binary bias must be finite")
    if args.chunk_size < 1 or args.batch_size < 1:
        raise ValueError("Chunk and batch sizes must be positive")
    if any(mode != "thirdgen" for mode in modes):
        sources = bool(args.second_probabilities) + bool(args.sequence_column) + bool(args.context_fasta)
        if args.second_probabilities and args.sequence_column:
            sources -= 1  # Optional sequence verification in the matched probability table.
        if sources > 1:
            raise ValueError("Use one second-generation source: a matched table, inline sequence, or context FASTA")
        if args.context_fasta:
            if not args.fasta_id_column or not args.center_column or args.coordinate_base is None:
                raise ValueError("FASTA input requires --fasta-id-column, --center-column, and --coordinate-base")
    return modes


def run(args):
    modes = validate_options(args)
    source = Path(args.input).expanduser().resolve()
    source_hash = sha256(source)
    head = read_table(source, nrows=min(args.chunk_size, 1000))
    if head.empty:
        raise ValueError("Input has no observations")
    if args.truth_column and args.truth_column not in head:
        raise ValueError(f"Input lacks truth annotation column: {args.truth_column}")
    chemistry = detect_chemistry(head, args.chemistry)
    selected = selection(args, chemistry, modes)
    output_labels = labels_by_mode(selected, modes)
    engines = Engines(args, chemistry)
    destination = Path(args.output).expanduser().resolve() if args.output else None
    if destination is None and args.command == "predict":
        raise ValueError("predict requires --output")
    if destination and destination.exists():
        raise ValueError(f"Output already exists; choose a new directory: {destination}")
    if destination and args.model_root and destination.is_relative_to(Path(args.model_root).expanduser().resolve()):
        raise ValueError("Output must not be inside the model directory")
    if destination:
        destination.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=".nanoremix-", dir=destination.parent if destination else None))
    reports, seen, offset = {}, set(), 0
    try:
        for chunk in read_table(source, chunksize=args.chunk_size):
            if args.data_type and "data_type" in chunk and not chunk.data_type.str.strip().str.lower().eq(args.data_type).all():
                raise ValueError("Input data_type conflicts with the explicitly selected --data-type")
            preparation_kind = "probabilities" if modes == ("secondgen",) else args.input_kind
            frame = prepare_frame(chunk, chemistry, offset, source_hash, args.task, selected, preparation_kind, args.condition)
            keys = set(frame.record_key)
            if keys & seen:
                raise ValueError("Duplicate record_key across chunks")
            seen.update(keys)
            offset += len(frame)
            for mode, strategy, output in produce(frame, args, engines, selected, modes, output_labels):
                key = f"{mode}/{strategy}"
                metrics = reports.setdefault(key, Metrics(output_labels[mode]))
                metrics.update(output, args.truth_column)
                if args.command != "predict":
                    continue
                folder = staging / mode / strategy
                folder.mkdir(parents=True, exist_ok=True)
                audit_cols = [c for c in output if c not in {f"prob_{x}" for x in output_labels[mode]}]
                audit_path = folder / "selection_audit.tsv"
                output[audit_cols].to_csv(audit_path, sep="\t", mode="a", header=not audit_path.exists(), index=False)
                keep = np.ones(len(output), dtype=bool) if args.output_filter == "all" else output.status.eq("ok")
                if args.output_filter == "modified":
                    keep &= ~output.predicted_label.isin(list("ACGU"))
                path = folder / "predictions.tsv"
                output.loc[keep].to_csv(path, sep="\t", mode="a", header=not path.exists(), index=False)
        if sha256(source) != source_hash:
            raise ValueError("Input changed while inference was running; refusing an inconsistent result")
        report = {
            "version": __version__, "created_utc": datetime.now(timezone.utc).isoformat(),
            "input": str(source), "input_sha256": source_hash, "input_rows": offset,
            "chemistry": chemistry, "task": args.task, "selected_labels": selected,
            "data_type": args.data_type or "unspecified",
            "output_policy": "data_type" if args.data_type else "explicit_legacy_mode",
            "primary_mode": modes[0], "output_modes": modes,
            "selected_labels_by_mode": output_labels,
            "excluded_labels_by_mode": {
                mode: [c for c in selected if c not in labels] for mode, labels in output_labels.items()
            },
            "unsupported_secondgen_modifications": [c for c in selected if c not in SECOND] if "secondgen" in modes else [],
            "selected_modifications": [c for c in selected if c not in "ACGU"],
            "paper_to_code_strategy": {"s1": "s1", "s2": "s3"},
            "mode": args.mode, "strategy": args.strategy, "fusion_scope": args.fusion_scope,
            "class_support": {
                c: {
                    "second_generation": c in SECOND,
                    "third_generation": c in (set("ACGU") | set(BINARY) if args.task == "binary" else NATIVE[chemistry]),
                } for c in selected
            },
            "weight_status": "fixed_user_parameters_not_fitted_here",
            "second_weight": args.second_weight, "third_weight": args.third_weight, "binary_bias": args.binary_bias,
            "thresholds": {k: getattr(args, k) for k in ("min_confidence", "min_selected_mass", "min_second_mass")},
            "class_selection": "select requested labels plus specified canonical controls, family-condition, then normalize and classify",
            "output_filter": args.output_filter,
            "both_semantics": "separate outputs; no implicit S1/S2 consensus",
            "model_assets_sha256": engines.assets, "model_checksums_verified": bool(args.model_checksums),
            "second_source_used": any(mode != "thirdgen" for mode in modes),
            "model_retrained": False,
            "software": software_versions(), "configuration": vars(args).copy(),
            "reports": {k: v.report() for k, v in reports.items()},
            "validation_limit": "Software validation is not independent biological validation; original trained assets must be supplied.",
        }
        if report["second_source_used"]:
            for field in ("second_probabilities", "context_fasta"):
                if getattr(args, field):
                    report[field + "_sha256"] = sha256(getattr(args, field))
        (staging / "run_manifest.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
        if args.command == "predict":
            if destination.exists():
                raise ValueError("Output was created concurrently; refusing to replace it")
            os.rename(staging, destination)
        return report
    finally:
        if staging.exists():
            shutil.rmtree(staging)
