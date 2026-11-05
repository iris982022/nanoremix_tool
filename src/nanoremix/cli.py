import argparse
import importlib
import json
import sys
from pathlib import Path

from . import __version__
from .labels import BASE, BINARY, NATIVE, ORDER, SECOND, normalize_strategy


def parser():
    root = argparse.ArgumentParser(
        prog="nanoremix",
        description="RNA002/RNA004 Binary/Multi inference with explicit class selection and generation provenance.",
    )
    root.add_argument("--version", action="version", version=__version__)
    commands = root.add_subparsers(dest="command", required=True)
    for command in ("predict", "validate"):
        p = commands.add_parser(command, help="Run inference" if command == "predict" else "Run complete input/model validation without saving predictions")
        p.add_argument("--input", required=True, help="CSV/TSV, optionally gzip-compressed")
        p.add_argument("--output", help="New output directory (required for predict)")
        p.add_argument("--task", choices=["binary", "multi"], required=True)
        p.add_argument("--chemistry", choices=["auto", "RNA002", "RNA004"], default="auto")
        select = p.add_mutually_exclusive_group(required=True)
        select.add_argument("--mods", help="Comma-separated modifications, or all; matching canonical controls are added")
        select.add_argument("--targets", help="Multi-only exact competing class list, including canonical controls when desired")
        p.add_argument("--canonical", choices=["matching", "all"], default="matching")
        routing = p.add_mutually_exclusive_group(required=True)
        routing.add_argument("--data-type", choices=["biological", "synthetic"],
                             help="Biological data defaults to fused output; synthetic data uses third-generation-only output")
        routing.add_argument("--mode", choices=["thirdgen", "fused", "both", "secondgen"],
                             help="Advanced compatibility override; use --data-type for the normal workflow")
        p.add_argument("--also-thirdgen", action="store_true", help="Biological data: also save pure third-generation predictions")
        p.add_argument("--also-secondgen", action="store_true", help="Biological data: also save pure second-generation predictions for supported selected labels")
        p.add_argument("--strategy", type=normalize_strategy, choices=["s1", "s2", "both"], default="s2",
                       help="Paper strategy; historical s3 is accepted as an alias for s2. both writes separate predictions.")
        p.add_argument("--input-kind", choices=["features", "probabilities"], default="features")
        p.add_argument("--condition", choices=["full", "mean", "std", "length"], default="full")
        p.add_argument("--model-root", help="Directory with RNA002/, RNA004/, and optional secondgen/")
        p.add_argument("--trust-model-files", action="store_true", help="Explicit consent to load verified local pickle/checkpoint/Python model files")
        p.add_argument("--model-checksums", help="JSON mapping model-root-relative file names to verified SHA256 values")
        p.add_argument("--second-model-root", help="Override secondgen directory containing scripts/train.py and models/")
        p.add_argument("--second-probabilities", help="Original 15-class CSV/TSV matched by explicit identity columns")
        p.add_argument("--join-on", help="Comma-separated identity keys; never a positional row-order join")
        p.add_argument("--sequence-column", help="81-nt sequence column, or optional validation column in the second probability table")
        p.add_argument("--context-fasta", help="Reference or same-run event-sequence FASTA, never inferred from short k-mers")
        p.add_argument("--fasta-id-column", help="Input column matching the first FASTA header token")
        p.add_argument("--center-column", help="Target nucleotide coordinate within that FASTA sequence")
        p.add_argument("--coordinate-base", type=int, choices=[0, 1])
        p.add_argument("--strand-column", help="+/- orientation; absent means already transcript-oriented")
        p.add_argument("--fusion-scope", choices=["common", "extended"], default="extended")
        p.add_argument("--second-weight", type=float, default=1.0)
        p.add_argument("--third-weight", type=float, default=1.0)
        p.add_argument("--binary-bias", type=float, default=0.0, help="Frozen Binary fusion log-odds bias; never fitted on prediction labels")
        p.add_argument("--min-confidence", type=float, default=0.0)
        p.add_argument("--min-selected-mass", type=float, default=0.0)
        p.add_argument("--min-second-mass", type=float, default=0.0)
        p.add_argument("--output-filter", choices=["all", "retained", "modified"], default="all")
        p.add_argument("--on-invalid", choices=["error", "mark"], default="error")
        p.add_argument("--truth-column", help="Optional evaluation labels; never used for predictions or parameter fitting")
        p.add_argument("--chunk-size", type=int, default=50000)
        p.add_argument("--batch-size", type=int, default=256)
        p.add_argument("--device", choices=["cpu", "cuda", "mps"], default="cpu")
        p.add_argument("--keep-features", action="store_true", help="Keep input feature columns, excluding all input probability columns")
    catalog = commands.add_parser("catalog", help="List supported modifications by chemistry, task and generation")
    catalog.add_argument("--chemistry", choices=["RNA002", "RNA004"])
    doctor = commands.add_parser("doctor", help="Check required model files, manifests and hashes without executing models")
    doctor.add_argument("--model-root", default="models")
    doctor.add_argument("--chemistry", choices=["RNA002", "RNA004"], help="Omit to check both chemistries")
    doctor.add_argument("--task", choices=["binary", "multi"], help="Omit to check both tasks")
    doctor.add_argument("--mods", default="all")
    doctor.add_argument("--strategy", type=normalize_strategy, choices=["s1", "s2", "both"], default="both")
    doctor.add_argument("--data-type", choices=["biological", "synthetic"], default="biological")
    doctor.add_argument("--second-source", choices=["sequence", "probabilities"], default="sequence",
                        help="probabilities: the user will supply precomputed second-generation scores")
    doctor.add_argument("--condition", choices=["full", "mean", "std", "length"], default="full")
    doctor.add_argument("--model-checksums")
    doctor.add_argument("--strict", action="store_true", help="Exit with status 2 if files or dependencies are missing/invalid")
    prepare = commands.add_parser("prepare", help="Convert legacy extracted features into a prepared context table")
    prepare.add_argument("--input", required=True)
    prepare.add_argument("--output", required=True)
    prepare.add_argument("--chemistry", required=True, choices=["RNA002", "RNA004"])
    prepare.add_argument("--format", required=True, choices=["step4", "extracted", "step6"])
    return root


def catalog(chemistry=None):
    result = {}
    for chem in (chemistry,) if chemistry else NATIVE:
        result[chem] = {
            "context_nt": 9 if chem == "RNA002" else 17,
            "signal_descriptors": 20 if chem == "RNA002" else 36,
            "binary_modifications": BINARY,
            "multi_thirdgen_labels": NATIVE[chem],
            "multi_fused_common_labels": [c for c in ORDER if c in NATIVE[chem] and c in SECOND],
            "multi_fused_third_only_labels": [c for c in NATIVE[chem] if c not in SECOND],
            "secondgen_labels": SECOND,
        }
    return result


def doctor(model_root=None, **scope):
    from .assets import inspect_models

    return inspect_models(model_root, check_imports=True, **scope)


def prepare(args):
    from .tables import separator

    output = Path(args.output)
    if output.exists():
        raise ValueError("Prepared output exists; choose a new filename")
    k = 5 if args.chemistry == "RNA002" else 9
    features = importlib.import_module(f"nanoremix._vendor.binary{k}.features")
    io = importlib.import_module(f"nanoremix._vendor.binary{k}.io")
    if k == 5 and args.format in {"step4", "extracted"}:
        raw = io.read_step4_tsv(args.input) if args.format == "step4" else io.read_extracted_csv(args.input)
        table = features.generate_9mer_candidates(raw)
    elif k == 9 and args.format == "step6":
        table = features.generate_17mer_candidates(io.read_step6_csv(args.input))
    else:
        raise ValueError("RNA002 accepts step4/extracted; RNA004 accepts step6")
    output.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(output, index=False, sep=separator(output))
    return {"prepared_rows": len(table), "output": str(output),
            "note": "Attach explicit 81-nt sequence identity/coordinates separately for fused inference."}


def main(argv=None):
    p = parser()
    args = p.parse_args(argv)
    try:
        if args.command == "catalog":
            result = catalog(args.chemistry)
        elif args.command == "doctor":
            result = doctor(**{k: v for k, v in vars(args).items() if k not in ("command", "strict")})
        elif args.command == "prepare":
            result = prepare(args)
        else:
            from .pipeline import run
            result = run(args)
            result = {k: result[k] for k in (
                "input_rows", "chemistry", "task", "data_type", "selected_labels",
                "selected_labels_by_mode", "mode", "primary_mode", "output_modes", "strategy",
            )}
            if args.command == "predict":
                result["output"] = str(Path(args.output).resolve())
            else:
                result["validation"] = "complete; no prediction files saved"
        print(json.dumps(result, indent=2))
        if args.command == "doctor" and args.strict and not result["preflight_passed"]:
            return 2
        return 0
    except (ValueError, KeyError, FileNotFoundError, ImportError, OSError) as error:
        print(f"nanoremix: {error}", file=sys.stderr)
        return 2
