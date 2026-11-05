"""Adapters retain the original trained feature schemas and ensemble manifests."""

import importlib
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

from .labels import BASE, BINARY, NATIVE, SECOND, SECOND_ORIGINAL, code_strategy
from .probability import binary_product, row_validity
from .tables import probability_columns, read_fasta, read_table, sequence81, sha256


class Engines:
    def __init__(self, args, chemistry):
        self.args = args
        self.chemistry = chemistry
        self.k = 5 if chemistry == "RNA002" else 9
        self.root = Path(args.model_root).expanduser().resolve() if args.model_root else None
        self.loaded = {}
        self.assets = {}
        self.second_table = None
        self.fasta = None
        self.checksums = json.loads(Path(args.model_checksums).read_text()) if args.model_checksums else None

    def trusted(self, path):
        if not self.args.trust_model_files:
            raise ValueError("Model loading executes trusted pickle/Python assets; pass --trust-model-files only for verified local models")
        path = Path(path).resolve()
        if not path.is_file():
            raise FileNotFoundError(f"Model asset is missing: {path}")
        digest = sha256(path)
        if self.checksums is not None:
            if self.root is None:
                raise ValueError("--model-checksums requires --model-root")
            try:
                key = path.relative_to(self.root).as_posix()
            except ValueError as error:
                raise ValueError("Checked model asset lies outside --model-root") from error
            if self.checksums.get(key) != digest:
                raise ValueError(f"Missing or mismatched SHA256 for model asset: {key}")
        self.assets[str(path)] = digest
        return path

    def require_third_dependencies(self):
        missing = [name for name in ("sklearn", "scipy", "xgboost", "joblib")
                   if importlib.util.find_spec(name) is None]
        if missing:
            raise ImportError(f"Feature inference needs the models extra; missing {missing}. Install nanoremix[models].")

    def binary(self, frame, strategy, selected):
        result = {}
        for mod in BINARY:
            if mod not in selected:
                continue
            output = np.full(len(frame), np.nan)
            rows = frame.center_ref_base.eq(BASE[mod]).to_numpy()
            if not rows.any():
                result[mod] = output
                continue
            part = frame.loc[rows].reset_index(drop=True)
            if self.args.input_kind == "probabilities":
                prefixes = [f"{strategy}_prob_"]
                if strategy == "s2":
                    prefixes.append("s3_prob_")
                output[rows] = probability_columns(part, prefixes, [mod])[:, 0]
            else:
                if self.root is None:
                    raise ValueError("Feature inference requires --model-root")
                self.require_third_dependencies()
                features = importlib.import_module(f"nanoremix._vendor.binary{self.k}.features")
                models = importlib.import_module(f"nanoremix._vendor.binary{self.k}.models")
                strategy_name = f"strategy{code_strategy(strategy)[-1]}"
                directory = self.root / self.chemistry / "binary" / "binary_classify" / f"{BASE[mod]}_{mod}" / strategy_name
                if self.args.condition != "full":
                    directory /= f"only_{self.args.condition}"
                part = (features.add_nine_mer_positions(part) if self.k == 5
                        else features.add_seventeen_mer_positions(part))
                probabilities = []
                for spec in features.binary_model_specs(strategy_name, self.args.condition):
                    path = directory / f"{spec.model_name}_calibrated.pkl"
                    if str(path) not in self.loaded:
                        self.trusted(path)
                        self.loaded[str(path)] = models.load_model(path)
                    model = self.loaded[str(path)]
                    classes = list(getattr(model, "classes_", []))
                    if classes not in ([0, 1], [1, 0], ["0", "1"], ["1", "0"]):
                        raise ValueError(f"Unrecognized binary class orientation: {classes}")
                    matrix = models.align_feature_matrix(features.build_feature_matrix(part, spec), model)
                    probabilities.append(models.binary_positive_probability(model, matrix))
                output[rows] = binary_product(probabilities) if len(probabilities) > 1 else probabilities[0]
            result[mod] = output
        return result

    def multi(self, frame, strategy):
        if self.args.input_kind == "probabilities":
            prefixes = [f"{strategy}_prob_"]
            if strategy == "s2":
                prefixes.append("s3_prob_")
            return probability_columns(frame, prefixes, NATIVE[self.chemistry])
        if self.root is None:
            raise ValueError("Feature inference requires --model-root")
        self.require_third_dependencies()
        key = f"multi:{strategy}"
        if key not in self.loaded:
            directory = self.root / self.chemistry / "multi" / strategy
            manifest_path = self.trusted(directory / "manifest.json")
            manifest = json.loads(manifest_path.read_text())
            name = f"train_5mer_{code_strategy(strategy)}" if self.k == 5 else f"train_9mer_limit_{code_strategy(strategy)}"
            module = importlib.import_module(f"nanoremix._vendor.multi{self.k}.{name}")
            schemas = module.schemas(self.k, code_strategy(strategy))
            if manifest.get("schemas") != schemas or len(manifest["components"]) != len(schemas):
                raise ValueError("Multi manifest feature schemas/component count differ from the original trained representation")
            fusion = manifest["fusion"]
            if not isinstance(fusion.get("gate"), bool) or set(fusion.get("families", {})) != set("ACGU"):
                raise ValueError("Multi manifest requires a valid family-fusion configuration")
            for spec in fusion["families"].values():
                weights = np.asarray(spec["weights"], dtype=float)
                temperature = float(spec["temperature"])
                if (weights.shape != (len(schemas),) or not np.isfinite(weights).all()
                        or (weights < 0).any() or not np.isclose(weights.sum(), 1)
                        or spec["method"] not in {"arithmetic", "geometric"}
                        or not np.isfinite(temperature) or temperature <= 0):
                    raise ValueError("Invalid frozen Multi fusion weights/method/temperature")
            for component in manifest["components"]:
                path = (directory / component["pickle"]).resolve()
                if not path.is_relative_to(directory.resolve()):
                    raise ValueError("Model component escapes its model directory")
                self.trusted(path)
            model = module.load_model(directory)
            if tuple(model.classes_) != NATIVE[self.chemistry]:
                raise ValueError("Multi model class order does not match the released chemistry-specific label space")
            if model.manifest.get("strategy") != code_strategy(strategy):
                raise ValueError("Multi manifest strategy does not match the requested paper strategy")
            self.loaded[key] = model
        return np.asarray(self.loaded[key].predict_proba(frame), dtype=float)

    def second(self, frame):
        args = self.args
        aliases = dict(zip(SECOND, SECOND_ORIGINAL))
        if args.second_probabilities:
            keys = [x.strip() for x in args.join_on.split(",")] if args.join_on else []
            if not keys or any(not x for x in keys):
                raise ValueError("--second-probabilities requires explicit --join-on identity columns")
            if self.second_table is None:
                table = read_table(args.second_probabilities)
                if not set(keys) <= set(table) or table[keys].eq("").any().any() or table.duplicated(keys).any():
                    raise ValueError("Second probability join keys must exist, be nonempty, and be unique")
                self.second_table = table
            if not set(keys) <= set(frame) or frame[keys].eq("").any().any():
                raise ValueError("Input is missing valid second-probability identity keys")
            payload = self.second_table.copy()
            payload["_second_present"] = "yes"
            joined = frame[keys].merge(payload, on=keys, how="left", sort=False, validate="many_to_one")
            values = probability_columns(joined, ["second_prob_", "probability_"], SECOND, aliases)
            status = np.where(joined["_second_present"].eq("yes"), "ok", "missing_second_match").astype(object)
            status[(status == "ok") & ~row_validity(values)] = "invalid_second_probability"
            if args.sequence_column:
                if args.sequence_column not in joined:
                    raise ValueError("Second probability table lacks --sequence-column")
                check = frame.copy()
                check[args.sequence_column] = joined[args.sequence_column].to_numpy()
                _, seq_status = sequence81(check, sequence_column=args.sequence_column)
                status[(status == "ok") & (seq_status != "ok")] = seq_status[(status == "ok") & (seq_status != "ok")]
            return values, status
        has_probabilities = any(c.startswith(("second_prob_", "probability_")) for c in frame)
        if has_probabilities and (args.sequence_column or args.context_fasta):
            raise ValueError("Both inline second probabilities and a sequence source were supplied; choose one explicitly")
        if has_probabilities:
            values = probability_columns(frame, ["second_prob_", "probability_"], SECOND, aliases)
            return values, np.where(row_validity(values), "ok", "invalid_second_probability")
        if args.sequence_column:
            if args.sequence_column not in frame:
                raise ValueError(f"Missing sequence column: {args.sequence_column}")
            seq, status = sequence81(frame, sequence_column=args.sequence_column)
        elif args.context_fasta:
            if not args.fasta_id_column or not args.center_column or args.coordinate_base is None:
                raise ValueError("FASTA extraction requires ID column, target-center column, and explicit coordinate base")
            if self.fasta is None:
                self.fasta = read_fasta(args.context_fasta)
            seq, status = sequence81(
                frame, fasta=self.fasta, id_column=args.fasta_id_column,
                center_column=args.center_column, strand_column=args.strand_column,
                coordinate_base=args.coordinate_base,
            )
        else:
            raise ValueError("Fused/secondgen mode requires second probabilities or an explicit 81-nt sequence source")
        values = np.full((len(frame), len(SECOND)), np.nan)
        good = status == "ok"
        if good.any():
            values[good] = self._second_network(seq[good])
        return values, status

    def _second_network(self, sequences):
        import torch

        root = Path(self.args.second_model_root).expanduser().resolve() if self.args.second_model_root else (
            self.root / "secondgen" if self.root else None
        )
        if root is None:
            raise ValueError("Sequence inference requires --model-root or --second-model-root")
        source = self.trusted(root / "scripts/train.py")
        if "second_source" not in self.loaded:
            name = "nanoremix_verified_second_source"
            spec = importlib.util.spec_from_file_location(name, source)
            module = importlib.util.module_from_spec(spec)
            sys.modules[name] = module
            spec.loader.exec_module(module)
            if tuple(module.CLASSES) != SECOND_ORIGINAL:
                raise ValueError("Second-generation source class axis mismatch")
            self.loaded["second_source"] = module
        module = self.loaded["second_source"]
        model_dir = root / "models/final_25fold_ensemble_81bp_cnn_gru_attention"
        paths = [model_dir / f"repeat_{r:02d}/fold_{f:02d}/best_model.pt" for r in range(1, 6) for f in range(1, 6)]
        if len(list(model_dir.rglob("best_model.pt"))) != 25:
            raise ValueError("Second-generation ensemble requires exactly 25 original checkpoints")
        unique, inverse = np.unique(sequences, return_inverse=True)
        encoded = module.encode_sequences(unique.tolist())
        result = np.zeros((len(unique), len(SECOND)), dtype=float)
        for path in paths:
            self.trusted(path)
            checkpoint = torch.load(path, map_location="cpu", weights_only=False)
            if (tuple(checkpoint["classes"]) != SECOND_ORIGINAL or checkpoint["sequence_length"] != 81
                    or checkpoint["architecture"] != "cnn_gru_attention"):
                raise ValueError("Second-generation checkpoint schema mismatch")
            temperature = float(checkpoint["temperature"])
            if not np.isfinite(temperature) or temperature <= 0:
                raise ValueError("Invalid second-generation temperature")
            network = module.CNGRUAttention()
            network.load_state_dict(checkpoint["state_dict"], strict=True)
            network.to(self.args.device).eval()
            with torch.inference_mode():
                for start in range(0, len(unique), self.args.batch_size):
                    batch = encoded[start:start + self.args.batch_size]
                    logits = network(module.to_onehot(batch, self.args.device))
                    result[start:start + len(batch)] += torch.softmax(logits / temperature, 1).cpu().numpy()
            del network, checkpoint
        result /= 25
        if not row_validity(result).all():
            raise ValueError("Second-generation ensemble produced an invalid distribution")
        return result[inverse]
