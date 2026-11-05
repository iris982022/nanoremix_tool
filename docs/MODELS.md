# Model Assets

**Pretrained models are distributed separately from the source repository.**
Download the model asset ZIPs required for your route from the GitHub Release.
Four third-generation assets separate chemistry and task; a fifth asset contains
all 25 shared second-generation checkpoints. This allows a smaller per-file
download while retaining the original inference weights. Model/checkpoint bytes
were verified against the supplied model checksums. Do not substitute a different
chemistry, demo estimators, or historical development models.

Extract each selected model asset ZIP into the repository root to create or
extend `models/`. Each contains `models/SHA256.json`; copies of this index in the
different assets are identical. The complete model set contains 94 required
assets plus the SHA256 index. Normal full-feature S1/S2 inference is included;
optional Binary ablation weights are not included.

## Included Models

| GitHub Release asset | Contents | Required for |
| --- | --- | --- |
| `NanoRemix_v0.3.0_models_RNA002_Binary_20261104.zip` | 18 Binary `.pkl` files | RNA002 Binary S1/S2 |
| `NanoRemix_v0.3.0_models_RNA002_Multi_20261104.zip` | 6 Multi `.pkl` files and both manifests | RNA002 Multi S1/S2 |
| `NanoRemix_v0.3.0_models_RNA004_Binary_20261104.zip` | 30 Binary `.pkl` files | RNA004 Binary S1/S2 |
| `NanoRemix_v0.3.0_models_RNA004_Multi_20261104.zip` | 10 Multi `.pkl` files and both manifests | RNA004 Multi S1/S2 |
| `NanoRemix_v0.3.0_models_SecondGen_20261104.zip` | Network definition and all 25 `.pt` checkpoints | Sequence-based biological fusion or second-generation-only prediction |

Only the route you use is needed locally. A release advertising all routes needs
all five assets. Binary has three pairs: A/m6A, C/m5C and U/pseudoU.
S1 requires five positional models per pair for RNA002 or nine for RNA004.
S2 requires one joint model per pair. Multi S1 likewise needs all positional
components, even when the user selects a single modification.
Second-generation inference requires all 25 checkpoints, not a single fold.

Existing third-generation probability input bypasses third-generation models.
Existing second-generation probability input bypasses the second-generation
network. This does not turn precomputed probabilities into pretrained weights.

## Weights Are Not Sufficient Alone

Distribute the trained objects **with their original inference metadata**:

- Binary: calibrated classifier objects, feature names/order and class orientation.
- Multi: `manifest.json`, all referenced component files, native class order,
  feature schemas, internal aggregation weights, method and temperature.
- Second generation: compatible `scripts/train.py` defining the network and
  sequence encoding; 25 checkpoints with `state_dict`, class order,
  `sequence_length=81`, `architecture=cnn_gru_attention` and positive temperature.
- A trusted SHA256 index and compatible dependency versions.
- Documented upstream feature preprocessing/normalization and a small,
  shareable verification input with frozen expected predictions.

Do not replace the calibrated objects with bare XGBoost estimators or infer
unknown temperatures from filenames. **Model weights** and **second/third fusion
coefficients** are different things. The CLI defaults to second/third coefficients
`1/1` and Binary bias `0`; these are operational defaults, not claimed fitted
paper parameters. Any separately validated operating coefficients must be
distributed and passed explicitly through the [documented options](OPTIONS.md).

## Directory Layout

```text
models/
  SHA256.json
  RNA002/
    binary/binary_classify/
      A_m6A/
        strategy1/m1_calibrated.pkl ... m5_calibrated.pkl
        strategy3/m_all_calibrated.pkl
      C_m5C/...
      U_pseudoU/...
    multi/
      s1/manifest.json
      s1/<all components named in manifest>
      s2/manifest.json
      s2/<all components named in manifest>
  RNA004/
    binary/binary_classify/
      A_m6A/
        strategy1/m1_calibrated.pkl ... m9_calibrated.pkl
        strategy3/m_all_calibrated.pkl
      C_m5C/...
      U_pseudoU/...
    multi/s1/manifest.json
    multi/s1/<all components named in manifest>
    multi/s2/manifest.json
    multi/s2/<all components named in manifest>
  secondgen/
    scripts/train.py
    models/final_25fold_ensemble_81bp_cnn_gru_attention/
      repeat_01/fold_01/best_model.pt
      ... (five repeats x five folds)
      repeat_05/fold_05/best_model.pt
```

The `multi/s2` directory holds the original **s3** model. Keep `strategy: "s3"`
inside its original manifest. Do not rewrite the class axis or feature schema.
Binary ablation models are optional and are not needed for normal prediction.

`SHA256.json` maps model-root-relative filenames to lowercase SHA256 strings,
covering every loaded manifest, component, checkpoint and Python source file.
It accompanies the model bundle. Original weight/checkpoint/network-source files
are unchanged. Multi manifests retain the inference fields, frozen parameters,
component hashes and original manifest hash; training logs and local absolute
paths are omitted. Hashes alone do not establish provenance.

## Install and Verify

Use Python 3.10/3.11 for the original scikit-learn 1.3.0 / XGBoost 2.1.4 objects:

```bash
python -m pip install -e '.[models,second]' -c constraints-models.txt
nanoremix doctor --model-root models --model-checksums models/SHA256.json --strict
```

On macOS, XGBoost also needs an OpenMP runtime. If it reports a missing
`libomp.dylib`, install `libomp` with your package manager (for Homebrew:
`brew install libomp`) before rerunning the check. On other platforms, follow
the native-library error reported by XGBoost.

Extracting a model ZIP places its paths under `models/`. For just one
synthetic-data route, download/extract the matching chemistry/task asset, then:

```bash
nanoremix doctor --model-root models \
  --task binary --chemistry RNA004 --mods m6A \
  --strategy s2 --data-type synthetic --strict
```

For biological data with precomputed second-generation probabilities, the
second-generation asset is unnecessary. Add `--second-source probabilities`
to `doctor`; that inventory excludes network checkpoints. The supplied
probability table is validated during prediction.

The checker detects missing/empty files, unresolved Git LFS pointers, incompatible
Multi metadata and component checksums. `--model-checksums` adds full index checks.
`--strict` returns exit code 2 if the selected assets or dependencies fail checks.
The report distinguishes `files_complete` from `preflight_passed`.
When files and Python packages are present, `doctor` also imports the dependency
libraries in a separate process to detect missing native runtimes such as OpenMP.
It does not import the model bundle's Python source or deserialize weights.

This static check does not load pickle/checkpoint files, validate their internal
tensors, or prove numeric equivalence. `nanoremix validate` accepts prediction
options and runs the actual input/model inference without writing predictions.
Before release, compare predictions to frozen expected outputs using real models.
For a real-weight execution check on artificial features:

```bash
python tools/smoke_models.py --model-root models \
  --output results/model-check --trust-model-files
```

This exercises both chemistries, both tasks, S1/S2, pure third-generation
prediction, 25-checkpoint sequence inference, fusion and component outputs.
It is not a biological accuracy test or a reproduction of paper benchmarks.

Pass `--trust-model-files` only for verified files. Add
`--model-checksums models/SHA256.json` to prediction commands to enforce the index.
See [Security](../SECURITY.md).

## Release Packaging

Keep source code in the repository and upload the five small-by-route model ZIPs
as Release assets. Users download one chemistry/task asset for pure third-generation
inference and add the shared second-generation asset for sequence-based fusion.
Upload every model asset to advertise every route. Publish their exact filenames
and external SHA256 checksums in the Release notes. The source ZIP itself contains
no trained model files; users should never have to retrain advertised models.
