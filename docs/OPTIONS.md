# Options

Normal use needs `--input`, `--output`, `--task`, `--chemistry`,
`--data-type` and `--mods`. The defaults are prepared feature input and strategy S2.
Use `nanoremix predict --help` for the full command reference.

## Supported Modifications

Binary supports **m6A, m5C, pseudoU** for both chemistries.

| Reference family | RNA002 Multi | RNA004 Multi |
| --- | --- | --- |
| A | m6A, m1A, Inosine, Am | m6A, m1A, Am |
| C | m5C, hm5C, ac4C, FormylC | m5C, hm5C, ac4C |
| G | m7G, o8G | o8G, Gm |
| U | pseudoU, m1psi, m5U, Um | m5U, pseudoU, Um, m1psi |

With A/C/G/U controls, the native Multi axes contain 18 and 16 labels respectively.
The second-generation axis is
`A,C,G,U,Am,m1A,m6A,Inosine,ac4C,m5C,Gm,m7G,m5U,pseudoU,Um`.
`nanoremix catalog` returns exact ordered native class lists.
Aliases include `RNA-editing`, `PseudoU`, `formylC`, `m1Psi` and `T`.
Unsupported selections raise an error.

## Output and Selection

| Option | Default | Effect |
| --- | --- | --- |
| `--also-thirdgen` | Off | Biological data: save pure third-generation output |
| `--also-secondgen` | Off | Biological data: save pure second-generation output |
| `--strategy s1/s2/both` | `s2` | Separate third-generation strategies; historical `s3` aliases `s2` |
| `--canonical matching/all` | `matching` | Retain matching or all canonical controls |
| `--targets` | Unset | Multi-only exact competing class list, instead of `--mods` |
| `--fusion-scope extended/common` | `extended` | Preserve third-only labels, or restrict to shared labels |
| `--output-filter all/retained/modified` | `all` | Select saved prediction rows; audit remains complete |
| `--on-invalid error/mark` | `error` | Abort invalid observations, or mark supported per-row errors |
| `--keep-features` | Off | Preserve additional input columns except probability columns |

Selection happens before normalization and the final decision, not as a filter
on a full-catalog winning label. Probabilities remain conditional on that selection.
A modified class cannot be the sole competing class in its family.

Second-only output retains the requested labels supported by that model. For
example, m5C/hm5C produces C/m5C/hm5C in extended fusion, but only C/m5C in
second-only output. Unsupported comparisons are marked, not forced to canonical.
Classes with no informative second-generation comparison remain third-generation
evidence in extended fusion.

## Thresholds and Fusion Parameters

| Option | Default | Meaning |
| --- | --- | --- |
| `--min-confidence` | `0` | Minimum winning conditional probability |
| `--min-selected-mass` | `0` | Minimum retained mass within the reference-base family |
| `--min-second-mass` | `0` | Minimum analogous second-generation retained mass |
| `--second-weight` | `1` | Second-generation evidence coefficient |
| `--third-weight` | `1` | Third-generation evidence coefficient |
| `--binary-bias` | `0` | Fixed Binary fusion log-odds bias |

Multi third/fused output reports retained third-generation family mass;
second-only output reports retained second-generation family mass.
Binary `selected_family_mass` is 1 for its pair; second-generation retained pair
mass is reported separately. These are not biological modification fractions.
Rejected observations may retain candidate scores but no accepted label.

Coefficients are nonnegative with positive total weight. Binary fusion combines
weighted log odds plus bias; Multi combines shared-class log probabilities.
These defaults are not fitted from user data, and are not a claim to reproduce
paper-specific calibration. Test labels are never used to fit parameters.

## Model and Runtime Options

| Option | Default / purpose |
| --- | --- |
| `--model-root` | Root containing chemistry-specific weights and optional `secondgen/` |
| `--second-model-root` | Override the second-generation directory |
| `--model-checksums` | Trusted JSON index of model-root-relative SHA256 values |
| `--trust-model-files` | Explicit permission to load trusted serialized/Python assets |
| `--condition` | `full`; Binary ablations `mean/std/length` need their own weights |
| `--chunk-size` | `50000` input rows; some identity/FASTA data stays memory-resident |
| `--batch-size` | `256` sequences for second-generation inference |
| `--device` | `cpu`; optional `cuda/mps` require a compatible local environment |
| `--truth-column` | Optional evaluation only, never training |

When overriding the second model directory with checksum checks enabled, it must
remain under `--model-root`. Device options are not a claim of tested accelerator
support.

Optional evaluation reports confusion counts and selected-scope
accuracy/precision/recall/F1 in the run manifest, before output-row filtering.
AUROC/AP and paper-specific benchmark aggregation are outside this CLI.

## Utilities

- `catalog`: supported labels.
- `doctor`: static file/dependency checks; `--strict` fails incomplete setups.
- `validate`: full input/model inference checks, without saving predictions.
- `prepare`: supported legacy feature-table conversion, not raw-signal extraction.

Advanced compatibility option `--mode thirdgen/fused/both/secondgen` may replace,
not accompany, `--data-type`. Here `both` means fused plus third-generation output.
Additional-output flags require `--data-type biological`.
