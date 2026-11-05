# Input Contract

## Data type and output selection

Use `--data-type biological` for fused inference, optionally adding
`--also-thirdgen` and/or `--also-secondgen`. Use `--data-type synthetic` for
third-generation-only inference; second input is not read in this workflow.
An explicit legacy `--mode` can be used instead, never together with data type.
No filename-based data-type inference is performed.

## Prepared feature input

UTF-8 CSV (`.csv`) or TSV (`.tsv`/`.tab`), optionally gzip-compressed.
One row is one candidate sequence-signal observation. The CLI does not aggregate
observations into sites or read-level top-10% benchmark scores.

Required:

- `9mer_final`: RNA002 has nine tokens; RNA004 has 17. Use canonical characters
  or underscore-delimited tokens, e.g. `A_C_G_U_A_C_G_U_A`. Modified tokens are
  canonicalized before feature construction, never used as a prediction label.
- `mean_output1` ... `mean_outputK`
- `median_output1` ... `median_outputK`
- `std_output1` ... `std_outputK`
- `length_output1` ... `length_outputK`, where K is 5 or 9.

Signal channels must be finite; dispersion and length must be nonnegative.
Column order is aligned to the original model schema. Values must already use
the original preprocessing and normalization; the package does not infer a scale.
Binary `--condition mean|std|length` uses only that signal channel and requires
the matching ablation weights, not full-feature weights.

Optional:

- `chemistry`: consistently `RNA002` or `RNA004`.
- `data_type`: if supplied with `--data-type`, every row must match that declaration.
- `center_ref_base`: canonical A/C/G/U (T aliases U); must agree with the context.
- `record_key`: globally unique, nonempty. Otherwise a source-hash/row key is
  generated. Do not use a nonunique site ID as the record key for multiple reads.
- `sample_id`, `read_id`, `site_id`, `reference_name`, `chrom`, `pos`, `strand`.
- User-named `--truth-column`, for evaluation only.

Chemistry auto-detection uses a consistent chemistry column or context token
length, not filenames or the truth labels. Mixed-chemistry input is rejected.
For probability-only input without short context, supply `center_ref_base` plus a
chemistry column or explicit `--chemistry`.

## Legacy extraction adapters

```bash
nanoremix prepare --input extracted002.csv --format extracted \
  --chemistry RNA002 --output prepared002.csv
nanoremix prepare --input step4.tsv --format step4 \
  --chemistry RNA002 --output prepared002-step4.csv
nanoremix prepare --input step6.csv --format step6 \
  --chemistry RNA004 --output prepared004.csv
```

These reuse the original candidate constructors. They are not FAST5/POD5
basecalling/event-alignment pipelines. Review the original extraction schema and
normalization before calling them. Output is a prepared CSV; explicit 81-nt
context identity must be attached separately for sequence-based fusion.

## Precomputed third-generation probabilities

Choose `--input-kind probabilities`. Values must be finite in [0,1].

Binary input has one pair probability per requested modification and strategy:
`s1_prob_m6A`, `s2_prob_m6A`, etc. The value is P(modification | its Binary pair),
not an unconditioned Multi score.

Multi requires the **complete chemistry-native class axis**, with columns
`s1_prob_A`, ..., `s1_prob_o8G`, and/or the equivalent `s2_prob_*` columns.
Each native row must sum to one within absolute tolerance 1e-6.
Use `nanoremix catalog` for class names. Supplying a selected subset alone loses
the original probability mass and is deliberately rejected. Predictions only
contain the selected output classes.

Legacy `s3_prob_*` columns are accepted for paper S2. Supplying both s2 and s3
columns for the same class is ambiguous and rejected.

## Second-generation input

Three alternatives are supported. Fused mode does not substitute missing values.

1. Inline probabilities: one column `second_prob_<label>` per canonical label in
   the 15-class second-generation axis.
2. A separate `--second-probabilities` table with an explicit `--join-on` identity.
   `probability_background_A`, ..., `probability_PseudoU` from the original
   exporter are also accepted. Exactly one alias column per class is allowed.
3. An 81-nt sequence source plus the original model source and 25 checkpoints.

Separate probability tables are joined by identity, never by row order.
Right-hand keys must be nonempty and unique. The user must include every needed
dimension (sample, read, strand/site) to avoid reusing the wrong observation.
There is no assumption that a site ID is unique across samples.

With `--sequence-column sequence81`, sequences must contain only A/C/G/T/U and
have length 81. The center must agree with the canonical reference base, and the
central 9/17-nt region must agree with the short context if present. If used
together with a separate probability table, this option verifies that table's
sequence column; it does not rerun the network.

```bash
nanoremix predict --input prepared004.csv --task multi --chemistry RNA004 \
  --mods m6A,m5C --data-type biological --strategy s2 \
  --sequence-column sequence81 --model-root models --trust-model-files \
  --output results/sequence-fused
```

FASTA extraction requires explicit coordinates:

```bash
nanoremix predict --input prepared004.csv --task multi --chemistry RNA004 \
  --mods m6A,m5C --data-type biological --strategy s2 \
  --context-fasta transcripts.fa --fasta-id-column reference_name \
  --center-column target_coordinate --coordinate-base 0 --strand-column strand \
  --model-root models --trust-model-files --output results/fasta-fused
```

The coordinate is within the named FASTA sequence, not an implicitly converted
genomic/transcript coordinate. `-` takes the reverse complement of the extracted
window. For already transcript-oriented FASTA, omit `--strand-column`; do not
reverse-complement twice. Incomplete 40-nt flanks and ambiguous bases are rejected.
Record whether sequences came from a reference or an event-derived sequence.

## Invalid rows

`--on-invalid error` aborts the run without a final output directory.
`--on-invalid mark` retains explicit per-row errors such as invalid signal,
missing second match, context mismatch or insufficient flanks. Schema errors,
unknown classes, duplicate record keys and untrusted/incompatible models always
abort. Rows outside requested reference families are ordinary audited exclusions.
