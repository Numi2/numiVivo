# Imported CellTypist inference

Native frozen gene-space inference now supports the complete Allen Parse L1
model: 888 features and seven classes. All model features match the Parse axis.
The implementation normalizes using the complete sparse source row, applies
log1p and training scaling, caps standardized values above 10, and returns
independent sigmoid scores. It does not fit models or make labels authoritative.

## Verified scope

- 32 numerical software rows match NumPy/SciPy.
- 14 malformed model, axis and count cases are rejected.
- All 267 cells in the first frozen parent PBS window (Donor10_PBS, run 1)
  match scikit-learn 1.9.0 labels exactly. Maximum decision error is 9.60e-14;
  maximum independent sigmoid error is 2.95e-15. The model records training
  with scikit-learn 1.7.2.
- Product source compiled with Swift optimization in a standalone harness.
  The actual scoped Omics CLI also built and reproduced all 267 retained real-cell
  results exactly across 529,404 count records. Output receipt hashes were verified;
  incorrect stream hashes, truncated records, incorrect row totals/cardinalities
  and existing destinations were rejected without publishing partial output.
  Full-package qualification and H5AD orchestration, whole-cohort comparison, measured
  memory scaling and held-out annotation accuracy remain unfinished.

The source window was chosen by control condition and original order, without
label or score selection. The model was trained on Parse PBS controls: this is
interoperability evidence, not independent biological validation. It does not
change the frozen B-cell cohort or establish perturbation prediction.

## Evidence and reproduction

### Streaming command

Build the actual scoped CLI with `Tools/Omics/H5AD/build.sh <build-directory> --with-cli`.
Run `numivivo-omics singlecell-celltypist-stream --model model.json --plan count-plan.json --stream-sha256 <sha256> --output new-bundle < counts.bin`.
The binary stream contains little-endian UInt64 pairs: row in the low 32 bits
and feature in the high 32 bits, followed by the positive integer count.
Rows and feature indices must be ordered; normalization includes every source
feature. The embedded count plan binds cell and feature identities, row
cardinalities and optional row totals. The receipt binds model, plan, results,
input stream and implementation identity. Predictions are independent sigmoid
scores, not calibrated probabilities or authoritative labels. The cell axis
remains resident; million-cell memory scaling is unqualified.

After unpacking the original archive below into a new directory, run:

```sh
python3 check_stream_cli.py <build-directory>/numivivo-omics <unpacked-evidence> <new-check-directory>
```

The check uses the retained real inputs and the adjacent compressed window
plan. It rebases original prepared rows 42:309 to 0:267 without filtering cells
or features. It verifies exact retained results, receipt hashes and five
transactional rejection cases. This test harness materializes the small window;
the product command streams count records. The original archive retains the
standalone harness version; the CLI extension lives in the repository sources.

### Original numerical evidence

[Archive](evidence.tar.gz) and [manifest](manifest.json) retain 28 files,
including model, JSON export, source, binaries, sparse real inputs, outputs,
checks and source-range receipts. Every archived member was independently
read and hash-verified. The source model SHA256 is
`959938f78eb6e0394a9da0d31a65d57316fcf1be0df9cdde3222cea3842f9a63`.

Unpack into a new directory, with Swift and NumPy/SciPy/scikit-learn available:

```sh
swiftc -O VivoCellTypistReference.swift main.swift -o native
python check.py
./native model.json real-input.json real-output.json
python check_sklearn.py
swiftc -O VivoCellTypistReference.swift rejections.swift -o rejections
./rejections model.json
```

`check_real.py` preserves original paths to the published Parse reader and
axes/chunks for fresh network execution. The sparse real inputs are embedded
for offline replay. `runtime.json` records the compiler and source hash.
Extraction is restricted to the pinned model hash and specified reconstruction
types; estimator objects become inert state holders. This is not a general
pickle importer. The initial extraction failed on `Scaler`; the observed
upstream key `Scaler_` corrected it.

Sources: [Allen downloads](https://apps.allenimmunology.org/aifi/resources/parse-10m-cytokines/downloads/),
[pinned CellTypist](https://github.com/Teichlab/celltypist/tree/fe357564a6625d3b1732a022fd39f18e55696e80).
Parse data/model derivatives retain Parse Biosciences and Allen Institute
attribution, CC BY-NC 4.0. Upstream code retains the included license.
