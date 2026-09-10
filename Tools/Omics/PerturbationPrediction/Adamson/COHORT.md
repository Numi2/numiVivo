# Original-author Adamson assignment cohort

The [original author library](https://github.com/thomasmaxwellnorman/perturbseq_demo/blob/4328ac9ed29a6afe7162e1822c2f7632448667ba/perturbseq/cell_population.py#L203)
defines an unambiguous cell as `number_of_cells == 1`, `good_coverage == True`,
and a guide identity other than `*`. We apply that exact rule to the restored
original GEO records, excluding missing records explicitly. The public author
code was retrieved at the pinned revision and matches SHA-256
`556cd498d006ac5c1553f3ef4aaa74cee4bad6444142b0c3c7553dff4ea1d7c1`.
This reusable loader rule is evidence for assignment selection; its accompanying
demo experiment does not establish this UPR experiment's control identities.

The [paper's methods](https://escholarship.org/content/qt9961h6m4/qt9961h6m4_noSplash_f242fbdc7b9f06f490340b23dde1a1ae.pdf)
describe coverage as guide reads per guide UMI. Good coverage requires an upper
coverage mode plus minimum read and UMI support. Multiple qualifying identities
are classified as multiples; none are unidentifiable. We use the deposited
original flags, without re-estimating the coverage threshold. The source's
`number of cells` field participates in this assignment rule; it is not an
independent replicate or donor count.

| Original GEO assignment cohort | Cells | RNA UMIs |
| --- | ---: | ---: |
| Selected | 50,440 | 781,977,660 |
| Excluded | 14,897 | 257,880,138 |
| Complete source | 65,337 | 1,039,857,798 |

The exclusion ledger retains original row, full barcode, guide and every failed
criterion. Reasons overlap: 14,644 do not have one assigned identity; 6,038 lack
good coverage; 112 carry `*`; 80 have no original guide record. Exclusive reason
combinations also partition the excluded cohort exactly. No RNA outcome,
expression threshold, GO support or model score determines selection.

The cohort has 94 guide groups and retains 82 of the previously captured 90
candidate gene prefixes. ATF4, ATF6, C7orf26, EIF2AK3, ERN1, PSMA1, PSMD12 and
XBP1 have no selected cells. The corresponding guides remain in the exclusion
ledger. IER3IP1 and TIMM23 remain unsupported by the frozen GO capture, leaving
80 supported candidates in the selected cohort. Annotations were not recaptured
or selected after scoring; no predictions have been fitted or scored.

The labels `63(mod)_pBA580` (4,595 selected cells) and `Gal4-4(mod)_pBA582`
(646) remain unresolved experimental control candidates. `62(mod)_pBA581` has
zero selected cells. Their names and these counts are not sufficient primary
experimental definitions, so no label is yet assigned the authoritative control
role. The [experimental-role audit](EXPERIMENTAL_ROLES.md) independently reproduces
this inventory and identifies 94 selected guide groups versus 93 guides in the
paper's experiment summary. Primary controls and the complete guide roster must
be reconciled before fitting. The
[frozen prediction protocol](PROTOCOL.md) is unchanged.

## Native execution and independent verification

The [native source-bound selection](../../H5AD/CELL_SELECTION.md) plan references
the restored H5AD hash
`a33bab097da97dc46414473d1761675940c59df8a53ee8517f60f98f950c334e`.
Swift/HDF5 scans all source counts and directly aggregates the selected rows.
It retains the complete source snapshot and explicit original-row mapping.
All 32,738 features remain; selected canonical nonzeros total 181,795,274.

`prepare_cohort.py` computes an independent integer SciPy reference using sparse
membership multiplication over 64-gene CSC blocks. `check_cohort.py` independently
reconstructs assignment decisions with Python's CSV reader, then verifies every
selected native count, cell identity, quality value, original-row mapping and
technical GEM-group membership. Selected and excluded UMI totals reconstruct the
complete source total. Numeric GEM groups still do not establish biological
replication or donor-held-out prediction.

The full native release product rebuilt on the physical M4 Pro. Executable
SHA-256 is `4eee4da0f3cdfef0ea47d341a6ce75d27042d2ff8e017b57be3f5e9602bb1884`.
Native aggregation took 35.70 seconds, peak RSS 489,668,608 bytes, and exact
native replay passed. Report SHA-256 is
`fd23394c260a1be7daecfc692a444cdef6b858b8bc7f55b89192bf559fae9a78`.
The new selection suite passed 26 checks; the existing streaming regression
suite passed 25. Timings describe this CPU/HDF5 run and do not establish a GPU
or end-to-end advantage over scverse.

Initial test-attempt failures remain recorded: one fixture run started before
executable transfer completed; the rejection harness initially expected lowercase
`unknown` while the CLI emitted `Unknown`; and the independent CSV checker
initially expected `True` while original records use `TRUE`. The corrected
checks pass without changing the native algorithm or cohort selection.

## Reproduce

```sh
python prepare_cohort.py --source restored.h5ad --barcodes GSM2406681_10X010_barcodes.tsv.gz --guides GSM2406681_10X010_cell_identities.csv.gz --author-code cell_population.py --mapping restored/pseudobulk-plan.json --out cohort
numivivo singlecell-h5ad-pseudobulk restored.h5ad --plan cohort/plan.json --output cohort/native
numivivo singlecell-h5ad-pseudobulk-verify cohort/native
python check_cohort.py --cohort cohort --barcodes GSM2406681_10X010_barcodes.tsv.gz --guides GSM2406681_10X010_cell_identities.csv.gz --full-report restored/pseudobulk/report.json --full-reference prepared/reference.npz --descriptors descriptors --out cohort/check.json
```

Native commands use the configured native HDF5 library. Python prepares evidence
and independently verifies it; the native CLI does not call Python. The
[evidence manifest](evidence/2026-09-10-author-cohort/manifest.json) preserves the
plans, exclusions, counts, receipts, checks and actual failed attempts. The large
source H5AD and executable are retained externally with hashes. A redundant
remote copy of the original deposited source was removed only after exact hash
and open-handle checks, retaining the verified original locally and the complete
historical and restored annotated bundles remotely.
