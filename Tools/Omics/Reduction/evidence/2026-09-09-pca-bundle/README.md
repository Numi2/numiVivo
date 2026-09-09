# Standalone PCA qualification, 2026-09-09

Complete Norman filtered release: 111,445 cells, 33,694 genes and 361,582,621
source nonzeros, selecting 2,000 HVGs and 20 principal components. Source is
Zenodo record 13350497, `NormanWeissman2019_filtered.h5ad`; source SHA256 and
exact compiled source/tool hashes are in `checks.json`. The source is retained
in the remote execution bundle, rather than duplicated in Git.

Native Apple M4 Pro runs use the isolated HDF5 runtime described in `host.txt`.
`publish.log.gz`, `verify.log.gz` and `legacy.log.gz` retain command output and
macOS time measurements. Maximum resident bytes are distinct from the smaller
peak footprint also reported by that tool.

| Native operation | Seconds | Maximum resident bytes |
| --- | ---: | ---: |
| Standalone publish | 136.11 | 443203584 |
| Standalone reconstruction | 137.46 | 443432960 |
| Current-product legacy pseudobulk/PCA | 158.85 | 1386496000 |

Reproduce from the source with the qualified executable and HDF5 runtime:

```sh
numivivo singlecell-h5ad-pca original.h5ad --plan input-plan.json --output new-pca
numivivo singlecell-h5ad-pca-verify new-pca
```

`bundle/` contains every generated artifact except the source H5AD. Decompress
`.gz` files to their original names before using the reference adapter; hashes
of decompressed bytes match `receipt.json`. The canonical plan explicitly retains
projection centers. `input-plan.json` is the original submitted plan.

The reference adapter checks every binary coordinate/value and compares all
scores, loadings, feature statistics, identities and QC with the previously
qualified full Norman report. The entire current legacy report remains byte
identical to that previous report. Run `pca_bundle_reference.py` followed by
`check_reference.py` as documented in `../../PCA_BUNDLE.md`; derived reference
JSON is reconstructible and omitted here. Independent Scanpy checks include
projection centers, HVGs, eigenspectrum and aligned score error (1.0971e-12).
Reference Python ran on a different host; no CPU/scverse speed comparison is made.

`tests.log.gz` records 69 passing Swift tests in 19 suites. `release.log.gz`
retains the successful release build and warnings. `regressions/` records the
standalone CLI (12 commands, 7 expected rejections), streamed H5AD and legacy
CSR/CSC checks. Expected rejection outputs include rehashed score/model tampering,
no-overwrite and work/cache limits; they are controls, not unexplained failures.

Metadata, QC, PCA basis and fitted arrays remain resident. This is complete-cohort
storage and numerical qualification, not million-cell execution, parallel/Metal
performance or held-out biological prediction.
