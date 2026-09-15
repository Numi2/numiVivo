# Current public Kang Metal NB run

This record captures a matched native run on the public Kang et al. 2018
IFNB-stimulated B-cell study ([doi:10.1038/nbt.4042](https://doi.org/10.1038/nbt.4042)).
The prepared H5AD contains 2,651 B cells, 15,706 genes and 16 paired
donor-condition pseudobulks. The mapping keeps the source batch metadata
explicitly unreported and uses the declared paired-donor design; no batch
correction is inferred.

The run used the source revision in commit `9df690091f486ee0e6b256e494999731b0369056`
and the release executable SHA-256
`99d72481ca29315e2c3ec2e8b5cf1bca1e31e5cc6fcb4192242d32fd625caad5` on a
physical Apple M4 with macOS 26.6 and Metal 4. The prepared H5AD SHA-256 is
`6e5609e5bd0635ac1bd2a8f9a2bf98ad4fa9a012e8e4e2ad926436d073d745f0`.
The mapping SHA-256 is
`297a523f383f62fd6f8ae8800f8ac395e701c57149b46e8c58f14c343c76b84f`.
The CPU and Metal plan hashes are, respectively,
`a50236019992ce19193f8ed47edb02848e58f7d9ad003bac65cc97b3e36cb812` and
`2780b301a15c08e13aaa09b4098c13e82c7692afa534b8b52ba82e70aa6f25c9`.

## Matched results

| Run | Wall time | Tested features | Successful final fits | Final-fit backend | Explicit numerical failures |
| --- | ---: | ---: | ---: | --- | ---: |
| Current CPU control | 30.499 s | 5,400 / 5,400 | 5,400 | CPU FP64 | 0 |
| Current Metal plan | 330.366 s | 3,803 / 5,400 | 3,803 | 3,766 `metalFP32`, 37 CPU outlier profiles | 1,597 |

Both runs retained 6,812 low-expression and 3,494 positive-support rank-deficient
features. The Metal plan did not silently fall back: its 1,597 failures remain
`numericalFailure` records, and the 37 CPU fits are the declared outlier path
that reuses an already-computed CPU profile. The Metal analysis receipt replayed
successfully with `singlecell-analysis-verify` in 340.791 s; the wrapper observed
a child maximum RSS of 782,123,008 bytes. The stored result fingerprint is
`248b7d03416d3071a12373f7d5eef1efa5fe9671369ead11bb5d1573fba5be61`.

Among the 3,803 features tested by both plans, the maximum absolute log2 effect
difference was `3.022550222303577e-6`, the mean absolute difference was
`2.389010076016236e-9`, and effect signs agreed for all features. The five
predeclared IFNB marker genes (ISG15, IFIT1, IFIT3, MX1 and OAS1) had positive
effects in both reports. Adjusted-p-value differences remain part of the
FP32 numerical profile and are not treated as a calibration result.

## Commands and retained failure

The executable was run with the following sequence, using new output paths for
every receipt and store:

```text
numivivo singlecell-h5ad-import prepared.h5ad --plan mapping.json --output imported
numivivo singlecell-run imported/manifest.json --store store --output count-receipt-current.json
numivivo singlecell-analyze count-receipt-current.json --plan analysis-cpu.json --store store --output analysis-cpu-current2-receipt.json
numivivo singlecell-analyze count-receipt-current.json --plan analysis-metal.json --store store --output analysis-metal-receipt.json
numivivo singlecell-analysis-verify analysis-metal-receipt.json --store store
```

The first full Metal attempt, before commit `9df6900`, was killed with
`SIGKILL` after 139.201 s while each gene fit still owned a command queue and
three large shared buffers. Commit `9df6900` pools those resources behind the
pipeline cache; the full run then completed without a kill and with stable
resident memory. This repairs lifecycle pressure, but does not convert the
remaining per-feature FP32 failures into successful fits.

This is execution and same-count numerical-comparison evidence for one public
paired-donor study. It does not establish GPU speedup, complete Metal model
fitting, FDR or interval calibration, independent-study replication, causal
mechanism, or prediction of a phenotype, tissue state or treatment outcome.
