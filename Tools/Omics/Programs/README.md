# Deterministic expression programs

Add `programs` to a `singlecell-analyze` plan or to the existing
`singlecell-h5ad-pseudobulk` plan. This uses the same analysis/publication/replay
owners as counts and PCA. Plans without programs retain their earlier behavior.

```json
"programs": {
  "definitions": [{
    "id": "declared-signature",
    "organism": "NCBITaxon:9606",
    "featureNamespace": "HUMAN_GENE_SYMBOL",
    "sourceURI": "urn:example:signature",
    "sourceVersion": "1",
    "sourceDescription": "Authored definition with explicit provenance",
    "members": [{"featureID": "GENE_A", "weight": 1}, {"featureID": "GENE_B", "weight": -1}],
    "minimumWeightCoverage": 1
  }],
  "maximumScoreValues": 5000000,
  "maximumUpdates": 100000000
}
```

## Definition and interpretation

For each cell and program, the score is the sum of signed member weights times
log1p library-normalized expression, divided by the sum of absolute matched
weights. All source genes contribute to library size. A positive equal-weight
program gives mean log-normalized expression; negative weights encode an
explicit opposing signature. This is not Scanpy's sampled-control `score_genes`,
GSEA, a fitted classifier, or a calibrated measurement of pathway activity.
There is no implicit control-gene sampling, automatic label, or maximum-score
cell assignment. Single-member programs are supported as marker measurements.

The resident route uses the analysis plan's normalization target and retained
cells. The H5AD route uses its reduction normalization target when present,
otherwise 10,000. Results record the actual target. Scoring uses original RNA,
independent of corrected PCA or embeddings. It does not modify source counts,
labels or PCA.

Matching is exact against feature IDs, in source feature order. Duplicate
members, zero/nonfinite weights and mismatching sample organisms are rejected.
There is no automatic alias, case, Ensembl-version or ortholog conversion.
`featureNamespace` declares the author's namespace; the current count model has
no independent namespace field with which to verify that declaration.

Every definition includes source URI, version and description. A SHA256 of the
sorted-key JSON definition binds those fields, ordered members, weights and
coverage policy to the result. This fingerprints the supplied definition; it
does not fetch or independently authenticate the external source. Public
benchmark snapshots below retain their downloaded source bytes and SHA256.

Missing genes are separate from measured zero expression. By default every
member must exist. An explicit lower `minimumWeightCoverage` permits a partial
program, with missing IDs, matched indices, effective weights and retained
absolute-weight fraction recorded. Remaining weights are renormalized; a
partial program is therefore a different observable. All-missing programs
always fail. A zero-count library has a null score, while a nonempty library
with no detected members has score zero. Per-cell detected-member counts remain
visible alongside scores.

## Sparse execution and limits

`VivoSingleCellProgramAccumulator` owns resolution and arithmetic for both
resident CSR and streamed canonical CSR/CSC. H5AD scoring adds one scan after
QC, with no dense cells-by-genes matrix. With PCA enabled, the combined workflow
uses four source scans: QC, programs, HVG moments and selected-cache writing.
The existing `reductionStorage.sourcePasses` counts the three scans used by the
reduction component, including shared QC, rather than all enabled analyses.

Scores and detection counts are resident cells-by-programs arrays. Allocation
is bounded before construction; sparse membership updates are counted and
bounded before each update. Defaults are five million score values and 100
million updates; configurable hard caps are 20 million and one billion.
Definitions are capped at 256 programs, 10,000 members each and 100,000 total
memberships. Existing source/count bounds remain. This does not qualify
million-cell annotation or GPU execution.

## Experimental checks

The fixed benchmark definitions are `kang-programs.json` and
`baron-programs.json`. The source snapshots in `sources/` come from MSigDB's
[interferon alpha response](https://www.gsea-msigdb.org/gsea/msigdb/cards/HALLMARK_INTERFERON_ALPHA_RESPONSE)
and [pancreas beta-cell](https://www.gsea-msigdb.org/gsea/msigdb/cards/HALLMARK_PANCREAS_BETA_CELLS)
Hallmark sets, attributed to Liberzon et al., PMID 26771021. All downloaded
`geneSymbols` are retained with unit weights. The explicit coverage threshold
was set to 0.8 before scoring. No members are selected based on the observed
benchmark response. A separate five-gene IFN definition reproduces the prior
integration program; it is not the full Hallmark set.

Kang uses all 2,651 B cells from eight paired donors. Its full 97-gene IFN
program has 95 exact matches; TENT5A and WARS1 are reported missing, with no
implicit renaming. Baron uses all 8,569 human cells and all 40 beta-cell program
genes. The independent reference normalizes the same sparse counts using
Scanpy, computes weighted sparse products and detection counts, and maps cell
identities explicitly.

Predeclared descriptive observations are increased mean IFN scores under
stimulation within each Kang donor, and increased mean beta-program scores in
author-labelled beta cells within each Baron donor. These are association
checks, not independently validated cell labels or held-out prediction.
Reference results retain the observed direction even if it contradicts the
expectation. The numerical status and biological observations are separate.

```sh
python check_cli.py --binary /path/to/numivivo --prepared /data/kang --programs kang-programs.json --imported /data/kang/imported --out /new/kang
python check_reference.py --prepared /data/kang --bundle /new/kang/bundle --benchmark kang --out /new/kang/reference.json
python check_cli.py --binary /path/to/numivivo --prepared /data/baron --programs baron-programs.json --out /new/baron
python check_reference.py --prepared /data/baron --bundle /new/baron/bundle --benchmark baron --out /new/baron/reference.json
python check_formats.py --binary /path/to/numivivo --out /new/formats
```

Native tests separately exercise signed arithmetic, missing/default coverage,
provenance fingerprint changes, empty libraries, species and budget rejection,
weight scaling and strict decoding. Numerical fixtures do not supply biological
evidence. Learned reference mapping, calibrated annotation confidence, marker
specificity across independent tissues, held-out perturbation prediction and
cross-study biological qualification remain open.

## Recorded qualification, 2026-09-09

The production release passes 41 Swift tests in 11 suites, 18 general CLI
assertions and both program workflows (11 Kang assertions, 10 Baron assertions).
Kang resident and streamed scores/detection match exactly by cell identity.
Baron's complete previous QC/pseudobulk/PCA report is exactly preserved after
excluding the newly added programs field. Native CSR/CSC fixture results are
exactly equal, including signed scores and missing/zero handling.

Independent maximum absolute score errors are 1.333e-15 (Kang) and 8.882e-16
(Baron). Both IFN programs increase in every one of eight paired donors; the
full Hallmark mean differences range from 0.6048 to 0.6622. Beta-program mean
differences versus other author-labelled cells range from 0.4011 to 0.5291
across four donors, with descriptive author-label AUCs 0.8629–0.9483. These
within-study observations do not establish prospective annotation accuracy.

[Evidence](evidence/2026-09-09/source-state.json) records source and binary hashes,
commands, native receipts, reference versions, rejections and the initial Swift
test-fixture compilation failure. The staging-cleanup check was corrected to
inspect the owner's actual `.numivivo-stream-*` prefix; both completed output
directories were then checked directly and had no residual staging.
