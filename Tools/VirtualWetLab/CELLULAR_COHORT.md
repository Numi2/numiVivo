# Source-bound cellular onboarding

This delivery admits measured Replogle–Nadig data into the **existing native count
stores and composite corpus**, then exposes bounded measured reads in NumiLab.
It does not train or promote a predictor. The frozen cellular and spatial
experiment adapters remain selectable, with their historical replay semantics.

## Scientific boundary

The engineering target vocabulary is selected by source metadata: at least 32
assigned cells in at least two contexts, at least 512 source controls, then a
fixed SHA-256 ordering capped at 256 targets. Each supported source contributes
32 deterministic cells per target and 512 distinct controls. Counts are an
engineering subset size, not a biological sample size. Exclusions remain in the
registration. The same preparation accepts `targetLimit: null` and cell limits
`null` for full admission; that larger import must be separately registered and
have sufficient storage. It is not implied by engineering-subset qualification.

Original source identifiers, barcodes, GEM captures, guide pairs, target labels,
CRISPRi modality and endpoint day are retained. GEM captures are **technical
blocks**, not independent animals or biological replicates. Biological units are
unresolved. Native historical sample metadata requires an identifier; its
`unresolved:SOURCE` value is an explicit sentinel, not an identified replicate.
The historical corpus `guideID` field requires the target label; exact guide
pairs remain in the hash-bound row sidecar and must be used for guide audits.

Three distinct development tasks are prepared: known-target capture transfer,
held-target transfer across all sources, and held-context transfer to RPE1.
Their assignments precede expression processing. Browser access is stricter than
any individual task: only treated rows assigned to training in **all** tasks
are accessible. Controls are explicitly open reference inputs. Held responses
remain reserved-development; this is not a newly independent biological cohort.

Raw integer counts are preserved. The union of source feature IDs is the corpus
axis; native source masks mark absent features. The original deposits themselves
filter genes by source-wide expression. We retain that panel limitation; no new
HVG selection, representation fitting or checkpoint selection happens here.
Read-time normalization is `log1p(count / retained-panel cell total * 10000)`.
It must not be labeled complete-transcriptome CPM or confused with raw UMI counts.

## Reproduce admission

Use `examples/replogle_nadig_admission.json` with
`prepare_cellular_cohort.py metadata --manifest MANIFEST --output METADATA`.
Metadata can also be extracted from an already downloaded source using
`--source-id ID --source FILE`. Metadata inspection does not read expression.

Freeze the exact metadata and rules before processing counts:

```
python prepare_cellular_cohort.py freeze --manifest MANIFEST --metadata METADATA --output REGISTRATION
python prepare_cellular_cohort.py fetch --manifest MANIFEST --source-id ID --output SOURCE.h5ad
python prepare_cellular_cohort.py admit --registration REGISTRATION --metadata METADATA --source-id ID --source SOURCE.h5ad --output COHORT --binary NUMIVIVO
```

Repeat fetch/admit for each source. Downloaded staging files can be removed after
the verified admission; source checksums, metadata, exact projected raw rows and
native receipts remain. Fetch verifies full-file published digests. The Arc
mirror has misleading `normalized` filenames for the Replogle files: accept
those only after matching the original raw-file MD5 and checking integer count
semantics. A filename is not normalization evidence.

The importer holds at most one dense source row and 256 sparse rows while writing
an H5AD projection, then invokes `singlecell-h5ad-store` and
`cell-response-prepare`. It does not construct a dense cohort tensor. The native
composite reader keeps separate source stores and source-qualified identities.
The constant descriptor used for admission is explicitly **not biological target
information**; replace it through the next training registration before fitting.

```
python package_cellular_cohort.py --registration REGISTRATION --cohort COHORT --binary NUMIVIVO --metadata METADATA
python qualify_cellular_cohort.py --config COHORT/assay.json --output QUALIFICATION.json
```

The native `*-sources.json` manifests are invocation location maps, not sealed
scientific identity. Regenerate them with `prepare_cellular_cohort.py sources --cohort COHORT --task
known-target --output NEW-LOCATION-MAP.json` after relocation; do not rewrite
native receipts or sealed experiments.

## Shared laboratory access

Register `COHORT/assay.json` in the existing service. Presentation is population,
geometry is absent, measured exploration is enabled, and prediction is disabled
until a compatible model is registered. A proposal retains the user's objective
and reports that missing model; it cannot silently execute an older artifact.

`numi wet-lab readout` uses the acknowledged shared selection. Gene distributions,
feature search, interventions, objective coverage and expression tiles are
bounded. Tiles allow at most 64 cells by 32 genes. Missing features return masks
and null values; reserved outcomes return no measured values. No exploration
operation authorizes observation reveal. First catalog metadata does not load
expression; readouts use cached native count-store indices and bounded rows.

Sources: [original Replogle deposit](https://plus.figshare.com/articles/dataset/20029387),
[Nadig GSE264667](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE264667),
[official State source-format example](https://github.com/ArcInstitute/state/tree/9bbfe78a434a55205e4de834e1ea99f85f7a3add).
No State model, code or weights run in this delivery. State Transition permission
and a matched learned comparator remain prerequisites for the later model study.
