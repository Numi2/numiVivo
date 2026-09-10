# Adamson experimental roles: metadata audit

**The independent prediction experiment is prepared but cannot yet be scored
under its frozen protocol.** Original GEO metadata independently reconstructs
all 50,440 selected cells and 94 guide groups. Neither expression counts nor
plausible guide names establish the missing experimental control assignments.
The audit also identifies a roster discrepancy requiring primary evidence.

## What was verified

`audit_experimental_roles.py` reads the previously archived original GEO barcodes
and guide records, verifies stored and decoded hashes, and reconstructs the
author's assignment rule. Every selected original row, barcode, excluded row,
guide group and per-guide cell count matches the frozen cohort. It reads no
expression matrix, fits no model and changes no cohort or annotation.

| Inventory from restored GEO records | Complete source cells | Selected cells |
| --- | ---: | ---: |
| `63(mod)_pBA580`, role unverified | 6,284 | 4,595 |
| `Gal4-4(mod)_pBA582`, role unverified | 1,343 | 646 |
| `62(mod)_pBA581`, role unverified | 2 | 0 |

These source counts use restored full-barcode assignments. They differ from the
historical scPerturb-label counts in the original ingestion report; the
[restoration audit](GEO_RESTORATION.md) explains that discrepancy.

The selected cohort contains **92 gene-like guide groups across 82 candidate
target prefixes**, plus the two non-gene labels above. The gene-like groups
contain 45,199 selected cells. Eighty candidate targets have frozen GO support.
All guide-to-target assignments remain candidate mappings until the experimental
roster is verified; an exact match to an expression feature is not that proof.
Current GO annotations also remain potentially informed by the study itself.

## Roster discrepancy

[Figure S1E of the original paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC5315571/)
describes 93 sgRNAs including two controls in the UPR experiment. The observed
selected inventory has 94 guide groups. Matching the two non-gene labels to the
reported number of controls would still leave this difference unresolved.

This is a discrepancy between a published experiment summary and deposited
selected labels, not proof that a particular guide is erroneous. No guide is
removed, merged or assigned a role to force agreement. The frozen cohort,
annotations and [prediction protocol](PROTOCOL.md) remain unchanged.

The next required source is Table S1 (protospacer sequences, including its
experiment-specific annotations), or original-author UPR analysis records
explicitly mapping these construct IDs to controls and targets. Reconcile each
selected guide with that roster and retain any unresolved or excluded identities
with reasons. The separate author demo's three-guide epistasis controls do not
establish controls for this pooled UPR screen.

## Retrieval status and next execution step

On 2026-09-10, publisher supplement requests returned HTTP 403; the Princeton
record returned HTTP 401; browser retrieval of PMC and eScholarship required
human verification. No access challenge was bypassed. The previously retrieved
164-byte `supplementary.zip` is an XML error saying supplementary access is
unavailable, and `author-manuscript.pdf` is empty. Neither is a scientific source
file. The indexed primary paper supports the experiment-level count above, but
its guide table was not retrieved or inspected.

Public GEO supplement directory listings succeeded and contain barcodes, guide
assignments, gene identities and the count matrix; they expose no separate
control/construct roster. The existing author repository inventory leads to the
already inspected demo and the later Norman study, without supplying the needed
UPR mapping. The [retrieval record and complete guide inventory](evidence/2026-09-10-experimental-roles/manifest.json)
preserve these results so subsequent work can seek missing evidence directly.

After verification, freeze the exact roles and all target folds, then execute
the existing native target kernel with lambda 1 and its unchanged baselines.
Freeze predictions before reading held-out outcomes. This will test the fixed
algorithm in an independent study; it will not test transfer of Norman-trained
coefficients, unseen tissues or clinical outcomes. No Adamson prediction result
is claimed by this audit.

## Reproduce the metadata audit

From the repository root, using Python's standard library:

```sh
python3 Tools/Omics/PerturbationPrediction/Adamson/audit_experimental_roles.py --out /tmp/adamson-roles.json
```

The destination must be new. The archived report retains every original guide,
including guides with zero selected cells, and explicitly leaves
`controlsVerified`, `experimentalTargetsVerified` and `readyForFrozenPrediction`
false. A passing metadata reconstruction does not clear those scientific gates.
