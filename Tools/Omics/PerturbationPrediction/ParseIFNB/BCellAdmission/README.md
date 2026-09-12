# Parse B-cell metadata and model-feature admission

**72,446** source rows carry the literal labels `B Intermediate/Memory` or
`B Naive` within the original 725,031-cell PBS/IFN-beta cohort. All twelve
donors have both conditions. Selection was frozen before reading cell-type codes;
it uses metadata, not response magnitude. Plasmablast remains a separate label.
Source labels are not authoritative biological identities.

The newer context model has **11,600 exact correspondences and 200 absent
features** among its 11,800 genes. The full coordinate mapping and missing list
are retained. No alias mapping, zero filling, merging or panel change was applied.
The earlier 409 missing duration-model features describe a different, unchanged
panel. This metadata audit did not aggregate counts or fit/score a model.
The subsequent [complete count workflow](../BCellCounts/README.md) now has all
twelve ingestions and replays verified; predictive admission remains unresolved.

The audit preserves label codes for all 725,031 original rows, all 432
donor/condition/label counts and exact original B-cell source row coordinates.
The smallest donor/condition B-cell group has 344 cells. This does not replace
the complete PBMC cohort or establish independent biological replication from
cell counts.

## Remaining experimental and model contract

The [context-panel correspondence audit](context-feature-correspondence.json)
uses the retained HGNC approved/previous-symbol table to investigate all 200
missing names. It finds **192 unique candidate correspondences** after checking
collisions against all 11,600 exact matches and the other candidates. A separate
direct table scan verifies every candidate identity and source coordinate.
These are nomenclature candidates, not proven sequence/annotation equivalence;
the existing exact-name admission gate and model panel remain unchanged.

Eight entries remain unresolved: CARD17, CTAG2, HBQ1, NPPC and ZNF781 have no
corresponding source symbol; QARS has ambiguous HGNC identity; TIAF1 conflicts
with another panel coordinate; VARS has ambiguous source correspondence.
Resolving names therefore cannot by itself admit the complete current model.
No generic aliases, fuzzy matches, count values or prediction outcomes were used.
[The script](context_feature_correspondence.py) retains original local paths and
input hashes; reproduction requires the existing metadata and HGNC dependencies.

The [primary preprint](https://pmc.ncbi.nlm.nih.gov/articles/PMC12724453/) links
its [supplementary workbook](https://pmc.ncbi.nlm.nih.gov/articles/instance/12724453/bin/media-2.xlsx).
The workbook download presented a browser challenge and was not bypassed.
Research Square returned 403. The IFN-beta dose/reagent remains unresolved;
24-hour exposure alone does not supply those details.

Resolve the experimental definition, establish the supported feature contract
before fitting and verify counts for the declared population before scoring.
Do not silently pass the current model an incomplete axis. Parse counts and DE
were previously inspected; a later prediction test must disclose that history.
No success or failure of the proposed prediction test is claimed here.

## Evidence and reproduction

The metadata reader checks HTTP 206, Content-Range and the original ETag on
every request, retaining SHA256 digests. No X dataset was read or count value interpreted. Prior donor/
condition arrays and roster match the published preparation manifest; the
model feature metadata match the published context-kernel archive.

An independent scalar checker reconstructs every original row decision, all
432 group/label counts and all 11,800 feature decisions. Run
`python verify_archive.py` to check archive integrity, then unpack into a fresh
directory and run `verify.py` with NumPy for the complete scalar metadata check.
The archive includes source roster/codes, selected labels, model metadata,
freeze, mapping, request hashes and scripts. Fresh online execution with
`admit.py` requires adjusting its original study/repo paths and using the
retained category file and version-pinned reader. Raw HTTP bytes are hash-bound,
not embedded; the source H5AD/count matrix is not included.

The first scalar checker repeatedly decompressed NPZ members inside the row
loop and was interrupted. Loading the two arrays once resolved the inefficiency;
the comparisons and row set are unchanged, and the failed source is retained.
Packaging initially failed with ENOSPC; verified published duplicate model
inputs were reclaimed before retrying. No research evidence was deleted.

Original Parse metadata: **Parse Biosciences, CC BY-NC 4.0**.
The software license does not relicense source data.
