# Source-bound HIRISA graph-community export

The native annotation executable exported the retained 30 graph communities to
`obs/numivivo_graph_community` for all 1,612,594 HIRISA cells. AnnData reopened the
result with a disk-backed CSR matrix (18,082 genes). Every label, both dataframe
indices and all 29 original observation/feature columns passed comparison.
Six sampled count rows (14,931 stored entries) matched exactly; this run does
not claim a full count-payload comparison.

The preparer verified the frozen result and receipt chain through graph and PCA
to the exact original H5AD SHA256. Every barcode and sample pair matched in order.
The PCA plan specifies `geo_accession` as its sample column. An initial checker
incorrectly assumed `sampleID`, rejected the comparison and wrote no annotation;
that failed assumption is retained in the evidence note.

Native export took 5.56 seconds with 141,770,752 bytes maximum RSS on the local
Apple M4. This is one observed run, not a comparative speedup qualification.
The existing qualified annotation executable was used; no model was retrained.
These are numerical graph communities, not authoritative cell types or evidence
of general biological outcome prediction.

Output retained locally:
`/Users/home/numivivo-cluster-export-20260912/hirisa-clustered.h5ad`.
The new file uses separate-inode APFS clone publication. The original dataset,
frozen clustering and earlier annotation output remain available.

## Reproduction and evidence

Use the frozen manifest at
`Tools/Omics/Benchmarks/HIRISA/evidence/2026-09-10-full-clustering/manifest.json`
to retrieve the five bundle gzip files and `run/pca-plan.json.gz` into a study
folder. The retained `prepare.py` validates their hashes and creates `plan.json`.
Run the qualified Atlas annotation binary with `annotate SOURCE PLAN NEW_OUTPUT`,
then `check_backed.py --source SOURCE --annotated NEW_OUTPUT --report REPORT`.
The scripts retain this study's absolute source path; adapt that path explicitly
when reproducing elsewhere. Keep the frozen result under `bundle/result.json.gz`.

[Receipts and verification](evidence/2026-09-12) retain the checker sources,
source manifest, native receipt, measured resource log and AnnData comparisons.
The multi-gigabyte H5AD files are not copied into Git.
