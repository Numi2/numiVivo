# GSE181897 donor and condition join

All 62 frozen context-model donors equal the complete set of retained source donors
with both B and C observations. The previously qualified primary author-code receipt
maps B to IFNB and C to control. The join includes 4,938 B cells across these paired
conditions, from the 15,272-cell retained B-cell source containing other conditions.
No new donor exclusion or treatment-role inference is introduced.

`mapping.json` retains every model ID, original experiment ID, condition, native sample
ID, batch ID, pooling code and cell count. It hashes the H5AD, author-code role receipt
and frozen model donor metadata. `map.py` recomputes the complete eligibility set and
requires exact equality to the 62 model donors before writing any mapping.

Pooling codes are preserved literally. Their mapping to primary GEO pool accessions
has not been established by this check, so accessions remain null and exposure
covariates are not admitted. The earlier 500 IU/mL / 9-hour protocol reference is not
silently assigned to each pool. This is metadata linkage, not a new raw-count replay,
reverification of donor expression means, participant-independence check or prediction
accuracy result. The model and its failed outcomes are unchanged.

Reproduction uses backed AnnData reads with the original paths in `map.py` and a
fresh output directory. Runtime evidence remains at
`/Users/n/numivivo-gse181897-exposure-mapping-20260912`.
