# Native tensor annotation writes

Native annotation now writes dense arrays through rank eight, matching the existing projection rank boundary. Numeric `obsm` and `varm` tensors preserve the first cell/feature axis and permit at most 256 components per row across all trailing dimensions. `uns` supports numeric, boolean and string tensors. Sparse matrices remain rank two; dataframe columns remain vectors. The two-million-element plan allowance, dimension limits, source fingerprint, atomic publication and typed-value checks remain enforced.

Previously, the annotation writer rejected the same input with `H5AD array rank or element allowance`, even though projection could preserve tensors. The retained `before` log demonstrates this failure. This does not add ragged arrays or arbitrary nested compound records.

The actual Omics library and product CLI were built on the physical Mac mini with `Tools/Omics/H5AD/build.sh --with-cli`. The binary SHA-256 is `d2684279f5dd741cf06bfaa947fb700ea34762ab6ffd071df340b419c1a99e6d`. This is the scoped Omics/product-router build, not a full package qualification.

The [checker](../check_tensor_annotations.py) passed eight command checks on both a software fixture and the complete Kang dataset: 24,673 cells by 15,706 genes. It verifies numeric types and tensor values with AnnData 0.13.3.post0, preservation of expression and obs/var, repeated-axis projection, reconstruction, empty dimensions and five invalid-input rejection cases. The existing annotation suite additionally passed all 26 checks, including unchanged fields, precision, provenance, atomic failure and unsupported storage dependencies.

Kang is read with AnnData and written as a modern source before native annotation. The upstream source hash is `1d48c1ff10bcfad7c1bc75863ce702a491c6c71b08e921a5e9754a986a9d19f3`; this rewrite is explicit and does not claim byte equality to the upstream file. Native edits leave the rewritten input bytes unchanged. Tensor values are deterministic software controls, not biological predictions. No dense cell-by-gene matrix is constructed.

`evidence.tar.gz` contains reports, plans, logs, checker, build source hashes and a manifest binding retained artifacts. Full H5AD artifacts and the executable remain under `/Users/n/numivivo-h5ad-tensors-20260912`. The original rejection and build warnings are retained.
