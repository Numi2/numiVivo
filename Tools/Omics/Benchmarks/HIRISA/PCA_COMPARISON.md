# Complete HIRISA PCA comparison gates

Frozen before native full-source PCA publication or score inspection. The
independent full-source fit has completed and its maximum relative residual is
1.339e-13. These numerical gates retain the existing smaller-cohort reference
thresholds; they are not selected from native HIRISA agreement outcomes.

Require exact source and artifact hashes, all 1,612,594 original row identities,
18,082 gene identities in source order, all QC library sizes/detected genes/mito
counts, and exact selected genes and mean bins. Compare all-feature means,
variances and log dispersions with rtol/atol 1e-9; normalized dispersions with
rtol/atol 2e-6 (the established Scanpy compatibility tolerance). Compare centers
with rtol/atol 1e-12. Undefined dispersion masks must match.

Require explained variances rtol 1e-7/atol 1e-9 and variance ratios rtol
1e-7/atol 1e-10; minimum loading subspace cosine greater than 1-1e-8;
Procrustes-aligned relative score Frobenius error below 1e-5 over every row.
Require score means at most 1e-10 in absolute value and score covariance equal
to the diagonal reported variance with rtol 1e-7/atol 1e-8. Independently compute
native loading residuals using the complete sparse reference operator and
require each at most the frozen native 1e-6 tolerance. Retain the native engine's
existing orthogonality gate and report the independently measured error.

Read score records in bounded row blocks, verify every stored row/component
coordinate and finite value, and compare all projections. Do not load the raw
cells-by-genes array or convert the native score matrix to a JSON list. Reference
and native execution occur on different hosts; their timings are observations,
not a controlled CPU/scverse speed comparison.

The existing complete native QC report contains exactly 345,933,053 quality
bytes and 149,506,727 metadata bytes. Its source-bound array hashes provide an
additional exact check on standalone publication. The 256 MiB standalone QC cap
must be repaired to 512 MiB across publisher, verifier and downstream snapshot
readers before this dataset can publish. Metadata and QC remain resident; this
admission repair does not establish fully streamed metadata execution.
