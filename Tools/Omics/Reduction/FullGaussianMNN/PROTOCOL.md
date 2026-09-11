# Eligible-level HNSW matching with full Gaussian correction

Declared before this trial's fitting and biological evaluation on 2026-09-11.
This isolates approximate anchor matching from the earlier combined HNSW/local64
failure. It is development on inspected cohorts, not independent validation.

Keep every original Hagai, Kang and Ding cell and original twenty PCs, donor/batch
eligibility, k=20, sigma=15, alignment threshold 0.1, global median norm, tie/order,
panorama/visit rules, merge multiplicity and duplicate anchors. No cell labels,
RNA programs or outcomes enter fit. No subsampling or parameter search.

Reuse the exact eligible-prefix/suffix native HNSW implementation and fixed
hnswlib 0.8.0 M=16, efConstruction=200, efSearch=128, seed=7, FP64 Manhattan
matching on normalized original PCs. Each index has the existing 500-million
metric-evaluation budget. Record every selected neighbor/distance and mutual pair.

Replace only local64 smoothing with the qualified production all-anchor Gaussian
C ABI from c2ce802bf8a77080339a6b46c3b3b12915199b60: direct FP64 squared distances,
32-query/256-anchor tiles, all anchor records retained, vForce exp/BLAS weighted
sum, original scalar fallback below 1e-280 total. Admit at most twice the exact
anchors*queries*dimensions distance work per correction, charging fallback terms.
Retain every full step source/bias/query/delta/total snapshot. No dense cells by
genes array, no truncated Gaussian, no approximate-kernel claim.

Freeze source/library/original input identities before fitting and all outputs
before biological evaluation. Replay every complete output/witness exactly.
Independently recompute every matching distance and all-anchor Gaussian update
from complete step snapshots with NumPy direct distances, exponentials and sums.
Require per-coordinate delta error <= 1e-10*(1+abs(reference)); reconstruct full
final coordinates and report all step and original-coordinate differences.
Compare every pair with the exact original anchors; require global precision and
recall >=0.95, preserve/report changed assembly order. Also compare candidate
matching and anchors with the previous eligible-level run to verify the isolated
intervention. A failure is retained, not repaired by changing parameters.

Use the unchanged complete-cohort 30-neighbor classifiers, every-type recall,
condition/program/within-stratum program margins and fixed original evaluation
inputs. Keep four missing Kang strata and partial Ding source labels. Record
fit/replay time separately from witness I/O and full process wall/RSS. The
physical M4 Pro original mapped-reader baseline remains active. This research
assembly driver is not native Swift public API, a prospective transform,
end-to-end speed comparison, million-cell integration or biological qualification.
Favorable results justify a separately qualified native matching implementation;
no production promotion from this trial alone.
