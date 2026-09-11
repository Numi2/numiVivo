# Scalar buffer-access repair

The isolated callback helper did not restore the intermediate scalar timing.
Extract the original scalar Gaussian sum into a non-inlined helper over borrowed,
already-validated buffers. Preserve coordinate order, anchor order, exp, numerator
and denominator arithmetic. Do not substitute a different default or call tiled
arithmetic from the scalar branch. Require every original scalar score/anchor byte
and report to remain exact on all three full cohorts. Repeat the final tiled,
replay, default-thread and Hagai lifecycle checks and retain intermediate evidence.


The Swift borrowed-buffer helper still increased scalar timing on Hagai. Move
only that unchanged scalar summation into the native C++ owner, with contraction
and reassociation disabled, finite-value validation and per-256-anchor cancellation.
Keep the existing Swift division, coordinate update, work charge and report.
The fixed full-cohort exact-byte gate remains mandatory; retain the unsuccessful
Swift extraction. This is implementation repair, not a changed numerical model.


The production-restored library provides the authoritative scalar timing baseline.
Amortize source/bias finite validation once per complete scalar assembly phase,
rather than repeating it for every query. Borrow the selected-row list and mutable
local score buffer; use one 64-double stack accumulator. Keep scalar summation,
query order, division, coordinate update and RMS accumulation unchanged. Capacity,
selected-index, work and finite checks remain explicit; failures never publish
partial results. Repeat the full protocol against the preserved production build.
