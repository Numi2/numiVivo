# All-anchor Gaussian kernel qualification protocol

Declared before execution on 2026-09-11, based on 9895348d5f14266a068b1626b04ece86e44c859d.

Preserve the exact exhaustive Manhattan matcher, grouped source ties, original
20 PCs, all cells, metadata, normalization, anchor duplicates, alignment order,
visit limits and Gaussian bandwidth of the original Hagai/Kang/Ding runs.
Replace only all-anchor Gaussian evaluation with fixed 32-query/256-anchor tiles.
Direct squared differences accumulate in coordinate order, with FP contraction
and reassociation disabled. Accelerate vvexp evaluates weights; BLAS multiplies
all weights by all corresponding biases. No distance cutoff or anchor truncation.
Underflowed zero total weight preserves the query row. Temporary storage is bounded
and charged to admission; cancellation is checked at each anchor tile.

Keep the old scalar default and encoding for historical artifact reproduction.
The new explicit kernel option identifies the changed floating point reductions.
No tuning of tiles, sigma, anchors or biological margins after inspecting results.

First qualify native small cases against independent scalar/NumPy sums: dimension
1/2/20/64, tile boundaries, duplicate anchors, large shared offsets, underflow,
nonfinite input, capacity/workspace rejection and cancellation. Scalar sum versus
tiled maximum delta error must be <=1e-10*(1+maximum absolute reference delta).
Run native tests and ASan/UBSan on the bridge.

Then execute one scalar and two tiled runs on all original 13,863 Hagai, 24,673
Kang and 44,031 Ding cells with unchanged inputs. Every anchor byte, alignment,
assembly order and panorama must match. Every final coordinate must differ by
<=1e-10*(1+maximum absolute original score); tiled scores/reports replay exactly.
Compare all step metrics with the same 1e-10 relative-plus-absolute tolerance and
require identical zero-weight counts. Check original artifact identities before
and after execution. Record actual owner timing, end-to-end command time, process
peak RSS, machine/compiler/runtime hashes and competing baseline ownership.

Freeze complete output hashes before applying the unchanged complete-cohort
biological evaluators and original margins. Retain every failure and unavailable
stratum. Record differences versus original exact MNN, not just gates versus
uncorrected PCA. This inspected-cohort test is not independent biological or
prospective outcome validation. No million-cell, out-of-core, Metal or competitive
end-to-end speed claim follows from these bounded native-owner timings.
