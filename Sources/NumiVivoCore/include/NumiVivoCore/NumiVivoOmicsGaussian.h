#ifndef NUMIVIVO_OMICS_GAUSSIAN_H
#define NUMIVIVO_OMICS_GAUSSIAN_H
#include "NumiVivoCore/NumiVivoCore.h"
#ifdef __cplusplus
extern "C" {
#endif
/* All-anchor FP64 Gaussian bias average, fixed 32-query / 256-anchor tiles.
 * Input/output/workspace buffers must be disjoint, contiguous and remain valid
 * for the synchronous call. Source/bias are anchor-major, queries/delta row-major.
 * Duplicate anchors are retained. All arrays are borrowed; no full matrix copy.
 * Workspace requirement is (256*dimensions + 256*query_rows) doubles.
 * Status 0 success; 1 arguments/capacity; 2 work; 4 nonfinite; 5 cancelled.
 * Outputs are not a result on failure. Zero-weight rows receive zero delta.
 * Deterministic on the pinned executable/OS; Accelerate reduction differs from
 * scalar order. Totals below 1e-280 use an explicitly budgeted scalar fallback.
 * evaluated_scalar_terms reports direct distance terms, including fallback. Framework internal workspace is outside the explicit buffer size.
 */
NVIVO_EXPORT int32_t nvivo_omics_gaussian_bias(
    const double *source, const double *bias, uint32_t anchors,
    const double *queries, uint32_t query_rows, uint32_t dimensions, double sigma,
    uint64_t source_capacity, uint64_t query_capacity, uint64_t maximum_scalar_terms,
    double *delta, uint64_t delta_capacity, double *totals, uint64_t totals_capacity,
    double *workspace, uint64_t workspace_capacity, uint64_t *evaluated_scalar_terms,
    int32_t (*cancel)(void *), void *context);
/* Legacy scalar summation order, one query, unnormalized numerator and total.
 * Same disjoint borrowed-buffer and failure-output contract. Status 0 success,
 * 1 arguments/capacity, 4 nonfinite, 5 cancelled. The Swift owner admits anchors*d
 * terms before entry. No workspace allocation; cancellation every 256 anchors. */
NVIVO_EXPORT int32_t nvivo_omics_gaussian_scalar(
    const double *source, const double *bias, uint32_t anchors, const double *query,
    uint32_t dimensions, double sigma, uint64_t source_capacity, uint64_t query_capacity,
    double *numerator, uint64_t output_capacity, double *total,
    int32_t (*cancel)(void *), void *context);
typedef struct NVivoGaussianScalarReport {
    uint64_t zero_weight_rows;
    double minimum_weight, maximum_weight, sum_squared;
} NVivoGaussianScalarReport;
/* Complete selected-row scalar phase. Immutable anchor snapshots are validated
 * once. Scores update in selected order, using the same per-query scalar sums.
 * All spans must be disjoint except the intentionally mutable scores span.
 * Maximum work counts anchors*dimensions*selected_rows. No heap allocation;
 * one 64-double stack accumulator. Partial scores/report are not a result on
 * failure. Status 0 success, 1 arguments, 2 work, 4 nonfinite, 5 cancellation. */
NVIVO_EXPORT int32_t nvivo_omics_gaussian_scalar_apply(
    const double *source, const double *bias, uint32_t anchors, uint32_t dimensions,
    double sigma, uint64_t source_capacity, double *scores, uint32_t rows, uint64_t score_capacity,
    const intptr_t *selected, uint32_t selected_rows, uint64_t maximum_scalar_terms,
    NVivoGaussianScalarReport *report, int32_t (*cancel)(void *), void *context);
#ifdef __cplusplus
}
#endif
#endif
