#ifndef NUMIVIVO_OMICS_MNN_H
#define NUMIVIVO_OMICS_MNN_H
#include "NumiVivoCore/NumiVivoCore.h"
#ifdef __cplusplus
extern "C" {
#endif
/* Exact Manhattan lower/upper-level neighbor heaps on row-major normalized PCA.
 * Each cross-level pair is evaluated once. Ties use (level, original row).
 * All borrowed spans must be disjoint and valid for the synchronous call.
 * Capacities count elements; neighbor spans each require rows*neighbors elements,
 * count spans each require rows. No internal heap allocation. Outputs and partial
 * counters are not usable as a result on failure. Status: 0 success, 1 invalid
 * arguments/capacity/design, 2 scalar-term budget, 4 nonfinite, 5 cancellation.
 * Boundary levels have zero neighbors in the absent direction.
 */
NVIVO_EXPORT int32_t nvivo_omics_mnn_exact(
    const double *scores, const intptr_t *levels, uint32_t rows, uint32_t dimensions,
    uint32_t neighbors, uint32_t level_count, uint64_t score_capacity,
    uint64_t level_capacity, uint64_t maximum_scalar_terms,
    intptr_t *lower_ids, double *lower_distances, intptr_t *upper_ids,
    double *upper_distances, uint64_t neighbor_capacity,
    intptr_t *lower_counts, intptr_t *upper_counts, uint64_t count_capacity,
    uint64_t *evaluated_scalar_terms, int32_t (*cancel)(void *), void *context);
#ifdef __cplusplus
}
#endif
#endif
