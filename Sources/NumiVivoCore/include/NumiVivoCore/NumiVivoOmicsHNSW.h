#ifndef NUMIVIVO_OMICS_HNSW_H
#define NUMIVIVO_OMICS_HNSW_H
#include "NumiVivoCore/NumiVivoCore.h"
#ifdef __cplusplus
extern "C" {
#endif

typedef struct NVivoHNSWOptions {
    uint32_t struct_size, abi_version;
    uint32_t rows, dimensions, neighbors, connections;
    uint32_t ef_construction, ef_search;
    uint64_t seed, maximum_distance_evaluations, score_cache_bytes;
} NVivoHNSWOptions;
typedef struct NVivoHNSWReport {
    uint64_t construction_distances, query_distances;
    uint64_t score_read_bytes, score_read_calls, score_cache_hits;
    uint64_t peak_cached_score_bytes, index_storage_bytes;
} NVivoHNSWReport;
typedef int32_t (*NVivoHNSWCancel)(void *context);
/* Serial, deterministic within the pinned executable/OS. Input is immutable
 * complete row-major u32-row/u32-component/f64-LE records. Output is self first,
 * then approximate neighbors in distance/index order, with Euclidean distances.
 * Status: 0 success, 1 arguments, 2 work budget, 3 score I/O/format, 4 nonfinite
 * distance, 5 cancelled, 6 internal index failure, 7 allocation failure.
 * Capacity is in entries, not bytes. Arrays are not a result on failure. */
NVIVO_EXPORT int32_t nvivo_omics_hnsw_neighbors(
    const char *score_path, const NVivoHNSWOptions *options,
    uint32_t *indices, double *distances, uint64_t output_capacity,
    NVivoHNSWReport *report, NVivoHNSWCancel cancel, void *cancel_context);
/* Row pointers are valid only during the callback. A nonzero callback result
 * stops the operation with status 8. No partial output is a published result. */
typedef int32_t (*NVivoHNSWRow)(uint32_t row, const uint32_t *indices,
    const double *distances, uint32_t count, void *context);
NVIVO_EXPORT int32_t nvivo_omics_hnsw_neighbors_stream(
    const char *score_path, const NVivoHNSWOptions *options, NVivoHNSWRow sink,
    NVivoHNSWReport *report, NVivoHNSWCancel cancel, void *context);
#ifdef __cplusplus
}
#endif
#endif
