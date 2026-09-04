#ifndef NUMIVIVO_MECHANISM_SYNTHESIS_H
#define NUMIVIVO_MECHANISM_SYNTHESIS_H

#include "NumiVivoCore/NumiVivoCore.h"

#ifdef __cplusplus
extern "C" {
#endif

#define NVIVO_MECHANISM_SYNTHESIS_ABI_VERSION 1

typedef enum NVivoMechanismSynthesisFlag {
    NVIVO_SYNTHESIS_REQUIRE_INDEPENDENT_SHUTDOWN = 1u << 0u,
    NVIVO_SYNTHESIS_REQUIRE_MONITOR = 1u << 1u,
    NVIVO_SYNTHESIS_REQUIRE_CONTEXT_INSULATION = 1u << 2u,
    NVIVO_SYNTHESIS_REQUIRE_RESOURCE_BUFFERING = 1u << 3u,
    NVIVO_SYNTHESIS_REQUIRE_DISTINCT_ORTHOGONALITY = 1u << 4u
} NVivoMechanismSynthesisFlag;

typedef struct NVivoMechanismSynthesisOptions {
    uint32_t struct_size;
    uint32_t abi_version;
    uint32_t flags;
    uint32_t maximum_solutions;
    uint64_t maximum_visited_nodes;
    uint32_t reserved[16];
} NVivoMechanismSynthesisOptions;

NVIVO_EXPORT void nvivo_default_mechanism_synthesis_options(
    NVivoMechanismSynthesisOptions *options);

/*
 * Synthesize abstract molecular mechanisms from two bounded UTF-8 JSON
 * documents. The problem document defines slots and constraints. The library
 * document defines sequence-free mechanism records and evidence envelopes.
 *
 * result_json contains valid Pareto candidates even when the search budget is
 * exhausted. diagnostics_json always contains parse, validation and search
 * findings. Both buffers are released with nvivo_buffer_release.
 */
NVIVO_EXPORT NVivoStatus nvivo_synthesize_mechanisms_json(
    const uint8_t *problem_json,
    size_t problem_size,
    const uint8_t *library_json,
    size_t library_size,
    const NVivoMechanismSynthesisOptions *options,
    NVivoByteBuffer *result_json,
    NVivoByteBuffer *diagnostics_json);

#ifdef __cplusplus
}
#endif

#endif
