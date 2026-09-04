#ifndef NUMIVIVO_NETWORK_ANALYSIS_H
#define NUMIVIVO_NETWORK_ANALYSIS_H

#include "NumiVivoCore/NumiVivoCore.h"

#ifdef __cplusplus
extern "C" {
#endif

#define NVIVO_NETWORK_ANALYSIS_ABI_VERSION 1

typedef enum NVivoNetworkAnalysisFlag {
    NVIVO_NETWORK_INFER_CONSERVATION = 1u << 0u,
    NVIVO_NETWORK_COMPUTE_COMPLEX_DEFICIENCY = 1u << 1u
} NVivoNetworkAnalysisFlag;

typedef struct NVivoNetworkAnalysisOptions {
    uint32_t struct_size;
    uint32_t abi_version;
    uint32_t flags;
    uint32_t maximum_dense_species;
    uint32_t maximum_dense_reactions;
    uint32_t maximum_conservation_laws;
    uint32_t reserved0;
    uint64_t maximum_dense_elements;
    double rank_tolerance;
    double conservation_tolerance;
    uint32_t reserved[16];
} NVivoNetworkAnalysisOptions;

NVIVO_EXPORT void nvivo_default_network_analysis_options(
    NVivoNetworkAnalysisOptions *options);

/*
 * Compiles a VivoProgram JSON document to the same typed reaction IR used by
 * production execution, then performs structural network analysis. Analysis is
 * sequence-free and does not execute the model.
 */
NVIVO_EXPORT NVivoStatus nvivo_analyze_program_json(
    const uint8_t *source,
    size_t source_size,
    const NVivoCompileOptions *compile_options,
    const NVivoNetworkAnalysisOptions *analysis_options,
    NVivoByteBuffer *analysis_json,
    NVivoByteBuffer *diagnostics_json);

#ifdef __cplusplus
}
#endif

#endif
