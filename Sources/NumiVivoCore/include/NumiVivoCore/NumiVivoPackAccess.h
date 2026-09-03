#ifndef NUMIVIVO_PACK_ACCESS_H
#define NUMIVIVO_PACK_ACCESS_H

#include "NumiVivoCore.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum NVivoPackSectionType {
    NVIVO_PACK_SECTION_STRINGS = 1,
    NVIVO_PACK_SECTION_SPECIES = 2,
    NVIVO_PACK_SECTION_PARAMETERS = 3,
    NVIVO_PACK_SECTION_REACTION_PARAMETER_INDICES = 4,
    NVIVO_PACK_SECTION_STOICHIOMETRY = 5,
    NVIVO_PACK_SECTION_REACTIONS = 6,
    NVIVO_PACK_SECTION_EXPRESSIONS = 7,
    NVIVO_PACK_SECTION_ACTIONS = 8,
    NVIVO_PACK_SECTION_RULES = 9,
    NVIVO_PACK_SECTION_MONITORS = 10,
    NVIVO_PACK_SECTION_COHORTS = 11,
    NVIVO_PACK_SECTION_SPECIES_INCIDENCE_OFFSETS = 12,
    NVIVO_PACK_SECTION_SPECIES_INCIDENCE = 13,
    NVIVO_PACK_SECTION_RUNTIME_CONTRACT = 14
} NVivoPackSectionType;

typedef struct NVivoPackSectionView {
    const uint8_t *data;
    size_t size;
    uint32_t stride;
    uint32_t count;
    uint32_t flags;
    uint32_t alignment;
} NVivoPackSectionView;

typedef struct NVivoGPUParameterRecord {
    float value;
    float minimum;
    float maximum;
    uint32_t flags;
} NVivoGPUParameterRecord;

/*
 * Returns a borrowed section view into `pack`. The caller must keep the pack
 * bytes alive and immutable while the view is used. The complete pack is
 * structurally and cryptographically validated before a view is returned.
 */
NVIVO_EXPORT NVivoStatus nvivo_program_pack_section(
    const uint8_t *pack,
    size_t pack_size,
    uint32_t section_type,
    NVivoPackSectionView *view);

/*
 * Converts authoritative FP64 parameter records into the exact FP32 hot-table
 * layout consumed by Metal kernels. The returned buffer is owned by the caller
 * and released with nvivo_buffer_release.
 */
NVIVO_EXPORT NVivoStatus nvivo_materialize_gpu_parameters(
    const uint8_t *pack,
    size_t pack_size,
    NVivoByteBuffer *gpu_parameter_records);

#ifdef __cplusplus
}
#endif

#endif
