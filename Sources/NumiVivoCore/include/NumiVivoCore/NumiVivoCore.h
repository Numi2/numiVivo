#ifndef NUMIVIVO_CORE_H
#define NUMIVIVO_CORE_H

#include <stddef.h>
#include <stdint.h>

#if defined(_WIN32)
  #if defined(NVIVO_BUILDING_CORE)
    #define NVIVO_EXPORT __declspec(dllexport)
  #else
    #define NVIVO_EXPORT __declspec(dllimport)
  #endif
#else
  #define NVIVO_EXPORT __attribute__((visibility("default")))
#endif

#ifdef __cplusplus
extern "C" {
#endif

#define NVIVO_VERSION_MAJOR 0
#define NVIVO_VERSION_MINOR 1
#define NVIVO_VERSION_PATCH 0
#define NVIVO_COMPILER_ABI_VERSION 1
#define NVIVO_PROGRAM_PACK_VERSION 1
#define NVIVO_FINGERPRINT_BYTES 32

typedef enum NVivoStatus {
    NVIVO_STATUS_OK = 0,
    NVIVO_STATUS_INVALID_ARGUMENT = 1,
    NVIVO_STATUS_PARSE_ERROR = 2,
    NVIVO_STATUS_VALIDATION_ERROR = 3,
    NVIVO_STATUS_SAFETY_REJECTED = 4,
    NVIVO_STATUS_COMPILE_ERROR = 5,
    NVIVO_STATUS_INVALID_PACK = 6,
    NVIVO_STATUS_RESOURCE_LIMIT = 7,
    NVIVO_STATUS_OUT_OF_MEMORY = 8,
    NVIVO_STATUS_INTERNAL_ERROR = 255
} NVivoStatus;

typedef enum NVivoFidelityLevel {
    NVIVO_FIDELITY_F0_LOGIC = 0,
    NVIVO_FIDELITY_F1_DETERMINISTIC = 1,
    NVIVO_FIDELITY_F2_STOCHASTIC = 2,
    NVIVO_FIDELITY_F3_SPATIAL = 3,
    NVIVO_FIDELITY_F4_TISSUE = 4
} NVivoFidelityLevel;

typedef enum NVivoCompileFlag {
    NVIVO_COMPILE_STRICT_UNITS = 1u << 0u,
    NVIVO_COMPILE_STRICT_SAFETY = 1u << 1u,
    NVIVO_COMPILE_DETERMINISTIC_PACK = 1u << 2u,
    NVIVO_COMPILE_PERMIT_HYPOTHETICAL_PARAMETERS = 1u << 3u,
    NVIVO_COMPILE_REQUIRE_TERMINATION = 1u << 4u
} NVivoCompileFlag;

typedef struct NVivoByteBuffer {
    uint8_t *data;
    size_t size;
} NVivoByteBuffer;

typedef struct NVivoResourceLimits {
    uint32_t maximum_species;
    uint32_t maximum_parameters;
    uint32_t maximum_reactions;
    uint32_t maximum_rules;
    uint32_t maximum_constraints;
    uint32_t maximum_expression_instructions;
    uint32_t maximum_stoichiometry_terms;
    uint32_t maximum_temporal_states;
} NVivoResourceLimits;

typedef struct NVivoCompileOptions {
    uint32_t struct_size;
    uint32_t abi_version;
    uint32_t flags;
    uint32_t requested_fidelity;
    NVivoResourceLimits limits;
    uint32_t reserved[16];
} NVivoCompileOptions;

typedef struct NVivoPackSummary {
    uint32_t struct_size;
    uint32_t pack_major;
    uint32_t pack_minor;
    uint32_t compiler_abi;
    uint32_t fidelity;
    uint32_t section_count;
    uint64_t total_bytes;
    uint8_t source_fingerprint[NVIVO_FINGERPRINT_BYTES];
    uint8_t content_fingerprint[NVIVO_FINGERPRINT_BYTES];
    uint32_t reserved[16];
} NVivoPackSummary;

NVIVO_EXPORT const char *nvivo_version_string(void);
NVIVO_EXPORT uint32_t nvivo_program_pack_version(void);
NVIVO_EXPORT void nvivo_default_compile_options(NVivoCompileOptions *options);

/*
 * Compiles one UTF-8 VivoProgram JSON document into an immutable ProgramPack.
 * Both output buffers are owned by the caller and must be released with
 * nvivo_buffer_release. diagnostics_json is populated for warnings as well as
 * errors, allowing successful compilation to retain all findings.
 */
NVIVO_EXPORT NVivoStatus nvivo_compile_program_json(
    const uint8_t *source,
    size_t source_size,
    const NVivoCompileOptions *options,
    NVivoByteBuffer *program_pack,
    NVivoByteBuffer *diagnostics_json);

/* Parses, decodes, validates, and performs static safety analysis without
 * serializing a ProgramPack. */
NVIVO_EXPORT NVivoStatus nvivo_validate_program_json(
    const uint8_t *source,
    size_t source_size,
    const NVivoCompileOptions *options,
    NVivoByteBuffer *diagnostics_json);

/* Verifies a ProgramPack header, section bounds, alignment, overlap, and
 * fingerprints. summary may be NULL. */
NVIVO_EXPORT NVivoStatus nvivo_inspect_program_pack(
    const uint8_t *pack,
    size_t pack_size,
    uint32_t verify_section_hashes,
    NVivoPackSummary *summary,
    NVivoByteBuffer *inspection_json);

/* Releases buffers returned by this ABI. Passing NULL is valid. */
NVIVO_EXPORT void nvivo_buffer_release(NVivoByteBuffer *buffer);

#ifdef __cplusplus
}
#endif

#endif
