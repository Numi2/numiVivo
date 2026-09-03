#ifndef NUMIVIVO_CORE_H
#define NUMIVIVO_CORE_H

#include <stddef.h>
#include <stdint.h>

#if defined(_WIN32)
#define NVIVO_EXPORT __declspec(dllexport)
#else
#define NVIVO_EXPORT __attribute__((visibility("default")))
#endif

#ifdef __cplusplus
extern "C" {
#endif

#define NVIVO_VERSION_MAJOR 0
#define NVIVO_VERSION_MINOR 1
#define NVIVO_VERSION_PATCH 0
#define NVIVO_PROGRAM_PACK_VERSION 1

NVIVO_EXPORT const char *nvivo_version_string(void);
NVIVO_EXPORT uint32_t nvivo_program_pack_version(void);

#ifdef __cplusplus
}
#endif

#endif
