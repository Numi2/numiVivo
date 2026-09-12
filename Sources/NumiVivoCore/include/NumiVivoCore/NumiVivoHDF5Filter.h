#ifndef NUMIVIVO_HDF5_FILTER_H
#define NUMIVIVO_HDF5_FILTER_H
#include <stddef.h>
#include <stdint.h>

/* H5Z_class2_t ABI for the dynamically loaded HDF5 filter registry. This
 * declaration requires no HDF5 headers or link-time HDF5 dependency. */
typedef struct NVivoHDF5FilterDescriptor {
    int version;
    int id;
    unsigned encoder_present;
    unsigned decoder_present;
    const char *name;
    int (*can_apply)(int64_t, int64_t, int64_t);
    int (*set_local)(int64_t, int64_t, int64_t);
    size_t (*filter)(unsigned, size_t, const unsigned *, size_t, size_t *, void **);
} NVivoHDF5FilterDescriptor;
#endif
