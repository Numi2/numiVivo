#include <metal_stdlib>
using namespace metal;

constant uint NVIVO_MIGRATION_INVALID_NUMBER = 1u << 0u;
constant uint NVIVO_MIGRATION_NEGATIVE_COUNT = 1u << 1u;
constant uint NVIVO_MIGRATION_COUNT_OVERFLOW = 1u << 2u;
constant uint NVIVO_MIGRATION_FRACTIONAL_COUNT = 1u << 3u;
constant uint NVIVO_MIGRATION_FLOAT_PRECISION_LOSS = 1u << 4u;
constant uint NVIVO_MIGRATION_INVALID_INDEX = 0xffffffffu;

struct NVivoMigrationCommand {
    uint elementCount;
    uint reserved0;
    float integerTolerance;
    float reserved1;
};

struct NVivoMigrationStatus {
    atomic_uint flags;
    atomic_uint invalidElementCount;
    atomic_uint firstInvalidElement;
    atomic_uint maximumErrorBits;
};

kernel void nvivo_migration_clear_status(
    device NVivoMigrationStatus* status [[buffer(0)]],
    uint index [[thread_position_in_grid]]
) {
    if (index != 0) return;
    atomic_store_explicit(&status->flags, 0u, memory_order_relaxed);
    atomic_store_explicit(&status->invalidElementCount, 0u, memory_order_relaxed);
    atomic_store_explicit(&status->firstInvalidElement, NVIVO_MIGRATION_INVALID_INDEX, memory_order_relaxed);
    atomic_store_explicit(&status->maximumErrorBits, 0u, memory_order_relaxed);
}

inline void nvivo_record_migration_failure(
    device NVivoMigrationStatus* status,
    uint index,
    uint failure,
    float error
) {
    atomic_fetch_or_explicit(&status->flags, failure, memory_order_relaxed);
    atomic_fetch_add_explicit(&status->invalidElementCount, 1u, memory_order_relaxed);
    atomic_fetch_min_explicit(&status->firstInvalidElement, index, memory_order_relaxed);
    if (isfinite(error) && error >= 0.0f) {
        atomic_fetch_max_explicit(&status->maximumErrorBits, as_type<uint>(error), memory_order_relaxed);
    }
}

kernel void nvivo_migration_f32_to_u32(
    const device float* source [[buffer(0)]],
    device uint* destination [[buffer(1)]],
    constant NVivoMigrationCommand& command [[buffer(2)]],
    device NVivoMigrationStatus* status [[buffer(3)]],
    uint index [[thread_position_in_grid]]
) {
    if (index >= command.elementCount) return;
    const float value = source[index];
    if (!isfinite(value)) {
        destination[index] = 0u;
        nvivo_record_migration_failure(status, index, NVIVO_MIGRATION_INVALID_NUMBER, INFINITY);
        return;
    }
    if (value < 0.0f) {
        destination[index] = 0u;
        nvivo_record_migration_failure(status, index, NVIVO_MIGRATION_NEGATIVE_COUNT, -value);
        return;
    }
    if (value > 4294967040.0f) {
        destination[index] = 0xffffffffu;
        nvivo_record_migration_failure(
            status,
            index,
            NVIVO_MIGRATION_COUNT_OVERFLOW,
            value - 4294967040.0f
        );
        return;
    }

    const float rounded = rint(value);
    const float error = fabs(value - rounded);
    destination[index] = uint(rounded);
    if (error > command.integerTolerance) {
        nvivo_record_migration_failure(
            status,
            index,
            NVIVO_MIGRATION_FRACTIONAL_COUNT,
            error
        );
    }
}

kernel void nvivo_migration_u32_to_f32(
    const device uint* source [[buffer(0)]],
    device float* destination [[buffer(1)]],
    constant NVivoMigrationCommand& command [[buffer(2)]],
    device NVivoMigrationStatus* status [[buffer(3)]],
    uint index [[thread_position_in_grid]]
) {
    if (index >= command.elementCount) return;
    const uint value = source[index];
    destination[index] = float(value);
    // Every integer through 2^24 is represented exactly by binary32. Larger
    // counts require an explicit authority decision rather than silent rounding.
    if (value > 16777216u) {
        const float roundTrip = float(value);
        const float error = fabs(float(value) - roundTrip);
        nvivo_record_migration_failure(
            status,
            index,
            NVIVO_MIGRATION_FLOAT_PRECISION_LOSS,
            max(error, 1.0f)
        );
    }
}

kernel void nvivo_migration_validate_f32(
    const device float* state [[buffer(0)]],
    constant NVivoMigrationCommand& command [[buffer(1)]],
    device NVivoMigrationStatus* status [[buffer(2)]],
    uint index [[thread_position_in_grid]]
) {
    if (index >= command.elementCount) return;
    if (!isfinite(state[index])) {
        nvivo_record_migration_failure(
            status,
            index,
            NVIVO_MIGRATION_INVALID_NUMBER,
            INFINITY
        );
    }
}
