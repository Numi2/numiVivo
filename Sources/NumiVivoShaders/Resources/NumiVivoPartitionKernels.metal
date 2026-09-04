#include <metal_stdlib>
using namespace metal;

namespace numivivo_partition {

constant uint kStatusNeedsSmallerStep = 1u << 0u;
constant uint kStatusNonFinite = 1u << 1u;
constant uint kStatusNegative = 1u << 2u;
constant uint kStatusAboveBound = 1u << 3u;
constant uint kInvalidIndex = 0xffffffffu;

struct CompartmentRecord {
    float volume;
    float inverseVolume;
    float reserved0;
    float reserved1;
};

struct AnalyteRecord {
    float minimumConcentration;
    float maximumConcentration;
    uint reserved0;
    uint reserved1;
};

struct EdgeRecord {
    uint analyteIndex;
    uint sourceCompartment;
    uint targetCompartment;
    uint flags;
    float partitionCoefficient;
    float clearance;
    float sourceUnboundFraction;
    float targetUnboundFraction;
};

struct RuntimeCommand {
    uint compartmentCount;
    uint analyteCount;
    uint edgeCount;
    uint stepIndex;
    float deltaTime;
    float absoluteTime;
    float negativeTolerance;
    uint attemptIndex;
};

struct RuntimeStatus {
    atomic_uint flags;
    atomic_uint invalidElementCount;
    atomic_uint firstInvalidElement;
    atomic_uint maximumConcentrationBits;
    atomic_uint maximumNegativeMagnitudeBits;
    atomic_uint reserved0;
    atomic_uint reserved1;
    atomic_uint reserved2;
};

inline float derivative(
    const device float* state,
    const device CompartmentRecord* compartments,
    const device EdgeRecord* edges,
    uint analyte,
    uint compartment,
    constant RuntimeCommand& command
) {
    float rate = 0.0f;
    for (uint edgeIndex = 0; edgeIndex < command.edgeCount; ++edgeIndex) {
        const device EdgeRecord& edge = edges[edgeIndex];
        if (edge.analyteIndex != analyte) continue;
        if (edge.sourceCompartment >= command.compartmentCount ||
            edge.targetCompartment >= command.compartmentCount ||
            !(edge.partitionCoefficient > 0.0f) ||
            edge.clearance < 0.0f) {
            return NAN;
        }

        const float source = state[
            analyte * command.compartmentCount + edge.sourceCompartment
        ];
        const float target = state[
            analyte * command.compartmentCount + edge.targetCompartment
        ];
        const float driving = edge.sourceUnboundFraction * source -
            edge.targetUnboundFraction * target / edge.partitionCoefficient;
        const float amountRate = edge.clearance * driving;
        if (compartment == edge.sourceCompartment) {
            rate -= amountRate * compartments[edge.sourceCompartment].inverseVolume;
        }
        if (compartment == edge.targetCompartment) {
            rate += amountRate * compartments[edge.targetCompartment].inverseVolume;
        }
    }
    return rate;
}

kernel void nvivo_partition_clear_status(
    device RuntimeStatus* status [[buffer(0)]],
    uint index [[thread_position_in_grid]]
) {
    if (index != 0) return;
    atomic_store_explicit(&status->flags, 0u, memory_order_relaxed);
    atomic_store_explicit(&status->invalidElementCount, 0u, memory_order_relaxed);
    atomic_store_explicit(&status->firstInvalidElement, kInvalidIndex, memory_order_relaxed);
    atomic_store_explicit(&status->maximumConcentrationBits, 0u, memory_order_relaxed);
    atomic_store_explicit(&status->maximumNegativeMagnitudeBits, 0u, memory_order_relaxed);
    atomic_store_explicit(&status->reserved0, 0u, memory_order_relaxed);
    atomic_store_explicit(&status->reserved1, 0u, memory_order_relaxed);
    atomic_store_explicit(&status->reserved2, 0u, memory_order_relaxed);
}

kernel void nvivo_partition_stage(
    const device float* baseState [[buffer(0)]],
    device float* stageState [[buffer(1)]],
    device float* firstDerivative [[buffer(2)]],
    const device CompartmentRecord* compartments [[buffer(3)]],
    const device EdgeRecord* edges [[buffer(4)]],
    constant RuntimeCommand& command [[buffer(5)]],
    uint element [[thread_position_in_grid]]
) {
    const uint elementCount = command.compartmentCount * command.analyteCount;
    if (element >= elementCount) return;
    const uint analyte = element / command.compartmentCount;
    const uint compartment = element - analyte * command.compartmentCount;
    const float rate = derivative(baseState, compartments, edges, analyte, compartment, command);
    firstDerivative[element] = rate;
    stageState[element] = baseState[element] + command.deltaTime * rate;
}

kernel void nvivo_partition_finalize(
    const device float* baseState [[buffer(0)]],
    const device float* stageState [[buffer(1)]],
    const device float* firstDerivative [[buffer(2)]],
    device float* candidateState [[buffer(3)]],
    const device CompartmentRecord* compartments [[buffer(4)]],
    const device EdgeRecord* edges [[buffer(5)]],
    constant RuntimeCommand& command [[buffer(6)]],
    uint element [[thread_position_in_grid]]
) {
    const uint elementCount = command.compartmentCount * command.analyteCount;
    if (element >= elementCount) return;
    const uint analyte = element / command.compartmentCount;
    const uint compartment = element - analyte * command.compartmentCount;
    const float secondDerivative = derivative(stageState, compartments, edges, analyte, compartment, command);
    candidateState[element] = baseState[element] +
        0.5f * command.deltaTime * (firstDerivative[element] + secondDerivative);
}

kernel void nvivo_partition_validate(
    device float* candidateState [[buffer(0)]],
    const device AnalyteRecord* analytes [[buffer(1)]],
    constant RuntimeCommand& command [[buffer(2)]],
    device RuntimeStatus* status [[buffer(3)]],
    uint element [[thread_position_in_grid]]
) {
    const uint elementCount = command.compartmentCount * command.analyteCount;
    if (element >= elementCount) return;
    const uint analyte = element / command.compartmentCount;
    const device AnalyteRecord& bounds = analytes[analyte];
    float value = candidateState[element];
    uint failure = 0u;

    if (!isfinite(value)) {
        failure = kStatusNonFinite;
    } else if (value < bounds.minimumConcentration) {
        if (value >= -command.negativeTolerance && bounds.minimumConcentration <= 0.0f) {
            value = 0.0f;
            candidateState[element] = 0.0f;
        } else {
            failure = kStatusNegative | kStatusNeedsSmallerStep;
            atomic_fetch_max_explicit(
                &status->maximumNegativeMagnitudeBits,
                as_type<uint>(fabs(value)),
                memory_order_relaxed
            );
        }
    } else if (value > bounds.maximumConcentration) {
        failure = kStatusAboveBound | kStatusNeedsSmallerStep;
    }

    if (value >= 0.0f && isfinite(value)) {
        atomic_fetch_max_explicit(
            &status->maximumConcentrationBits,
            as_type<uint>(value),
            memory_order_relaxed
        );
    }
    if (failure != 0u) {
        atomic_fetch_or_explicit(&status->flags, failure, memory_order_relaxed);
        atomic_fetch_add_explicit(&status->invalidElementCount, 1u, memory_order_relaxed);
        atomic_fetch_min_explicit(&status->firstInvalidElement, element, memory_order_relaxed);
    }
}

} // namespace numivivo_partition
