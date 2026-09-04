#include <metal_stdlib>
using namespace metal;

struct NVivoPhysiologyIncidence {
    uint sourcePairIndex;
    float coefficientPerSecond;
    uint flags;
    uint reserved;
};

struct NVivoPhysiologyClearance {
    float firstOrderRate;
    float maximumRate;
    float halfSaturation;
    uint flags;
};

struct NVivoPhysiologyTransform {
    uint pairIndex;
    uint environmentIndex;
    uint flags;
    uint reserved;
    float replacement;
    float additiveDelta;
    float minimum;
    float maximum;
};

struct NVivoPhysiologyPublicationRequest {
    uint pairIndex;
    uint environmentIndex;
    uint outputIndex;
    uint flags;
};

struct NVivoPhysiologyCommand {
    uint pairCount;
    uint environmentCount;
    uint stateElementCount;
    uint incidenceCount;

    uint preTransformCount;
    uint postTransformCount;
    uint publicationCount;
    uint stepIndex;

    uint appliedDoseCount;
    uint runtimeFlags;
    uint transactionWord0;
    uint transactionWord1;

    uint transactionWord2;
    uint transactionWord3;
    uint reservedWord0;
    uint reservedWord1;

    float dt;
    float absoluteTime;
    float minimumTimeStep;
    float maximumDerivative;
};

struct NVivoPhysiologyStatus {
    atomic_uint flags;
    atomic_uint violationCount;
    atomic_uint firstViolationPair;
    atomic_uint firstViolationEnvironment;

    atomic_uint maximumViolationBits;
    atomic_uint maximumDerivativeBits;
    atomic_uint suggestedTimeStepBits;
    atomic_uint requestedResponse;

    atomic_uint appliedDoseCount;
    atomic_uint appliedPreTransformCount;
    atomic_uint appliedPostTransformCount;
    atomic_uint publicationCount;

    atomic_uint stepIndex;
    atomic_uint transactionWord0;
    atomic_uint transactionWord1;
    atomic_uint reserved;
};

constant uint NVIVO_PHYS_REPLACE = 1u << 0u;
constant uint NVIVO_PHYS_MINIMUM = 1u << 1u;
constant uint NVIVO_PHYS_MAXIMUM = 1u << 2u;
constant uint NVIVO_PHYS_SATURABLE = 1u << 0u;

constant uint NVIVO_PHYS_NONFINITE = 1u << 0u;
constant uint NVIVO_PHYS_BELOW_MINIMUM = 1u << 1u;
constant uint NVIVO_PHYS_ABOVE_MAXIMUM = 1u << 2u;
constant uint NVIVO_PHYS_EXCESSIVE_DERIVATIVE = 1u << 3u;
constant uint NVIVO_PHYS_INVALID_TRANSFORM = 1u << 4u;
constant uint NVIVO_PHYS_INVALID_INCIDENCE = 1u << 5u;
constant uint NVIVO_PHYS_REQUIRES_SUBSTEP = 1u << 6u;
constant uint NVIVO_PHYS_REJECTED = 1u << 7u;

inline void nvivo_phys_atomic_max_positive(device atomic_uint* target, float value) {
    if (isfinite(value) && value > 0.0f) {
        atomic_fetch_max_explicit(target, as_type<uint>(value), memory_order_relaxed);
    }
}

inline void nvivo_phys_request_smaller_step(
    device NVivoPhysiologyStatus* status,
    constant NVivoPhysiologyCommand& command
) {
    atomic_fetch_or_explicit(
        &status->flags,
        NVIVO_PHYS_REQUIRES_SUBSTEP,
        memory_order_relaxed
    );
    const float proposed = max(command.minimumTimeStep, command.dt * 0.5f);
    atomic_fetch_min_explicit(
        &status->suggestedTimeStepBits,
        as_type<uint>(proposed),
        memory_order_relaxed
    );
    atomic_fetch_max_explicit(&status->requestedResponse, 1u, memory_order_relaxed);
}

inline void nvivo_phys_record_violation(
    device NVivoPhysiologyStatus* status,
    uint flag,
    uint pairIndex,
    uint environmentIndex,
    float magnitude
) {
    atomic_fetch_or_explicit(&status->flags, flag, memory_order_relaxed);
    atomic_fetch_add_explicit(&status->violationCount, 1u, memory_order_relaxed);
    uint expected = 0xffffffffu;
    if (atomic_compare_exchange_weak_explicit(
            &status->firstViolationPair,
            &expected,
            pairIndex,
            memory_order_relaxed,
            memory_order_relaxed)) {
        atomic_store_explicit(
            &status->firstViolationEnvironment,
            environmentIndex,
            memory_order_relaxed
        );
    }
    nvivo_phys_atomic_max_positive(&status->maximumViolationBits, magnitude);
}

inline float nvivo_phys_derivative(
    uint pairIndex,
    uint environmentIndex,
    device const float* state,
    device const uint* incidenceOffsets,
    device const NVivoPhysiologyIncidence* incidence,
    device const NVivoPhysiologyClearance* clearances,
    constant NVivoPhysiologyCommand& command,
    device NVivoPhysiologyStatus* status
) {
    const uint stateIndex = pairIndex * command.environmentCount + environmentIndex;
    float derivative = 0.0f;
    const uint begin = incidenceOffsets[pairIndex];
    const uint end = incidenceOffsets[pairIndex + 1u];
    if (begin > end || end > command.incidenceCount) {
        nvivo_phys_record_violation(
            status,
            NVIVO_PHYS_INVALID_INCIDENCE | NVIVO_PHYS_REJECTED,
            pairIndex,
            environmentIndex,
            1.0f
        );
        return 0.0f;
    }

    for (uint index = begin; index < end; ++index) {
        const NVivoPhysiologyIncidence entry = incidence[index];
        if (entry.sourcePairIndex >= command.pairCount || !isfinite(entry.coefficientPerSecond)) {
            nvivo_phys_record_violation(
                status,
                NVIVO_PHYS_INVALID_INCIDENCE | NVIVO_PHYS_REJECTED,
                pairIndex,
                environmentIndex,
                1.0f
            );
            continue;
        }
        const uint sourceIndex = entry.sourcePairIndex * command.environmentCount + environmentIndex;
        derivative = fma(entry.coefficientPerSecond, state[sourceIndex], derivative);
    }

    const float concentration = state[stateIndex];
    const NVivoPhysiologyClearance clearance = clearances[pairIndex];
    derivative = fma(-clearance.firstOrderRate, concentration, derivative);
    if ((clearance.flags & NVIVO_PHYS_SATURABLE) != 0u) {
        const float nonnegative = max(concentration, 0.0f);
        const float denominator = clearance.halfSaturation + nonnegative;
        if (clearance.maximumRate < 0.0f || clearance.halfSaturation <= 0.0f ||
            !isfinite(clearance.maximumRate) || !isfinite(clearance.halfSaturation) ||
            denominator <= 0.0f) {
            nvivo_phys_record_violation(
                status,
                NVIVO_PHYS_INVALID_INCIDENCE | NVIVO_PHYS_REJECTED,
                pairIndex,
                environmentIndex,
                1.0f
            );
        } else {
            derivative -= clearance.maximumRate * nonnegative / denominator;
        }
    }

    if (!isfinite(derivative)) {
        nvivo_phys_record_violation(
            status,
            NVIVO_PHYS_NONFINITE | NVIVO_PHYS_REJECTED,
            pairIndex,
            environmentIndex,
            1.0f
        );
        return 0.0f;
    }
    const float magnitude = abs(derivative);
    nvivo_phys_atomic_max_positive(&status->maximumDerivativeBits, magnitude);
    if (command.maximumDerivative > 0.0f && magnitude > command.maximumDerivative) {
        atomic_fetch_or_explicit(
            &status->flags,
            NVIVO_PHYS_EXCESSIVE_DERIVATIVE,
            memory_order_relaxed
        );
        nvivo_phys_request_smaller_step(status, command);
    }
    return derivative;
}

kernel void nvivo_phys_clear_status(
    constant NVivoPhysiologyCommand& command [[buffer(0)]],
    device NVivoPhysiologyStatus* status [[buffer(1)]],
    uint gid [[thread_position_in_grid]]
) {
    if (gid != 0u) return;
    atomic_store_explicit(&status->flags, 0u, memory_order_relaxed);
    atomic_store_explicit(&status->violationCount, 0u, memory_order_relaxed);
    atomic_store_explicit(&status->firstViolationPair, 0xffffffffu, memory_order_relaxed);
    atomic_store_explicit(&status->firstViolationEnvironment, 0xffffffffu, memory_order_relaxed);
    atomic_store_explicit(&status->maximumViolationBits, 0u, memory_order_relaxed);
    atomic_store_explicit(&status->maximumDerivativeBits, 0u, memory_order_relaxed);
    atomic_store_explicit(&status->suggestedTimeStepBits, as_type<uint>(command.dt), memory_order_relaxed);
    atomic_store_explicit(&status->requestedResponse, 0u, memory_order_relaxed);
    atomic_store_explicit(&status->appliedDoseCount, command.appliedDoseCount, memory_order_relaxed);
    atomic_store_explicit(&status->appliedPreTransformCount, command.preTransformCount, memory_order_relaxed);
    atomic_store_explicit(&status->appliedPostTransformCount, command.postTransformCount, memory_order_relaxed);
    atomic_store_explicit(&status->publicationCount, command.publicationCount, memory_order_relaxed);
    atomic_store_explicit(&status->stepIndex, command.stepIndex, memory_order_relaxed);
    atomic_store_explicit(&status->transactionWord0, command.transactionWord0, memory_order_relaxed);
    atomic_store_explicit(&status->transactionWord1, command.transactionWord1, memory_order_relaxed);
    atomic_store_explicit(&status->reserved, 0u, memory_order_relaxed);
}

kernel void nvivo_phys_prepare_transaction(
    device const float* currentState [[buffer(0)]],
    device float* baseState [[buffer(1)]],
    device float* candidateState [[buffer(2)]],
    device float* stageState [[buffer(3)]],
    device float* derivativeK1 [[buffer(4)]],
    constant NVivoPhysiologyCommand& command [[buffer(5)]],
    uint gid [[thread_position_in_grid]]
) {
    if (gid >= command.stateElementCount) return;
    const float value = currentState[gid];
    baseState[gid] = value;
    candidateState[gid] = value;
    stageState[gid] = value;
    derivativeK1[gid] = 0.0f;
}

kernel void nvivo_phys_apply_transforms(
    device const NVivoPhysiologyTransform* transforms [[buffer(0)]],
    device float* state [[buffer(1)]],
    constant NVivoPhysiologyCommand& command [[buffer(2)]],
    device NVivoPhysiologyStatus* status [[buffer(3)]],
    constant uint& transformCount [[buffer(4)]],
    uint gid [[thread_position_in_grid]]
) {
    if (gid >= transformCount) return;
    const NVivoPhysiologyTransform transform = transforms[gid];
    if (transform.pairIndex >= command.pairCount ||
        transform.environmentIndex >= command.environmentCount ||
        !isfinite(transform.replacement) ||
        !isfinite(transform.additiveDelta) ||
        !isfinite(transform.minimum) ||
        !isfinite(transform.maximum) ||
        transform.minimum > transform.maximum) {
        nvivo_phys_record_violation(
            status,
            NVIVO_PHYS_INVALID_TRANSFORM | NVIVO_PHYS_REJECTED,
            transform.pairIndex,
            transform.environmentIndex,
            1.0f
        );
        return;
    }
    const uint index = transform.pairIndex * command.environmentCount + transform.environmentIndex;
    float value = (transform.flags & NVIVO_PHYS_REPLACE) != 0u
        ? transform.replacement
        : state[index];
    value += transform.additiveDelta;
    if ((transform.flags & NVIVO_PHYS_MINIMUM) != 0u) value = max(value, transform.minimum);
    if ((transform.flags & NVIVO_PHYS_MAXIMUM) != 0u) value = min(value, transform.maximum);
    state[index] = value;
}

kernel void nvivo_phys_heun_predict(
    device const float* baseState [[buffer(0)]],
    device float* stageState [[buffer(1)]],
    device float* derivativeK1 [[buffer(2)]],
    device const uint* incidenceOffsets [[buffer(3)]],
    device const NVivoPhysiologyIncidence* incidence [[buffer(4)]],
    device const NVivoPhysiologyClearance* clearances [[buffer(5)]],
    constant NVivoPhysiologyCommand& command [[buffer(6)]],
    device NVivoPhysiologyStatus* status [[buffer(7)]],
    uint gid [[thread_position_in_grid]]
) {
    if (gid >= command.stateElementCount) return;
    const uint pairIndex = gid / command.environmentCount;
    const uint environmentIndex = gid - pairIndex * command.environmentCount;
    const float derivative = nvivo_phys_derivative(
        pairIndex,
        environmentIndex,
        baseState,
        incidenceOffsets,
        incidence,
        clearances,
        command,
        status
    );
    derivativeK1[gid] = derivative;
    stageState[gid] = fma(command.dt, derivative, baseState[gid]);
}

kernel void nvivo_phys_heun_correct(
    device const float* baseState [[buffer(0)]],
    device const float* stageState [[buffer(1)]],
    device const float* derivativeK1 [[buffer(2)]],
    device float* candidateState [[buffer(3)]],
    device const uint* incidenceOffsets [[buffer(4)]],
    device const NVivoPhysiologyIncidence* incidence [[buffer(5)]],
    device const NVivoPhysiologyClearance* clearances [[buffer(6)]],
    constant NVivoPhysiologyCommand& command [[buffer(7)]],
    device NVivoPhysiologyStatus* status [[buffer(8)]],
    uint gid [[thread_position_in_grid]]
) {
    if (gid >= command.stateElementCount) return;
    const uint pairIndex = gid / command.environmentCount;
    const uint environmentIndex = gid - pairIndex * command.environmentCount;
    const float derivativeK2 = nvivo_phys_derivative(
        pairIndex,
        environmentIndex,
        stageState,
        incidenceOffsets,
        incidence,
        clearances,
        command,
        status
    );
    candidateState[gid] = baseState[gid] + 0.5f * command.dt * (derivativeK1[gid] + derivativeK2);
}

kernel void nvivo_phys_validate_candidate(
    device const float* candidateState [[buffer(0)]],
    device const float2* bounds [[buffer(1)]],
    constant NVivoPhysiologyCommand& command [[buffer(2)]],
    device NVivoPhysiologyStatus* status [[buffer(3)]],
    uint gid [[thread_position_in_grid]]
) {
    if (gid >= command.stateElementCount) return;
    const uint pairIndex = gid / command.environmentCount;
    const uint environmentIndex = gid - pairIndex * command.environmentCount;
    const float value = candidateState[gid];
    const float2 interval = bounds[pairIndex];
    const float tolerance = as_type<float>(command.reservedWord0);
    if (!isfinite(value)) {
        nvivo_phys_record_violation(
            status,
            NVIVO_PHYS_NONFINITE | NVIVO_PHYS_REJECTED,
            pairIndex,
            environmentIndex,
            1.0f
        );
        return;
    }
    if (value < interval.x - tolerance) {
        nvivo_phys_record_violation(
            status,
            NVIVO_PHYS_BELOW_MINIMUM,
            pairIndex,
            environmentIndex,
            interval.x - value
        );
        nvivo_phys_request_smaller_step(status, command);
    }
    if (value > interval.y + tolerance) {
        nvivo_phys_record_violation(
            status,
            NVIVO_PHYS_ABOVE_MAXIMUM,
            pairIndex,
            environmentIndex,
            value - interval.y
        );
        nvivo_phys_request_smaller_step(status, command);
    }
}

kernel void nvivo_phys_publish(
    device const float* candidateState [[buffer(0)]],
    device const NVivoPhysiologyPublicationRequest* requests [[buffer(1)]],
    device float* output [[buffer(2)]],
    constant NVivoPhysiologyCommand& command [[buffer(3)]],
    device NVivoPhysiologyStatus* status [[buffer(4)]],
    uint gid [[thread_position_in_grid]]
) {
    if (gid >= command.publicationCount) return;
    const NVivoPhysiologyPublicationRequest request = requests[gid];
    if (request.pairIndex >= command.pairCount ||
        request.environmentIndex >= command.environmentCount ||
        request.outputIndex >= command.publicationCount) {
        nvivo_phys_record_violation(
            status,
            NVIVO_PHYS_INVALID_TRANSFORM | NVIVO_PHYS_REJECTED,
            request.pairIndex,
            request.environmentIndex,
            1.0f
        );
        return;
    }
    const uint stateIndex = request.pairIndex * command.environmentCount + request.environmentIndex;
    output[request.outputIndex] = candidateState[stateIndex];
}
