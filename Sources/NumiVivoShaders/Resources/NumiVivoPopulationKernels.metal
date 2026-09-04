#include <metal_stdlib>
using namespace metal;

namespace numivivo_population {

constant uint kBoundaryNoFlux = 0u;
constant uint kBoundaryPeriodic = 1u;
constant uint kBoundaryAbsorbing = 2u;
constant uint kInvalidIndex = 0xffffffffu;

constant uint kTransitionConstitutive = 0u;
constant uint kTransitionActivated = 1u;
constant uint kTransitionRepressed = 2u;

constant uint kStatusNeedsSmallerStep = 1u << 0u;
constant uint kStatusNonFinite = 1u << 1u;
constant uint kStatusNegativeDensity = 1u << 2u;
constant uint kStatusDensityOverflow = 1u << 3u;
constant uint kStatusInvalidTopology = 1u << 4u;

struct PhenotypeRecord {
    float birthRate;
    float deathRate;
    float carryingCapacity;
    float diffusionCoefficient;
    float minimumDensity;
    float maximumDensity;
    uint flags;
    uint reserved;
};

struct TransitionRecord {
    uint sourcePhenotype;
    uint destinationPhenotype;
    uint regulatorField;
    uint mode;
    float baseRate;
    float maximumRegulatedRate;
    float threshold;
    float hillCoefficient;
};

struct RuntimeCommand {
    uint phenotypeCount;
    uint voxelCount;
    uint regulatorFieldCount;
    uint transitionCount;
    uint width;
    uint height;
    uint depth;
    uint boundaryMode;
    float spacingX;
    float spacingY;
    float spacingZ;
    float deltaTime;
    float negativeTolerance;
    float absoluteTime;
    uint stepIndex;
    uint attemptIndex;
};

struct RuntimeStatus {
    atomic_uint flags;
    atomic_uint invalidElementCount;
    atomic_uint firstInvalidElement;
    atomic_uint maximumDensityBits;
    atomic_uint minimumDensityMagnitudeBits;
    atomic_uint reserved0;
    atomic_uint reserved1;
    atomic_uint reserved2;
};

inline uint flatten(uint3 coordinate, constant RuntimeCommand& command) {
    return coordinate.x + command.width * (coordinate.y + command.height * coordinate.z);
}

inline uint3 unflatten(uint voxel, constant RuntimeCommand& command) {
    const uint plane = command.width * command.height;
    const uint z = voxel / plane;
    const uint remainder = voxel - z * plane;
    const uint y = remainder / command.width;
    const uint x = remainder - y * command.width;
    return uint3(x, y, z);
}

inline float load_boundary(
    const device float* state,
    uint phenotype,
    int3 coordinate,
    constant RuntimeCommand& command
) {
    bool outside = coordinate.x < 0 || coordinate.y < 0 || coordinate.z < 0 ||
                   coordinate.x >= int(command.width) ||
                   coordinate.y >= int(command.height) ||
                   coordinate.z >= int(command.depth);
    if (!outside) {
        const uint voxel = flatten(uint3(coordinate), command);
        return state[phenotype * command.voxelCount + voxel];
    }

    if (command.boundaryMode == kBoundaryAbsorbing) return 0.0f;
    if (command.boundaryMode == kBoundaryPeriodic) {
        int3 wrapped = coordinate;
        wrapped.x = (wrapped.x % int(command.width) + int(command.width)) % int(command.width);
        wrapped.y = (wrapped.y % int(command.height) + int(command.height)) % int(command.height);
        wrapped.z = (wrapped.z % int(command.depth) + int(command.depth)) % int(command.depth);
        const uint voxel = flatten(uint3(wrapped), command);
        return state[phenotype * command.voxelCount + voxel];
    }

    const int3 clamped = int3(
        clamp(coordinate.x, 0, int(command.width) - 1),
        clamp(coordinate.y, 0, int(command.height) - 1),
        clamp(coordinate.z, 0, int(command.depth) - 1)
    );
    const uint voxel = flatten(uint3(clamped), command);
    return state[phenotype * command.voxelCount + voxel];
}

inline float regulated_rate(
    const device TransitionRecord& transition,
    const device float* regulatorFields,
    uint voxel,
    constant RuntimeCommand& command
) {
    float rate = max(transition.baseRate, 0.0f);
    if (transition.mode == kTransitionConstitutive || transition.regulatorField == kInvalidIndex) {
        return rate;
    }
    if (transition.regulatorField >= command.regulatorFieldCount) return NAN;

    const float signal = max(
        regulatorFields[transition.regulatorField * command.voxelCount + voxel],
        0.0f
    );
    const float threshold = max(transition.threshold, 1.0e-20f);
    const float coefficient = max(transition.hillCoefficient, 1.0f);
    const float poweredSignal = metal::powr(signal, coefficient);
    const float poweredThreshold = metal::powr(threshold, coefficient);
    const float denominator = poweredSignal + poweredThreshold;
    if (!(denominator > 0.0f) || !isfinite(denominator)) return NAN;

    if (transition.mode == kTransitionActivated) {
        rate += max(transition.maximumRegulatedRate, 0.0f) * poweredSignal / denominator;
    } else if (transition.mode == kTransitionRepressed) {
        rate += max(transition.maximumRegulatedRate, 0.0f) * poweredThreshold / denominator;
    } else {
        return NAN;
    }
    return rate;
}

inline float derivative(
    const device float* state,
    const device PhenotypeRecord* phenotypes,
    const device TransitionRecord* transitions,
    const device float* interactionMatrix,
    const device float* regulatorFields,
    const device float4* velocity,
    uint phenotype,
    uint voxel,
    constant RuntimeCommand& command
) {
    const device PhenotypeRecord& phenotypeRecord = phenotypes[phenotype];
    const float density = state[phenotype * command.voxelCount + voxel];
    if (!isfinite(density)) return NAN;

    float totalDensity = 0.0f;
    for (uint other = 0; other < command.phenotypeCount; ++other) {
        totalDensity += max(state[other * command.voxelCount + voxel], 0.0f);
    }

    const float carryingCapacity = max(phenotypeRecord.carryingCapacity, 1.0e-20f);
    const float logistic = phenotypeRecord.birthRate * density *
        (1.0f - totalDensity / carryingCapacity);
    float result = logistic - phenotypeRecord.deathRate * density;

    float interactionRate = 0.0f;
    const uint interactionBase = phenotype * command.phenotypeCount;
    for (uint actor = 0; actor < command.phenotypeCount; ++actor) {
        interactionRate += interactionMatrix[interactionBase + actor] *
            max(state[actor * command.voxelCount + voxel], 0.0f);
    }
    result += density * interactionRate;

    for (uint transitionIndex = 0; transitionIndex < command.transitionCount; ++transitionIndex) {
        const device TransitionRecord& transition = transitions[transitionIndex];
        const float rate = regulated_rate(transition, regulatorFields, voxel, command);
        if (!isfinite(rate) || rate < 0.0f) return NAN;
        if (transition.sourcePhenotype == phenotype) {
            result -= rate * density;
        }
        if (transition.destinationPhenotype == phenotype) {
            const float sourceDensity = state[
                transition.sourcePhenotype * command.voxelCount + voxel
            ];
            result += rate * max(sourceDensity, 0.0f);
        }
    }

    const uint3 coordinate = unflatten(voxel, command);
    const int3 location = int3(coordinate);
    const float left = load_boundary(state, phenotype, location + int3(-1, 0, 0), command);
    const float right = load_boundary(state, phenotype, location + int3(1, 0, 0), command);
    const float down = load_boundary(state, phenotype, location + int3(0, -1, 0), command);
    const float up = load_boundary(state, phenotype, location + int3(0, 1, 0), command);
    const float back = load_boundary(state, phenotype, location + int3(0, 0, -1), command);
    const float front = load_boundary(state, phenotype, location + int3(0, 0, 1), command);

    const float inverseDX2 = 1.0f / max(command.spacingX * command.spacingX, 1.0e-30f);
    const float inverseDY2 = 1.0f / max(command.spacingY * command.spacingY, 1.0e-30f);
    const float inverseDZ2 = 1.0f / max(command.spacingZ * command.spacingZ, 1.0e-30f);
    const float laplacian = (left - 2.0f * density + right) * inverseDX2 +
                            (down - 2.0f * density + up) * inverseDY2 +
                            (back - 2.0f * density + front) * inverseDZ2;
    result += phenotypeRecord.diffusionCoefficient * laplacian;

    const float3 localVelocity = velocity[voxel].xyz;
    const float inverseDX = 1.0f / max(command.spacingX, 1.0e-20f);
    const float inverseDY = 1.0f / max(command.spacingY, 1.0e-20f);
    const float inverseDZ = 1.0f / max(command.spacingZ, 1.0e-20f);
    const float gradientX = localVelocity.x >= 0.0f
        ? (density - left) * inverseDX
        : (right - density) * inverseDX;
    const float gradientY = localVelocity.y >= 0.0f
        ? (density - down) * inverseDY
        : (up - density) * inverseDY;
    const float gradientZ = localVelocity.z >= 0.0f
        ? (density - back) * inverseDZ
        : (front - density) * inverseDZ;
    result -= dot(localVelocity, float3(gradientX, gradientY, gradientZ));
    return result;
}

kernel void nvivo_population_clear_status(
    device RuntimeStatus* status [[buffer(0)]],
    uint index [[thread_position_in_grid]]
) {
    if (index != 0) return;
    atomic_store_explicit(&status->flags, 0u, memory_order_relaxed);
    atomic_store_explicit(&status->invalidElementCount, 0u, memory_order_relaxed);
    atomic_store_explicit(&status->firstInvalidElement, kInvalidIndex, memory_order_relaxed);
    atomic_store_explicit(&status->maximumDensityBits, 0u, memory_order_relaxed);
    atomic_store_explicit(&status->minimumDensityMagnitudeBits, 0u, memory_order_relaxed);
    atomic_store_explicit(&status->reserved0, 0u, memory_order_relaxed);
    atomic_store_explicit(&status->reserved1, 0u, memory_order_relaxed);
    atomic_store_explicit(&status->reserved2, 0u, memory_order_relaxed);
}

kernel void nvivo_population_stage(
    const device float* baseState [[buffer(0)]],
    device float* stageState [[buffer(1)]],
    device float* firstDerivative [[buffer(2)]],
    const device PhenotypeRecord* phenotypes [[buffer(3)]],
    const device TransitionRecord* transitions [[buffer(4)]],
    const device float* interactionMatrix [[buffer(5)]],
    const device float* regulatorFields [[buffer(6)]],
    const device float4* velocity [[buffer(7)]],
    constant RuntimeCommand& command [[buffer(8)]],
    uint element [[thread_position_in_grid]]
) {
    const uint elementCount = command.phenotypeCount * command.voxelCount;
    if (element >= elementCount) return;
    const uint phenotype = element / command.voxelCount;
    const uint voxel = element - phenotype * command.voxelCount;
    const float rate = derivative(
        baseState,
        phenotypes,
        transitions,
        interactionMatrix,
        regulatorFields,
        velocity,
        phenotype,
        voxel,
        command
    );
    firstDerivative[element] = rate;
    stageState[element] = baseState[element] + command.deltaTime * rate;
}

kernel void nvivo_population_finalize(
    const device float* baseState [[buffer(0)]],
    const device float* stageState [[buffer(1)]],
    const device float* firstDerivative [[buffer(2)]],
    device float* candidateState [[buffer(3)]],
    const device PhenotypeRecord* phenotypes [[buffer(4)]],
    const device TransitionRecord* transitions [[buffer(5)]],
    const device float* interactionMatrix [[buffer(6)]],
    const device float* regulatorFields [[buffer(7)]],
    const device float4* velocity [[buffer(8)]],
    constant RuntimeCommand& command [[buffer(9)]],
    uint element [[thread_position_in_grid]]
) {
    const uint elementCount = command.phenotypeCount * command.voxelCount;
    if (element >= elementCount) return;
    const uint phenotype = element / command.voxelCount;
    const uint voxel = element - phenotype * command.voxelCount;
    const float secondDerivative = derivative(
        stageState,
        phenotypes,
        transitions,
        interactionMatrix,
        regulatorFields,
        velocity,
        phenotype,
        voxel,
        command
    );
    candidateState[element] = baseState[element] +
        0.5f * command.deltaTime * (firstDerivative[element] + secondDerivative);
}

kernel void nvivo_population_validate(
    device float* candidateState [[buffer(0)]],
    const device PhenotypeRecord* phenotypes [[buffer(1)]],
    constant RuntimeCommand& command [[buffer(2)]],
    device RuntimeStatus* status [[buffer(3)]],
    uint element [[thread_position_in_grid]]
) {
    const uint elementCount = command.phenotypeCount * command.voxelCount;
    if (element >= elementCount) return;
    const uint phenotype = element / command.voxelCount;
    const device PhenotypeRecord& definition = phenotypes[phenotype];
    float value = candidateState[element];

    uint failure = 0u;
    if (!isfinite(value)) {
        failure = kStatusNonFinite;
    } else if (value < definition.minimumDensity) {
        if (value >= -command.negativeTolerance && definition.minimumDensity <= 0.0f) {
            value = 0.0f;
            candidateState[element] = 0.0f;
        } else {
            failure = kStatusNegativeDensity | kStatusNeedsSmallerStep;
            atomic_fetch_max_explicit(
                &status->minimumDensityMagnitudeBits,
                as_type<uint>(fabs(value)),
                memory_order_relaxed
            );
        }
    } else if (value > definition.maximumDensity) {
        failure = kStatusDensityOverflow | kStatusNeedsSmallerStep;
    }

    if (value >= 0.0f && isfinite(value)) {
        atomic_fetch_max_explicit(
            &status->maximumDensityBits,
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

} // namespace numivivo_population
