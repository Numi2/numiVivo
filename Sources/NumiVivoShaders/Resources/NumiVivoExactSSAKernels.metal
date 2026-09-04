#include <metal_stdlib>
using namespace metal;

namespace numivivo_exact_ssa {

constant uint kFlagNeedsSmallerStep = 1u << 0u;
constant uint kFlagInvalidPropensity = 1u << 1u;
constant uint kFlagNegativeCount = 1u << 2u;
constant uint kFlagCountOverflow = 1u << 3u;
constant uint kFlagInvalidModel = 1u << 4u;

constant uint kLawZeroOrder = 0u;
constant uint kLawFirstOrder = 1u;
constant uint kLawSecondOrderDistinct = 2u;
constant uint kLawSecondOrderSame = 3u;
constant uint kLawHillActivation = 4u;
constant uint kLawHillRepression = 5u;
constant uint kInvalidIndex = 0xffffffffu;

struct ReactionRecord {
    uint law;
    uint reactantA;
    uint reactantB;
    uint stoichiometryOffset;
    uint stoichiometryCount;
    uint flags;
    uint reserved0;
    uint reserved1;
    float rate;
    float parameter1;
    float parameter2;
    float parameter3;
};

struct StoichiometryRecord {
    uint speciesIndex;
    int delta;
};

struct RuntimeCommand {
    uint laneCount;
    uint speciesCount;
    uint reactionCount;
    uint maximumEventsPerLane;
    uint stepIndex;
    uint attemptIndex;
    uint reserved0;
    uint reserved1;
    float deltaTime;
    float absoluteTime;
    ulong seed;
};

struct RuntimeStatus {
    atomic_uint flags;
    atomic_uint lanesNeedingSmallerStep;
    atomic_uint invalidLaneCount;
    atomic_uint totalEventCount;
    atomic_uint maximumEventsObserved;
    atomic_uint firstInvalidLane;
    atomic_uint firstInvalidReaction;
    atomic_uint reserved;
};

inline uint multiply_high(uint left, uint right) {
    return uint((ulong(left) * ulong(right)) >> 32u);
}

inline uint4 philox_round(uint4 counter, uint2 key) {
    constexpr uint multiplier0 = 0xD2511F53u;
    constexpr uint multiplier1 = 0xCD9E8D57u;
    const uint high0 = multiply_high(multiplier0, counter.x);
    const uint low0 = multiplier0 * counter.x;
    const uint high1 = multiply_high(multiplier1, counter.z);
    const uint low1 = multiplier1 * counter.z;
    return uint4(high1 ^ counter.y ^ key.x, low1, high0 ^ counter.w ^ key.y, low0);
}

inline uint4 philox(uint4 counter, uint2 key) {
    constexpr uint keyIncrement0 = 0x9E3779B9u;
    constexpr uint keyIncrement1 = 0xBB67AE85u;
    for (uint round = 0; round < 10; ++round) {
        counter = philox_round(counter, key);
        key += uint2(keyIncrement0, keyIncrement1);
    }
    return counter;
}

inline float uniform_open(uint value) {
    // The half-unit offset prevents exact zero and one.
    return (float(value) + 0.5f) * 0x1.0p-32f;
}

inline float reaction_propensity(
    const device ReactionRecord& reaction,
    const device uint* counts,
    uint lane,
    constant RuntimeCommand& command
) {
    const uint laneCount = command.laneCount;
    const float a = reaction.reactantA == kInvalidIndex
        ? 0.0f
        : float(counts[reaction.reactantA * laneCount + lane]);
    const float b = reaction.reactantB == kInvalidIndex
        ? 0.0f
        : float(counts[reaction.reactantB * laneCount + lane]);

    switch (reaction.law) {
        case kLawZeroOrder:
            return reaction.rate;
        case kLawFirstOrder:
            return reaction.rate * a;
        case kLawSecondOrderDistinct:
            return reaction.rate * a * b;
        case kLawSecondOrderSame:
            return reaction.rate * 0.5f * a * max(a - 1.0f, 0.0f);
        case kLawHillActivation: {
            const float coefficient = max(reaction.parameter2, 1.0f);
            const float threshold = max(reaction.parameter1, 1.0e-20f);
            const float poweredSignal = metal::powr(max(a, 0.0f), coefficient);
            const float poweredThreshold = metal::powr(threshold, coefficient);
            const float denominator = poweredThreshold + poweredSignal;
            return denominator > 0.0f ? reaction.rate * poweredSignal / denominator : 0.0f;
        }
        case kLawHillRepression: {
            const float coefficient = max(reaction.parameter2, 1.0f);
            const float threshold = max(reaction.parameter1, 1.0e-20f);
            const float poweredSignal = metal::powr(max(a, 0.0f), coefficient);
            const float poweredThreshold = metal::powr(threshold, coefficient);
            const float denominator = poweredThreshold + poweredSignal;
            return denominator > 0.0f ? reaction.rate * poweredThreshold / denominator : 0.0f;
        }
        default:
            return NAN;
    }
}

inline bool can_apply_reaction(
    const device StoichiometryRecord* stoichiometry,
    const device uint* counts,
    const device ReactionRecord& reaction,
    uint lane,
    uint laneCount,
    thread uint& failureFlag
) {
    for (uint index = 0; index < reaction.stoichiometryCount; ++index) {
        const device StoichiometryRecord& term = stoichiometry[reaction.stoichiometryOffset + index];
        const uint current = counts[term.speciesIndex * laneCount + lane];
        if (term.delta < 0) {
            const uint consumed = uint(-term.delta);
            if (current < consumed) {
                failureFlag = kFlagNegativeCount;
                return false;
            }
        } else if (term.delta > 0) {
            const uint produced = uint(term.delta);
            if (current > 0xffffffffu - produced) {
                failureFlag = kFlagCountOverflow;
                return false;
            }
        }
    }
    return true;
}

inline void apply_reaction(
    const device StoichiometryRecord* stoichiometry,
    device uint* counts,
    const device ReactionRecord& reaction,
    uint lane,
    uint laneCount
) {
    for (uint index = 0; index < reaction.stoichiometryCount; ++index) {
        const device StoichiometryRecord& term = stoichiometry[reaction.stoichiometryOffset + index];
        const uint offset = term.speciesIndex * laneCount + lane;
        if (term.delta < 0) {
            counts[offset] -= uint(-term.delta);
        } else {
            counts[offset] += uint(term.delta);
        }
    }
}

kernel void nvivo_exact_ssa_clear_status(
    device RuntimeStatus* status [[buffer(0)]],
    uint index [[thread_position_in_grid]]
) {
    if (index != 0) return;
    atomic_store_explicit(&status->flags, 0u, memory_order_relaxed);
    atomic_store_explicit(&status->lanesNeedingSmallerStep, 0u, memory_order_relaxed);
    atomic_store_explicit(&status->invalidLaneCount, 0u, memory_order_relaxed);
    atomic_store_explicit(&status->totalEventCount, 0u, memory_order_relaxed);
    atomic_store_explicit(&status->maximumEventsObserved, 0u, memory_order_relaxed);
    atomic_store_explicit(&status->firstInvalidLane, kInvalidIndex, memory_order_relaxed);
    atomic_store_explicit(&status->firstInvalidReaction, kInvalidIndex, memory_order_relaxed);
    atomic_store_explicit(&status->reserved, 0u, memory_order_relaxed);
}

kernel void nvivo_exact_ssa_advance(
    device uint* counts [[buffer(0)]],
    const device ReactionRecord* reactions [[buffer(1)]],
    const device StoichiometryRecord* stoichiometry [[buffer(2)]],
    constant RuntimeCommand& command [[buffer(3)]],
    device RuntimeStatus* status [[buffer(4)]],
    device uint* eventCounts [[buffer(5)]],
    uint lane [[thread_position_in_grid]]
) {
    if (lane >= command.laneCount) return;

    float localTime = 0.0f;
    uint eventCount = 0u;
    uint failureFlag = 0u;
    uint invalidReaction = kInvalidIndex;

    while (localTime < command.deltaTime && eventCount < command.maximumEventsPerLane) {
        float totalPropensity = 0.0f;
        for (uint reactionIndex = 0; reactionIndex < command.reactionCount; ++reactionIndex) {
            const float propensity = reaction_propensity(
                reactions[reactionIndex],
                counts,
                lane,
                command
            );
            if (!isfinite(propensity) || propensity < 0.0f) {
                failureFlag = kFlagInvalidPropensity;
                invalidReaction = reactionIndex;
                break;
            }
            totalPropensity += propensity;
        }
        if (failureFlag != 0u || !isfinite(totalPropensity)) {
            if (failureFlag == 0u) failureFlag = kFlagInvalidPropensity;
            break;
        }
        if (totalPropensity <= 0.0f) break;

        const uint4 random = philox(
            uint4(eventCount, lane, command.stepIndex, command.attemptIndex),
            uint2(uint(command.seed), uint(command.seed >> 32u))
        );
        const float timeUniform = uniform_open(random.x);
        const float reactionUniform = uniform_open(random.y);
        const float waitingTime = -metal::fast::log(timeUniform) / totalPropensity;
        if (!isfinite(waitingTime) || waitingTime < 0.0f) {
            failureFlag = kFlagInvalidPropensity;
            break;
        }
        if (localTime + waitingTime > command.deltaTime) break;

        const float threshold = reactionUniform * totalPropensity;
        float cumulative = 0.0f;
        uint selected = command.reactionCount - 1u;
        for (uint reactionIndex = 0; reactionIndex < command.reactionCount; ++reactionIndex) {
            const float propensity = reaction_propensity(
                reactions[reactionIndex],
                counts,
                lane,
                command
            );
            cumulative += max(propensity, 0.0f);
            if (threshold <= cumulative) {
                selected = reactionIndex;
                break;
            }
        }

        if (!can_apply_reaction(
                stoichiometry,
                counts,
                reactions[selected],
                lane,
                command.laneCount,
                failureFlag)) {
            invalidReaction = selected;
            break;
        }
        apply_reaction(
            stoichiometry,
            counts,
            reactions[selected],
            lane,
            command.laneCount
        );
        localTime += waitingTime;
        ++eventCount;
    }

    eventCounts[lane] = eventCount;
    atomic_fetch_add_explicit(&status->totalEventCount, eventCount, memory_order_relaxed);
    atomic_fetch_max_explicit(&status->maximumEventsObserved, eventCount, memory_order_relaxed);

    if (eventCount == command.maximumEventsPerLane && localTime < command.deltaTime) {
        atomic_fetch_or_explicit(&status->flags, kFlagNeedsSmallerStep, memory_order_relaxed);
        atomic_fetch_add_explicit(&status->lanesNeedingSmallerStep, 1u, memory_order_relaxed);
    }
    if (failureFlag != 0u) {
        atomic_fetch_or_explicit(&status->flags, failureFlag, memory_order_relaxed);
        atomic_fetch_add_explicit(&status->invalidLaneCount, 1u, memory_order_relaxed);
        atomic_fetch_min_explicit(&status->firstInvalidLane, lane, memory_order_relaxed);
        if (invalidReaction != kInvalidIndex) {
            atomic_fetch_min_explicit(&status->firstInvalidReaction, invalidReaction, memory_order_relaxed);
        }
    }
}

kernel void nvivo_exact_ssa_validate_counts(
    const device uint* counts [[buffer(0)]],
    constant RuntimeCommand& command [[buffer(1)]],
    device RuntimeStatus* status [[buffer(2)]],
    uint index [[thread_position_in_grid]]
) {
    const ulong elementCount = ulong(command.speciesCount) * ulong(command.laneCount);
    if (ulong(index) >= elementCount) return;
    // UInt32 state is finite by construction. This pass deliberately remains
    // separate so future ProgramPack bounds can be attached without changing
    // the stochastic event kernel or its random stream.
    (void)counts[index];
    (void)status;
}

} // namespace numivivo_exact_ssa
