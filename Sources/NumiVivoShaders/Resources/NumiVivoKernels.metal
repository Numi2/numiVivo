#ifndef NUMIVIVO_METAL_ABI_H
#include "NumiVivoMetalABI.h"
#endif

namespace numivivo {

constant float kBooleanThreshold = 0.5f;
constant float kSmallValue = 1.0e-20f;
constant float kTwoPi = 6.2831853071795864769f;

inline ulong state_index(uint species, uint cell, constant NVivoStepUniforms& uniforms) {
    return ulong(species) * ulong(uniforms.cellCapacity) + ulong(cell);
}

inline ulong temporal_index(uint temporalState, uint cell, constant NVivoStepUniforms& uniforms) {
    return ulong(temporalState) * ulong(uniforms.cellCapacity) + ulong(cell);
}

inline ulong reaction_index(uint reaction, uint cell, constant NVivoStepUniforms& uniforms) {
    return ulong(reaction) * ulong(uniforms.cellCapacity) + ulong(cell);
}

inline bool cell_is_active(constant NVivoProgramArguments& program,
                           uint cell,
                           constant NVivoStepUniforms& uniforms) {
    return cell < uniforms.activeCellCount && program.cellActiveMask[cell] != 0u;
}

inline void record_first(device atomic_uint* target, uint value) {
    uint expected = NVIVO_INVALID_INDEX;
    while (expected == NVIVO_INVALID_INDEX &&
           !atomic_compare_exchange_weak_explicit(
               target,
               &expected,
               value,
               memory_order_relaxed,
               memory_order_relaxed)) {
    }
}

inline void record_expression_fault(constant NVivoProgramArguments& program,
                                    uint cell,
                                    uint subject) {
    atomic_fetch_or_explicit(
        &program.diagnostics->flags,
        NVIVO_DIAGNOSTIC_EXPRESSION_FAULT | NVIVO_DIAGNOSTIC_REJECT,
        memory_order_relaxed
    );
    atomic_fetch_add_explicit(&program.diagnostics->expressionFaultCount, 1u, memory_order_relaxed);
    record_first(&program.diagnostics->firstCell, cell);
    record_first(&program.diagnostics->firstSubject, subject);
}

inline void record_flux_truncation(constant NVivoProgramArguments& program,
                                   uint cell,
                                   uint reaction) {
    atomic_fetch_or_explicit(
        &program.diagnostics->flags,
        NVIVO_DIAGNOSTIC_FLUX_TRUNCATION,
        memory_order_relaxed
    );
    atomic_fetch_add_explicit(&program.diagnostics->fluxTruncationCount, 1u, memory_order_relaxed);
    record_first(&program.diagnostics->firstCell, cell);
    record_first(&program.diagnostics->firstSubject, reaction);
}

inline uint4 philox_round(uint4 counter, uint2 key) {
    constexpr uint m0 = 0xD2511F53u;
    constexpr uint m1 = 0xCD9E8D57u;
    const uint hi0 = mulhi(m0, counter.x);
    const uint lo0 = m0 * counter.x;
    const uint hi1 = mulhi(m1, counter.z);
    const uint lo1 = m1 * counter.z;
    return uint4(
        hi1 ^ counter.y ^ key.x,
        lo1,
        hi0 ^ counter.w ^ key.y,
        lo0
    );
}

inline uint4 philox4x32_10(uint4 counter, uint2 key) {
    constexpr uint w0 = 0x9E3779B9u;
    constexpr uint w1 = 0xBB67AE85u;
    for (uint round = 0u; round < 10u; ++round) {
        counter = philox_round(counter, key);
        key += uint2(w0, w1);
    }
    return counter;
}

struct RandomStream {
    uint4 counter;
    uint2 key;
    uint4 block;
    uint lane;
};

inline RandomStream make_random_stream(uint cell,
                                       uint reaction,
                                       constant NVivoStepUniforms& uniforms) {
    RandomStream stream;
    stream.counter = uint4(
        cell,
        reaction,
        uniforms.logicalStepLow,
        uniforms.logicalStepHigh ^ (uniforms.substepIndex * 0x9E3779B9u)
    );
    stream.key = uint2(uniforms.seedLow, uniforms.seedHigh);
    stream.block = uint4(0u);
    stream.lane = 4u;
    return stream;
}

inline uint next_random_uint(thread RandomStream& stream) {
    if (stream.lane >= 4u) {
        stream.block = philox4x32_10(stream.counter, stream.key);
        stream.counter.w += 1u;
        stream.lane = 0u;
    }
    const uint value = stream.block[stream.lane];
    stream.lane += 1u;
    return value;
}

inline float next_random_uniform(thread RandomStream& stream) {
    const uint value = next_random_uint(stream);
    return (float(value) + 0.5f) * 2.3283064365386963e-10f;
}

inline float next_random_normal(thread RandomStream& stream) {
    const float u1 = max(next_random_uniform(stream), 1.0e-7f);
    const float u2 = next_random_uniform(stream);
    return sqrt(-2.0f * log(u1)) * cos(kTwoPi * u2);
}

inline float log_factorial(uint value) {
    switch (value) {
        case 0u: case 1u: return 0.0f;
        case 2u: return 0.6931471805599453f;
        case 3u: return 1.7917594692280550f;
        case 4u: return 3.1780538303479458f;
        case 5u: return 4.7874917427820460f;
        case 6u: return 6.5792512120101010f;
        case 7u: return 8.5251613610654150f;
        case 8u: return 10.604602902745251f;
        case 9u: return 12.801827480081469f;
        case 10u: return 15.104412573075516f;
        case 11u: return 17.502307845873887f;
        case 12u: return 19.987214495661885f;
        default: {
            const float x = float(value);
            const float inverse = 1.0f / x;
            const float inverse3 = inverse * inverse * inverse;
            return (x + 0.5f) * log(x) - x + 0.9189385332046727f +
                   inverse / 12.0f - inverse3 / 360.0f;
        }
    }
}

struct PoissonSample {
    uint count;
    bool fallback;
};

inline uint bounded_count_from_float(float value) {
    if (!(value > 0.0f)) return 0u;
    return value >= 4294967040.0f ? 0xffffffffu : uint(value);
}

inline PoissonSample sample_poisson(float lambda, thread RandomStream& stream) {
    PoissonSample result{0u, false};
    if (!(lambda > 0.0f) || !isfinite(lambda)) return result;

    if (lambda < 12.0f) {
        const float threshold = exp(-lambda);
        float product = 1.0f;
        uint count = 0u;
        for (uint attempt = 0u; attempt < 64u; ++attempt) {
            product *= next_random_uniform(stream);
            if (product <= threshold) {
                result.count = count;
                return result;
            }
            count += 1u;
        }
        result.fallback = true;
        result.count = bounded_count_from_float(max(0.0f, rint(lambda + sqrt(lambda) * next_random_normal(stream))));
        return result;
    }

    const float root = sqrt(lambda);
    const float b = 0.931f + 2.53f * root;
    const float a = -0.059f + 0.02483f * b;
    const float inverseAlpha = 1.1239f + 1.1328f / (b - 3.4f);
    const float vR = 0.9277f - 3.6224f / (b - 2.0f);

    for (uint attempt = 0u; attempt < NVIVO_MAX_POISSON_ATTEMPTS; ++attempt) {
        const float u = next_random_uniform(stream) - 0.5f;
        const float v = next_random_uniform(stream);
        const float us = 0.5f - fabs(u);
        if (us <= 0.0f) continue;
        const int candidate = int(floor((2.0f * a / us + b) * u + lambda + 0.43f));
        if (candidate < 0) continue;
        if (us >= 0.07f && v <= vR) {
            result.count = uint(candidate);
            return result;
        }
        if (us < 0.013f && v > us) continue;
        const float left = log(v * inverseAlpha / (a / (us * us) + b));
        const float right = -lambda + float(candidate) * log(lambda) - log_factorial(uint(candidate));
        if (left <= right) {
            result.count = uint(candidate);
            return result;
        }
    }

    result.fallback = true;
    result.count = bounded_count_from_float(max(0.0f, rint(lambda + root * next_random_normal(stream))));
    return result;
}

struct ExpressionResult {
    float value;
    bool fault;
};

inline ExpressionResult evaluate_expression(constant NVivoProgramArguments& program,
                                            constant NVivoStepUniforms& uniforms,
                                            uint expressionOffset,
                                            uint expressionCount,
                                            uint cell,
                                            uint subject) {
    thread float stack[64];
    uint depth = 0u;
    bool terminated = false;
    const uint maximumInstructions = expressionCount == 0u ? 4096u : expressionCount;

    for (uint index = 0u; index < maximumInstructions; ++index) {
        const NVivoExpressionInstruction instruction = program.expressions[expressionOffset + index];
        const ushort opcode = instruction.opcode;

        if (opcode == 255u) {
            terminated = true;
            break;
        }
        if (opcode == 0u) {
            if (depth >= 64u) return {0.0f, true};
            stack[depth++] = instruction.immediate;
            continue;
        }
        if (opcode == 1u) {
            if (depth >= 64u || instruction.operand >= uniforms.speciesCount) return {0.0f, true};
            stack[depth++] = program.shadowState[state_index(instruction.operand, cell, uniforms)];
            continue;
        }
        if (opcode == 2u) {
            if (depth >= 64u) return {0.0f, true};
            if ((instruction.flags & NVIVO_EXPRESSION_REFERENCE_IS_TIME) != 0u) {
                stack[depth++] = uniforms.absoluteTime;
            } else {
                if (instruction.operand >= uniforms.parameterCount) return {0.0f, true};
                stack[depth++] = program.parameters[instruction.operand].value;
            }
            continue;
        }

        if (depth == 0u) return {0.0f, true};
        if (opcode == 3u) {
            stack[depth - 1u] = stack[depth - 1u] > kBooleanThreshold ? 0.0f : 1.0f;
            continue;
        }
        if (opcode == 19u || opcode == 20u || opcode == 21u || opcode == 22u) {
            if (instruction.auxiliary >= uniforms.temporalStateCount) return {0.0f, true};
            const ulong temporal = temporal_index(instruction.auxiliary, cell, uniforms);
            const bool condition = stack[depth - 1u] > kBooleanThreshold;
            float state = program.shadowTemporalState[temporal];
            if (opcode == 19u) {
                state = condition ? state + uniforms.deltaTime : 0.0f;
                stack[depth - 1u] = state >= instruction.immediate ? 1.0f : 0.0f;
            } else if (opcode == 20u) {
                state = condition ? instruction.immediate : max(0.0f, state - uniforms.deltaTime);
                stack[depth - 1u] = state > 0.0f ? 1.0f : 0.0f;
            } else if (opcode == 21u) {
                const bool previous = state > kBooleanThreshold;
                stack[depth - 1u] = condition && !previous ? 1.0f : 0.0f;
                state = condition ? 1.0f : 0.0f;
            } else {
                const bool previous = state > kBooleanThreshold;
                stack[depth - 1u] = !condition && previous ? 1.0f : 0.0f;
                state = condition ? 1.0f : 0.0f;
            }
            program.shadowTemporalState[temporal] = state;
            continue;
        }

        if (opcode == 18u) {
            if (depth < 3u) return {0.0f, true};
            const float maximum = stack[--depth];
            const float minimum = stack[--depth];
            stack[depth - 1u] = clamp(stack[depth - 1u], minimum, maximum);
            continue;
        }

        if (depth < 2u) return {0.0f, true};
        const float right = stack[--depth];
        float& left = stack[depth - 1u];
        switch (opcode) {
            case 4u: left = (left > kBooleanThreshold && right > kBooleanThreshold) ? 1.0f : 0.0f; break;
            case 5u: left = (left > kBooleanThreshold || right > kBooleanThreshold) ? 1.0f : 0.0f; break;
            case 6u: left = left > right ? 1.0f : 0.0f; break;
            case 7u: left = left >= right ? 1.0f : 0.0f; break;
            case 8u: left = left < right ? 1.0f : 0.0f; break;
            case 9u: left = left <= right ? 1.0f : 0.0f; break;
            case 10u: left = left == right ? 1.0f : 0.0f; break;
            case 11u: left = left != right ? 1.0f : 0.0f; break;
            case 12u: left += right; break;
            case 13u: left -= right; break;
            case 14u: left *= right; break;
            case 15u:
                if (fabs(right) <= kSmallValue) return {0.0f, true};
                left /= right;
                break;
            case 16u: left = min(left, right); break;
            case 17u: left = max(left, right); break;
            default: return {0.0f, true};
        }
        if (!isfinite(left)) return {0.0f, true};
    }

    if (!terminated || depth != 1u || !isfinite(stack[0])) {
        record_expression_fault(program, cell, subject);
        return {0.0f, true};
    }
    return {stack[0], false};
}

inline float parameter_value(constant NVivoProgramArguments& program,
                             constant NVivoStepUniforms& uniforms,
                             const NVivoReactionRecord& reaction,
                             uint position,
                             thread bool& fault) {
    if (position >= reaction.parameterCount) {
        fault = true;
        return 0.0f;
    }
    const uint parameterIndex = program.reactionParameterIndices[reaction.parameterOffset + position];
    if (parameterIndex >= uniforms.parameterCount) {
        fault = true;
        return 0.0f;
    }
    return program.parameters[parameterIndex].value;
}

inline float species_value(constant NVivoProgramArguments& program,
                           constant NVivoStepUniforms& uniforms,
                           uint species,
                           uint cell,
                           thread bool& fault) {
    if (species >= uniforms.speciesCount) {
        fault = true;
        return 0.0f;
    }
    const float value = program.shadowState[state_index(species, cell, uniforms)];
    if (!isfinite(value)) {
        fault = true;
        return 0.0f;
    }
    return max(value, 0.0f);
}

inline float mass_action_product(constant NVivoProgramArguments& program,
                                 constant NVivoStepUniforms& uniforms,
                                 uint termOffset,
                                 uint termCount,
                                 uint cell,
                                 thread bool& fault) {
    float product = 1.0f;
    for (uint index = 0u; index < termCount; ++index) {
        const NVivoStoichiometryRecord term = program.stoichiometry[termOffset + index];
        const float value = species_value(program, uniforms, term.speciesIndex, cell, fault);
        if (fault) return 0.0f;
        const int coefficient = max(1, int(term.coefficient));
        product *= pow(value, float(coefficient));
        if (!isfinite(product)) {
            fault = true;
            return 0.0f;
        }
    }
    return product;
}

inline float reaction_rate(constant NVivoProgramArguments& program,
                           constant NVivoStepUniforms& uniforms,
                           const NVivoReactionRecord& reaction,
                           uint reactionIndexValue,
                           uint cell,
                           thread bool& fault) {
    if ((reaction.flags & NVIVO_REACTION_HAS_GATE) != 0u) {
        if (reaction.gateExpressionOffset == NVIVO_INVALID_INDEX) {
            fault = true;
            return 0.0f;
        }
        const auto gate = evaluate_expression(
            program,
            uniforms,
            reaction.gateExpressionOffset,
            0u,
            cell,
            reactionIndexValue
        );
        if (gate.fault) {
            fault = true;
            return 0.0f;
        }
        if (gate.value <= kBooleanThreshold) return 0.0f;
    }

    if (reaction.rateLaw == 255u) {
        if (reaction.expressionCount == 0u) {
            fault = true;
            return 0.0f;
        }
        const auto custom = evaluate_expression(
            program,
            uniforms,
            reaction.expressionOffset,
            reaction.expressionCount,
            cell,
            reactionIndexValue
        );
        fault = custom.fault;
        return custom.value;
    }

    const float reactants = mass_action_product(
        program,
        uniforms,
        reaction.reactantOffset,
        reaction.reactantCount,
        cell,
        fault
    );
    if (fault) return 0.0f;

    switch (reaction.rateLaw) {
        case 0u:
            return parameter_value(program, uniforms, reaction, 0u, fault);
        case 1u:
            return parameter_value(program, uniforms, reaction, 0u, fault) * reactants;
        case 2u: {
            if (reaction.reactantCount == 0u) { fault = true; return 0.0f; }
            const auto inputTerm = program.stoichiometry[reaction.reactantOffset];
            const float input = species_value(program, uniforms, inputTerm.speciesIndex, cell, fault);
            const float maximum = parameter_value(program, uniforms, reaction, 0u, fault);
            const float half = max(parameter_value(program, uniforms, reaction, 1u, fault), kSmallValue);
            const float exponent = max(parameter_value(program, uniforms, reaction, 2u, fault), kSmallValue);
            const float numerator = pow(input, exponent);
            return maximum * numerator / (pow(half, exponent) + numerator + kSmallValue);
        }
        case 3u: {
            if (reaction.reactantCount == 0u) { fault = true; return 0.0f; }
            const auto inputTerm = program.stoichiometry[reaction.reactantOffset];
            const float input = species_value(program, uniforms, inputTerm.speciesIndex, cell, fault);
            const float maximum = parameter_value(program, uniforms, reaction, 0u, fault);
            const float half = max(parameter_value(program, uniforms, reaction, 1u, fault), kSmallValue);
            const float exponent = max(parameter_value(program, uniforms, reaction, 2u, fault), kSmallValue);
            const float halfPower = pow(half, exponent);
            return maximum * halfPower / (halfPower + pow(input, exponent) + kSmallValue);
        }
        case 4u: {
            if (reaction.reactantCount == 0u) { fault = true; return 0.0f; }
            const auto inputTerm = program.stoichiometry[reaction.reactantOffset];
            const float input = species_value(program, uniforms, inputTerm.speciesIndex, cell, fault);
            const float maximum = parameter_value(program, uniforms, reaction, 0u, fault);
            const float half = max(parameter_value(program, uniforms, reaction, 1u, fault), kSmallValue);
            return maximum * input / (half + input);
        }
        case 5u: {
            const float forward = parameter_value(program, uniforms, reaction, 0u, fault) * reactants;
            const float products = mass_action_product(
                program,
                uniforms,
                reaction.productOffset,
                reaction.productCount,
                cell,
                fault
            );
            const float reverse = parameter_value(program, uniforms, reaction, 1u, fault) * products;
            return forward - reverse;
        }
        case 6u: {
            if (reaction.reactantCount == 0u || reaction.productCount == 0u) { fault = true; return 0.0f; }
            const uint source = program.stoichiometry[reaction.reactantOffset].speciesIndex;
            const uint destination = program.stoichiometry[reaction.productOffset].speciesIndex;
            const float sourceValue = species_value(program, uniforms, source, cell, fault);
            const float destinationValue = species_value(program, uniforms, destination, cell, fault);
            return parameter_value(program, uniforms, reaction, 0u, fault) * (sourceValue - destinationValue);
        }
        case 7u: {
            if (reaction.reactantCount == 0u || reaction.productCount == 0u) { fault = true; return 0.0f; }
            const uint source = program.stoichiometry[reaction.reactantOffset].speciesIndex;
            const uint destination = program.stoichiometry[reaction.productOffset].speciesIndex;
            const float sourceValue = species_value(program, uniforms, source, cell, fault);
            const float destinationValue = species_value(program, uniforms, destination, cell, fault);
            const float maximum = parameter_value(program, uniforms, reaction, 0u, fault);
            const float half = max(parameter_value(program, uniforms, reaction, 1u, fault), kSmallValue);
            return maximum * (sourceValue / (half + sourceValue) - destinationValue / (half + destinationValue));
        }
        case 8u:
            return parameter_value(program, uniforms, reaction, 0u, fault) * reactants;
        default:
            fault = true;
            return 0.0f;
    }
}

inline float cap_extent(constant NVivoProgramArguments& program,
                        constant NVivoStepUniforms& uniforms,
                        const NVivoReactionRecord& reaction,
                        uint cell,
                        float extent,
                        thread bool& truncated) {
    const bool forward = extent >= 0.0f;
    const uint offset = forward ? reaction.reactantOffset : reaction.productOffset;
    const uint count = forward ? reaction.reactantCount : reaction.productCount;
    float maximumExtent = fabs(extent);

    for (uint index = 0u; index < count; ++index) {
        const auto term = program.stoichiometry[offset + index];
        const auto species = program.species[term.speciesIndex];
        if ((species.flags & NVIVO_SPECIES_EXTERNALLY_OWNED) != 0u) continue;
        const float available = max(0.0f, program.shadowState[state_index(term.speciesIndex, cell, uniforms)]);
        maximumExtent = min(maximumExtent, available / float(max(1, int(term.coefficient))));
    }

    if (maximumExtent + 1.0e-6f < fabs(extent)) truncated = true;
    return forward ? maximumExtent : -maximumExtent;
}

inline void schedule_delay(constant NVivoProgramArguments& program,
                           constant NVivoStepUniforms& uniforms,
                           const NVivoReactionRecord& reaction,
                           uint reactionIndexValue,
                           uint cell,
                           thread float& extent,
                           thread bool& fault) {
    if ((reaction.flags & NVIVO_REACTION_DELAYED) == 0u) return;
    if (uniforms.delaySlotCount < 2u || !(uniforms.deltaTime > 0.0f)) {
        fault = true;
        return;
    }

    const ulong cellCapacity = ulong(uniforms.cellCapacity);
    const ulong slotStride = ulong(uniforms.reactionCount) * cellCapacity;
    const ulong reactionBase = ulong(reactionIndexValue) * cellCapacity + ulong(cell);
    const ulong dueIndex = ulong(uniforms.delayWriteSlot) * slotStride + reactionBase;
    const float due = program.shadowDelayedFlux[dueIndex];
    program.shadowDelayedFlux[dueIndex] = 0.0f;

    const uint futureSteps = max(1u, uint(ceil(reaction.delaySeconds / uniforms.deltaTime)));
    if (futureSteps >= uniforms.delaySlotCount) {
        fault = true;
        return;
    }
    const uint futureSlot = (uniforms.delayWriteSlot + futureSteps) % uniforms.delaySlotCount;
    const ulong futureIndex = ulong(futureSlot) * slotStride + reactionBase;
    program.shadowDelayedFlux[futureIndex] += extent;
    extent = due;
}

[[host_name("nvivo_initialize_program")]] kernel void nvivo_initialize_program(
    constant NVivoProgramArguments& program [[buffer(0)]],
    constant NVivoStepUniforms& uniforms [[buffer(1)]],
    uint gid [[thread_position_in_grid]]) {

    const ulong stateCount = ulong(uniforms.speciesCount) * ulong(uniforms.cellCapacity);
    const ulong temporalCount = ulong(uniforms.temporalStateCount) * ulong(uniforms.cellCapacity);
    const ulong fluxCount = ulong(uniforms.reactionCount) * ulong(uniforms.cellCapacity);
    const ulong delayCount = fluxCount * ulong(uniforms.delaySlotCount);
    const ulong refractoryCount = ulong(uniforms.ruleCount) * ulong(uniforms.cellCapacity);
    const ulong index = ulong(gid);

    if (index < stateCount) {
        const uint species = uint(index / ulong(uniforms.cellCapacity));
        const float value = program.species[species].initialValue;
        program.committedState[index] = value;
        program.shadowState[index] = value;
    }
    if (index < temporalCount) {
        program.committedTemporalState[index] = 0.0f;
        program.shadowTemporalState[index] = 0.0f;
    }
    if (index < fluxCount) program.reactionFlux[index] = 0.0f;
    if (index < delayCount) {
        program.committedDelayedFlux[index] = 0.0f;
        program.shadowDelayedFlux[index] = 0.0f;
    }
    if (index < refractoryCount) {
        program.committedRuleRefractory[index] = 0.0f;
        program.shadowRuleRefractory[index] = 0.0f;
    }

    if (gid == 0u) {
        atomic_store_explicit(program.eventCount, 0u, memory_order_relaxed);
        program.publication->committedStepLow = 0u;
        program.publication->committedStepHigh = 0u;
        program.publication->stateVersion = 0u;
        program.publication->status = NVIVO_PUBLICATION_COMMITTED;
        program.publication->diagnosticFlags = 0u;
        program.publication->shutdownState = 0u;
        program.publication->activeCellCount = uniforms.activeCellCount;
        program.publication->eventCount = 0u;
    }
}

[[host_name("nvivo_prepare_step")]] kernel void nvivo_prepare_step(
    constant NVivoProgramArguments& program [[buffer(0)]],
    constant NVivoStepUniforms& uniforms [[buffer(1)]],
    uint gid [[thread_position_in_grid]]) {

    const ulong stateCount = ulong(uniforms.speciesCount) * ulong(uniforms.cellCapacity);
    const ulong temporalCount = ulong(uniforms.temporalStateCount) * ulong(uniforms.cellCapacity);
    const ulong fluxCount = ulong(uniforms.reactionCount) * ulong(uniforms.cellCapacity);
    const ulong delayCount = fluxCount * ulong(uniforms.delaySlotCount);
    const ulong refractoryCount = ulong(uniforms.ruleCount) * ulong(uniforms.cellCapacity);
    const ulong index = ulong(gid);

    if (index < stateCount) program.shadowState[index] = program.committedState[index];
    if (index < temporalCount) program.shadowTemporalState[index] = program.committedTemporalState[index];
    if (index < fluxCount) program.reactionFlux[index] = 0.0f;
    if (index < delayCount) program.shadowDelayedFlux[index] = program.committedDelayedFlux[index];
    if (index < refractoryCount) program.shadowRuleRefractory[index] = program.committedRuleRefractory[index];

    if (gid == 0u) {
        atomic_store_explicit(&program.diagnostics->flags, 0u, memory_order_relaxed);
        atomic_store_explicit(&program.diagnostics->nonFiniteCount, 0u, memory_order_relaxed);
        atomic_store_explicit(&program.diagnostics->boundViolationCount, 0u, memory_order_relaxed);
        atomic_store_explicit(&program.diagnostics->monitorViolationCount, 0u, memory_order_relaxed);
        atomic_store_explicit(&program.diagnostics->shutdownCount, 0u, memory_order_relaxed);
        atomic_store_explicit(&program.diagnostics->expressionFaultCount, 0u, memory_order_relaxed);
        atomic_store_explicit(&program.diagnostics->stochasticFallbackCount, 0u, memory_order_relaxed);
        atomic_store_explicit(&program.diagnostics->fluxTruncationCount, 0u, memory_order_relaxed);
        atomic_store_explicit(&program.diagnostics->eventOverflowCount, 0u, memory_order_relaxed);
        atomic_store_explicit(&program.diagnostics->firstCell, NVIVO_INVALID_INDEX, memory_order_relaxed);
        atomic_store_explicit(&program.diagnostics->firstSubject, NVIVO_INVALID_INDEX, memory_order_relaxed);
        atomic_store_explicit(&program.diagnostics->requestedResponse, 0u, memory_order_relaxed);
        atomic_store_explicit(&program.diagnostics->maximumSeverity, 0u, memory_order_relaxed);
        atomic_store_explicit(program.eventCount, 0u, memory_order_relaxed);
    }
}

struct NVivoInputUpdate {
    uint speciesIndex;
    uint cellIndex;
    float value;
    uint mode;
};

[[host_name("nvivo_stage_input_updates")]] kernel void nvivo_stage_input_updates(
    constant NVivoProgramArguments& program [[buffer(0)]],
    device const NVivoInputUpdate* updates [[buffer(1)]],
    constant uint& updateCount [[buffer(2)]],
    constant NVivoStepUniforms& uniforms [[buffer(3)]],
    uint gid [[thread_position_in_grid]]) {

    if (gid >= updateCount) return;
    const NVivoInputUpdate update = updates[gid];
    if (update.speciesIndex >= uniforms.speciesCount || update.cellIndex >= uniforms.activeCellCount) {
        record_expression_fault(program, update.cellIndex, update.speciesIndex);
        return;
    }
    const NVivoSpeciesRecord species = program.species[update.speciesIndex];
    if ((species.flags & NVIVO_SPECIES_EXTERNALLY_OWNED) == 0u) {
        record_expression_fault(program, update.cellIndex, update.speciesIndex);
        return;
    }
    const ulong index = state_index(update.speciesIndex, update.cellIndex, uniforms);
    const float next = update.mode == 1u ? program.shadowState[index] + update.value : update.value;
    program.shadowState[index] = next;
}

[[host_name("nvivo_evaluate_reaction_cohort")]] kernel void nvivo_evaluate_reaction_cohort(
    constant NVivoProgramArguments& program [[buffer(0)]],
    constant NVivoStepUniforms& uniforms [[buffer(1)]],
    constant NVivoCohortUniforms& cohort [[buffer(2)]],
    uint gid [[thread_position_in_grid]]) {

    if (cohort.dispatchCellCount == 0u) return;
    const uint localReaction = gid / cohort.dispatchCellCount;
    const uint cell = gid - localReaction * cohort.dispatchCellCount;
    if (localReaction >= cohort.reactionCount || cell >= uniforms.activeCellCount) return;
    const uint reactionIndexValue = cohort.reactionOffset + localReaction;
    if (reactionIndexValue >= uniforms.reactionCount) return;
    const ulong outputIndex = reaction_index(reactionIndexValue, cell, uniforms);
    if (!cell_is_active(program, cell, uniforms)) {
        program.reactionFlux[outputIndex] = 0.0f;
        return;
    }

    const NVivoReactionRecord reaction = program.reactions[reactionIndexValue];
    bool fault = false;
    float rate = reaction_rate(program, uniforms, reaction, reactionIndexValue, cell, fault);
    if (fault || !isfinite(rate)) {
        record_expression_fault(program, cell, reactionIndexValue);
        program.reactionFlux[outputIndex] = 0.0f;
        return;
    }

    float extent = 0.0f;
    if (uniforms.mode == NVIVO_MODE_F2_STOCHASTIC &&
        (reaction.flags & NVIVO_REACTION_STOCHASTIC_ELIGIBLE) != 0u) {
        RandomStream random = make_random_stream(cell, reactionIndexValue, uniforms);
        const float mean = fabs(rate) * uniforms.deltaTime;
        const PoissonSample sample = sample_poisson(mean, random);
        extent = rate < 0.0f ? -float(sample.count) : float(sample.count);
        if (sample.fallback) {
            atomic_fetch_or_explicit(
                &program.diagnostics->flags,
                NVIVO_DIAGNOSTIC_STOCHASTIC_FALLBACK,
                memory_order_relaxed
            );
            atomic_fetch_add_explicit(
                &program.diagnostics->stochasticFallbackCount,
                1u,
                memory_order_relaxed
            );
        }
    } else {
        extent = rate * uniforms.deltaTime;
    }

    bool truncated = false;
    extent = cap_extent(program, uniforms, reaction, cell, extent, truncated);
    if (truncated) record_flux_truncation(program, cell, reactionIndexValue);

    schedule_delay(program, uniforms, reaction, reactionIndexValue, cell, extent, fault);
    if (fault) {
        record_expression_fault(program, cell, reactionIndexValue);
        extent = 0.0f;
    } else if ((reaction.flags & NVIVO_REACTION_DELAYED) != 0u) {
        truncated = false;
        extent = cap_extent(program, uniforms, reaction, cell, extent, truncated);
        if (truncated) record_flux_truncation(program, cell, reactionIndexValue);
    }
    program.reactionFlux[outputIndex] = extent;
}

[[host_name("nvivo_apply_incidence")]] kernel void nvivo_apply_incidence(
    constant NVivoProgramArguments& program [[buffer(0)]],
    constant NVivoStepUniforms& uniforms [[buffer(1)]],
    uint gid [[thread_position_in_grid]]) {

    if (uniforms.activeCellCount == 0u) return;
    const uint speciesIndexValue = gid / uniforms.activeCellCount;
    const uint cell = gid - speciesIndexValue * uniforms.activeCellCount;
    if (speciesIndexValue >= uniforms.speciesCount || cell >= uniforms.activeCellCount) return;
    if (!cell_is_active(program, cell, uniforms)) return;

    const NVivoSpeciesRecord species = program.species[speciesIndexValue];
    if ((species.flags & NVIVO_SPECIES_EXTERNALLY_OWNED) != 0u) return;

    float delta = 0.0f;
    const uint begin = program.speciesIncidenceOffsets[speciesIndexValue];
    const uint end = program.speciesIncidenceOffsets[speciesIndexValue + 1u];
    for (uint index = begin; index < end; ++index) {
        const NVivoIncidenceRecord incidence = program.speciesIncidence[index];
        if (incidence.reactionIndex >= uniforms.reactionCount) {
            record_expression_fault(program, cell, speciesIndexValue);
            return;
        }
        delta += float(incidence.netCoefficient) *
                 program.reactionFlux[reaction_index(incidence.reactionIndex, cell, uniforms)];
    }

    const ulong index = state_index(speciesIndexValue, cell, uniforms);
    float value = program.shadowState[index] + delta;
    if ((species.flags & NVIVO_SPECIES_COUNT_VALUED) != 0u && uniforms.mode == NVIVO_MODE_F2_STOCHASTIC) {
        value = rint(value);
    }
    program.shadowState[index] = value;
}

inline void emit_event(constant NVivoProgramArguments& program,
                       constant NVivoStepUniforms& uniforms,
                       uint cell,
                       uint kind,
                       uint subject,
                       float value0,
                       float value1,
                       uint flags) {
    const uint index = atomic_fetch_add_explicit(program.eventCount, 1u, memory_order_relaxed);
    if (index >= uniforms.eventCapacity) {
        atomic_fetch_or_explicit(
            &program.diagnostics->flags,
            NVIVO_DIAGNOSTIC_EVENT_OVERFLOW,
            memory_order_relaxed
        );
        atomic_fetch_add_explicit(&program.diagnostics->eventOverflowCount, 1u, memory_order_relaxed);
        return;
    }
    program.events[index] = {
        cell,
        kind,
        subject,
        uniforms.logicalStepLow,
        value0,
        value1,
        flags,
        uniforms.logicalStepHigh
    };
}

inline float action_value(constant NVivoProgramArguments& program,
                          constant NVivoStepUniforms& uniforms,
                          const NVivoActionRecord& action,
                          uint cell,
                          uint subject,
                          thread bool& fault) {
    if (action.expressionCount == 0u) return action.constantValue;
    const auto result = evaluate_expression(
        program,
        uniforms,
        action.expressionOffset,
        action.expressionCount,
        cell,
        subject
    );
    fault = result.fault;
    return result.value;
}

[[host_name("nvivo_apply_rules")]] kernel void nvivo_apply_rules(
    constant NVivoProgramArguments& program [[buffer(0)]],
    constant NVivoStepUniforms& uniforms [[buffer(1)]],
    uint cell [[thread_position_in_grid]]) {

    if (!cell_is_active(program, cell, uniforms)) return;

    for (uint ruleIndex = 0u; ruleIndex < uniforms.ruleCount; ++ruleIndex) {
        const NVivoRuleRecord rule = program.rules[ruleIndex];
        const ulong refractoryIndex = ulong(ruleIndex) * ulong(uniforms.cellCapacity) + ulong(cell);
        float refractory = max(0.0f, program.shadowRuleRefractory[refractoryIndex] - uniforms.deltaTime);
        program.shadowRuleRefractory[refractoryIndex] = refractory;
        if (refractory > 0.0f) continue;

        const auto condition = evaluate_expression(
            program,
            uniforms,
            rule.conditionOffset,
            rule.conditionCount,
            cell,
            ruleIndex
        );
        if (condition.fault || condition.value <= kBooleanThreshold) continue;

        for (uint actionOffset = 0u; actionOffset < rule.actionCount; ++actionOffset) {
            const uint actionIndexValue = rule.actionOffset + actionOffset;
            const NVivoActionRecord action = program.actions[actionIndexValue];
            bool fault = false;
            const float value = action_value(program, uniforms, action, cell, actionIndexValue, fault);
            if (fault || !isfinite(value)) {
                record_expression_fault(program, cell, actionIndexValue);
                continue;
            }

            if (action.kind == 10u || action.kind == 11u) {
                const uint flag = action.kind == 11u
                    ? NVIVO_DIAGNOSTIC_PERMANENT_SHUTDOWN
                    : NVIVO_DIAGNOSTIC_REVERSIBLE_SHUTDOWN;
                atomic_fetch_or_explicit(&program.diagnostics->flags, flag, memory_order_relaxed);
                atomic_fetch_add_explicit(&program.diagnostics->shutdownCount, 1u, memory_order_relaxed);
                continue;
            }

            if ((action.flags & NVIVO_ACTION_TARGET_IS_STRING) != 0u) {
                emit_event(program, uniforms, cell, action.kind, action.targetIndex, value, 0.0f, action.flags);
                continue;
            }
            if (action.targetIndex >= uniforms.speciesCount) {
                record_expression_fault(program, cell, actionIndexValue);
                continue;
            }

            const ulong target = state_index(action.targetIndex, cell, uniforms);
            const float current = program.shadowState[target];
            float next = current;
            const float limitedRate = action.maximumRate > 0.0f
                ? min(fabs(value), action.maximumRate)
                : fabs(value);

            switch (action.kind) {
                case 0u:
                case 5u: {
                    next = value;
                    if (action.maximumRate > 0.0f) {
                        const float maximumChange = action.maximumRate * uniforms.deltaTime;
                        next = current + clamp(next - current, -maximumChange, maximumChange);
                    }
                    break;
                }
                case 1u:
                case 6u:
                    next = current + value;
                    break;
                case 2u:
                    next = current + copysign(limitedRate * uniforms.deltaTime, value);
                    break;
                case 3u:
                case 4u:
                    next = current - limitedRate * uniforms.deltaTime;
                    break;
                case 7u:
                case 8u:
                case 9u:
                    emit_event(program, uniforms, cell, action.kind, action.targetIndex, value, current, action.flags);
                    break;
                default:
                    record_expression_fault(program, cell, actionIndexValue);
                    continue;
            }
            program.shadowState[target] = next;
        }
        program.shadowRuleRefractory[refractoryIndex] = max(0.0f, rule.refractorySeconds);
    }
}

[[host_name("nvivo_validate_state")]] kernel void nvivo_validate_state(
    constant NVivoProgramArguments& program [[buffer(0)]],
    constant NVivoStepUniforms& uniforms [[buffer(1)]],
    uint gid [[thread_position_in_grid]]) {

    if (uniforms.activeCellCount == 0u) return;
    const uint speciesIndexValue = gid / uniforms.activeCellCount;
    const uint cell = gid - speciesIndexValue * uniforms.activeCellCount;
    if (speciesIndexValue >= uniforms.speciesCount || cell >= uniforms.activeCellCount) return;
    if (!cell_is_active(program, cell, uniforms)) return;

    const ulong index = state_index(speciesIndexValue, cell, uniforms);
    const float value = program.shadowState[index];
    const NVivoSpeciesRecord species = program.species[speciesIndexValue];

    if (!isfinite(value)) {
        atomic_fetch_or_explicit(
            &program.diagnostics->flags,
            NVIVO_DIAGNOSTIC_NONFINITE | NVIVO_DIAGNOSTIC_REJECT,
            memory_order_relaxed
        );
        atomic_fetch_add_explicit(&program.diagnostics->nonFiniteCount, 1u, memory_order_relaxed);
        record_first(&program.diagnostics->firstCell, cell);
        record_first(&program.diagnostics->firstSubject, speciesIndexValue);
        return;
    }

    const float tolerance = 8.0f * max(1.0f, fabs(value)) * 1.1920929e-7f;
    if (value < species.minimum - tolerance || value > species.maximum + tolerance) {
        atomic_fetch_or_explicit(
            &program.diagnostics->flags,
            NVIVO_DIAGNOSTIC_BOUND_VIOLATION | NVIVO_DIAGNOSTIC_REJECT,
            memory_order_relaxed
        );
        atomic_fetch_add_explicit(&program.diagnostics->boundViolationCount, 1u, memory_order_relaxed);
        record_first(&program.diagnostics->firstCell, cell);
        record_first(&program.diagnostics->firstSubject, speciesIndexValue);
    }
}

[[host_name("nvivo_evaluate_monitors")]] kernel void nvivo_evaluate_monitors(
    constant NVivoProgramArguments& program [[buffer(0)]],
    constant NVivoStepUniforms& uniforms [[buffer(1)]],
    uint cell [[thread_position_in_grid]]) {

    if (!cell_is_active(program, cell, uniforms)) return;

    for (uint monitorIndex = 0u; monitorIndex < uniforms.monitorCount; ++monitorIndex) {
        const NVivoMonitorRecord monitor = program.monitors[monitorIndex];
        const auto result = evaluate_expression(
            program,
            uniforms,
            monitor.expressionOffset,
            monitor.expressionCount,
            cell,
            monitorIndex
        );
        if (result.fault || result.value <= kBooleanThreshold) continue;

        atomic_fetch_or_explicit(&program.diagnostics->flags, NVIVO_DIAGNOSTIC_MONITOR, memory_order_relaxed);
        atomic_fetch_add_explicit(&program.diagnostics->monitorViolationCount, 1u, memory_order_relaxed);
        atomic_fetch_max_explicit(&program.diagnostics->maximumSeverity, monitor.severity, memory_order_relaxed);
        atomic_fetch_max_explicit(&program.diagnostics->requestedResponse, monitor.response, memory_order_relaxed);
        record_first(&program.diagnostics->firstCell, cell);
        record_first(&program.diagnostics->firstSubject, monitorIndex);

        switch (monitor.response) {
            case 0u:
            case 1u:
                break;
            case 2u:
                atomic_fetch_or_explicit(&program.diagnostics->flags, NVIVO_DIAGNOSTIC_REJECT, memory_order_relaxed);
                break;
            case 3u:
                atomic_fetch_or_explicit(&program.diagnostics->flags, NVIVO_DIAGNOSTIC_SUBSTEP, memory_order_relaxed);
                break;
            case 4u:
                atomic_fetch_or_explicit(&program.diagnostics->flags, NVIVO_DIAGNOSTIC_REVERSIBLE_SHUTDOWN, memory_order_relaxed);
                atomic_fetch_add_explicit(&program.diagnostics->shutdownCount, 1u, memory_order_relaxed);
                break;
            case 5u:
                atomic_fetch_or_explicit(&program.diagnostics->flags, NVIVO_DIAGNOSTIC_PERMANENT_SHUTDOWN, memory_order_relaxed);
                atomic_fetch_add_explicit(&program.diagnostics->shutdownCount, 1u, memory_order_relaxed);
                break;
            default:
                record_expression_fault(program, cell, monitorIndex);
                break;
        }
    }
}

[[host_name("nvivo_commit_if_valid")]] kernel void nvivo_commit_if_valid(
    constant NVivoProgramArguments& program [[buffer(0)]],
    constant NVivoStepUniforms& uniforms [[buffer(1)]],
    uint gid [[thread_position_in_grid]]) {

    const uint flags = atomic_load_explicit(&program.diagnostics->flags, memory_order_relaxed);
    const uint noCommitMask = NVIVO_DIAGNOSTIC_NONFINITE |
                              NVIVO_DIAGNOSTIC_BOUND_VIOLATION |
                              NVIVO_DIAGNOSTIC_REJECT |
                              NVIVO_DIAGNOSTIC_SUBSTEP |
                              NVIVO_DIAGNOSTIC_REVERSIBLE_SHUTDOWN |
                              NVIVO_DIAGNOSTIC_PERMANENT_SHUTDOWN |
                              NVIVO_DIAGNOSTIC_EXPRESSION_FAULT;
    const bool commit = (flags & noCommitMask) == 0u;

    const ulong stateCount = ulong(uniforms.speciesCount) * ulong(uniforms.cellCapacity);
    const ulong temporalCount = ulong(uniforms.temporalStateCount) * ulong(uniforms.cellCapacity);
    const ulong delayCount = ulong(uniforms.reactionCount) * ulong(uniforms.cellCapacity) * ulong(uniforms.delaySlotCount);
    const ulong refractoryCount = ulong(uniforms.ruleCount) * ulong(uniforms.cellCapacity);
    const ulong index = ulong(gid);

    if (commit) {
        if (index < stateCount) program.committedState[index] = program.shadowState[index];
        if (index < temporalCount) program.committedTemporalState[index] = program.shadowTemporalState[index];
        if (index < delayCount) program.committedDelayedFlux[index] = program.shadowDelayedFlux[index];
        if (index < refractoryCount) program.committedRuleRefractory[index] = program.shadowRuleRefractory[index];
    }

    if (gid == 0u) {
        uint status = NVIVO_PUBLICATION_COMMITTED;
        uint shutdownState = program.publication->shutdownState;
        if ((flags & NVIVO_DIAGNOSTIC_PERMANENT_SHUTDOWN) != 0u) {
            status = NVIVO_PUBLICATION_PERMANENT_SHUTDOWN;
            shutdownState = 2u;
        } else if ((flags & NVIVO_DIAGNOSTIC_REVERSIBLE_SHUTDOWN) != 0u) {
            status = NVIVO_PUBLICATION_REVERSIBLE_SHUTDOWN;
            shutdownState = max(shutdownState, 1u);
        } else if ((flags & NVIVO_DIAGNOSTIC_SUBSTEP) != 0u) {
            status = NVIVO_PUBLICATION_SUBSTEP_REQUIRED;
        } else if (!commit) {
            status = NVIVO_PUBLICATION_REJECTED;
        }

        if (commit) {
            program.publication->committedStepLow = uniforms.logicalStepLow;
            program.publication->committedStepHigh = uniforms.logicalStepHigh;
            program.publication->stateVersion += 1u;
        } else {
            atomic_store_explicit(program.eventCount, 0u, memory_order_relaxed);
        }
        program.publication->status = status;
        program.publication->diagnosticFlags = flags;
        program.publication->shutdownState = shutdownState;
        program.publication->activeCellCount = uniforms.activeCellCount;
        program.publication->eventCount = commit
            ? min(atomic_load_explicit(program.eventCount, memory_order_relaxed), uniforms.eventCapacity)
            : 0u;
    }
}

[[host_name("nvivo_diffuse_csr")]] kernel void nvivo_diffuse_csr(
    device const uint* neighborOffsets [[buffer(0)]],
    device const NVivoSpatialNeighbor* neighbors [[buffer(1)]],
    device const float* inputState [[buffer(2)]],
    device float* outputState [[buffer(3)]],
    constant NVivoSpatialUniforms& uniforms [[buffer(4)]],
    uint gid [[thread_position_in_grid]]) {

    if (uniforms.nodeCount == 0u) return;
    const uint species = gid / uniforms.nodeCount;
    const uint node = gid - species * uniforms.nodeCount;
    if (species >= uniforms.speciesCount || node >= uniforms.nodeCount) return;

    const ulong index = ulong(species) * ulong(uniforms.nodeCapacity) + ulong(node);
    const float current = inputState[index];
    float derivative = 0.0f;
    const uint begin = neighborOffsets[node];
    const uint end = neighborOffsets[node + 1u];
    for (uint neighborIndex = begin; neighborIndex < end; ++neighborIndex) {
        const NVivoSpatialNeighbor neighbor = neighbors[neighborIndex];
        if (neighbor.nodeIndex >= uniforms.nodeCount || neighbor.speciesIndex != species) continue;
        const ulong otherIndex = ulong(species) * ulong(uniforms.nodeCapacity) + ulong(neighbor.nodeIndex);
        const float other = inputState[otherIndex];
        derivative += neighbor.conductance * (neighbor.partitionCoefficient * other - current);
    }
    outputState[index] = max(uniforms.minimumValue, current + uniforms.deltaTime * derivative);
}

[[host_name("nvivo_reduce_sum_stage1")]] kernel void nvivo_reduce_sum_stage1(
    device const float* values [[buffer(0)]],
    device float* partials [[buffer(1)]],
    constant NVivoReductionUniforms& uniforms [[buffer(2)]],
    uint gid [[thread_position_in_grid]],
    uint lane [[thread_index_in_threadgroup]],
    uint group [[threadgroup_position_in_grid]],
    uint threadsPerGroup [[threads_per_threadgroup]]) {

    threadgroup float shared[256];
    float value = 0.0f;
    if (gid < uniforms.valueCount) value = values[ulong(gid) * ulong(uniforms.stride)];
    shared[lane] = value;
    threadgroup_barrier(mem_flags::mem_threadgroup);

    for (uint width = threadsPerGroup / 2u; width > 0u; width >>= 1u) {
        if (lane < width) shared[lane] += shared[lane + width];
        threadgroup_barrier(mem_flags::mem_threadgroup);
    }
    if (lane == 0u) partials[uniforms.outputOffset + group] = shared[0];
}

} // namespace numivivo
