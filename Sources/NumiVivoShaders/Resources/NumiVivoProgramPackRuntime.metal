#include <metal_stdlib>
using namespace metal;

// ProgramPack v1 / VivoRuntimeCommandABI v1. These buffer bindings match
// VivoTransactionalMolecularRuntime, not the historical argument-buffer ABI.
namespace vivo_pack_runtime {
struct Command {
    uint abiVersion, fidelity, laneCount, speciesCount;
    uint parameterCount, parameterStride, reactionCount, expressionCount;
    uint ruleCount, monitorCount, temporalStateCount, stepIndex;
    uint substepIndex, flags, voxelCount, eventCapacity;
    float dt, absoluteTime, minimumDt, maximumDt;
    uint gridX, gridY, gridZ, boundaryMode;
    float spacingX, spacingY, spacingZ, reservedFloat;
    ulong seedLow, seedHigh;
};
struct Status {
    atomic_uint flags, violationCount, firstLane, firstSpecies;
    atomic_uint firstMonitor, response, requiredSubsteps, eventCount;
    atomic_uint eventDropped, maximumViolation, maximumRate, reservedAtomic;
    uint reserved0, reserved1, reserved2, reserved3;
};
struct Species { uint name, compartment, unit, flags; float initial, minimum, maximum, reserved; };
struct Term { uint species; short coefficient; ushort role; };
struct Reaction {
    uint name, compartment, reactantOffset, productOffset;
    uint reactantCount, productCount, parameterOffset, parameterCount;
    uint expressionOffset, expressionCount, law, flags;
    float delaySeconds, characteristicRate; uint cohort, gateOffset;
};
struct Instruction { ushort opcode, flags; uint operand; float immediate; uint auxiliary; };
struct Action { uint target, expressionOffset, expressionCount, kind; float value, maximumRate; uint unit, flags; };
struct Rule { uint name, conditionOffset, conditionCount, actionOffset; uint actionCount; int priority; float refractory; uint temporalOffset; };
struct Monitor { uint name, expressionOffset, expressionCount, message; uint severity, response, temporalOffset, flags; };
struct Incidence { uint reaction; short coefficient; ushort reserved16; uint reserved0, reserved1; };
struct Update { uint species, lane, mode; float value; };
struct Publication { uint species, lane, output, flags; };
struct Transport { float diffusion, permeability, decay; uint flags; };
struct Event { uint lane, step, kind, subject; float value, time; uint flags, reserved; };
static_assert(sizeof(Command) == 128, "command ABI");
static_assert(sizeof(Status) == 64, "status ABI");
static_assert(sizeof(Species) == 32 && sizeof(Reaction) == 64, "ProgramPack records");
static_assert(sizeof(Term) == 8 && sizeof(Instruction) == 16 && sizeof(Incidence) == 16, "sparse records");
static_assert(sizeof(Action) == 32 && sizeof(Rule) == 32 && sizeof(Monitor) == 32, "control records");
static_assert(sizeof(Event) == 32, "event ABI");
constant uint invalid = 0xffffffffu;
constant uint nonfinite = 1u, bounds = 2u, monitorFailure = 4u, reject = 8u;
constant uint substep = 16u, reversible = 32u, permanent = 64u;
constant uint randomExhausted = 128u, eventOverflow = 256u, invalidRate = 512u, abiMismatch = 1024u;

inline uint u32(device const uchar* bytes, ulong offset) {
    return uint(bytes[offset]) | (uint(bytes[offset + 1]) << 8) |
           (uint(bytes[offset + 2]) << 16) | (uint(bytes[offset + 3]) << 24);
}
inline ulong u64(device const uchar* bytes, ulong offset) {
    return ulong(u32(bytes, offset)) | (ulong(u32(bytes, offset + 4)) << 32);
}
inline device const uchar* section(device const uchar* pack, uint kind) {
    uint count = u32(pack, 28);
    for (uint i = 0; i < count; ++i) {
        ulong descriptor = 128ul + ulong(i) * 72ul;
        if (u32(pack, descriptor) == kind) return pack + u64(pack, descriptor + 8);
    }
    return nullptr; // Host preflight requires all runtime sections.
}
struct Tables {
    device const Species* species;
    device const uint* parameterIndices;
    device const Term* terms;
    device const Reaction* reactions;
    device const Instruction* instructions;
    device const Action* actions;
    device const Rule* rules;
    device const Monitor* monitors;
    device const uint* offsets;
    device const Incidence* incidence;
};
inline Tables tables(device const uchar* pack) {
    return {reinterpret_cast<device const Species*>(section(pack, 2)),
            reinterpret_cast<device const uint*>(section(pack, 4)),
            reinterpret_cast<device const Term*>(section(pack, 5)),
            reinterpret_cast<device const Reaction*>(section(pack, 6)),
            reinterpret_cast<device const Instruction*>(section(pack, 7)),
            reinterpret_cast<device const Action*>(section(pack, 8)),
            reinterpret_cast<device const Rule*>(section(pack, 9)),
            reinterpret_cast<device const Monitor*>(section(pack, 10)),
            reinterpret_cast<device const uint*>(section(pack, 12)),
            reinterpret_cast<device const Incidence*>(section(pack, 13))};
}
inline ulong ix(uint species, uint lane, constant Command& cmd) {
    return ulong(species) * ulong(cmd.laneCount) + ulong(lane);
}
inline void fault(device Status& status, uint flags, uint lane, uint species = invalid, float amount = 0) {
    atomic_fetch_or_explicit(&status.flags, flags, memory_order_relaxed);
    atomic_fetch_add_explicit(&status.violationCount, 1u, memory_order_relaxed);
    atomic_fetch_min_explicit(&status.firstLane, lane, memory_order_relaxed);
    atomic_fetch_min_explicit(&status.firstSpecies, species, memory_order_relaxed);
    if (isfinite(amount)) atomic_fetch_max_explicit(&status.maximumViolation, as_type<uint>(fabs(amount)), memory_order_relaxed);
}
inline bool validCommand(constant Command& cmd, device Status& status) {
    bool valid = cmd.abiVersion == 1 && cmd.laneCount > 0 && cmd.voxelCount > 0 &&
                 cmd.parameterStride > 0 && isfinite(cmd.dt) && cmd.dt > 0;
    if (!valid) fault(status, abiMismatch | reject, 0);
    return valid;
}
inline float parameter(uint index, uint lane, device const float* values, constant Command& cmd) {
    return values[ulong(index) * cmd.parameterStride + lane / cmd.voxelCount];
}

inline float expression(Tables t, uint start, uint count, uint lane,
                        device const float* state, device const float* parameters,
                        device float2* temporal, constant Command& cmd,
                        float time, thread bool& valid) {
    thread float stack[256];
    uint depth = 0;
    if (start >= cmd.expressionCount) { valid = false; return 0; }
    uint available = cmd.expressionCount - start;
    uint limit = count == 0 ? min(available, 4096u) : min(available, count);
    for (uint i = 0; i < limit; ++i) {
        Instruction instruction = t.instructions[start + i];
        uint op = instruction.opcode;
        if (op == 255) {
            valid = valid && depth == 1 && isfinite(stack[0]);
            return valid ? stack[0] : 0;
        }
        if (op <= 2) {
            if (depth >= 256) { valid = false; return 0; }
            float value = 0;
            if (op == 0) value = instruction.immediate;
            else if (op == 1) {
                if (instruction.operand >= cmd.speciesCount) { valid = false; return 0; }
                value = state[ix(instruction.operand, lane, cmd)];
            } else if ((instruction.flags & 1u) != 0) value = time;
            else {
                if (instruction.operand >= cmd.parameterCount) { valid = false; return 0; }
                value = parameter(instruction.operand, lane, parameters, cmd);
            }
            if (!isfinite(value)) { valid = false; return 0; }
            stack[depth++] = value;
            continue;
        }
        if (depth == 0) { valid = false; return 0; }
        if (op == 3) { stack[depth - 1] = stack[depth - 1] > 0.5f ? 0 : 1; continue; }
        if (op >= 19 && op <= 22) {
            // Temporal expressions are legal in serial per-lane rules/monitors,
            // not in concurrently evaluated reaction-rate or gate expressions.
            if (temporal == nullptr || instruction.auxiliary >= cmd.temporalStateCount) { valid = false; return 0; }
            ulong position = ix(instruction.auxiliary, lane, cmd);
            float prior = temporal[position].x;
            bool condition = stack[depth - 1] > 0.5f;
            float next = prior, result = 0;
            if (op == 19) { next = condition ? prior + cmd.dt : 0; result = next >= instruction.immediate; }
            else if (op == 20) { next = condition ? instruction.immediate : max(0.0f, prior - cmd.dt); result = next > 0; }
            else if (op == 21) { result = condition && prior <= 0.5f; next = condition; }
            else { result = !condition && prior > 0.5f; next = condition; }
            temporal[position] = float2(next, 0);
            stack[depth - 1] = result;
            continue;
        }
        if (op == 18) {
            if (depth < 3) { valid = false; return 0; }
            float high = stack[--depth], low = stack[--depth];
            if (low > high) { valid = false; return 0; }
            stack[depth - 1] = clamp(stack[depth - 1], low, high);
            continue;
        }
        if (depth < 2) { valid = false; return 0; }
        float right = stack[--depth];
        thread float& left = stack[depth - 1];
        switch (op) {
            case 4: left = left > 0.5f && right > 0.5f; break;
            case 5: left = left > 0.5f || right > 0.5f; break;
            case 6: left = left > right; break;
            case 7: left = left >= right; break;
            case 8: left = left < right; break;
            case 9: left = left <= right; break;
            case 10: left = left == right; break;
            case 11: left = left != right; break;
            case 12: left += right; break;
            case 13: left -= right; break;
            case 14: left *= right; break;
            case 15: if (right == 0) { valid = false; return 0; } left /= right; break;
            case 16: left = min(left, right); break;
            case 17: left = max(left, right); break;
            default: valid = false; return 0;
        }
        if (!isfinite(left)) { valid = false; return 0; }
    }
    valid = false;
    return 0;
}
inline float massProduct(Tables t, uint offset, uint count, uint lane, device const float* state,
                         constant Command& cmd, thread bool& valid) {
    float result = 1;
    for (uint i = 0; i < count; ++i) {
        Term term = t.terms[offset + i];
        if (term.species >= cmd.speciesCount || term.coefficient <= 0) { valid = false; return 0; }
        float value = state[ix(term.species, lane, cmd)];
        if (!isfinite(value) || value < 0) { valid = false; return 0; }
        result *= pow(value, float(term.coefficient));
    }
    return result;
}
inline float kineticParameter(Tables t, Reaction r, uint position, uint lane,
                              device const float* parameters, constant Command& cmd, thread bool& valid) {
    if (position >= r.parameterCount) { valid = false; return 0; }
    uint index = t.parameterIndices[r.parameterOffset + position];
    if (index >= cmd.parameterCount) { valid = false; return 0; }
    float value = parameter(index, lane, parameters, cmd);
    if (!isfinite(value)) valid = false;
    return value;
}
inline float2 hazards(Tables t, Reaction r, uint lane, device const float* state,
                      device const float* parameters, constant Command& cmd,
                      float time, thread bool& valid) {
    if ((r.flags & 4u) != 0 || r.delaySeconds != 0) { valid = false; return 0; }
    if ((r.flags & 2u) != 0) {
        float gate = expression(t, r.gateOffset, 0, lane, state, parameters, nullptr, cmd, time, valid);
        if (!valid || gate <= 0.5f) return 0;
    }
    if (r.law == 255) {
        float value = expression(t, r.expressionOffset, r.expressionCount, lane, state, parameters, nullptr, cmd, time, valid);
        if (value < 0) valid = false;
        return float2(value, 0);
    }
    float k = kineticParameter(t, r, 0, lane, parameters, cmd, valid);
    if (k < 0) { valid = false; return 0; }
    if (r.law == 0) return float2(k, 0);
    if (r.law == 1 || r.law == 8) return float2(k * massProduct(t, r.reactantOffset, r.reactantCount, lane, state, cmd, valid), 0);
    if (r.law == 5) {
        float reverse = kineticParameter(t, r, 1, lane, parameters, cmd, valid);
        if (reverse < 0) valid = false;
        return float2(k * massProduct(t, r.reactantOffset, r.reactantCount, lane, state, cmd, valid),
                      reverse * massProduct(t, r.productOffset, r.productCount, lane, state, cmd, valid));
    }
    if (r.reactantCount == 0) { valid = false; return 0; }
    uint aIndex = t.terms[r.reactantOffset].species;
    float a = state[ix(aIndex, lane, cmd)];
    if (!isfinite(a) || a < 0) { valid = false; return 0; }
    if (r.law >= 2 && r.law <= 4) {
        float halfSaturation = kineticParameter(t, r, 1, lane, parameters, cmd, valid);
        if (!(halfSaturation > 0)) { valid = false; return 0; }
        if (r.law == 4) return float2(k * (a / (halfSaturation + a)), 0);
        float exponent = kineticParameter(t, r, 2, lane, parameters, cmd, valid);
        if (!(exponent > 0)) { valid = false; return 0; }
        float activation = 0;
        if (a > 0) {
            float logit = exponent * (log(a) - log(halfSaturation));
            activation = logit >= 0 ? 1 / (1 + exp(-logit)) : exp(logit) / (1 + exp(logit));
        }
        return float2(k * (r.law == 2 ? activation : 1 - activation), 0);
    }
    if (r.productCount == 0) { valid = false; return 0; }
    float b = state[ix(t.terms[r.productOffset].species, lane, cmd)];
    if (!isfinite(b) || b < 0) { valid = false; return 0; }
    if (r.law == 6) return float2(k * a, k * b);
    if (r.law == 7) {
        float halfSaturation = kineticParameter(t, r, 1, lane, parameters, cmd, valid);
        if (!(halfSaturation > 0)) { valid = false; return 0; }
        return float2(k * (a / (halfSaturation + a)), k * (b / (halfSaturation + b)));
    }
    valid = false;
    return 0;
}
inline float derivative(Tables t, uint species, uint lane, device const float* state,
                        device const float* parameters, constant Command& cmd,
                        float time, device Status& status) {
    float sum = 0, correction = 0;
    for (uint i = t.offsets[species]; i < t.offsets[species + 1]; ++i) {
        Incidence edge = t.incidence[i];
        if (edge.reaction >= cmd.reactionCount) { fault(status, abiMismatch | reject, lane, species); return 0; }
        bool valid = true;
        float2 pair = hazards(t, t.reactions[edge.reaction], lane, state, parameters, cmd, time, valid);
        if (!valid || !all(isfinite(pair))) { fault(status, invalidRate | reject, lane, species); return 0; }
        atomic_fetch_max_explicit(&status.maximumRate, as_type<uint>(max(pair.x, pair.y)), memory_order_relaxed);
        float term = float(edge.coefficient) * (pair.x - pair.y) - correction;
        float next = sum + term; correction = (next - sum) - term; sum = next;
    }
    return sum;
}

[[host_name("nvivo_clear_status")]] kernel void nvivo_clear_status(device Status& s [[buffer(0)]], uint gid [[thread_position_in_grid]]) {
    if (gid != 0) return;
    atomic_store_explicit(&s.flags, 0, memory_order_relaxed);
    atomic_store_explicit(&s.violationCount, 0, memory_order_relaxed);
    atomic_store_explicit(&s.firstLane, invalid, memory_order_relaxed);
    atomic_store_explicit(&s.firstSpecies, invalid, memory_order_relaxed);
    atomic_store_explicit(&s.firstMonitor, invalid, memory_order_relaxed);
    atomic_store_explicit(&s.response, 0, memory_order_relaxed);
    atomic_store_explicit(&s.requiredSubsteps, 1, memory_order_relaxed);
    atomic_store_explicit(&s.eventCount, 0, memory_order_relaxed);
    atomic_store_explicit(&s.eventDropped, 0, memory_order_relaxed);
    atomic_store_explicit(&s.maximumViolation, 0, memory_order_relaxed);
    atomic_store_explicit(&s.maximumRate, 0, memory_order_relaxed);
    atomic_store_explicit(&s.reservedAtomic, 0, memory_order_relaxed);
    s.reserved0 = s.reserved1 = s.reserved2 = s.reserved3 = 0;
}
[[host_name("nvivo_prepare_transaction")]] kernel void nvivo_prepare_transaction(
    device const float* current [[buffer(0)]], device float* base [[buffer(1)]],
    device float* candidate [[buffer(2)]], device const float2* temporalCurrent [[buffer(3)]],
    device float2* temporalCandidate [[buffer(4)]], device Status& status [[buffer(5)]],
    constant Command& cmd [[buffer(6)]], uint gid [[thread_position_in_grid]]) {
    if (!validCommand(cmd, status)) return;
    if (ulong(gid) < ulong(cmd.speciesCount) * cmd.laneCount) base[gid] = candidate[gid] = current[gid];
    if (ulong(gid) < ulong(cmd.temporalStateCount) * cmd.laneCount) temporalCandidate[gid] = temporalCurrent[gid];
}
[[host_name("nvivo_apply_coupling_updates")]] kernel void nvivo_apply_coupling_updates(
    device float* state [[buffer(0)]], device const Update* updates [[buffer(1)]],
    device Status& status [[buffer(2)]], constant Command& cmd [[buffer(3)]],
    constant uint& count [[buffer(4)]], uint gid [[thread_position_in_grid]]) {
    if (gid >= count) return;
    Update u = updates[gid];
    if (u.species >= cmd.speciesCount || u.lane >= cmd.laneCount || !isfinite(u.value) || u.mode > 2) {
        fault(status, abiMismatch | reject, u.lane, u.species); return;
    }
    ulong index = ix(u.species, u.lane, cmd);
    state[index] = u.mode == 0 ? u.value : (u.mode == 1 ? state[index] + u.value : 0.5f * state[index] + 0.5f * u.value);
}
[[host_name("nvivo_f1_heun_predict")]] kernel void nvivo_f1_heun_predict(
    device const float* state [[buffer(0)]], device float* stage [[buffer(1)]], device float* k1 [[buffer(2)]],
    device const float* parameters [[buffer(3)]], device const uchar* pack [[buffer(4)]],
    device Status& status [[buffer(5)]], constant Command& cmd [[buffer(6)]], uint gid [[thread_position_in_grid]]) {
    uint species = gid / cmd.laneCount, lane = gid % cmd.laneCount;
    if (species >= cmd.speciesCount) return;
    Tables t = tables(pack);
    if ((t.species[species].flags & 2u) != 0) { stage[gid] = state[gid]; k1[gid] = 0; return; }
    float d = derivative(t, species, lane, state, parameters, cmd, cmd.absoluteTime - cmd.dt, status);
    k1[gid] = d; stage[gid] = state[gid] + cmd.dt * d;
    if (!isfinite(stage[gid])) fault(status, nonfinite | reject, lane, species);
}
[[host_name("nvivo_f1_heun_correct")]] kernel void nvivo_f1_heun_correct(
    device const float* state [[buffer(0)]], device const float* stage [[buffer(1)]], device const float* k1 [[buffer(2)]],
    device float* candidate [[buffer(3)]], device const float* parameters [[buffer(4)]],
    device const uchar* pack [[buffer(5)]], device Status& status [[buffer(6)]],
    constant Command& cmd [[buffer(7)]], uint gid [[thread_position_in_grid]]) {
    uint species = gid / cmd.laneCount, lane = gid % cmd.laneCount;
    if (species >= cmd.speciesCount) return;
    Tables t = tables(pack);
    if ((t.species[species].flags & 2u) != 0) { candidate[gid] = state[gid]; return; }
    float d = derivative(t, species, lane, stage, parameters, cmd, cmd.absoluteTime, status);
    candidate[gid] = state[gid] + 0.5f * cmd.dt * (k1[gid] + d);
}

inline uint4 philox(uint4 c, uint2 key) {
    for (uint r = 0; r < 10; ++r) {
        uint hi0 = mulhi(0xD2511F53u, c.x), hi1 = mulhi(0xCD9E8D57u, c.z);
        c = uint4(hi1 ^ c.y ^ key.x, 0xCD9E8D57u * c.z, hi0 ^ c.w ^ key.y, 0xD2511F53u * c.x);
        key += uint2(0x9E3779B9u, 0xBB67AE85u);
    }
    return c;
}
inline float random(uint lane, uint reaction, thread uint& draw, constant Command& cmd) {
    uint2 key = uint2(uint(cmd.seedLow) ^ uint(cmd.seedHigh >> 32), uint(cmd.seedLow >> 32) ^ uint(cmd.seedHigh));
    key ^= uint2(cmd.substepIndex * 0x9E3779B9u, cmd.substepIndex * 0xBB67AE85u);
    uint4 values = philox(uint4(lane, reaction, cmd.stepIndex, draw >> 2), key);
    return (float(values[(draw++) & 3u] >> 9) + 0.5f) * 0x1p-23f;
}
inline uint poisson(float lambda, uint lane, uint reaction, thread uint& draw,
                    constant Command& cmd, thread bool& valid) {
    if (lambda == 0) return 0;
    if (!isfinite(lambda) || lambda < 0 || lambda > 1024) { valid = false; return 0; }
    uint pieces = max(1u, uint(ceil(lambda / 8.0f))), result = 0;
    float threshold = exp(-lambda / float(pieces));
    for (uint piece = 0; piece < pieces; ++piece) {
        float product = 1;
        bool completed = false;
        for (uint k = 0; k < 128; ++k) {
            product *= random(lane, reaction, draw, cmd);
            if (product <= threshold) { result += k; completed = true; break; }
        }
        if (!completed) { valid = false; return 0; }
    }
    return result;
}
[[host_name("nvivo_f2_sample_reactions")]] kernel void nvivo_f2_sample_reactions(
    device const float* state [[buffer(0)]], device const float* parameters [[buffer(1)]],
    device int* events [[buffer(2)]], device const uchar* pack [[buffer(3)]],
    device Status& status [[buffer(4)]], constant Command& cmd [[buffer(5)]], uint gid [[thread_position_in_grid]]) {
    uint reaction = gid / cmd.laneCount, lane = gid % cmd.laneCount;
    if (reaction >= cmd.reactionCount) return;
    Tables t = tables(pack); Reaction r = t.reactions[reaction];
    bool valid = r.law != 255;
    float2 pair = hazards(t, r, lane, state, parameters, cmd, cmd.absoluteTime - cmd.dt, valid);
    // ProgramPack's declared rate law defines the two nonnegative intensities.
    // Reversible reactions sample forward and reverse separately, not Poisson
    // of their signed difference. General combinatorial laws use the hybrid API.
    if (!valid || !all(isfinite(pair)) || any(pair < 0)) { events[gid] = 0; fault(status, invalidRate | reject, lane); return; }
    uint draw = 0;
    uint forward = poisson(pair.x * cmd.dt, lane, reaction, draw, cmd, valid);
    uint backward = poisson(pair.y * cmd.dt, lane, reaction, draw, cmd, valid);
    if (!valid) { events[gid] = 0; fault(status, randomExhausted | substep | reject, lane); return; }
    events[gid] = int(forward) - int(backward);
    atomic_fetch_max_explicit(&status.maximumRate, as_type<uint>(max(pair.x, pair.y)), memory_order_relaxed);
}
[[host_name("nvivo_f2_apply_reactions")]] kernel void nvivo_f2_apply_reactions(
    device const float* state [[buffer(0)]], device float* candidate [[buffer(1)]],
    device const int* events [[buffer(2)]], device const uchar* pack [[buffer(3)]],
    device Status& status [[buffer(4)]], constant Command& cmd [[buffer(5)]], uint gid [[thread_position_in_grid]]) {
    uint species = gid / cmd.laneCount, lane = gid % cmd.laneCount;
    if (species >= cmd.speciesCount) return;
    Tables t = tables(pack); Species metadata = t.species[species];
    candidate[gid] = state[gid];
    if ((metadata.flags & 2u) != 0 || t.offsets[species] == t.offsets[species + 1]) return;
    if ((metadata.flags & 32u) == 0 || state[gid] < 0 || state[gid] > 16777216.0f || state[gid] != rint(state[gid])) {
        fault(status, invalidRate | reject, lane, species); return;
    }
    long value = long(state[gid]);
    for (uint i = t.offsets[species]; i < t.offsets[species + 1]; ++i) {
        Incidence edge = t.incidence[i];
        value += long(edge.coefficient) * long(events[ix(edge.reaction, lane, cmd)]);
    }
    if (value < 0 || value > 16777216l) { fault(status, bounds | reject, lane, species); return; }
    candidate[gid] = float(value);
}

[[host_name("nvivo_f3_transport")]] kernel void nvivo_f3_transport(
    device const float* source [[buffer(0)]], device float* destination [[buffer(1)]],
    device const Transport* transports [[buffer(2)]], device const float4* velocities [[buffer(3)]],
    device const float* fractions [[buffer(4)]], device Status& status [[buffer(5)]],
    constant Command& cmd [[buffer(6)]], device const uchar* pack [[buffer(8)]], uint gid [[thread_position_in_grid]]) {
    uint species = gid / cmd.laneCount, lane = gid % cmd.laneCount;
    if (species >= cmd.speciesCount) return;
    Tables t = tables(pack);
    destination[gid] = source[gid];
    if ((t.species[species].flags & 2u) != 0) return;
    Transport transport = transports[species];
    if (transport.permeability != 0) { fault(status, abiMismatch | reject, lane, species); return; }
    float alpha = fractions[lane];
    if (alpha <= 0) return; // Inactive cut cell; all adjacent face fractions are zero.
    uint local = lane % cmd.voxelCount, environment = lane - local;
    uint3 dimensions = uint3(cmd.gridX, cmd.gridY, cmd.gridZ);
    uint3 position = uint3(local % cmd.gridX, (local / cmd.gridX) % cmd.gridY, local / (cmd.gridX * cmd.gridY));
    float3 spacing = float3(cmd.spacingX, cmd.spacingY, cmd.spacingZ);
    float value = source[gid], derivativeValue = -transport.decay * value;
    float outgoing = transport.decay;
    for (uint axis = 0; axis < 3; ++axis) {
        for (int direction = -1; direction <= 1; direction += 2) {
            int coordinate = int(position[axis]) + direction;
            bool outside = coordinate < 0 || coordinate >= int(dimensions[axis]);
            if (outside && cmd.boundaryMode == 0) continue;
            uint3 neighbor = position;
            if (outside && cmd.boundaryMode == 1) coordinate = coordinate < 0 ? int(dimensions[axis]) - 1 : 0;
            bool absorbing = outside && cmd.boundaryMode == 2;
            uint otherLane = lane;
            if (!absorbing) {
                neighbor[axis] = uint(coordinate);
                otherLane = environment + neighbor.x + cmd.gridX * (neighbor.y + cmd.gridY * neighbor.z);
            }
            float face = absorbing ? alpha : min(alpha, fractions[otherLane]);
            float other = absorbing ? 0 : source[ix(species, otherLane, cmd)];
            float velocity = absorbing ? velocities[lane][axis] : 0.5f * (velocities[lane][axis] + velocities[otherLane][axis]);
            float normalVelocity = float(direction) * velocity;
            float h = spacing[axis];
            float diffusion = face * transport.diffusion / (alpha * h * h);
            float advective = face * normalVelocity / (alpha * h);
            derivativeValue += diffusion * (other - value) - advective * (normalVelocity >= 0 ? value : other);
            outgoing += diffusion + max(advective, 0.0f);
        }
    }
    if (!isfinite(outgoing) || cmd.dt * outgoing > 0.9f) {
        fault(status, substep | reject, lane, species);
        atomic_fetch_max_explicit(&status.requiredSubsteps, 2u, memory_order_relaxed);
        return;
    }
    float next = value + cmd.dt * derivativeValue;
    if (!isfinite(next)) { fault(status, nonfinite | reject, lane, species); return; }
    // Transport of discrete molecular counts needs a stochastic hopping model,
    // not fractional FP32 concentrations. Host preflight rejects that pairing.
    destination[gid] = next;
}

inline void event(device Event* events, device Status& status, uint lane, uint kind,
                  uint subject, float value, uint flags, constant Command& cmd) {
    uint slot = atomic_fetch_add_explicit(&status.eventCount, 1u, memory_order_relaxed);
    if (slot >= cmd.eventCapacity) {
        atomic_fetch_add_explicit(&status.eventDropped, 1u, memory_order_relaxed);
        fault(status, eventOverflow | reject, lane); return;
    }
    events[slot] = {lane, cmd.stepIndex, kind, subject, value, cmd.absoluteTime, flags, 0};
}
[[host_name("nvivo_execute_rules")]] kernel void nvivo_execute_rules(
    device float* state [[buffer(0)]], device float2* temporal [[buffer(1)]],
    device const float* parameters [[buffer(2)]], device const uchar* pack [[buffer(3)]],
    device Event* events [[buffer(4)]], device Status& status [[buffer(5)]],
    constant Command& cmd [[buffer(6)]], uint lane [[thread_position_in_grid]]) {
    if (lane >= cmd.laneCount) return;
    Tables t = tables(pack);
    // One thread owns the complete ordered rule sequence for its lane. Parallel
    // rule writers would race on shared state and ignore compiled priority.
    for (uint index = 0; index < cmd.ruleCount; ++index) {
        Rule rule = t.rules[index];
        if (rule.refractory != 0) { fault(status, abiMismatch | reject, lane); return; }
        bool valid = true;
        float condition = expression(t, rule.conditionOffset, rule.conditionCount, lane, state, parameters, temporal, cmd, cmd.absoluteTime, valid);
        if (!valid) { fault(status, reject, lane); return; }
        if (condition <= 0.5f) continue;
        for (uint j = 0; j < rule.actionCount; ++j) {
            Action action = t.actions[rule.actionOffset + j];
            float value = action.expressionCount == 0 ? action.value : expression(t, action.expressionOffset, action.expressionCount, lane, state, parameters, temporal, cmd, cmd.absoluteTime, valid);
            if (!valid || !isfinite(value)) { fault(status, reject | nonfinite, lane); return; }
            if (action.kind == 10 || action.kind == 11) {
                atomic_fetch_or_explicit(&status.flags, action.kind == 11 ? permanent : reversible, memory_order_relaxed);
                atomic_fetch_max_explicit(&status.response, action.kind == 11 ? 5u : 4u, memory_order_relaxed);
                return;
            }
            if ((action.flags & 1u) != 0) { event(events, status, lane, action.kind, action.target, value, action.flags, cmd); continue; }
            if (action.target >= cmd.speciesCount || (t.species[action.target].flags & 2u) != 0) {
                fault(status, abiMismatch | reject, lane); return;
            }
            ulong target = ix(action.target, lane, cmd);
            float current = state[target], next = current;
            float rate = action.maximumRate > 0 ? min(fabs(value), action.maximumRate) : fabs(value);
            switch (action.kind) {
                case 0: case 5:
                    next = action.maximumRate > 0 ? current + clamp(value - current, -action.maximumRate * cmd.dt, action.maximumRate * cmd.dt) : value;
                    break;
                case 1: case 6: next = current + value; break;
                case 2: next = current + copysign(rate * cmd.dt, value); break;
                case 3: case 4: next = current - rate * cmd.dt; break;
                case 7: case 8: case 9: event(events, status, lane, action.kind, action.target, value, action.flags, cmd); break;
                default: fault(status, abiMismatch | reject, lane); return;
            }
            state[target] = next;
        }
    }
}
[[host_name("nvivo_evaluate_monitors")]] kernel void nvivo_evaluate_monitors(
    device const float* state [[buffer(0)]], device float2* temporal [[buffer(1)]],
    device const float* parameters [[buffer(2)]], device const uchar* pack [[buffer(3)]],
    device Event* events [[buffer(4)]], device Status& status [[buffer(5)]],
    constant Command& cmd [[buffer(6)]], uint lane [[thread_position_in_grid]]) {
    if (lane >= cmd.laneCount) return;
    Tables t = tables(pack);
    for (uint i = 0; i < cmd.monitorCount; ++i) {
        Monitor monitor = t.monitors[i];
        bool valid = true;
        float value = expression(t, monitor.expressionOffset, monitor.expressionCount, lane, state, parameters, temporal, cmd, cmd.absoluteTime, valid);
        if (!valid) { fault(status, reject, lane); return; }
        bool triggered = (monitor.flags & 1u) != 0 ? value > 0.5f : value <= 0.5f;
        if (!triggered) continue;
        uint flags = monitorFailure;
        if (monitor.response == 1 || monitor.response == 2) flags |= reject;
        if (monitor.response == 3) flags |= substep | reject;
        if (monitor.response == 4) flags |= reversible;
        if (monitor.response == 5) flags |= permanent;
        fault(status, flags, lane);
        atomic_fetch_min_explicit(&status.firstMonitor, i, memory_order_relaxed);
        atomic_fetch_max_explicit(&status.response, monitor.response, memory_order_relaxed);
        event(events, status, lane, 0x100u, i, value, monitor.flags, cmd);
    }
}
[[host_name("nvivo_validate_shadow")]] kernel void nvivo_validate_shadow(
    device const float* state [[buffer(0)]], device const uchar* pack [[buffer(1)]],
    device Status& status [[buffer(2)]], constant Command& cmd [[buffer(3)]], uint gid [[thread_position_in_grid]]) {
    uint species = gid / cmd.laneCount, lane = gid % cmd.laneCount;
    if (species >= cmd.speciesCount) return;
    Species metadata = tables(pack).species[species];
    float value = state[gid];
    if (!isfinite(value)) { fault(status, nonfinite | reject, lane, species); return; }
    if (value < metadata.minimum || value > metadata.maximum) {
        fault(status, bounds | reject, lane, species, max(metadata.minimum - value, value - metadata.maximum)); return;
    }
    if ((metadata.flags & 32u) != 0 && (value < 0 || value > 16777216.0f || value != rint(value))) {
        fault(status, bounds | reject, lane, species);
    }
}
[[host_name("nvivo_publish")]] kernel void nvivo_publish(
    device const float* state [[buffer(0)]], device const Publication* requests [[buffer(1)]],
    device float* output [[buffer(2)]], device Status& status [[buffer(3)]],
    constant Command& cmd [[buffer(4)]], constant uint& count [[buffer(5)]], uint gid [[thread_position_in_grid]]) {
    if (gid >= count) return;
    Publication request = requests[gid];
    if (request.species >= cmd.speciesCount || request.lane >= cmd.laneCount || request.output >= count) {
        fault(status, abiMismatch | reject, request.lane, request.species); return;
    }
    output[request.output] = state[ix(request.species, request.lane, cmd)];
}
} // namespace vivo_pack_runtime
