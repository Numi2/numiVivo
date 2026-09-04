#include <metal_stdlib>
using namespace metal;

// Executable hybrid ABI v1. Independent of the legacy ProgramPack shader ABI.
namespace vivo_hybrid {
struct Reaction {
    uint law, reactantA, reactantB, changeOffset;
    uint changeCount, reserved, reserved1, reserved2;
    float rate, parameter1, parameter2, parameter3;
};
struct Change { uint index; int delta; };
struct Cohort { uint reactionOffset; uint reactionCount; };
struct Command {
    uint laneCount, speciesCount, reactionCount, exactCohortCount;
    uint stepLow, stepHigh, stage, eventsPerDispatch;
    float dt; uint reserved; ulong seed;
};
struct Progress { float elapsed; uint draw; uint events; uint done; };
struct Status {
    atomic_uint flags, firstLane, firstReaction, unfinishedExactLanes;
    atomic_uint maximumExactEvents, reserved0, reserved1, reserved2;
};
struct Publication { float continuous; uint count; uint mode; uint reserved; };
static_assert(sizeof(Reaction) == 48, "hybrid reaction ABI");
static_assert(sizeof(Command) == 48, "hybrid command ABI");
static_assert(sizeof(Progress) == 16, "hybrid progress ABI");
static_assert(sizeof(Status) == 32, "hybrid status ABI");

inline void fail(device Status& status, uint flag, uint lane, uint reaction = 0xffffffffu) {
    atomic_fetch_or_explicit(&status.flags, flag, memory_order_relaxed);
    atomic_fetch_min_explicit(&status.firstLane, lane, memory_order_relaxed);
    atomic_fetch_min_explicit(&status.firstReaction, reaction, memory_order_relaxed);
}
inline ulong ix(uint species, uint lane, constant Command& cmd) {
    return ulong(species) * ulong(cmd.laneCount) + ulong(lane);
}
inline uint4 philox(uint4 c, uint2 key) {
    for (uint round = 0; round < 10; ++round) {
        uint hi0 = mulhi(0xD2511F53u, c.x), lo0 = 0xD2511F53u * c.x;
        uint hi1 = mulhi(0xCD9E8D57u, c.z), lo1 = 0xCD9E8D57u * c.z;
        c = uint4(hi1 ^ c.y ^ key.x, lo1, hi0 ^ c.w ^ key.y, lo0);
        key += uint2(0x9E3779B9u, 0xBB67AE85u);
    }
    return c;
}
inline float uniform(uint lane, uint stream, uint domain, thread uint& draw, constant Command& cmd) {
    uint2 key = uint2(uint(cmd.seed), uint(cmd.seed >> 32));
    key ^= uint2(domain ^ cmd.stepHigh, (domain * 0x9E3779B9u) ^ (cmd.stepHigh * 0xBB67AE85u));
    uint4 block = philox(uint4(lane, stream, draw >> 2, cmd.stepLow), key);
    uint bits = block[draw & 3u];
    ++draw;
    // FP32 midpoint grid strictly inside (0,1), even for the largest UInt32.
    return (float(bits >> 9) + 0.5f) * 0x1p-23f;
}
inline float hill(float a, float half, float exponent) {
    if (a <= 0) return 0;
    float z = exponent * (log(a) - log(half));
    return z >= 0 ? 1.0f / (1.0f + exp(-z)) : exp(z) / (1.0f + exp(z));
}
inline float law(Reaction r, float a, float b, bool discrete) {
    switch (r.law) {
        case 0: return r.rate;
        case 1: return r.rate * a;
        case 2: return r.rate * a * b;
        case 3: return 0.5f * r.rate * a * (discrete ? max(a - 1.0f, 0.0f) : a);
        case 4: return r.rate * hill(a, r.parameter1, r.parameter2);
        case 5: return r.rate * (1.0f - hill(a, r.parameter1, r.parameter2));
        default: return NAN;
    }
}
inline float propensity(Reaction r, device const uint* counts, device const Change* changes,
                        uint lane, constant Command& cmd) {
    for (uint j = 0; j < r.changeCount; ++j) {
        Change c = changes[r.changeOffset + j];
        if (c.delta < 0 && ulong(counts[ix(c.index, lane, cmd)]) < ulong(-long(c.delta))) return 0;
    }
    float a = r.reactantA == 0xffffffffu ? 0 : float(counts[ix(r.reactantA, lane, cmd)]);
    float b = r.reactantB == 0xffffffffu ? 0 : float(counts[ix(r.reactantB, lane, cmd)]);
    return law(r, a, b, true);
}

kernel void nvivo_hybrid_clear(device Status& status [[buffer(17)]], uint gid [[thread_position_in_grid]]) {
    if (gid != 0) return;
    atomic_store_explicit(&status.flags, 0u, memory_order_relaxed);
    atomic_store_explicit(&status.firstLane, 0xffffffffu, memory_order_relaxed);
    atomic_store_explicit(&status.firstReaction, 0xffffffffu, memory_order_relaxed);
    atomic_store_explicit(&status.unfinishedExactLanes, 0u, memory_order_relaxed);
    atomic_store_explicit(&status.maximumExactEvents, 0u, memory_order_relaxed);
    atomic_store_explicit(&status.reserved0, 0u, memory_order_relaxed);
    atomic_store_explicit(&status.reserved1, 0u, memory_order_relaxed);
    atomic_store_explicit(&status.reserved2, 0u, memory_order_relaxed);
}
kernel void nvivo_hybrid_reset_continuation(device Status& status [[buffer(17)]], uint gid [[thread_position_in_grid]]) {
    if (gid == 0) atomic_store_explicit(&status.unfinishedExactLanes, 0u, memory_order_relaxed);
}

kernel void nvivo_hybrid_rk_rates(
    device const float* current [[buffer(0)]], device const float* stage [[buffer(1)]],
    device const Reaction* reactions [[buffer(6)]], device const uint* modes [[buffer(11)]],
    device float* flux [[buffer(12)]], device Status& status [[buffer(17)]],
    constant Command& cmd [[buffer(18)]], uint gid [[thread_position_in_grid]]) {
    uint lane = gid % cmd.laneCount, reaction = gid / cmd.laneCount;
    if (reaction >= cmd.reactionCount || modes[reaction] != 3u) return;
    device const float* source = cmd.stage == 0 ? current : stage;
    Reaction r = reactions[reaction];
    float a = r.reactantA == 0xffffffffu ? 0 : source[ix(r.reactantA, lane, cmd)];
    float b = r.reactantB == 0xffffffffu ? 0 : source[ix(r.reactantB, lane, cmd)];
    float value = law(r, a, b, false);
    if (!isfinite(a) || !isfinite(b) || a < 0 || b < 0 || !isfinite(value) || value < 0) {
        fail(status, 1u, lane, reaction); value = 0;
    }
    flux[gid] = value;
}
inline float derivative(uint species, uint lane, device const uint* offsets,
                        device const Change* incidence, device const float* flux,
                        constant Command& cmd) {
    float sum = 0, compensation = 0;
    for (uint p = offsets[species]; p < offsets[species + 1]; ++p) {
        Change edge = incidence[p];
        float term = float(edge.delta) * flux[ix(edge.index, lane, cmd)] - compensation;
        float next = sum + term;
        compensation = (next - sum) - term;
        sum = next;
    }
    return sum;
}
kernel void nvivo_hybrid_rk_predict(
    device const float* current [[buffer(0)]], device float* stage [[buffer(1)]],
    device float* k1 [[buffer(3)]], device const uint* offsets [[buffer(8)]],
    device const Change* incidence [[buffer(9)]], device const uint* modes [[buffer(10)]],
    device const float* flux [[buffer(12)]], device Status& status [[buffer(17)]],
    constant Command& cmd [[buffer(18)]], uint gid [[thread_position_in_grid]]) {
    uint lane = gid % cmd.laneCount, species = gid / cmd.laneCount;
    if (species >= cmd.speciesCount || modes[species] != 3u) return;
    float d = derivative(species, lane, offsets, incidence, flux, cmd);
    float value = current[gid] + cmd.dt * d;
    k1[gid] = d; stage[gid] = value;
    if (!isfinite(value) || value < 0) fail(status, 8u, lane);
}
kernel void nvivo_hybrid_rk_correct(
    device const float* current [[buffer(0)]], device float* candidate [[buffer(2)]],
    device const float* k1 [[buffer(3)]], device const uint* offsets [[buffer(8)]],
    device const Change* incidence [[buffer(9)]], device const uint* modes [[buffer(10)]],
    device const float* flux [[buffer(12)]], device Status& status [[buffer(17)]],
    constant Command& cmd [[buffer(18)]], uint gid [[thread_position_in_grid]]) {
    uint lane = gid % cmd.laneCount, species = gid / cmd.laneCount;
    if (species >= cmd.speciesCount || modes[species] != 3u) return;
    float d = derivative(species, lane, offsets, incidence, flux, cmd);
    float value = current[gid] + (0.5f * cmd.dt) * (k1[gid] + d);
    candidate[gid] = value;
    if (!isfinite(value) || value < 0) fail(status, 8u, lane);
}

inline float small_log_factorial(uint k) {
    float result = 0;
    for (uint j = 2; j <= k; ++j) result += log(float(j));
    return result;
}
// Bounded inversion/PTRS. No Gaussian substitution or clipping on exhaustion.
inline uint poisson(float lambda, uint lane, uint reaction, thread uint& draw,
                    constant Command& cmd, thread bool& valid) {
    if (lambda == 0) return 0;
    if (!isfinite(lambda) || lambda < 0 || lambda > 1.0e6f) { valid = false; return 0; }
    if (lambda < 10) {
        float threshold = exp(-lambda), product = 1;
        for (uint k = 0; k < 256; ++k) {
            product *= uniform(lane, reaction, 0x54415531u, draw, cmd);
            if (product <= threshold) return k;
        }
        valid = false; return 0;
    }
    float root = sqrt(lambda), b = 0.931f + 2.53f * root;
    float a = -0.059f + 0.02483f * b;
    float inverseAlpha = 1.1239f + 1.1328f / (b - 3.4f);
    float vr = 0.9277f - 3.6224f / (b - 2.0f);
    for (uint attempt = 0; attempt < 256; ++attempt) {
        float u = uniform(lane, reaction, 0x54415531u, draw, cmd) - 0.5f;
        float v = uniform(lane, reaction, 0x54415531u, draw, cmd);
        float us = 0.5f - fabs(u);
        if (us <= 0) continue;
        float proposal = floor((2 * a / us + b) * u + lambda + 0.43f);
        if (proposal < 0 || proposal > 2147483520.0f) continue;
        uint k = uint(proposal);
        if (us >= 0.07f && v <= vr) return k;
        if (us < 0.013f && v > us) continue;
        float lhs = log(v * inverseAlpha / (a / (us * us) + b));
        float logPMF;
        if (k == 0) logPMF = -lambda;
        else if (k < 32) logPMF = -lambda + float(k) * log(lambda) - small_log_factorial(k);
        else {
            // Stable deviance form near the mode instead of cancelling two
            // O(lambda log lambda) values in the Poisson log probability.
            float x = float(k), vdev = (x - lambda) / (x + lambda), dev;
            if (fabs(x - lambda) < 0.1f * (x + lambda)) {
                float term = 2 * x * vdev, series = (x - lambda) * vdev;
                float v2 = vdev * vdev;
                for (uint j = 1; j <= 12; ++j) { term *= v2; series += term / float(2 * j + 1); }
                dev = series;
            } else dev = x * log(x / lambda) + lambda - x;
            float inv = 1 / x;
            float stirling = inv / 12 - inv * inv * inv / 360;
            logPMF = -dev - stirling - 0.5f * log(6.283185307179586f * x);
        }
        if (lhs <= logPMF) return k;
    }
    valid = false; return 0;
}
kernel void nvivo_hybrid_tau_sample(
    device const uint* counts [[buffer(4)]], device const Reaction* reactions [[buffer(6)]],
    device const Change* changes [[buffer(7)]], device const uint* modes [[buffer(11)]],
    device uint* firings [[buffer(13)]], device Status& status [[buffer(17)]],
    constant Command& cmd [[buffer(18)]], uint gid [[thread_position_in_grid]]) {
    uint lane = gid % cmd.laneCount, reaction = gid / cmd.laneCount;
    if (reaction >= cmd.reactionCount || modes[reaction] != 2u) return;
    float rate = propensity(reactions[reaction], counts, changes, lane, cmd);
    if (!isfinite(rate) || rate < 0) { fail(status, 1u, lane, reaction); firings[gid] = 0; return; }
    uint draw = 0; bool valid = true;
    firings[gid] = poisson(rate * cmd.dt, lane, reaction, draw, cmd, valid);
    if (!valid) fail(status, 16u, lane, reaction);
}
kernel void nvivo_hybrid_tau_apply(
    device const uint* current [[buffer(4)]], device uint* candidate [[buffer(5)]],
    device const uint* offsets [[buffer(8)]], device const Change* incidence [[buffer(9)]],
    device const uint* modes [[buffer(10)]], device const uint* firings [[buffer(13)]],
    device Status& status [[buffer(17)]], constant Command& cmd [[buffer(18)]],
    uint gid [[thread_position_in_grid]]) {
    uint lane = gid % cmd.laneCount, species = gid / cmd.laneCount;
    if (species >= cmd.speciesCount || modes[species] != 2u) return;
    long value = long(current[gid]);
    for (uint p = offsets[species]; p < offsets[species + 1]; ++p) {
        Change edge = incidence[p];
        long term = long(edge.delta) * long(firings[ix(edge.index, lane, cmd)]);
        if ((term > 0 && value > 9223372036854775807l - term) ||
            (term < 0 && value < (-9223372036854775807l - 1l) - term)) {
            fail(status, 4u, lane, edge.index); return;
        }
        value += term;
    }
    if (value < 0) { fail(status, 2u, lane); return; }
    if (ulong(value) > 0xfffffffful) { fail(status, 4u, lane); return; }
    candidate[gid] = uint(value);
}

kernel void nvivo_hybrid_exact_advance(
    device uint* counts [[buffer(5)]], device const Reaction* reactions [[buffer(6)]],
    device const Change* changes [[buffer(7)]], device const Cohort* cohorts [[buffer(14)]],
    device const uint* indices [[buffer(15)]], device Progress* progress [[buffer(16)]],
    device Status& status [[buffer(17)]], constant Command& cmd [[buffer(18)]],
    uint gid [[thread_position_in_grid]]) {
    uint lane = gid % cmd.laneCount, group = gid / cmd.laneCount;
    if (group >= cmd.exactCohortCount) return;
    Progress p = progress[gid];
    if (p.done != 0) return;
    Cohort cohort = cohorts[group];
    for (uint iteration = 0; iteration < cmd.eventsPerDispatch; ++iteration) {
        float total = 0, compensation = 0;
        for (uint j = 0; j < cohort.reactionCount; ++j) {
            uint r = indices[cohort.reactionOffset + j];
            float rate = propensity(reactions[r], counts, changes, lane, cmd);
            if (!isfinite(rate) || rate < 0) { fail(status, 1u, lane, r); p.done = 1; break; }
            float term = rate - compensation, next = total + term;
            compensation = (next - total) - term; total = next;
        }
        if (p.done != 0) break;
        if (!isfinite(total)) { fail(status, 1u, lane); p.done = 1; break; }
        if (total <= 0) { p.elapsed = cmd.dt; p.done = 1; break; }
        if (p.draw > 0xfffffffbu) { fail(status, 32u, lane); p.done = 1; break; }
        float wait = -log(uniform(lane, group, 0x53534131u, p.draw, cmd)) / total;
        if (wait >= cmd.dt - p.elapsed) { p.elapsed = cmd.dt; p.done = 1; break; }
        if (!(wait > 0) || !(p.elapsed + wait > p.elapsed)) { fail(status, 32u, lane); p.done = 1; break; }
        float threshold = uniform(lane, group, 0x53534131u, p.draw, cmd) * total;
        uint selected = 0xffffffffu, lastPositive = 0xffffffffu;
        float cumulative = 0, correction = 0;
        for (uint j = 0; j < cohort.reactionCount; ++j) {
            uint r = indices[cohort.reactionOffset + j];
            float rate = propensity(reactions[r], counts, changes, lane, cmd);
            if (rate > 0) lastPositive = r;
            float term = rate - correction, next = cumulative + term;
            correction = (next - cumulative) - term; cumulative = next;
            if (threshold < cumulative) { selected = r; break; }
        }
        if (selected == 0xffffffffu) selected = lastPositive;
        if (selected == 0xffffffffu) { fail(status, 1u, lane); p.done = 1; break; }
        Reaction r = reactions[selected];
        bool valid = true;
        for (uint j = 0; j < r.changeCount; ++j) {
            Change c = changes[r.changeOffset + j];
            long next = long(counts[ix(c.index, lane, cmd)]) + long(c.delta);
            if (next < 0) { fail(status, 2u, lane, selected); valid = false; }
            if (next > 4294967295l) { fail(status, 4u, lane, selected); valid = false; }
        }
        if (!valid) { p.done = 1; break; }
        for (uint j = 0; j < r.changeCount; ++j) {
            Change c = changes[r.changeOffset + j];
            ulong position = ix(c.index, lane, cmd);
            counts[position] = uint(long(counts[position]) + long(c.delta));
        }
        p.elapsed += wait;
        ++p.events;
    }
    progress[gid] = p;
    atomic_fetch_max_explicit(&status.maximumExactEvents, p.events, memory_order_relaxed);
    if (p.done == 0) atomic_fetch_add_explicit(&status.unfinishedExactLanes, 1u, memory_order_relaxed);
}

kernel void nvivo_hybrid_validate(
    device const float* values [[buffer(2)]], device const uint* modes [[buffer(10)]],
    device Status& status [[buffer(17)]], constant Command& cmd [[buffer(18)]],
    uint gid [[thread_position_in_grid]]) {
    uint species = gid / cmd.laneCount;
    if (species >= cmd.speciesCount) return;
    if (modes[species] == 3u && (!isfinite(values[gid]) || values[gid] < 0)) fail(status, 8u, gid % cmd.laneCount);
}
kernel void nvivo_hybrid_publish(
    device const float* values [[buffer(2)]], device const uint* counts [[buffer(5)]],
    device const uint* modes [[buffer(10)]], constant Command& cmd [[buffer(18)]],
    device const uint* requests [[buffer(19)]], device Publication* output [[buffer(20)]],
    constant uint& count [[buffer(21)]], uint gid [[thread_position_in_grid]]) {
    if (gid >= count) return;
    uint index = requests[gid], mode = modes[index / cmd.laneCount];
    output[gid] = {mode == 3u ? values[index] : 0.0f, mode == 3u ? 0u : counts[index], mode, 0u};
}
} // namespace vivo_hybrid
