#include <metal_stdlib>
using namespace metal;

// This file contains observation/input operations only. Reaction evolution uses
// the existing nvivo_f1_heun_predict/correct ProgramPack kernels, not new kinetics.
struct TargetInputCommand {
    uint lanes;
    uint drugSpecies;
    uint competitorSpecies; // UINT_MAX means no competitor in the compiled graph.
    float baseExposureM;
};
struct TargetObservationCommand {
    uint4 targetSpecies; // free, reversible, covalent, competitor
    float4 data;        // observed fraction, reported SD (0 when unavailable), tolerance, reserved
    uint4 shape;        // observable (0..3), lane capacity, reserved, reserved
};

kernel void nvivo_target_likelihood_expose(
    device float* state [[buffer(0)]],
    device const float2* inputs [[buffer(1)]], // exposure multiplier, competitor concentration
    constant TargetInputCommand& cmd [[buffer(2)]],
    uint lane [[thread_position_in_grid]]) {
    if (lane >= cmd.lanes) return;
    state[cmd.drugSpecies * cmd.lanes + lane] = cmd.baseExposureM * inputs[lane].x;
    if (cmd.competitorSpecies != 0xffffffffu)
        state[cmd.competitorSpecies * cmd.lanes + lane] = inputs[lane].y;
}

kernel void nvivo_target_likelihood_observe(
    device const float* state [[buffer(0)]],
    device const float4* noise [[buffer(1)]], // scale, extra SD, bias, reserved
    device float2* accumulator [[buffer(2)]], // compensated sum, compensation
    device uint* failures [[buffer(3)]],
    constant TargetObservationCommand& cmd [[buffer(4)]],
    uint lane [[thread_position_in_grid]]) {
    const uint lanes = cmd.shape.y;
    if (lane >= lanes || failures[lane] != 0u) return;
    const uint4 idx = cmd.targetSpecies * lanes + lane;
    const float4 f = float4(state[idx.x], state[idx.y], state[idx.z], state[idx.w]);
    const float total = (f.x + f.y) + (f.z + f.w);
    if (!all(isfinite(f)) || any(f < 0.0f) || !isfinite(total) || total <= 0.0f ||
        abs(total - 1.0f) > cmd.data.z || cmd.shape.x > 3u) {
        failures[lane] |= 1u;
        return;
    }
    float predicted;
    switch (cmd.shape.x) {
        case 0u: predicted = (f.y + f.z) / total; break;
        case 1u: predicted = f.z / total; break;
        case 2u: predicted = f.x / total; break;
        default: predicted = f.w / total; break;
    }
    const float4 model = noise[lane];
    const float sd = hypot(cmd.data.y * model.x, model.y);
    if (!all(isfinite(model)) || model.x <= 0.0f || model.y < 0.0f ||
        !isfinite(sd) || sd < 1.1754943508222875e-38f) {
        failures[lane] |= 2u;
        return;
    }
    const float z = (predicted + model.z - cmd.data.x) / sd;
    const float term = -0.5f * z * z - log(sd) - 0.91893853320467274178f;
    const float2 old = accumulator[lane];
    const float adjusted = term - old.y;
    const float sum = old.x + adjusted;
    const float correction = (sum - old.x) - adjusted;
    if (!isfinite(term) || !isfinite(sum) || !isfinite(correction)) {
        failures[lane] |= 4u;
        return;
    }
    accumulator[lane] = float2(sum, correction);
}
