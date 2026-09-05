#ifndef NUMIVIVO_METAL_ABI_H
#include "NumiVivoMetalABI.h"
#endif

namespace numivivo {

[[host_name("nvivo_commit_prepared_authoritative")]] kernel void nvivo_commit_prepared_authoritative(
    constant NVivoProgramArguments& program [[buffer(0)]],
    constant NVivoStepUniforms& uniforms [[buffer(1)]],
    uint gid [[thread_position_in_grid]]) {

    const ulong index = ulong(gid);
    const ulong stateCount = ulong(uniforms.speciesCount) * ulong(uniforms.cellCapacity);
    const ulong temporalCount = ulong(uniforms.temporalStateCount) * ulong(uniforms.cellCapacity);
    const ulong delayCount = ulong(uniforms.reactionCount) * ulong(uniforms.cellCapacity) *
                             ulong(uniforms.delaySlotCount);
    const ulong refractoryCount = ulong(uniforms.ruleCount) * ulong(uniforms.cellCapacity);

    if (index < stateCount) program.committedState[index] = program.shadowState[index];
    if (index < temporalCount) program.committedTemporalState[index] = program.shadowTemporalState[index];
    if (index < delayCount) program.committedDelayedFlux[index] = program.shadowDelayedFlux[index];
    if (index < refractoryCount) program.committedRuleRefractory[index] = program.shadowRuleRefractory[index];

    // No thread consults publication.status. The actor has exclusive ownership
    // of the prepared state, and the CPU observes publication only after this
    // command buffer completes. This avoids a global status race inside a
    // dispatch, for which Metal provides no grid-wide barrier.
    if (gid == 0u) {
        program.publication->committedStepLow = uniforms.logicalStepLow;
        program.publication->committedStepHigh = uniforms.logicalStepHigh;
        program.publication->stateVersion += 1u;
        program.publication->status = NVIVO_PUBLICATION_COMMITTED;
        program.publication->diagnosticFlags = 0u;
        program.publication->eventCount = min(
            atomic_load_explicit(program.eventCount, memory_order_relaxed),
            uniforms.eventCapacity
        );
    }
}

} // namespace numivivo
