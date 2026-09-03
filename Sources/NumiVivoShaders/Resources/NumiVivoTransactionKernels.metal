#ifndef NUMIVIVO_METAL_ABI_H
#include "NumiVivoMetalABI.h"
#endif

namespace numivivo {

constant uint NVIVO_PUBLICATION_PREPARED = 5u;

inline ulong transaction_state_count(constant NVivoStepUniforms& uniforms) {
    return ulong(uniforms.speciesCount) * ulong(uniforms.cellCapacity);
}

inline ulong transaction_temporal_count(constant NVivoStepUniforms& uniforms) {
    return ulong(uniforms.temporalStateCount) * ulong(uniforms.cellCapacity);
}

inline ulong transaction_delay_count(constant NVivoStepUniforms& uniforms) {
    return ulong(uniforms.reactionCount) * ulong(uniforms.cellCapacity) *
           ulong(uniforms.delaySlotCount);
}

inline ulong transaction_refractory_count(constant NVivoStepUniforms& uniforms) {
    return ulong(uniforms.ruleCount) * ulong(uniforms.cellCapacity);
}

kernel void nvivo_finalize_prepare(
    constant NVivoProgramArguments& program [[buffer(0)]],
    constant NVivoStepUniforms& uniforms [[buffer(1)]],
    uint gid [[thread_position_in_grid]]) {

    if (gid != 0u) return;
    const uint flags = atomic_load_explicit(&program.diagnostics->flags, memory_order_relaxed);
    uint status = NVIVO_PUBLICATION_PREPARED;
    if ((flags & NVIVO_DIAGNOSTIC_PERMANENT_SHUTDOWN) != 0u) {
        status = NVIVO_PUBLICATION_PERMANENT_SHUTDOWN;
    } else if ((flags & NVIVO_DIAGNOSTIC_REVERSIBLE_SHUTDOWN) != 0u) {
        status = NVIVO_PUBLICATION_REVERSIBLE_SHUTDOWN;
    } else if ((flags & NVIVO_DIAGNOSTIC_SUBSTEP) != 0u) {
        status = NVIVO_PUBLICATION_SUBSTEP_REQUIRED;
    } else if ((flags & (NVIVO_DIAGNOSTIC_NONFINITE |
                         NVIVO_DIAGNOSTIC_BOUND_VIOLATION |
                         NVIVO_DIAGNOSTIC_REJECT |
                         NVIVO_DIAGNOSTIC_EXPRESSION_FAULT)) != 0u) {
        status = NVIVO_PUBLICATION_REJECTED;
    }

    if (status != NVIVO_PUBLICATION_PREPARED) {
        atomic_store_explicit(program.eventCount, 0u, memory_order_relaxed);
    }
    program.publication->status = status;
    program.publication->diagnosticFlags = flags;
    program.publication->activeCellCount = uniforms.activeCellCount;
    program.publication->eventCount = status == NVIVO_PUBLICATION_PREPARED
        ? min(atomic_load_explicit(program.eventCount, memory_order_relaxed), uniforms.eventCapacity)
        : 0u;
}

kernel void nvivo_commit_prepared(
    constant NVivoProgramArguments& program [[buffer(0)]],
    constant NVivoStepUniforms& uniforms [[buffer(1)]],
    uint gid [[thread_position_in_grid]]) {

    if (program.publication->status != NVIVO_PUBLICATION_PREPARED) return;
    const ulong index = ulong(gid);
    const ulong stateCount = transaction_state_count(uniforms);
    const ulong temporalCount = transaction_temporal_count(uniforms);
    const ulong delayCount = transaction_delay_count(uniforms);
    const ulong refractoryCount = transaction_refractory_count(uniforms);

    if (index < stateCount) program.committedState[index] = program.shadowState[index];
    if (index < temporalCount) program.committedTemporalState[index] = program.shadowTemporalState[index];
    if (index < delayCount) program.committedDelayedFlux[index] = program.shadowDelayedFlux[index];
    if (index < refractoryCount) program.committedRuleRefractory[index] = program.shadowRuleRefractory[index];

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

kernel void nvivo_rollback_prepared(
    constant NVivoProgramArguments& program [[buffer(0)]],
    uint gid [[thread_position_in_grid]]) {

    if (gid != 0u) return;
    atomic_store_explicit(program.eventCount, 0u, memory_order_relaxed);
    program.publication->status = NVIVO_PUBLICATION_REJECTED;
    program.publication->diagnosticFlags = 0u;
    program.publication->eventCount = 0u;
}

struct NVivoHashUniforms {
    uint byteCountLow;
    uint byteCountHigh;
    uint chunkBytes;
    uint chunkCount;
    uint outputWordOffset;
    uint domain;
    uint logicalStepLow;
    uint logicalStepHigh;
};

constant uint nvivoSHA256Constants[64] = {
    0x428a2f98u, 0x71374491u, 0xb5c0fbcfu, 0xe9b5dba5u,
    0x3956c25bu, 0x59f111f1u, 0x923f82a4u, 0xab1c5ed5u,
    0xd807aa98u, 0x12835b01u, 0x243185beu, 0x550c7dc3u,
    0x72be5d74u, 0x80deb1feu, 0x9bdc06a7u, 0xc19bf174u,
    0xe49b69c1u, 0xefbe4786u, 0x0fc19dc6u, 0x240ca1ccu,
    0x2de92c6fu, 0x4a7484aau, 0x5cb0a9dcu, 0x76f988dau,
    0x983e5152u, 0xa831c66du, 0xb00327c8u, 0xbf597fc7u,
    0xc6e00bf3u, 0xd5a79147u, 0x06ca6351u, 0x14292967u,
    0x27b70a85u, 0x2e1b2138u, 0x4d2c6dfcu, 0x53380d13u,
    0x650a7354u, 0x766a0abbu, 0x81c2c92eu, 0x92722c85u,
    0xa2bfe8a1u, 0xa81a664bu, 0xc24b8b70u, 0xc76c51a3u,
    0xd192e819u, 0xd6990624u, 0xf40e3585u, 0x106aa070u,
    0x19a4c116u, 0x1e376c08u, 0x2748774cu, 0x34b0bcb5u,
    0x391c0cb3u, 0x4ed8aa4au, 0x5b9cca4fu, 0x682e6ff3u,
    0x748f82eeu, 0x78a5636fu, 0x84c87814u, 0x8cc70208u,
    0x90befffau, 0xa4506cebu, 0xbef9a3f7u, 0xc67178f2u
};

struct NVivoSHA256Context {
    uint state[8];
    uchar block[64];
    uint blockLength;
    ulong totalLength;
};

inline uint nvivo_rotr(uint value, uint amount) {
    return (value >> amount) | (value << (32u - amount));
}

inline void nvivo_sha256_transform(thread NVivoSHA256Context& context) {
    thread uint words[64];
    for (uint index = 0u; index < 16u; ++index) {
        const uint offset = index * 4u;
        words[index] = (uint(context.block[offset]) << 24u) |
                       (uint(context.block[offset + 1u]) << 16u) |
                       (uint(context.block[offset + 2u]) << 8u) |
                       uint(context.block[offset + 3u]);
    }
    for (uint index = 16u; index < 64u; ++index) {
        const uint s0 = nvivo_rotr(words[index - 15u], 7u) ^
                        nvivo_rotr(words[index - 15u], 18u) ^
                        (words[index - 15u] >> 3u);
        const uint s1 = nvivo_rotr(words[index - 2u], 17u) ^
                        nvivo_rotr(words[index - 2u], 19u) ^
                        (words[index - 2u] >> 10u);
        words[index] = words[index - 16u] + s0 + words[index - 7u] + s1;
    }

    uint a = context.state[0];
    uint b = context.state[1];
    uint c = context.state[2];
    uint d = context.state[3];
    uint e = context.state[4];
    uint f = context.state[5];
    uint g = context.state[6];
    uint h = context.state[7];

    for (uint index = 0u; index < 64u; ++index) {
        const uint sigma1 = nvivo_rotr(e, 6u) ^ nvivo_rotr(e, 11u) ^ nvivo_rotr(e, 25u);
        const uint choice = (e & f) ^ ((~e) & g);
        const uint temporary1 = h + sigma1 + choice + nvivoSHA256Constants[index] + words[index];
        const uint sigma0 = nvivo_rotr(a, 2u) ^ nvivo_rotr(a, 13u) ^ nvivo_rotr(a, 22u);
        const uint majority = (a & b) ^ (a & c) ^ (b & c);
        const uint temporary2 = sigma0 + majority;
        h = g;
        g = f;
        f = e;
        e = d + temporary1;
        d = c;
        c = b;
        b = a;
        a = temporary1 + temporary2;
    }

    context.state[0] += a;
    context.state[1] += b;
    context.state[2] += c;
    context.state[3] += d;
    context.state[4] += e;
    context.state[5] += f;
    context.state[6] += g;
    context.state[7] += h;
}

inline void nvivo_sha256_initialize(thread NVivoSHA256Context& context) {
    context.state[0] = 0x6a09e667u;
    context.state[1] = 0xbb67ae85u;
    context.state[2] = 0x3c6ef372u;
    context.state[3] = 0xa54ff53au;
    context.state[4] = 0x510e527fu;
    context.state[5] = 0x9b05688cu;
    context.state[6] = 0x1f83d9abu;
    context.state[7] = 0x5be0cd19u;
    context.blockLength = 0u;
    context.totalLength = 0ul;
}

inline void nvivo_sha256_update_byte(thread NVivoSHA256Context& context, uchar value) {
    context.block[context.blockLength++] = value;
    context.totalLength += 1ul;
    if (context.blockLength == 64u) {
        nvivo_sha256_transform(context);
        context.blockLength = 0u;
    }
}

inline void nvivo_sha256_update_uint(thread NVivoSHA256Context& context, uint value) {
    nvivo_sha256_update_byte(context, uchar(value & 0xffu));
    nvivo_sha256_update_byte(context, uchar((value >> 8u) & 0xffu));
    nvivo_sha256_update_byte(context, uchar((value >> 16u) & 0xffu));
    nvivo_sha256_update_byte(context, uchar((value >> 24u) & 0xffu));
}

inline void nvivo_sha256_finalize(thread NVivoSHA256Context& context,
                                  device uint* output) {
    const ulong bitLength = context.totalLength * 8ul;
    context.block[context.blockLength++] = 0x80u;
    if (context.blockLength > 56u) {
        while (context.blockLength < 64u) context.block[context.blockLength++] = 0u;
        nvivo_sha256_transform(context);
        context.blockLength = 0u;
    }
    while (context.blockLength < 56u) context.block[context.blockLength++] = 0u;
    for (uint byteIndex = 0u; byteIndex < 8u; ++byteIndex) {
        const uint shift = (7u - byteIndex) * 8u;
        context.block[context.blockLength++] = uchar((bitLength >> shift) & 0xfful);
    }
    nvivo_sha256_transform(context);
    for (uint index = 0u; index < 8u; ++index) output[index] = context.state[index];
}

kernel void nvivo_sha256_state_chunks(
    device const uchar* bytes [[buffer(0)]],
    device uint* outputWords [[buffer(1)]],
    constant NVivoHashUniforms& uniforms [[buffer(2)]],
    uint gid [[thread_position_in_grid]]) {

    if (gid >= uniforms.chunkCount || uniforms.chunkBytes == 0u) return;
    const ulong byteCount = ulong(uniforms.byteCountLow) |
                            (ulong(uniforms.byteCountHigh) << 32ul);
    const ulong start = ulong(gid) * ulong(uniforms.chunkBytes);
    const ulong end = min(byteCount, start + ulong(uniforms.chunkBytes));

    NVivoSHA256Context context;
    nvivo_sha256_initialize(context);
    nvivo_sha256_update_uint(context, 0x4e565348u);
    nvivo_sha256_update_uint(context, uniforms.domain);
    nvivo_sha256_update_uint(context, gid);
    nvivo_sha256_update_uint(context, uniforms.chunkCount);
    nvivo_sha256_update_uint(context, uniforms.logicalStepLow);
    nvivo_sha256_update_uint(context, uniforms.logicalStepHigh);
    nvivo_sha256_update_uint(context, uniforms.byteCountLow);
    nvivo_sha256_update_uint(context, uniforms.byteCountHigh);
    for (ulong index = start; index < end; ++index) {
        nvivo_sha256_update_byte(context, bytes[index]);
    }
    nvivo_sha256_finalize(
        context,
        outputWords + ulong(uniforms.outputWordOffset) + ulong(gid) * 8ul
    );
}

} // namespace numivivo
