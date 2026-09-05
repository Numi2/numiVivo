#include <metal_stdlib>
using namespace metal;

namespace nvivo_md_grid {
struct Command {
    uint particleCount;
    uint cellsX;
    uint cellsY;
    uint cellsZ;
    uint cellCount;
    uint cellCapacity;
    uint neighborCapacity;
    uint periodic;
    float neighborRadiusNM;
    float reserved0;
    float reserved1;
    float reserved2;
    float4 cellA;
    float4 cellB;
    float4 cellC;
    float4 reciprocalA;
    float4 reciprocalB;
    float4 reciprocalC;
};
struct Status { atomic_uint flags, firstParticle, violationCount, reserved; };
constant uint statusNeighborOverflow = 16u;
constant uint statusGridOverflow = 64u;

inline void fail(device Status& s, uint flag, uint particle) {
    atomic_fetch_or_explicit(&s.flags, flag, memory_order_relaxed);
    atomic_fetch_min_explicit(&s.firstParticle, particle, memory_order_relaxed);
    atomic_fetch_add_explicit(&s.violationCount, 1u, memory_order_relaxed);
}
inline float3 fractional(float3 p, constant Command& c) {
    float3 f = float3(dot(c.reciprocalA.xyz,p), dot(c.reciprocalB.xyz,p), dot(c.reciprocalC.xyz,p));
    return f - floor(f);
}
inline uint cellIndex(uint3 q, constant Command& c) {
    return (q.z*c.cellsY + q.y)*c.cellsX + q.x;
}
inline uint3 cellOf(float3 p, constant Command& c) {
    float3 f = fractional(p,c);
    return uint3(min(uint(f.x*float(c.cellsX)), c.cellsX-1),
                 min(uint(f.y*float(c.cellsY)), c.cellsY-1),
                 min(uint(f.z*float(c.cellsZ)), c.cellsZ-1));
}
inline int wrapCell(int value, uint n) {
    int m = value % int(n);
    return m < 0 ? m + int(n) : m;
}
inline float3 minimumImage(float3 d, constant Command& c) {
    float3 f = float3(dot(c.reciprocalA.xyz,d),dot(c.reciprocalB.xyz,d),dot(c.reciprocalC.xyz,d));
    f -= rint(f);
    return c.cellA.xyz*f.x + c.cellB.xyz*f.y + c.cellC.xyz*f.z;
}

[[host_name("nvivo_md_grid_clear")]] kernel void nvivo_md_grid_clear(device atomic_uint* counts [[buffer(0)]],
                                device Status& status [[buffer(1)]],
                                constant Command& c [[buffer(2)]],
                                uint gid [[thread_position_in_grid]]) {
    if (gid < c.cellCount) atomic_store_explicit(&counts[gid], 0u, memory_order_relaxed);
}

[[host_name("nvivo_md_grid_bin")]] kernel void nvivo_md_grid_bin(device const float4* positions [[buffer(0)]],
                              device atomic_uint* counts [[buffer(1)]],
                              device uint* particles [[buffer(2)]],
                              device Status& status [[buffer(3)]],
                              constant Command& c [[buffer(4)]],
                              uint gid [[thread_position_in_grid]]) {
    if (gid >= c.particleCount) return;
    uint cell = cellIndex(cellOf(positions[gid].xyz,c),c);
    uint slot = atomic_fetch_add_explicit(&counts[cell], 1u, memory_order_relaxed);
    if (slot >= c.cellCapacity) {
        fail(status,statusGridOverflow,gid);
        return;
    }
    particles[ulong(cell)*ulong(c.cellCapacity)+slot] = gid;
}

[[host_name("nvivo_md_grid_build_neighbors")]] kernel void nvivo_md_grid_build_neighbors(device const float4* positions [[buffer(0)]],
                                           device const atomic_uint* counts [[buffer(1)]],
                                           device const uint* particles [[buffer(2)]],
                                           device uint* neighborCounts [[buffer(3)]],
                                           device uint* neighborIndices [[buffer(4)]],
                                           device float4* reference [[buffer(5)]],
                                           device Status& status [[buffer(6)]],
                                           constant Command& c [[buffer(7)]],
                                           uint gid [[thread_position_in_grid]]) {
    if (gid >= c.particleCount) return;
    if ((atomic_load_explicit(&status.flags,memory_order_relaxed) & statusGridOverflow) != 0u) {
        neighborCounts[gid] = 0;
        reference[gid] = positions[gid];
        return;
    }
    uint3 own = cellOf(positions[gid].xyz,c);
    uint seen[27];
    uint seenCount = 0;
    uint n = 0;
    float radius2 = c.neighborRadiusNM*c.neighborRadiusNM;
    for (int dz=-1; dz<=1; ++dz) {
        for (int dy=-1; dy<=1; ++dy) {
            for (int dx=-1; dx<=1; ++dx) {
                uint3 q = uint3(uint(wrapCell(int(own.x)+dx,c.cellsX)),
                                uint(wrapCell(int(own.y)+dy,c.cellsY)),
                                uint(wrapCell(int(own.z)+dz,c.cellsZ)));
                uint cell = cellIndex(q,c);
                bool duplicate = false;
                for (uint s=0; s<seenCount; ++s) if (seen[s] == cell) { duplicate = true; break; }
                if (duplicate) continue;
                seen[seenCount++] = cell;
                uint occupancy = min(atomic_load_explicit(&counts[cell],memory_order_relaxed), c.cellCapacity);
                for (uint slot=0; slot<occupancy; ++slot) {
                    uint other = particles[ulong(cell)*ulong(c.cellCapacity)+slot];
                    if (other == gid) continue;
                    float3 d = minimumImage(positions[gid].xyz-positions[other].xyz,c);
                    if (dot(d,d) < radius2) {
                        if (n < c.neighborCapacity)
                            neighborIndices[ulong(gid)*ulong(c.neighborCapacity)+n] = other;
                        ++n;
                    }
                }
            }
        }
    }
    if (n > c.neighborCapacity) fail(status,statusNeighborOverflow,gid);
    neighborCounts[gid] = min(n,c.neighborCapacity);
    reference[gid] = positions[gid];
}
} // namespace nvivo_md_grid
