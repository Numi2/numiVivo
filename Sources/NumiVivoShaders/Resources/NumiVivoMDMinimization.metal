#include <metal_stdlib>
using namespace metal;

namespace nvivo_md_min {
struct MDCommand {
    uint particleCount,typeCount,electrostatics,periodic;
    uint stepLow,stepHigh,seedLow,seedHigh;
    float dtPS,cutoffNM,coulombPrefactor,reactionFieldK;
    float reactionFieldC,minimumDistanceNM,constraintTolerance,langevinA;
    float targetTemperatureK,boltzmannKJPerMolK;
    uint neighborCapacity; float neighborRadiusNM;
    float4 cellA,cellB,cellC,reciprocalA,reciprocalB,reciprocalC;
};
struct MinCommand { uint particleCount; float stepScale; float maxDisplacementNM; float reserved; };
struct Reduction { atomic_uint potentialBits; atomic_uint maximumForceBits; atomic_uint nonFiniteCount; atomic_uint reserved; };
inline float3 wrapPosition(float3 p,constant MDCommand&c){if(c.periodic==0)return p;float3 f=float3(dot(c.reciprocalA.xyz,p),dot(c.reciprocalB.xyz,p),dot(c.reciprocalC.xyz,p));f-=floor(f);return c.cellA.xyz*f.x+c.cellB.xyz*f.y+c.cellC.xyz*f.z;}
inline void atomicAddFloat(device atomic_uint*address,float value){uint old=atomic_load_explicit(address,memory_order_relaxed);while(true){float current=as_type<float>(old);uint desired=as_type<uint>(current+value);uint expected=old;if(atomic_compare_exchange_weak_explicit(address,&expected,desired,memory_order_relaxed,memory_order_relaxed))return;old=expected;}}

kernel void nvivo_md_minimize_clear(device Reduction&r[[buffer(0)]],uint gid[[thread_position_in_grid]]){if(gid==0){atomic_store_explicit(&r.potentialBits,as_type<uint>(0.0f),memory_order_relaxed);atomic_store_explicit(&r.maximumForceBits,as_type<uint>(0.0f),memory_order_relaxed);atomic_store_explicit(&r.nonFiniteCount,0,memory_order_relaxed);atomic_store_explicit(&r.reserved,0,memory_order_relaxed);}}

kernel void nvivo_md_minimize_reduce(device const float4*forceEnergy[[buffer(0)]],device Reduction&r[[buffer(1)]],constant MDCommand&c[[buffer(2)]],uint gid[[thread_position_in_grid]]){if(gid>=c.particleCount)return;float4 v=forceEnergy[gid];float fm=length(v.xyz);if(!all(isfinite(v))||!isfinite(fm)){atomic_fetch_add_explicit(&r.nonFiniteCount,1u,memory_order_relaxed);return;}atomicAddFloat(&r.potentialBits,v.w);atomic_fetch_max_explicit(&r.maximumForceBits,as_type<uint>(max(fm,0.0f)),memory_order_relaxed);}

kernel void nvivo_md_minimize_position(device float4*position[[buffer(0)]],device const float4*forceEnergy[[buffer(1)]],device const float4*dynamics[[buffer(2)]],constant MDCommand&md[[buffer(3)]],constant MinCommand&mc[[buffer(4)]],uint gid[[thread_position_in_grid]]){if(gid>=mc.particleCount||dynamics[gid].y==0)return;float3 delta=mc.stepScale*forceEnergy[gid].xyz;float lengthDelta=length(delta);if(lengthDelta>mc.maxDisplacementNM&&lengthDelta>0)delta*=mc.maxDisplacementNM/lengthDelta;position[gid]=float4(wrapPosition(position[gid].xyz+delta,md),0);}

kernel void nvivo_md_zero_velocity(device float4*velocity[[buffer(0)]],constant MDCommand&c[[buffer(1)]],uint gid[[thread_position_in_grid]]){if(gid<c.particleCount)velocity[gid]=0;}
} // namespace nvivo_md_min
