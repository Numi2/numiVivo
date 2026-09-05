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
struct ReduceCommand { uint sourceCount; uint reserved0; uint reserved1; uint reserved2; };
inline float3 wrapPosition(float3 p,constant MDCommand&c){if(c.periodic==0)return p;float3 f=float3(dot(c.reciprocalA.xyz,p),dot(c.reciprocalB.xyz,p),dot(c.reciprocalC.xyz,p));f-=floor(f);return c.cellA.xyz*f.x+c.cellB.xyz*f.y+c.cellC.xyz*f.z;}

/// x = potential energy contribution, y = force magnitude.
kernel void nvivo_md_minimize_terms(device const float4*forceEnergy[[buffer(0)]],
                                     device float2*terms[[buffer(1)]],
                                     constant MDCommand&c[[buffer(2)]],
                                     uint gid[[thread_position_in_grid]]){
    if(gid>=c.particleCount)return;float4 value=forceEnergy[gid];float fm=length(value.xyz);terms[gid]=all(isfinite(value))&&isfinite(fm)?float2(value.w,fm):float2(NAN,NAN);
}

/// Fixed pairwise tree: every output element always combines the same two input
/// indices. This avoids scheduling-dependent atomic summation.
kernel void nvivo_md_minimize_reduce_stage(device const float2*source[[buffer(0)]],
                                            device float2*destination[[buffer(1)]],
                                            constant ReduceCommand&c[[buffer(2)]],
                                            uint gid[[thread_position_in_grid]]){
    uint left=gid*2u;if(left>=c.sourceCount)return;float2 a=source[left];float2 b=left+1u<c.sourceCount?source[left+1u]:float2(0.0f,0.0f);destination[gid]=float2(a.x+b.x,max(a.y,b.y));
}

kernel void nvivo_md_minimize_position(device float4*position[[buffer(0)]],
                                       device const float4*forceEnergy[[buffer(1)]],
                                       device const float4*dynamics[[buffer(2)]],
                                       constant MDCommand&md[[buffer(3)]],
                                       constant MinCommand&mc[[buffer(4)]],
                                       uint gid[[thread_position_in_grid]]){
    if(gid>=mc.particleCount||dynamics[gid].y==0)return;float3 delta=mc.stepScale*forceEnergy[gid].xyz;float lengthDelta=length(delta);if(lengthDelta>mc.maxDisplacementNM&&lengthDelta>0)delta*=mc.maxDisplacementNM/lengthDelta;position[gid]=float4(wrapPosition(position[gid].xyz+delta,md),0);
}

kernel void nvivo_md_zero_velocity(device float4*velocity[[buffer(0)]],
                                   constant MDCommand&c[[buffer(1)]],
                                   uint gid[[thread_position_in_grid]]){if(gid<c.particleCount)velocity[gid]=0;}
} // namespace nvivo_md_min
