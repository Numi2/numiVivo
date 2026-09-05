#include <metal_stdlib>
#include "NumiVivoErrorFunctions.metalh"
using namespace metal;

namespace nvivo_pme_correction {
struct MDCommand {
    uint particleCount,typeCount,electrostatics,periodic;
    uint stepLow,stepHigh,seedLow,seedHigh;
    float dtPS,cutoffNM,coulombPrefactor,reactionFieldK;
    float reactionFieldC,minimumDistanceNM,constraintTolerance,langevinA;
    float targetTemperatureK,boltzmannKJPerMolK;
    uint neighborCapacity; float neighborRadiusNM;
    float4 cellA,cellB,cellC,reciprocalA,reciprocalB,reciprocalC;
};
struct Status { atomic_uint flags,firstParticle,violationCount,reserved; };
struct PairException { uint2 atoms; float2 scales; float2 overrideC12C6; uint flags; };
constant uint statusNonFinite=1u,statusInvalidGeometry=4u;
inline void fail(device Status&s,uint flag,uint particle){atomic_fetch_or_explicit(&s.flags,flag,memory_order_relaxed);atomic_fetch_min_explicit(&s.firstParticle,particle,memory_order_relaxed);atomic_fetch_add_explicit(&s.violationCount,1u,memory_order_relaxed);}
inline float3 minimumImage(float3 d,constant MDCommand&c){float3 f=float3(dot(c.reciprocalA.xyz,d),dot(c.reciprocalB.xyz,d),dot(c.reciprocalC.xyz,d));f-=rint(f);return c.cellA.xyz*f.x+c.cellB.xyz*f.y+c.cellC.xyz*f.z;}

[[host_name("nvivo_pme_exception_correction")]] kernel void nvivo_pme_exception_correction(device const float4*positions[[buffer(0)]],
                                            device const float4*dynamics[[buffer(1)]],
                                            device float4*forceEnergy[[buffer(2)]],
                                            device const PairException*exceptions[[buffer(3)]],
                                            device const uint*offsets[[buffer(4)]],
                                            device const uint*partners[[buffer(5)]],
                                            device const uint*indices[[buffer(6)]],
                                            device Status&status[[buffer(7)]],
                                            constant MDCommand&c[[buffer(8)]],
                                            uint gid[[thread_position_in_grid]]){
    if(gid>=c.particleCount)return;float3 force=0;float energy=0;float beta=c.reactionFieldK;
    for(uint k=offsets[gid];k<offsets[gid+1];++k){uint other=partners[k];PairException ex=exceptions[indices[k]];float deltaScale=ex.scales.x-1.0f;if(deltaScale==0.0f)continue;float3 d=minimumImage(positions[gid].xyz-positions[other].xyz,c);float r2=dot(d,d);if(!(r2>c.minimumDistanceNM*c.minimumDistanceNM)||!isfinite(r2)){fail(status,statusInvalidGeometry,gid);continue;}float ir=rsqrt(r2),r=r2*ir,ir2=ir*ir,br=beta*r;float erfv=nvivo_math::erf(br),gaussian=exp(-br*br);float qq=dynamics[gid].z*dynamics[other].z;if(qq==0)continue;float pref=c.coulombPrefactor*qq*deltaScale;energy+=0.5f*pref*erfv*ir;float scale=pref*(nvivo_math::erfMinusGaussian(br)*ir*ir2);force+=scale*d;}
    float4 value=forceEnergy[gid]+float4(force,energy);if(!all(isfinite(value))){fail(status,statusNonFinite,gid);return;}forceEnergy[gid]=value;
}

/// Tin-foil Ewald with a uniform neutralizing background for a non-neutral cell.
/// The correction is a scalar energy and therefore does not alter particle forces.
[[host_name("nvivo_pme_background_energy")]] kernel void nvivo_pme_background_energy(device float4*forceEnergy[[buffer(0)]],
                                         constant float&energy[[buffer(1)]],
                                         uint gid[[thread_position_in_grid]]){
    if(gid==0) forceEnergy[0].w += energy;
}
} // namespace nvivo_pme_correction
