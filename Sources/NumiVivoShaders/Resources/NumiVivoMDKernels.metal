#include <metal_stdlib>
using namespace metal;

namespace nvivo_md {
struct Command {
    uint particleCount;
    uint typeCount;
    uint electrostatics; // 0 cutoff, 1 reaction field
    uint periodic;
    float dtPS;
    float cutoffNM;
    float coulombPrefactor;
    float reactionFieldK;
    float reactionFieldC;
    float minimumDistanceNM;
    float reserved0;
    float reserved1;
    float4 cellA;
    float4 cellB;
    float4 cellC;
    float4 reciprocalA;
    float4 reciprocalB;
    float4 reciprocalC;
};
struct Status {
    atomic_uint flags;
    atomic_uint firstParticle;
    atomic_uint violationCount;
    atomic_uint reserved;
};
struct Bond { uint2 atoms; float2 parameters; };
struct Angle { uint4 atoms; float2 parameters; };
struct Torsion { uint4 atoms; float4 parameters; };
struct Incidence { uint termIndex; uint localIndex; };
struct PairException {
    uint2 atoms;
    float2 scales;
    float2 overrideC12C6;
    uint flags;
};
constant uint statusNonFinite = 1u;
constant uint statusOverlap = 2u;
constant uint statusInvalidGeometry = 4u;

inline void fail(device Status& status, uint flag, uint particle) {
    atomic_fetch_or_explicit(&status.flags, flag, memory_order_relaxed);
    atomic_fetch_min_explicit(&status.firstParticle, particle, memory_order_relaxed);
    atomic_fetch_add_explicit(&status.violationCount, 1u, memory_order_relaxed);
}
inline float3 minimumImage(float3 delta, constant Command& cmd) {
    if (cmd.periodic == 0) return delta;
    float3 fractional = float3(dot(cmd.reciprocalA.xyz, delta),
                               dot(cmd.reciprocalB.xyz, delta),
                               dot(cmd.reciprocalC.xyz, delta));
    fractional -= rint(fractional);
    return cmd.cellA.xyz * fractional.x + cmd.cellB.xyz * fractional.y + cmd.cellC.xyz * fractional.z;
}
inline float3 wrapPosition(float3 p, constant Command& cmd) {
    if (cmd.periodic == 0) return p;
    float3 fractional = float3(dot(cmd.reciprocalA.xyz, p),
                               dot(cmd.reciprocalB.xyz, p),
                               dot(cmd.reciprocalC.xyz, p));
    fractional -= floor(fractional);
    return cmd.cellA.xyz * fractional.x + cmd.cellB.xyz * fractional.y + cmd.cellC.xyz * fractional.z;
}
inline float4 bondContribution(uint local, Bond term, device const float4* position,
                               constant Command& cmd, device Status& status, uint owner) {
    float3 delta = minimumImage(position[term.atoms.x].xyz - position[term.atoms.y].xyz, cmd);
    float r2 = dot(delta, delta);
    if (!(r2 > cmd.minimumDistanceNM * cmd.minimumDistanceNM) || !isfinite(r2)) {
        fail(status, statusInvalidGeometry, owner); return 0;
    }
    float r = sqrt(r2);
    float dr = r - term.parameters.x;
    float energy = 0.5f * term.parameters.y * dr * dr;
    float3 forceA = -term.parameters.y * dr * delta / r;
    float3 force = local == 0 ? forceA : -forceA;
    return float4(force, 0.5f * energy);
}
inline float4 angleContribution(uint local, Angle term, device const float4* position,
                                constant Command& cmd, device Status& status, uint owner) {
    float3 u = minimumImage(position[term.atoms.x].xyz - position[term.atoms.y].xyz, cmd);
    float3 v = minimumImage(position[term.atoms.z].xyz - position[term.atoms.y].xyz, cmd);
    float ru2 = dot(u,u), rv2 = dot(v,v);
    if (!(ru2 > 1e-16f && rv2 > 1e-16f)) { fail(status, statusInvalidGeometry, owner); return 0; }
    float ru = sqrt(ru2), rv = sqrt(rv2);
    float cosine = clamp(dot(u,v)/(ru*rv), -1.0f, 1.0f);
    float theta = acos(cosine);
    float sine = sqrt(max(1.0f-cosine*cosine, 1e-12f));
    float dEdTheta = term.parameters.y * (theta-term.parameters.x);
    float3 gradA = (cosine*u/ru2 - v/(ru*rv))/sine;
    float3 gradC = (cosine*v/rv2 - u/(ru*rv))/sine;
    float3 forceA = -dEdTheta*gradA;
    float3 forceC = -dEdTheta*gradC;
    float3 forceB = -(forceA+forceC);
    float3 force = local == 0 ? forceA : (local == 1 ? forceB : forceC);
    float energy = 0.5f * term.parameters.y * (theta-term.parameters.x)*(theta-term.parameters.x);
    return float4(force, energy/3.0f);
}
inline float4 torsionContribution(uint local, Torsion term, device const float4* position,
                                  constant Command& cmd, device Status& status, uint owner) {
    float3 d0 = minimumImage(position[term.atoms.x].xyz-position[term.atoms.y].xyz, cmd);
    float3 d1 = minimumImage(position[term.atoms.z].xyz-position[term.atoms.y].xyz, cmd);
    float3 d2 = minimumImage(position[term.atoms.z].xyz-position[term.atoms.w].xyz, cmd);
    float3 c1 = cross(d0,d1), c2 = cross(d1,d2);
    float c1sq=dot(c1,c1), c2sq=dot(c2,c2), d1sq=dot(d1,d1);
    if (!(c1sq>1e-16f && c2sq>1e-16f && d1sq>1e-16f)) { fail(status,statusInvalidGeometry,owner); return 0; }
    float d1len=sqrt(d1sq);
    float phi=atan2(dot(d1,cross(c1,c2))/d1len, dot(c1,c2));
    float n=term.parameters.x, phase=term.parameters.y, k=term.parameters.z;
    float energy=k*(1.0f+cos(n*phi-phase));
    float dEdPhi=-k*n*sin(n*phi-phase);
    float3 f0=c1*(dEdPhi*d1len/c1sq);
    float3 f3=c2*(-dEdPhi*d1len/c2sq);
    float dot0=dot(d0,d1)/d1sq, dot2=dot(d2,d1)/d1sq;
    float3 f1=f0*(dot0-1.0f)-f3*dot2;
    float3 f2=-f0*dot0+f3*(dot2-1.0f);
    float3 force = local==0 ? f0 : (local==1 ? f1 : (local==2 ? f2 : f3));
    return float4(force, energy*0.25f);
}

kernel void nvivo_md_clear(device float4* forceEnergy [[buffer(0)]],
                           device Status& status [[buffer(1)]],
                           constant Command& cmd [[buffer(2)]],
                           uint gid [[thread_position_in_grid]]) {
    if (gid < cmd.particleCount) forceEnergy[gid]=0;
    if (gid==0) {
        atomic_store_explicit(&status.flags,0,memory_order_relaxed);
        atomic_store_explicit(&status.firstParticle,0xffffffffu,memory_order_relaxed);
        atomic_store_explicit(&status.violationCount,0,memory_order_relaxed);
        atomic_store_explicit(&status.reserved,0,memory_order_relaxed);
    }
}

kernel void nvivo_md_bonded(device const float4* position [[buffer(0)]],
                            device float4* forceEnergy [[buffer(1)]],
                            device const Bond* bonds [[buffer(2)]],
                            device const uint* bondOffsets [[buffer(3)]],
                            device const Incidence* bondIncidence [[buffer(4)]],
                            device const Angle* angles [[buffer(5)]],
                            device const uint* angleOffsets [[buffer(6)]],
                            device const Incidence* angleIncidence [[buffer(7)]],
                            device const Torsion* torsions [[buffer(8)]],
                            device const uint* torsionOffsets [[buffer(9)]],
                            device const Incidence* torsionIncidence [[buffer(10)]],
                            device Status& status [[buffer(11)]],
                            constant Command& cmd [[buffer(12)]],
                            uint gid [[thread_position_in_grid]]) {
    if (gid>=cmd.particleCount) return;
    float4 sum=0;
    for (uint i=bondOffsets[gid]; i<bondOffsets[gid+1]; ++i) {
        Incidence e=bondIncidence[i]; sum += bondContribution(e.localIndex,bonds[e.termIndex],position,cmd,status,gid);
    }
    for (uint i=angleOffsets[gid]; i<angleOffsets[gid+1]; ++i) {
        Incidence e=angleIncidence[i]; sum += angleContribution(e.localIndex,angles[e.termIndex],position,cmd,status,gid);
    }
    for (uint i=torsionOffsets[gid]; i<torsionOffsets[gid+1]; ++i) {
        Incidence e=torsionIncidence[i]; sum += torsionContribution(e.localIndex,torsions[e.termIndex],position,cmd,status,gid);
    }
    if (!all(isfinite(sum))) { fail(status,statusNonFinite,gid); return; }
    forceEnergy[gid] += sum;
}

inline int findException(uint owner, uint partner,
                         device const uint* offsets, device const uint* partners,
                         device const uint* indices) {
    uint lo=offsets[owner], hi=offsets[owner+1];
    while (lo<hi) {
        uint mid=lo+(hi-lo)/2, value=partners[mid];
        if (value<partner) lo=mid+1; else hi=mid;
    }
    return lo<offsets[owner+1] && partners[lo]==partner ? int(indices[lo]) : -1;
}

kernel void nvivo_md_nonbonded_direct(device const float4* position [[buffer(0)]],
                                      device float4* forceEnergy [[buffer(1)]],
                                      device const float4* dynamics [[buffer(2)]],
                                      device const uint* typeIndex [[buffer(3)]],
                                      device const float2* pairC12C6 [[buffer(4)]],
                                      device const PairException* exceptions [[buffer(5)]],
                                      device const uint* exceptionOffsets [[buffer(6)]],
                                      device const uint* exceptionPartners [[buffer(7)]],
                                      device const uint* exceptionIndices [[buffer(8)]],
                                      device Status& status [[buffer(9)]],
                                      constant Command& cmd [[buffer(10)]],
                                      uint gid [[thread_position_in_grid]]) {
    if (gid>=cmd.particleCount) return;
    float3 force=0; float energy=0;
    float cutoff2=cmd.cutoffNM*cmd.cutoffNM;
    for (uint j=0; j<cmd.particleCount; ++j) {
        if (j==gid) continue;
        float3 delta=minimumImage(position[gid].xyz-position[j].xyz,cmd);
        float r2=dot(delta,delta);
        if (r2>=cutoff2) continue;
        if (!(r2>cmd.minimumDistanceNM*cmd.minimumDistanceNM) || !isfinite(r2)) { fail(status,statusOverlap,gid); continue; }
        float invR=rsqrt(r2), invR2=invR*invR;
        float coulombScale=1, ljScale=1;
        uint ti=typeIndex[gid], tj=typeIndex[j];
        float2 coefficients=pairC12C6[ti*cmd.typeCount+tj];
        int ex=findException(gid,j,exceptionOffsets,exceptionPartners,exceptionIndices);
        if (ex>=0) {
            PairException value=exceptions[ex];
            coulombScale=value.scales.x; ljScale=value.scales.y;
            if ((value.flags&1u)!=0) coefficients=value.overrideC12C6;
        }
        float invR6=invR2*invR2*invR2, invR12=invR6*invR6;
        float lj=ljScale*(coefficients.x*invR12-coefficients.y*invR6);
        float ljForceScale=ljScale*(12.0f*coefficients.x*invR12-6.0f*coefficients.y*invR6)*invR2;
        float qq=dynamics[gid].z*dynamics[j].z;
        float electrostatic=0, electrostaticForceScale=0;
        if (coulombScale!=0 && qq!=0) {
            if (cmd.electrostatics==0) {
                electrostatic=cmd.coulombPrefactor*qq*invR;
                electrostaticForceScale=cmd.coulombPrefactor*qq*invR*invR2;
            } else {
                electrostatic=cmd.coulombPrefactor*qq*(invR+cmd.reactionFieldK*r2-cmd.reactionFieldC);
                electrostaticForceScale=cmd.coulombPrefactor*qq*(invR*invR2-2.0f*cmd.reactionFieldK);
            }
            electrostatic*=coulombScale; electrostaticForceScale*=coulombScale;
        }
        force += (ljForceScale+electrostaticForceScale)*delta;
        energy += 0.5f*(lj+electrostatic);
    }
    float4 value=forceEnergy[gid]+float4(force,energy);
    if (!all(isfinite(value))) { fail(status,statusNonFinite,gid); return; }
    forceEnergy[gid]=value;
}

kernel void nvivo_md_half_kick(device float4* velocity [[buffer(0)]],
                               device const float4* forceEnergy [[buffer(1)]],
                               device const float4* dynamics [[buffer(2)]],
                               constant Command& cmd [[buffer(3)]],
                               uint gid [[thread_position_in_grid]]) {
    if (gid>=cmd.particleCount) return;
    float invMass=dynamics[gid].y;
    if (invMass==0) return;
    velocity[gid].xyz += 0.5f*cmd.dtPS*invMass*forceEnergy[gid].xyz;
}

kernel void nvivo_md_drift(device float4* position [[buffer(0)]],
                           device const float4* velocity [[buffer(1)]],
                           device const float4* dynamics [[buffer(2)]],
                           constant Command& cmd [[buffer(3)]],
                           uint gid [[thread_position_in_grid]]) {
    if (gid>=cmd.particleCount || dynamics[gid].y==0) return;
    position[gid].xyz=wrapPosition(position[gid].xyz+cmd.dtPS*velocity[gid].xyz,cmd);
}

kernel void nvivo_md_kinetic(device const float4* velocity [[buffer(0)]],
                             device const float4* dynamics [[buffer(1)]],
                             device float* kinetic [[buffer(2)]],
                             constant Command& cmd [[buffer(3)]],
                             uint gid [[thread_position_in_grid]]) {
    if (gid>=cmd.particleCount) return;
    kinetic[gid]=0.5f*dynamics[gid].x*dot(velocity[gid].xyz,velocity[gid].xyz);
}

kernel void nvivo_md_validate(device const float4* position [[buffer(0)]],
                              device const float4* velocity [[buffer(1)]],
                              device const float4* forceEnergy [[buffer(2)]],
                              device Status& status [[buffer(3)]],
                              constant Command& cmd [[buffer(4)]],
                              uint gid [[thread_position_in_grid]]) {
    if (gid>=cmd.particleCount) return;
    if (!all(isfinite(position[gid])) || !all(isfinite(velocity[gid])) || !all(isfinite(forceEnergy[gid])))
        fail(status,statusNonFinite,gid);
}
} // namespace nvivo_md
