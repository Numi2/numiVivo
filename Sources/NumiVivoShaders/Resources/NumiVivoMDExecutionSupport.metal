#include <metal_stdlib>
using namespace metal;
namespace nvivo_md_support {
struct Command {
    uint particleCount,typeCount,electrostatics,periodic;
    uint stepLow,stepHigh,seedLow,seedHigh;
    float dtPS,cutoffNM,coulombPrefactor,reactionFieldK;
    float reactionFieldC,minimumDistanceNM,constraintTolerance,langevinA;
    float targetTemperatureK,boltzmannKJPerMolK;
    uint neighborCapacity; float neighborRadiusNM;
    float4 cellA,cellB,cellC,reciprocalA,reciprocalB,reciprocalC;
};
struct Constraint { uint2 atoms; float distanceNM; uint reserved; };
struct Incidence { uint termIndex, localIndex; };
struct Status { atomic_uint flags,firstParticle,violationCount,reserved; };
static_assert(sizeof(Command)==176, "MD support command ABI");
inline float3 image(float3 d,constant Command&c){
    if(c.periodic==0)return d;
    float3 f=float3(dot(c.reciprocalA.xyz,d),dot(c.reciprocalB.xyz,d),dot(c.reciprocalC.xyz,d));
    f-=rint(f);return c.cellA.xyz*f.x+c.cellB.xyz*f.y+c.cellC.xyz*f.z;
}
inline void fail(device Status&s,uint g,uint flag){
    atomic_fetch_or_explicit(&s.flags,flag,memory_order_relaxed);
    atomic_fetch_min_explicit(&s.firstParticle,g,memory_order_relaxed);
    atomic_fetch_add_explicit(&s.violationCount,1u,memory_order_relaxed);
}

// Position projection changes a constrained drift; retain its corresponding
// velocity increment instead of silently discarding the constraint impulse.
[[host_name("nvivo_md_drift_constraint_impulse")]] kernel void nvivo_md_drift_constraint_impulse(
    device const float4*position[[buffer(0)]],device const float4*unconstrained[[buffer(1)]],
    device float4*velocity[[buffer(2)]],device const float4*dynamics[[buffer(3)]],
    device Status&status[[buffer(4)]],constant Command&c[[buffer(5)]],uint g[[thread_position_in_grid]]){
    if(g>=c.particleCount||dynamics[g].y==0)return;
    float3 delta=image(position[g].xyz-unconstrained[g].xyz,c);
    float3 v=velocity[g].xyz+delta/c.dtPS;
    if(!all(isfinite(v))){fail(status,g,1u);return;}
    velocity[g]=float4(v,0);
}

[[host_name("nvivo_md_projected_force_seed")]] kernel void nvivo_md_projected_force_seed(
    device const float4*forceEnergy[[buffer(0)]],device const float4*dynamics[[buffer(1)]],
    device float4*direction[[buffer(2)]],constant Command&c[[buffer(3)]],uint g[[thread_position_in_grid]]){
    if(g<c.particleCount)direction[g]=dynamics[g].y>0?float4(forceEnergy[g].xyz,0):float4(0);
}

// The same sparse projection used for velocities runs with an identity metric
// on massive particles for minimization. Certify its tangent residual separately
// from the force magnitude; an unfinished projection is not convergence.
[[host_name("nvivo_md_validate_projected_force")]] kernel void nvivo_md_validate_projected_force(
    device const float4*position[[buffer(0)]],device const float4*direction[[buffer(1)]],
    device const Constraint*constraints[[buffer(2)]],device const uint*offsets[[buffer(3)]],
    device const Incidence*incidence[[buffer(4)]],device Status&status[[buffer(5)]],
    constant Command&c[[buffer(6)]],uint g[[thread_position_in_grid]]){
    if(g>=c.particleCount)return;
    for(uint i=offsets[g];i<offsets[g+1];++i){
        Incidence e=incidence[i];if(e.localIndex!=0)continue;
        Constraint q=constraints[e.termIndex];
        float3 r=image(position[q.atoms.x].xyz-position[q.atoms.y].xyz,c);
        float3 a=direction[q.atoms.x].xyz,b=direction[q.atoms.y].xyz;
        float denominator=length(r)*(length(a)+length(b));
        float residual=fabs(dot(r,a-b));
        if(!isfinite(residual)||!isfinite(denominator)||
           residual>c.constraintTolerance*max(denominator,1e-20f))fail(status,g,8u);
    }
}

[[host_name("nvivo_md_observation_terms")]] kernel void nvivo_md_observation_terms(
    device const float4*forceEnergy[[buffer(0)]],device const float4*velocity[[buffer(1)]],
    device const float4*dynamics[[buffer(2)]],device float2*terms[[buffer(3)]],
    constant Command&c[[buffer(4)]],uint g[[thread_position_in_grid]]){
    if(g<c.particleCount)terms[g]=float2(forceEnergy[g].w,
        dynamics[g].x>0?0.5f*dynamics[g].x*dot(velocity[g].xyz,velocity[g].xyz):0.0f);
}

[[host_name("nvivo_md_sum_pair_reduce")]] kernel void nvivo_md_sum_pair_reduce(
    device const float2*source[[buffer(0)]],device float2*destination[[buffer(1)]],
    constant uint&count[[buffer(2)]],uint g[[thread_position_in_grid]]){
    uint first=2u*g;if(first>=count)return;
    destination[g]=source[first]+(first+1u<count?source[first+1u]:float2(0));
}
} // namespace nvivo_md_support
