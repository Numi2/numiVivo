#include <metal_stdlib>
using namespace metal;
namespace nvivo_md_barostat {
struct Command {
    uint particleCount,componentCount;float scale,toleranceNM;
    float4 cellA,cellB,cellC,reciprocalA,reciprocalB,reciprocalC;
};
struct Status { atomic_uint flags,firstParticle,violationCount,reserved; };
static_assert(sizeof(Command)==112,"NPT command ABI");
inline float3 image(float3 d,constant Command&c){
    float3 f=float3(dot(c.reciprocalA.xyz,d),dot(c.reciprocalB.xyz,d),dot(c.reciprocalC.xyz,d));
    f-=rint(f);return c.cellA.xyz*f.x+c.cellB.xyz*f.y+c.cellC.xyz*f.z;
}
inline void fail(device Status&s,uint p){
    atomic_fetch_or_explicit(&s.flags,4u,memory_order_relaxed);
    atomic_fetch_min_explicit(&s.firstParticle,p,memory_order_relaxed);
    atomic_fetch_add_explicit(&s.violationCount,1u,memory_order_relaxed);
}

// One owner per molecule walks a breadth-first spanning tree. Unwrapping every
// atom independently against the first atom fails for molecules larger than
// half the box; every tree edge here is an explicit local topology relation.
kernel void nvivo_md_barostat_centers(
    device const float4*positions[[buffer(0)]],device const float4*dynamics[[buffer(1)]],
    device const uint*offsets[[buffer(2)]],device const uint*members[[buffer(3)]],
    device const uint*parents[[buffer(4)]],device const uint*edgeOffsets[[buffer(5)]],
    device const uint2*edges[[buffer(6)]],device float4*unwrapped[[buffer(7)]],
    device float4*centers[[buffer(8)]],device Status&status[[buffer(9)]],
    constant Command&c[[buffer(10)]],uint component[[thread_position_in_grid]]){
    if(component>=c.componentCount)return;
    float3 sum=0,correction=0;float mass=0,massCorrection=0;
    for(uint j=offsets[component];j<offsets[component+1];++j){
        uint particle=members[j],parent=parents[particle];
        float3 p=parent==0xffffffffu?positions[particle].xyz:
            unwrapped[parent].xyz+image(positions[particle].xyz-positions[parent].xyz,c);
        float m=dynamics[particle].x;
        if(!all(isfinite(p))||!isfinite(m)||m<=0){fail(status,particle);return;}
        unwrapped[particle]=float4(p,0);
        float3 term=m*p-correction,next=sum+term;correction=(next-sum)-term;sum=next;
        float dm=m-massCorrection,nm=mass+dm;massCorrection=(nm-mass)-dm;mass=nm;
    }
    if(!isfinite(mass)||mass<=0||!all(isfinite(sum))){fail(status,0);return;}
    // A nonzero winding cycle is not an isolated molecule whose center can be
    // scaled while preserving its internal coordinates. Reject it explicitly.
    for(uint j=edgeOffsets[component];j<edgeOffsets[component+1];++j){
        uint2 edge=edges[j];
        float3 expected=image(positions[edge.x].xyz-positions[edge.y].xyz,c);
        float error=length(unwrapped[edge.x].xyz-unwrapped[edge.y].xyz-expected);
        if(!isfinite(error)||error>c.toleranceNM){fail(status,edge.x);return;}
    }
    centers[component]=float4(sum/mass,mass);
}

kernel void nvivo_md_barostat_scale(
    device float4*positions[[buffer(0)]],device const float4*dynamics[[buffer(1)]],
    device const uint*componentIndices[[buffer(2)]],device const float4*centers[[buffer(3)]],
    device const float4*unwrapped[[buffer(4)]],device Status&status[[buffer(5)]],
    constant Command&c[[buffer(6)]],uint particle[[thread_position_in_grid]]){
    if(particle>=c.particleCount||dynamics[particle].y==0)return;
    if(atomic_load_explicit(&status.flags,memory_order_relaxed)!=0)return;
    float3 center=centers[componentIndices[particle]].xyz;
    float3 p=unwrapped[particle].xyz+(c.scale-1.0f)*center;
    float3 f=float3(dot(c.reciprocalA.xyz,p),dot(c.reciprocalB.xyz,p),dot(c.reciprocalC.xyz,p))/c.scale;
    f-=floor(f);
    float3 result=c.scale*(c.cellA.xyz*f.x+c.cellB.xyz*f.y+c.cellC.xyz*f.z);
    if(!all(isfinite(result))){fail(status,particle);return;}
    positions[particle]=float4(result,0);
    // Virtual positions are reconstructed from their parents by the next pass.
}
} // namespace nvivo_md_barostat
