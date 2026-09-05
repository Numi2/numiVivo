#include <metal_stdlib>
using namespace metal;
namespace nvivo_md_virtual {
struct Command{uint particleCount,typeCount,electrostatics,periodic;uint stepLow,stepHigh,seedLow,seedHigh;float dtPS,cutoffNM,coulombPrefactor,reactionFieldK;float reactionFieldC,minimumDistanceNM,constraintTolerance,langevinA;float targetTemperatureK,boltzmannKJPerMolK;uint neighborCapacity;float neighborRadiusNM;float4 cellA,cellB,cellC,reciprocalA,reciprocalB,reciprocalC;};
struct Status{atomic_uint flags,firstParticle,violationCount,reserved;};
struct LinearVirtualSite{uint2 siteAndCount;uint4 parents;float4 weights;};
struct ParentIncidence{uint siteIndex;float weight;};
static_assert(sizeof(Command)==176,"MD command ABI");
static_assert(sizeof(LinearVirtualSite)==48,"linear virtual site ABI");
static_assert(sizeof(ParentIncidence)==8,"virtual parent incidence ABI");
constant uint statusNonFinite=1u,statusInvalidGeometry=4u;
inline void fail(device Status&s,uint flag,uint p){atomic_fetch_or_explicit(&s.flags,flag,memory_order_relaxed);atomic_fetch_min_explicit(&s.firstParticle,p,memory_order_relaxed);atomic_fetch_add_explicit(&s.violationCount,1u,memory_order_relaxed);}
inline float3 minimumImage(float3 d,constant Command&c){if(c.periodic==0)return d;float3 f=float3(dot(c.reciprocalA.xyz,d),dot(c.reciprocalB.xyz,d),dot(c.reciprocalC.xyz,d));f-=rint(f);return c.cellA.xyz*f.x+c.cellB.xyz*f.y+c.cellC.xyz*f.z;}
inline float3 wrapPosition(float3 p,constant Command&c){if(c.periodic==0)return p;float3 f=float3(dot(c.reciprocalA.xyz,p),dot(c.reciprocalB.xyz,p),dot(c.reciprocalC.xyz,p));f-=floor(f);return c.cellA.xyz*f.x+c.cellB.xyz*f.y+c.cellC.xyz*f.z;}

kernel void nvivo_md_update_virtual_position(device float4*position[[buffer(0)]],device const LinearVirtualSite*sites[[buffer(1)]],device const uint*siteIndexByParticle[[buffer(2)]],device Status&s[[buffer(3)]],constant Command&c[[buffer(4)]],uint gid[[thread_position_in_grid]]){
    if(gid>=c.particleCount)return;uint index=siteIndexByParticle[gid];if(index==0xffffffffu)return;LinearVirtualSite site=sites[index];uint count=site.siteAndCount.y;if(site.siteAndCount.x!=gid||count<2||count>4){fail(s,statusInvalidGeometry,gid);return;}uint first=site.parents[0];if(first>=c.particleCount){fail(s,statusInvalidGeometry,gid);return;}float3 base=position[first].xyz,result=base;for(uint i=1;i<count;++i){uint parent=site.parents[i];if(parent>=c.particleCount){fail(s,statusInvalidGeometry,gid);return;}result+=site.weights[i]*minimumImage(position[parent].xyz-base,c);}if(!all(isfinite(result))){fail(s,statusNonFinite,gid);return;}position[gid]=float4(wrapPosition(result,c),0);
}
kernel void nvivo_md_update_virtual_velocity(device float4*velocity[[buffer(0)]],device const LinearVirtualSite*sites[[buffer(1)]],device const uint*siteIndexByParticle[[buffer(2)]],constant Command&c[[buffer(3)]],uint gid[[thread_position_in_grid]]){
    if(gid>=c.particleCount)return;uint index=siteIndexByParticle[gid];if(index==0xffffffffu)return;LinearVirtualSite site=sites[index];float3 value=0;for(uint i=0;i<site.siteAndCount.y;++i)value+=site.weights[i]*velocity[site.parents[i]].xyz;velocity[gid]=float4(value,0);
}
kernel void nvivo_md_redistribute_virtual_force(device float4*forceEnergy[[buffer(0)]],device const LinearVirtualSite*sites[[buffer(1)]],device const uint*parentOffsets[[buffer(2)]],device const ParentIncidence*incidence[[buffer(3)]],device Status&s[[buffer(4)]],constant Command&c[[buffer(5)]],uint gid[[thread_position_in_grid]]){
    if(gid>=c.particleCount)return;float3 add=0;for(uint i=parentOffsets[gid];i<parentOffsets[gid+1];++i){ParentIncidence edge=incidence[i];LinearVirtualSite site=sites[edge.siteIndex];add+=edge.weight*forceEnergy[site.siteAndCount.x].xyz;}if(!all(isfinite(add))){fail(s,statusNonFinite,gid);return;}forceEnergy[gid].xyz+=add;
}
}