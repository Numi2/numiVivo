#include <metal_stdlib>
#include "NumiVivoMDPeriodicGeometry.metalh"
using namespace metal;
namespace nvivo_md_virtual {
struct Command{uint particleCount,typeCount,electrostatics,periodic;uint stepLow,stepHigh,seedLow,seedHigh;float dtPS,cutoffNM,coulombPrefactor,reactionFieldK;float reactionFieldC,minimumDistanceNM,constraintTolerance,langevinA;float targetTemperatureK,boltzmannKJPerMolK;uint neighborCapacity;float neighborRadiusNM;float4 cellA,cellB,cellC,reciprocalA,reciprocalB,reciprocalC;};
struct Status{atomic_uint flags,firstParticle,violationCount,reserved;};
struct LinearVirtualSite{uint2 siteAndCount;uint4 parents;float4 weights;};
struct ParentIncidence{uint siteIndex;float weight;};
static_assert(sizeof(Command)==176,"MD command ABI");
static_assert(sizeof(LinearVirtualSite)==48,"linear virtual site ABI");
static_assert(sizeof(ParentIncidence)==8,"virtual parent incidence ABI");
inline void fail(device Status&s,uint flag,uint p){atomic_fetch_or_explicit(&s.flags,flag,memory_order_relaxed);atomic_fetch_min_explicit(&s.firstParticle,p,memory_order_relaxed);atomic_fetch_add_explicit(&s.violationCount,1u,memory_order_relaxed);}
inline float3 image(float3 d,constant Command&c){if(c.periodic==0)return d;return nvivo_md_periodic::minimumImage(d,c.cellA.xyz,c.cellB.xyz,c.cellC.xyz,c.reciprocalA.xyz,c.reciprocalB.xyz,c.reciprocalC.xyz);}
inline float3 wrap(float3 p,constant Command&c){if(c.periodic==0)return p;return nvivo_md_periodic::wrapPosition(p,c.cellA.xyz,c.cellB.xyz,c.cellC.xyz,c.reciprocalA.xyz,c.reciprocalB.xyz,c.reciprocalC.xyz);}

[[host_name("nvivo_md_update_virtual_position")]] kernel void nvivo_md_update_virtual_position(device float4*position[[buffer(0)]],device const LinearVirtualSite*sites[[buffer(1)]],device const uint*siteIndexByParticle[[buffer(2)]],device Status&s[[buffer(3)]],constant Command&c[[buffer(4)]],uint gid[[thread_position_in_grid]]){
    if(gid>=c.particleCount)return;uint index=siteIndexByParticle[gid];if(index==0xffffffffu)return;
    LinearVirtualSite site=sites[index];uint count=site.siteAndCount.y;
    if(site.siteAndCount.x!=gid||count<2||count>4){fail(s,4u,gid);return;}
    uint first=site.parents[0];if(first>=c.particleCount){fail(s,4u,gid);return;}
    float3 base=position[first].xyz,result=base;
    for(uint i=1;i<count;++i){uint parent=site.parents[i];if(parent>=c.particleCount){fail(s,4u,gid);return;}result+=site.weights[i]*image(position[parent].xyz-base,c);}
    if(!all(isfinite(result))){fail(s,1u,gid);return;}
    position[gid]=float4(wrap(result,c),0);
}
[[host_name("nvivo_md_update_virtual_velocity")]] kernel void nvivo_md_update_virtual_velocity(device float4*velocity[[buffer(0)]],device const LinearVirtualSite*sites[[buffer(1)]],device const uint*siteIndexByParticle[[buffer(2)]],constant Command&c[[buffer(3)]],uint gid[[thread_position_in_grid]]){
    if(gid>=c.particleCount)return;uint index=siteIndexByParticle[gid];if(index==0xffffffffu)return;
    LinearVirtualSite site=sites[index];float3 value=0;
    for(uint i=0;i<site.siteAndCount.y;++i)value+=site.weights[i]*velocity[site.parents[i]].xyz;
    velocity[gid]=float4(value,0);
}
[[host_name("nvivo_md_redistribute_virtual_force")]] kernel void nvivo_md_redistribute_virtual_force(device float4*forceEnergy[[buffer(0)]],device const LinearVirtualSite*sites[[buffer(1)]],device const uint*parentOffsets[[buffer(2)]],device const ParentIncidence*incidence[[buffer(3)]],device Status&s[[buffer(4)]],constant Command&c[[buffer(5)]],uint gid[[thread_position_in_grid]]){
    if(gid>=c.particleCount)return;
    uint begin=parentOffsets[gid],end=parentOffsets[gid+1];
    // A virtual site is never a physical parent. Do not issue even an add-zero
    // write to its force slot while physical-parent threads are reading it.
    if(begin==end)return;
    float3 add=0;
    for(uint i=begin;i<end;++i){ParentIncidence edge=incidence[i];LinearVirtualSite site=sites[edge.siteIndex];add+=edge.weight*forceEnergy[site.siteAndCount.x].xyz;}
    float3 value=forceEnergy[gid].xyz+add;
    if(!all(isfinite(value))){fail(s,1u,gid);return;}
    forceEnergy[gid].xyz=value;
}
} // namespace nvivo_md_virtual

// General dependent sites. Construction is forward by depth; force transport is
// reverse by depth. A pass only writes parents of a selected higher-depth child,
// so no force slot is read while that same slot is written in the same pass.
namespace nvivo_md_dependent {
using nvivo_md_virtual::Command;
using nvivo_md_virtual::Status;
using nvivo_md_virtual::image;
using nvivo_md_virtual::wrap;
using nvivo_md_virtual::fail;
struct Site { uint4 identity; uint4 parents; float4 originWeights,xWeights,yWeights,local; };
struct Pass { uint particleCount,siteCount,depth,reserved; };
struct Jacobian { float4 x,y,z; };
static_assert(sizeof(Site)==96,"dependent site ABI");
static_assert(sizeof(Jacobian)==48,"dependent-site Jacobian ABI");
inline float3 basis(uint i) { return i==0 ? float3(1,0,0) : (i==1 ? float3(0,1,0) : float3(0,0,1)); }
[[host_name("nvivo_md_construct_dependent_sites")]] kernel void construct(
    device float4* positions [[buffer(0)]],device const Site* sites [[buffer(1)]],device Jacobian* jacobians [[buffer(2)]],
    device Status& status [[buffer(3)]],constant Command& md [[buffer(4)]],constant Pass& pass [[buffer(5)]],uint gid [[thread_position_in_grid]]) {
    if(gid>=pass.siteCount) return;
    Site s=sites[gid]; if(s.identity.z!=pass.depth) return;
    uint count=s.identity.y,output=s.identity.x;
    if(count==0||count>4||output>=pass.particleCount) {fail(status,4u,output);return;}
    float3 p[4];
    for(uint i=0;i<count;++i) { if(s.parents[i]>=pass.particleCount){fail(status,4u,output);return;} }
    float3 base=positions[s.parents[0]].xyz;
    for(uint i=0;i<count;++i) p[i]=base+image(positions[s.parents[i]].xyz-base,md);
    float3 result=0,u=0,v=0,ex=0,ey=0,ez=0;float nx=0,nz=0;
    if(s.identity.w==0) {
        for(uint i=0;i<count;++i) result+=p[i]*s.originWeights[i];
    } else if(s.identity.w==1) {
        u=p[1]-p[0];v=p[2]-p[0];result=p[0]+s.local.x*u+s.local.y*v+s.local.z*cross(u,v);
    } else if(s.identity.w==2) {
        for(uint i=0;i<count;++i) {result+=p[i]*s.originWeights[i];u+=(p[i]-base)*s.xWeights[i];v+=(p[i]-base)*s.yWeights[i];}
        float3 z=cross(u,v);nx=length(u);nz=length(z);
        if(!(nx>1e-12f)||!(nz>1e-12f*max(1.0f,nx*length(v)))) {fail(status,1u,output);return;}
        ex=u/nx;ez=z/nz;ey=cross(ez,ex);result+=ex*s.local.x+ey*s.local.y+ez*s.local.z;
    } else {fail(status,4u,output);return;}
    if(!all(isfinite(result))) {fail(status,1u,output);return;}
    positions[output]=float4(wrap(result,md),0);
    for(uint i=0;i<count;++i) {
        float3 columns[3];
        for(uint axis=0;axis<3;++axis) {
            float3 e=basis(axis),column=0;
            if(s.identity.w==0) column=e*s.originWeights[i];
            else if(s.identity.w==1) {
                float3 d0=i==0?e:float3(0),d1=i==1?e:float3(0),d2=i==2?e:float3(0),du=d1-d0,dv=d2-d0;
                column=d0+s.local.x*du+s.local.y*dv+s.local.z*(cross(du,v)+cross(u,dv));
            } else {
                float3 dx=e*s.xWeights[i],dy=e*s.yWeights[i],dex=(dx-ex*dot(ex,dx))/nx;
                float3 dz=cross(dx,v)+cross(u,dy),dez=(dz-ez*dot(ez,dz))/nz;
                float3 dey=cross(dez,ex)+cross(ez,dex);
                column=e*s.originWeights[i]+dex*s.local.x+dey*s.local.y+dez*s.local.z;
            }
            if(!all(isfinite(column))) {fail(status,1u,output);return;}
            columns[axis]=column;
        }
        jacobians[gid*4+i]={float4(columns[0],0),float4(columns[1],0),float4(columns[2],0)};
    }
}
[[host_name("nvivo_md_dependent_site_velocities")]] kernel void velocities(
    device float4* velocity [[buffer(0)]],device const Site* sites [[buffer(1)]],device const Jacobian* jacobians [[buffer(2)]],
    device Status& status [[buffer(3)]],constant Pass& pass [[buffer(4)]],uint gid [[thread_position_in_grid]]) {
    if(gid>=pass.siteCount)return;Site s=sites[gid];if(s.identity.z!=pass.depth)return;
    float3 value=0;
    for(uint i=0;i<s.identity.y;++i){Jacobian j=jacobians[gid*4+i];float3 v=velocity[s.parents[i]].xyz;value+=j.x.xyz*v.x+j.y.xyz*v.y+j.z.xyz*v.z;}
    if(!all(isfinite(value))){fail(status,1u,s.identity.x);return;}velocity[s.identity.x]=float4(value,0);
}
[[host_name("nvivo_md_dependent_site_forces")]] kernel void forces(
    device float4* forceEnergy [[buffer(0)]],device const Site* sites [[buffer(1)]],device const Jacobian* jacobians [[buffer(2)]],
    device const uint* offsets [[buffer(3)]],device const uint2* incidences [[buffer(4)]],device Status& status [[buffer(5)]],
    constant Pass& pass [[buffer(6)]],uint gid [[thread_position_in_grid]]) {
    if(gid>=pass.particleCount)return;float3 added=0;bool written=false;
    for(uint edge=offsets[gid];edge<offsets[gid+1];++edge){uint2 relation=incidences[edge];Site s=sites[relation.x];if(s.identity.z!=pass.depth)continue;
        Jacobian j=jacobians[relation.x*4+relation.y];float3 f=forceEnergy[s.identity.x].xyz;
        added+=float3(dot(j.x.xyz,f),dot(j.y.xyz,f),dot(j.z.xyz,f));written=true;}
    if(!written)return;float3 value=forceEnergy[gid].xyz+added;if(!all(isfinite(value))){fail(status,1u,gid);return;}forceEnergy[gid].xyz=value;
}
[[host_name("nvivo_md_zero_dependent_site_forces")]] kernel void zero(
    device float4* forceEnergy [[buffer(0)]],device const Site* sites [[buffer(1)]],constant Pass& pass [[buffer(2)]],uint gid [[thread_position_in_grid]]) {
    if(gid<pass.siteCount)forceEnergy[sites[gid].identity.x].xyz=0; // preserve energy component
}
} // namespace nvivo_md_dependent
