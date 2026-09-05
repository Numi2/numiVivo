#include <metal_stdlib>
using namespace metal;
namespace nvivo_md_min {
struct Command {
    uint particleCount,typeCount,electrostatics,periodic;
    uint stepLow,stepHigh,seedLow,seedHigh;
    float dtPS,cutoffNM,coulombPrefactor,reactionFieldK;
    float reactionFieldC,minimumDistanceNM,constraintTolerance,langevinA;
    float targetTemperatureK,boltzmannKJPerMolK;
    uint neighborCapacity;float neighborRadiusNM;
    float4 cellA,cellB,cellC,reciprocalA,reciprocalB,reciprocalC;
};
struct MinCommand {uint particleCount;float stepScale,maxDisplacementNM,reserved;};
struct ReduceCommand {uint sourceCount,reserved0,reserved1,reserved2;};
static_assert(sizeof(Command)==176,"MD minimization command ABI");
inline float3 wrap(float3 p,constant Command&c){
    if(c.periodic==0)return p;
    float3 f=float3(dot(c.reciprocalA.xyz,p),dot(c.reciprocalB.xyz,p),dot(c.reciprocalC.xyz,p));
    f-=floor(f);return c.cellA.xyz*f.x+c.cellB.xyz*f.y+c.cellC.xyz*f.z;
}

// Energy is counted on every particle, including virtual sites. The stationarity
// residual is a projected force on independent massive particles only.
kernel void nvivo_md_minimize_terms(
    device const float4*forceEnergy[[buffer(0)]],device const float4*direction[[buffer(1)]],
    device const float4*dynamics[[buffer(2)]],device float2*terms[[buffer(3)]],
    constant Command&c[[buffer(4)]],uint g[[thread_position_in_grid]]){
    if(g>=c.particleCount)return;
    float residual=dynamics[g].y>0?length(direction[g].xyz):0.0f;
    terms[g]=float2(forceEnergy[g].w,residual);
}
kernel void nvivo_md_minimize_reduce_stage(
    device const float2*source[[buffer(0)]],device float2*destination[[buffer(1)]],
    constant ReduceCommand&c[[buffer(2)]],uint g[[thread_position_in_grid]]){
    uint first=2u*g;if(first>=c.sourceCount)return;
    float2 a=source[first],b=first+1u<c.sourceCount?source[first+1u]:float2(0);
    destination[g]=float2(a.x+b.x,max(a.y,b.y));
}
kernel void nvivo_md_minimize_position(
    device float4*position[[buffer(0)]],device const float4*direction[[buffer(1)]],
    device const float4*dynamics[[buffer(2)]],constant Command&md[[buffer(3)]],
    constant MinCommand&c[[buffer(4)]],uint g[[thread_position_in_grid]]){
    if(g>=c.particleCount||dynamics[g].y==0)return;
    float3 delta=c.stepScale*direction[g].xyz;
    float norm=length(delta);
    if(norm>c.maxDisplacementNM)delta*=c.maxDisplacementNM/norm;
    position[g]=float4(wrap(position[g].xyz+delta,md),0);
}
kernel void nvivo_md_zero_velocity(device float4*velocity[[buffer(0)]],
    constant Command&c[[buffer(1)]],uint g[[thread_position_in_grid]]){
    if(g<c.particleCount)velocity[g]=0;
}
} // namespace nvivo_md_min
