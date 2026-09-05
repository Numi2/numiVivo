#include <metal_stdlib>
using namespace metal;

namespace nvivo_md_barostat {
struct Command {
    uint particleCount;
    uint componentCount;
    float scale;
    float reserved;
    float4 cellA;
    float4 cellB;
    float4 cellC;
    float4 reciprocalA;
    float4 reciprocalB;
    float4 reciprocalC;
};
static_assert(sizeof(Command)==112,"barostat command ABI");
inline float3 fractional(float3 p,constant Command&c){float3 f=float3(dot(c.reciprocalA.xyz,p),dot(c.reciprocalB.xyz,p),dot(c.reciprocalC.xyz,p));return f-floor(f);}
inline float3 wrapWithScaledCell(float3 p,constant Command&c){float inverseScale=1.0f/c.scale;float3 f=float3(dot(c.reciprocalA.xyz,p),dot(c.reciprocalB.xyz,p),dot(c.reciprocalC.xyz,p))*inverseScale;f-=floor(f);return c.scale*(c.cellA.xyz*f.x+c.cellB.xyz*f.y+c.cellC.xyz*f.z);}

kernel void nvivo_md_barostat_centers(device const float4*positions[[buffer(0)]],
                                      device const float4*dynamics[[buffer(1)]],
                                      device const uint*offsets[[buffer(2)]],
                                      device const uint*particles[[buffer(3)]],
                                      device float4*centers[[buffer(4)]],
                                      constant Command&c[[buffer(5)]],
                                      uint gid[[thread_position_in_grid]]){
    if(gid>=c.componentCount)return;float3 sumCos=0,sumSin=0;float totalMass=0;
    for(uint k=offsets[gid];k<offsets[gid+1];++k){uint p=particles[k];float mass=dynamics[p].x;if(!(mass>0))continue;float3 f=fractional(positions[p].xyz,c);float3 angle=6.283185307179586f*f;sumCos+=mass*cos(angle);sumSin+=mass*sin(angle);totalMass+=mass;}
    if(!(totalMass>0)){centers[gid]=float4(0);return;}float3 angle=atan2(sumSin,sumCos);float3 f=angle/6.283185307179586f;f=select(f,f+1.0f,f<0.0f);float3 center=c.cellA.xyz*f.x+c.cellB.xyz*f.y+c.cellC.xyz*f.z;centers[gid]=float4(center,0);
}

kernel void nvivo_md_barostat_scale(device float4*positions[[buffer(0)]],
                                    device const uint*componentIndex[[buffer(1)]],
                                    device const float4*centers[[buffer(2)]],
                                    constant Command&c[[buffer(3)]],
                                    uint gid[[thread_position_in_grid]]){
    if(gid>=c.particleCount)return;uint component=componentIndex[gid];float3 shifted=positions[gid].xyz+(c.scale-1.0f)*centers[component].xyz;positions[gid]=float4(wrapWithScaledCell(shifted,c),0);
}
} // namespace nvivo_md_barostat
