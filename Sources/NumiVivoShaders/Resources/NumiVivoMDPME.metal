#include <metal_stdlib>
using namespace metal;

namespace nvivo_pme {
struct MDCommand {
    uint particleCount,typeCount,electrostatics,periodic;
    uint stepLow,stepHigh,seedLow,seedHigh;
    float dtPS,cutoffNM,coulombPrefactor,reactionFieldK;
    float reactionFieldC,minimumDistanceNM,constraintTolerance,langevinA;
    float targetTemperatureK,boltzmannKJPerMolK;
    uint neighborCapacity; float neighborRadiusNM;
    float4 cellA,cellB,cellC,reciprocalA,reciprocalB,reciprocalC;
};
struct PMECommand {
    uint particleCount;
    uint gridX;
    uint gridY;
    uint gridZ;
    uint gridPointCount;
    uint axis;
    uint stage;
    uint inverse;
    float betaPerNM;
    float volumeNM3;
    float coulombPrefactor;
    float inverseGridCount;
    float4 reciprocalA;
    float4 reciprocalB;
    float4 reciprocalC;
};
struct Status { atomic_uint flags,firstParticle,violationCount,reserved; };
struct PairException { uint2 atoms; float2 scales; float2 overrideC12C6; uint flags; };
constant uint statusNonFinite=1u,statusOverlap=2u;

inline void fail(device Status&s,uint flag,uint particle){
    atomic_fetch_or_explicit(&s.flags,flag,memory_order_relaxed);
    atomic_fetch_min_explicit(&s.firstParticle,particle,memory_order_relaxed);
    atomic_fetch_add_explicit(&s.violationCount,1u,memory_order_relaxed);
}
inline float3 minimumImage(float3 d,constant MDCommand&c){
    float3 f=float3(dot(c.reciprocalA.xyz,d),dot(c.reciprocalB.xyz,d),dot(c.reciprocalC.xyz,d));
    f-=rint(f);return c.cellA.xyz*f.x+c.cellB.xyz*f.y+c.cellC.xyz*f.z;
}
inline int findException(uint owner,uint partner,device const uint*o,device const uint*p,device const uint*i){
    uint lo=o[owner],hi=o[owner+1];while(lo<hi){uint m=lo+(hi-lo)/2,v=p[m];if(v<partner)lo=m+1;else hi=m;}
    return lo<o[owner+1]&&p[lo]==partner?int(i[lo]):-1;
}
inline float2 ljCoefficients(uint owner,uint other,device const uint*types,device const float2*pair,
                             device const PairException*exs,device const uint*eo,device const uint*ep,
                             device const uint*ei,constant MDCommand&c,thread float&cs,thread float&ls){
    cs=1.0f;ls=1.0f;float2 coeff=pair[types[owner]*c.typeCount+types[other]];
    int ex=findException(owner,other,eo,ep,ei);if(ex>=0){PairException v=exs[ex];cs=v.scales.x;ls=v.scales.y;if((v.flags&1u)!=0)coeff=v.overrideC12C6;}return coeff;
}

// Real-space half of Ewald. The reciprocal engine supplies the complementary erf term.
kernel void nvivo_pme_realspace_neighbor(device const float4*p[[buffer(0)]],
                                         device float4*fe[[buffer(1)]],
                                         device const float4*dyn[[buffer(2)]],
                                         device const uint*types[[buffer(3)]],
                                         device const float2*pair[[buffer(4)]],
                                         device const PairException*exs[[buffer(5)]],
                                         device const uint*eo[[buffer(6)]],
                                         device const uint*ep[[buffer(7)]],
                                         device const uint*ei[[buffer(8)]],
                                         device const uint*counts[[buffer(9)]],
                                         device const uint*neighbors[[buffer(10)]],
                                         device Status&s[[buffer(11)]],
                                         constant MDCommand&c[[buffer(12)]],
                                         uint g[[thread_position_in_grid]]){
    if(g>=c.particleCount)return;
    float4 sum=0;float beta=c.reactionFieldK;float cutoff2=c.cutoffNM*c.cutoffNM;
    for(uint n=0;n<counts[g];++n){uint j=neighbors[ulong(g)*c.neighborCapacity+n];
        float3 d=minimumImage(p[g].xyz-p[j].xyz,c);float r2=dot(d,d);if(r2>=cutoff2)continue;
        if(!(r2>c.minimumDistanceNM*c.minimumDistanceNM)||!isfinite(r2)){fail(s,statusOverlap,g);continue;}
        float ir=rsqrt(r2),r=r2*ir,ir2=ir*ir,cs,ls;float2 coeff=ljCoefficients(g,j,types,pair,exs,eo,ep,ei,c,cs,ls);
        float ir6=ir2*ir2*ir2,ir12=ir6*ir6;float lj=ls*(coeff.x*ir12-coeff.y*ir6);
        float ljScale=ls*(12.0f*coeff.x*ir12-6.0f*coeff.y*ir6)*ir2;
        float qq=dyn[g].z*dyn[j].z;float ce=0,cf=0;
        if(cs!=0&&qq!=0){float br=beta*r;float erfcv=erfc(br);float gaussian=exp(-br*br);
            ce=c.coulombPrefactor*qq*erfcv*ir*cs;
            cf=c.coulombPrefactor*qq*cs*(erfcv*ir*ir2+(2.0f*beta*0.5641895835477563f)*gaussian*ir2);
        }
        sum+=float4((ljScale+cf)*d,0.5f*(lj+ce));
    }
    float4 value=fe[g]+sum;if(!all(isfinite(value))){fail(s,statusNonFinite,g);return;}fe[g]=value;
}

inline uint gridIndex(uint x,uint y,uint z,constant PMECommand&c){return (z*c.gridY+y)*c.gridX+x;}
inline uint wrapIndex(int i,uint n){int m=i%int(n);return uint(m<0?m+int(n):m);}
inline float atomicAddFloat(device atomic_uint* address,float value){
    uint old=atomic_load_explicit(address,memory_order_relaxed);
    while(true){float current=as_type<float>(old);uint desired=as_type<uint>(current+value);uint expected=old;
        if(atomic_compare_exchange_weak_explicit(address,&expected,desired,memory_order_relaxed,memory_order_relaxed))return current;old=expected;}
}
inline void bspline4(float t,thread float w[4],thread float dw[4]){
    float omt=1.0f-t,t2=t*t,t3=t2*t,omt2=omt*omt;
    w[0]=omt2*omt/6.0f;w[1]=(3.0f*t3-6.0f*t2+4.0f)/6.0f;w[2]=(-3.0f*t3+3.0f*t2+3.0f*t+1.0f)/6.0f;w[3]=t3/6.0f;
    dw[0]=-0.5f*omt2;dw[1]=1.5f*t2-2.0f*t;dw[2]=-1.5f*t2+t+0.5f;dw[3]=0.5f*t2;
}
inline float3 fractional(float3 p,constant PMECommand&c){float3 f=float3(dot(c.reciprocalA.xyz,p),dot(c.reciprocalB.xyz,p),dot(c.reciprocalC.xyz,p));return f-floor(f);}

kernel void nvivo_pme_clear_grid(device float2*grid[[buffer(0)]],constant PMECommand&c[[buffer(1)]],uint gid[[thread_position_in_grid]]){if(gid<c.gridPointCount)grid[gid]=0;}

kernel void nvivo_pme_spread(device const float4*positions[[buffer(0)]],device const float4*dynamics[[buffer(1)]],device atomic_uint*gridRealBits[[buffer(2)]],device Status&s[[buffer(3)]],constant PMECommand&c[[buffer(4)]],uint gid[[thread_position_in_grid]]){
    if(gid>=c.particleCount)return;float q=dynamics[gid].z;if(q==0)return;float3 f=fractional(positions[gid].xyz,c);
    float3 u=f*float3(c.gridX,c.gridY,c.gridZ);int3 base=int3(floor(u))-1;float3 t=u-floor(u);
    float wx[4],wy[4],wz[4],d[4];bspline4(t.x,wx,d);bspline4(t.y,wy,d);bspline4(t.z,wz,d);
    for(uint iz=0;iz<4;++iz)for(uint iy=0;iy<4;++iy)for(uint ix=0;ix<4;++ix){uint x=wrapIndex(base.x+int(ix),c.gridX),y=wrapIndex(base.y+int(iy),c.gridY),z=wrapIndex(base.z+int(iz),c.gridZ);float value=q*wx[ix]*wy[iy]*wz[iz];atomicAddFloat(&gridRealBits[2ul*gridIndex(x,y,z,c)],value);}
}

inline uint reverseBitsN(uint v,uint bits){uint r=0;for(uint i=0;i<bits;++i){r=(r<<1)|(v&1u);v>>=1;}return r;}
inline uint log2Exact(uint n){return 31u-clz(n);}
inline uint3 decodeGrid(uint index,constant PMECommand&c){uint x=index%c.gridX;uint q=index/c.gridX;uint y=q%c.gridY;uint z=q/c.gridY;return uint3(x,y,z);}

kernel void nvivo_pme_bit_reverse(device const float2*src[[buffer(0)]],device float2*dst[[buffer(1)]],constant PMECommand&c[[buffer(2)]],uint gid[[thread_position_in_grid]]){
    if(gid>=c.gridPointCount)return;uint3 q=decodeGrid(gid,c);uint n=c.axis==0?c.gridX:(c.axis==1?c.gridY:c.gridZ);uint bits=log2Exact(n);uint coordinate=c.axis==0?q.x:(c.axis==1?q.y:q.z);uint r=reverseBitsN(coordinate,bits);if(c.axis==0)q.x=r;else if(c.axis==1)q.y=r;else q.z=r;dst[gridIndex(q.x,q.y,q.z,c)]=src[gid];
}

kernel void nvivo_pme_fft_stage(device const float2*src[[buffer(0)]],device float2*dst[[buffer(1)]],constant PMECommand&c[[buffer(2)]],uint gid[[thread_position_in_grid]]){
    uint n=c.axis==0?c.gridX:(c.axis==1?c.gridY:c.gridZ);uint lineCount=c.gridPointCount/n;uint butterflies=lineCount*(n>>1);if(gid>=butterflies)return;
    uint line=gid/(n>>1),b=gid%(n>>1);uint span=1u<<(c.stage+1u),half=span>>1u,group=b/half,j=b%half;uint i0=group*span+j,i1=i0+half;
    uint3 q;if(c.axis==0){uint y=line%c.gridY,z=line/c.gridY;q=uint3(i0,y,z);}else if(c.axis==1){uint x=line%c.gridX,z=line/c.gridX;q=uint3(x,i0,z);}else{uint x=line%c.gridX,y=line/c.gridX;q=uint3(x,y,i0);}
    uint idx0=gridIndex(q.x,q.y,q.z,c);if(c.axis==0)q.x=i1;else if(c.axis==1)q.y=i1;else q.z=i1;uint idx1=gridIndex(q.x,q.y,q.z,c);
    float sign=c.inverse!=0?1.0f:-1.0f;float angle=sign*6.283185307179586f*float(j)/float(span);float2 tw=float2(cos(angle),sin(angle));float2 v=src[idx1];float2 t=float2(tw.x*v.x-tw.y*v.y,tw.x*v.y+tw.y*v.x),u=src[idx0];dst[idx0]=u+t;dst[idx1]=u-t;
}

inline int signedMode(uint index,uint n){return index<=n/2?int(index):int(index)-int(n);}
inline float sincPi(float x){if(abs(x)<1e-7f)return 1.0f;float p=3.141592653589793f*x;return sin(p)/p;}
kernel void nvivo_pme_influence(device const float2*chargeK[[buffer(0)]],device float2*potentialK[[buffer(1)]],constant PMECommand&c[[buffer(2)]],uint gid[[thread_position_in_grid]]){
    if(gid>=c.gridPointCount)return;uint3 q=decodeGrid(gid,c);int mx=signedMode(q.x,c.gridX),my=signedMode(q.y,c.gridY),mz=signedMode(q.z,c.gridZ);if(mx==0&&my==0&&mz==0){potentialK[gid]=0;return;}
    float3 k=6.283185307179586f*(float(mx)*c.reciprocalA.xyz+float(my)*c.reciprocalB.xyz+float(mz)*c.reciprocalC.xyz);float k2=dot(k,k);float beta=c.betaPerNM;
    float sx=sincPi(float(mx)/float(c.gridX)),sy=sincPi(float(my)/float(c.gridY)),sz=sincPi(float(mz)/float(c.gridZ));float b=sx*sx*sx*sx*sy*sy*sy*sy*sz*sz*sz*sz;float deconv=max(b*b,1e-12f);
    float influence=float(c.gridPointCount)*(c.coulombPrefactor/c.volumeNM3)*12.566370614359172f*exp(-k2/(4.0f*beta*beta))/(k2*deconv);potentialK[gid]=chargeK[gid]*influence;
}

kernel void nvivo_pme_scale_inverse(device float2*grid[[buffer(0)]],constant PMECommand&c[[buffer(1)]],uint gid[[thread_position_in_grid]]){if(gid<c.gridPointCount)grid[gid]*=c.inverseGridCount;}

kernel void nvivo_pme_gather(device const float4*positions[[buffer(0)]],device const float4*dynamics[[buffer(1)]],device const float2*potential[[buffer(2)]],device float4*forceEnergy[[buffer(3)]],device Status&s[[buffer(4)]],constant PMECommand&c[[buffer(5)]],uint gid[[thread_position_in_grid]]){
    if(gid>=c.particleCount)return;float qcharge=dynamics[gid].z;if(qcharge==0)return;float3 f=fractional(positions[gid].xyz,c),u=f*float3(c.gridX,c.gridY,c.gridZ);int3 base=int3(floor(u))-1;float3 t=u-floor(u);
    float wx[4],wy[4],wz[4],dx[4],dy[4],dz[4];bspline4(t.x,wx,dx);bspline4(t.y,wy,dy);bspline4(t.z,wz,dz);float phi=0,dux=0,duy=0,duz=0;
    for(uint iz=0;iz<4;++iz)for(uint iy=0;iy<4;++iy)for(uint ix=0;ix<4;++ix){uint x=wrapIndex(base.x+int(ix),c.gridX),y=wrapIndex(base.y+int(iy),c.gridY),z=wrapIndex(base.z+int(iz),c.gridZ);float v=potential[gridIndex(x,y,z,c)].x;phi+=wx[ix]*wy[iy]*wz[iz]*v;dux+=dx[ix]*wy[iy]*wz[iz]*v;duy+=wx[ix]*dy[iy]*wz[iz]*v;duz+=wx[ix]*wy[iy]*dz[iz]*v;}
    float3 grad=float(c.gridX)*dux*c.reciprocalA.xyz+float(c.gridY)*duy*c.reciprocalB.xyz+float(c.gridZ)*duz*c.reciprocalC.xyz;float self=-c.coulombPrefactor*c.betaPerNM*0.5641895835477563f*qcharge*qcharge;float4 add=float4(-qcharge*grad,0.5f*qcharge*phi+self);float4 value=forceEnergy[gid]+add;if(!all(isfinite(value))){fail(s,statusNonFinite,gid);return;}forceEnergy[gid]=value;
}
} // namespace nvivo_pme
