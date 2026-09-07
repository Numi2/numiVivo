#include <metal_stdlib>
#include "NumiVivoErrorFunctions.metalh"
#include "NumiVivoMDPeriodicGeometry.metalh"
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
    return nvivo_md_periodic::minimumImage(d,c.cellA.xyz,c.cellB.xyz,c.cellC.xyz,
                                         c.reciprocalA.xyz,c.reciprocalB.xyz,c.reciprocalC.xyz);
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

// Real-space halfSpan of Ewald. The reciprocal engine supplies the complementary erf term.
[[host_name("nvivo_pme_realspace_neighbor")]] kernel void nvivo_pme_realspace_neighbor(device const float4*p[[buffer(0)]],
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
        if(c.cellA.w>0 && r>c.cellA.w){float t=(r-c.cellA.w)/(c.cutoffNM-c.cellA.w),t2=t*t,t3=t2*t;float sw=1-10*t3+15*t3*t-6*t3*t2,ds=(-30*t2+60*t3-30*t3*t)/(c.cutoffNM-c.cellA.w);ljScale=sw*ljScale-lj*ds*ir;lj*=sw;}
        float qq=dyn[g].z*dyn[j].z;float ce=0,cf=0;
        if(cs!=0&&qq!=0){float br=beta*r;float erfcv=nvivo_math::erfc(br);float gaussian=exp(-br*br);
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
// Shared sixth-order cardinal assignment. Components are the value and its
// first three derivatives with respect to the fractional grid coordinate.
// Charge gather differentiates the same assignment used to spread its energy;
// multipoles also consume the higher derivatives without a separate basis.
constant int m6coeff[6][6] = {{1,-5,10,-10,5,-1},{26,-50,20,20,-20,5},{66,0,-60,0,30,-10},{26,50,20,-20,-20,10},{1,5,10,10,5,-5},{0,0,0,0,0,1}};
inline void bspline6(float t,thread float4 out[6]) {
    for(uint j=0;j<6;++j){float v=float(m6coeff[j][5]),d=0,d2=0,d3=0;
        for(int power=4;power>=0;--power){d3=d3*t+3*d2;d2=d2*t+2*d;d=d*t+v;v=v*t+float(m6coeff[j][power]);}
        out[j]=float4(v,d,d2,d3)/120.0f;
    }
}
inline float3 fractional(float3 p,constant PMECommand&c){float3 f=float3(dot(c.reciprocalA.xyz,p),dot(c.reciprocalB.xyz,p),dot(c.reciprocalC.xyz,p));return f-floor(f);}

[[host_name("nvivo_pme_clear_grid")]] kernel void nvivo_pme_clear_grid(device float2*grid[[buffer(0)]],constant PMECommand&c[[buffer(1)]],uint gid[[thread_position_in_grid]]){if(gid<c.gridPointCount)grid[gid]=0;}

[[host_name("nvivo_pme_spread")]] kernel void nvivo_pme_spread(device const float4*positions[[buffer(0)]],device const float4*dynamics[[buffer(1)]],device atomic_uint*gridRealBits[[buffer(2)]],device Status&s[[buffer(3)]],constant PMECommand&c[[buffer(4)]],uint gid[[thread_position_in_grid]]){
    if(gid>=c.particleCount)return;float q=dynamics[gid].z;if(q==0)return;float3 f=fractional(positions[gid].xyz,c);
    float3 u=f*float3(c.gridX,c.gridY,c.gridZ);int3 base=int3(floor(u))-2;float3 t=u-floor(u);
    float4 wx[6],wy[6],wz[6];bspline6(t.x,wx);bspline6(t.y,wy);bspline6(t.z,wz);
    for(uint iz=0;iz<6;++iz)for(uint iy=0;iy<6;++iy)for(uint ix=0;ix<6;++ix){uint x=wrapIndex(base.x+int(ix),c.gridX),y=wrapIndex(base.y+int(iy),c.gridY),z=wrapIndex(base.z+int(iz),c.gridZ);float value=q*wx[ix].x*wy[iy].x*wz[iz].x;atomicAddFloat(&gridRealBits[2ul*gridIndex(x,y,z,c)],value);}
}

inline uint reverseBitsN(uint v,uint bits){uint r=0;for(uint i=0;i<bits;++i){r=(r<<1)|(v&1u);v>>=1;}return r;}
inline uint log2Exact(uint n){return 31u-clz(n);}
inline uint3 decodeGrid(uint index,constant PMECommand&c){uint x=index%c.gridX;uint q=index/c.gridX;uint y=q%c.gridY;uint z=q/c.gridY;return uint3(x,y,z);}

[[host_name("nvivo_pme_bit_reverse")]] kernel void nvivo_pme_bit_reverse(device const float2*src[[buffer(0)]],device float2*dst[[buffer(1)]],constant PMECommand&c[[buffer(2)]],uint gid[[thread_position_in_grid]]){
    if(gid>=c.gridPointCount)return;uint3 q=decodeGrid(gid,c);uint n=c.axis==0?c.gridX:(c.axis==1?c.gridY:c.gridZ);uint bits=log2Exact(n);uint coordinate=c.axis==0?q.x:(c.axis==1?q.y:q.z);uint r=reverseBitsN(coordinate,bits);if(c.axis==0)q.x=r;else if(c.axis==1)q.y=r;else q.z=r;dst[gridIndex(q.x,q.y,q.z,c)]=src[gid];
}

[[host_name("nvivo_pme_fft_stage")]] kernel void nvivo_pme_fft_stage(device const float2*src[[buffer(0)]],device float2*dst[[buffer(1)]],constant PMECommand&c[[buffer(2)]],uint gid[[thread_position_in_grid]]){
    uint n=c.axis==0?c.gridX:(c.axis==1?c.gridY:c.gridZ);uint lineCount=c.gridPointCount/n;uint butterflies=lineCount*(n>>1);if(gid>=butterflies)return;
    uint line=gid/(n>>1),b=gid%(n>>1);uint span=1u<<(c.stage+1u),halfSpan=span>>1u,group=b/halfSpan,j=b%halfSpan;uint i0=group*span+j,i1=i0+halfSpan;
    uint3 q;if(c.axis==0){uint y=line%c.gridY,z=line/c.gridY;q=uint3(i0,y,z);}else if(c.axis==1){uint x=line%c.gridX,z=line/c.gridX;q=uint3(x,i0,z);}else{uint x=line%c.gridX,y=line/c.gridX;q=uint3(x,y,i0);}
    uint idx0=gridIndex(q.x,q.y,q.z,c);if(c.axis==0)q.x=i1;else if(c.axis==1)q.y=i1;else q.z=i1;uint idx1=gridIndex(q.x,q.y,q.z,c);
    float sign=c.inverse!=0?1.0f:-1.0f;float angle=sign*6.283185307179586f*float(j)/float(span);float2 tw=float2(cos(angle),sin(angle));float2 v=src[idx1];float2 t=float2(tw.x*v.x-tw.y*v.y,tw.x*v.y+tw.y*v.x),u=src[idx0];dst[idx0]=u+t;dst[idx1]=u-t;
}

inline int signedMode(uint index,uint n){return index<=n/2?int(index):int(index)-int(n);}
inline float sincPi(float x){if(abs(x)<1e-7f)return 1.0f;float p=3.141592653589793f*x;return sin(p)/p;}
// The centered nodal weights of M6 are [1,26,66,26,1]/120.
// Their discrete Fourier modulus, not the continuum sinc transform, matches
// the cardinal interpolation. It is positive, including the Nyquist mode.
inline float cardinal6Modulus(int mode,uint count){
    float theta=6.283185307179586f*float(mode)/float(count);
    return (66.0f+52.0f*cos(theta)+2.0f*cos(2.0f*theta))/120.0f;
}
[[host_name("nvivo_pme_influence")]] kernel void nvivo_pme_influence(device const float2*chargeK[[buffer(0)]],device float2*potentialK[[buffer(1)]],constant PMECommand&c[[buffer(2)]],uint gid[[thread_position_in_grid]]){
    if(gid>=c.gridPointCount)return;uint3 q=decodeGrid(gid,c);int mx=signedMode(q.x,c.gridX),my=signedMode(q.y,c.gridY),mz=signedMode(q.z,c.gridZ);if(mx==0&&my==0&&mz==0){potentialK[gid]=0;return;}
    float3 k=6.283185307179586f*(float(mx)*c.reciprocalA.xyz+float(my)*c.reciprocalB.xyz+float(mz)*c.reciprocalC.xyz);float k2=dot(k,k);float beta=c.betaPerNM;
    float b=cardinal6Modulus(mx,c.gridX)*cardinal6Modulus(my,c.gridY)*cardinal6Modulus(mz,c.gridZ);
    float deconv=b*b;
    float influence=float(c.gridPointCount)*(c.coulombPrefactor/c.volumeNM3)*12.566370614359172f*exp(-k2/(4.0f*beta*beta))/(k2*deconv);potentialK[gid]=chargeK[gid]*influence;
}

[[host_name("nvivo_pme_scale_inverse")]] kernel void nvivo_pme_scale_inverse(device float2*grid[[buffer(0)]],constant PMECommand&c[[buffer(1)]],uint gid[[thread_position_in_grid]]){if(gid<c.gridPointCount)grid[gid]*=c.inverseGridCount;}

[[host_name("nvivo_pme_gather")]] kernel void nvivo_pme_gather(device const float4*positions[[buffer(0)]],device const float4*dynamics[[buffer(1)]],device const float2*potential[[buffer(2)]],device float4*forceEnergy[[buffer(3)]],device Status&s[[buffer(4)]],constant PMECommand&c[[buffer(5)]],uint gid[[thread_position_in_grid]]){
    if(gid>=c.particleCount)return;float qcharge=dynamics[gid].z;if(qcharge==0)return;float3 f=fractional(positions[gid].xyz,c),u=f*float3(c.gridX,c.gridY,c.gridZ);int3 base=int3(floor(u))-2;float3 t=u-floor(u);
    float4 wx[6],wy[6],wz[6];bspline6(t.x,wx);bspline6(t.y,wy);bspline6(t.z,wz);float phi=0,dux=0,duy=0,duz=0;
    for(uint iz=0;iz<6;++iz)for(uint iy=0;iy<6;++iy)for(uint ix=0;ix<6;++ix){uint x=wrapIndex(base.x+int(ix),c.gridX),y=wrapIndex(base.y+int(iy),c.gridY),z=wrapIndex(base.z+int(iz),c.gridZ);float v=potential[gridIndex(x,y,z,c)].x;phi+=wx[ix].x*wy[iy].x*wz[iz].x*v;dux+=wx[ix].y*wy[iy].x*wz[iz].x*v;duy+=wx[ix].x*wy[iy].y*wz[iz].x*v;duz+=wx[ix].x*wy[iy].x*wz[iz].y*v;}
    float3 grad=float(c.gridX)*dux*c.reciprocalA.xyz+float(c.gridY)*duy*c.reciprocalB.xyz+float(c.gridZ)*duz*c.reciprocalC.xyz;float self=-c.coulombPrefactor*c.betaPerNM*0.5641895835477563f*qcharge*qcharge;float4 add=float4(-qcharge*grad,0.5f*qcharge*phi+self);float4 value=forceEnergy[gid]+add;if(!all(isfinite(value))){fail(s,statusNonFinite,gid);return;}forceEnergy[gid]=value;
}
} // namespace nvivo_pme


// Sixth-order derivative-consistent multipolar assignment. Fourier transforms
// remain the same mdPMEBitReverse/mdPMEFFTStage kernels used by charge PME.
namespace nvivo_pme {
struct MultipoleSource { float4 positionCharge,dipole,q0,q1; };
struct MultipoleMeshExtra { uint4 halfWidths; };
inline float weightDerivative(uint3 order,float4 x,float4 y,float4 z){return x[order.x]*y[order.y]*z[order.z];}
constant uint3 momentPowers[10]={uint3(0),uint3(1,0,0),uint3(0,1,0),uint3(0,0,1),uint3(2,0,0),uint3(1,1,0),uint3(1,0,1),uint3(0,2,0),uint3(0,1,1),uint3(0,0,2)};
constant float momentWeights[10]={1,1,1,1,0.5f,1,1,0.5f,1,0.5f};
inline void momentsOf(MultipoleSource s,thread float m[10]){
    m[0]=s.positionCharge.w;m[1]=s.dipole.x;m[2]=s.dipole.y;m[3]=s.dipole.z;
    m[4]=s.q0.x;m[5]=s.q0.y;m[6]=s.q0.z;m[7]=s.q0.w;m[8]=s.q1.x;m[9]=s.q1.y;
}
[[host_name("nvivo_pme_multipole_spread")]] kernel void multipoleSpread(
    device const MultipoleSource*sources[[buffer(0)]],device atomic_uint*grid[[buffer(1)]],
    constant PMECommand&c[[buffer(2)]],uint gid[[thread_position_in_grid]]){
    if(gid>=c.particleCount)return;MultipoleSource s=sources[gid];float3 u=s.positionCharge.xyz;
    int3 base=int3(floor(u))-2;float3 t=u-floor(u);float4 wx[6],wy[6],wz[6];
    bspline6(t.x,wx);bspline6(t.y,wy);bspline6(t.z,wz);float m[10];momentsOf(s,m);
    for(uint z=0;z<6;++z)for(uint y=0;y<6;++y)for(uint x=0;x<6;++x){float value=0;
        for(uint a=0;a<10;++a)value+=m[a]*momentWeights[a]*weightDerivative(momentPowers[a],wx[x],wy[y],wz[z]);
        uint index=gridIndex(wrapIndex(base.x+int(x),c.gridX),wrapIndex(base.y+int(y),c.gridY),wrapIndex(base.z+int(z),c.gridZ),c);
        atomicAddFloat(&grid[2ul*index],value);
    }
}
[[host_name("nvivo_pme_multipole_influence")]] kernel void multipoleInfluence(
    device const float2*rho[[buffer(0)]],device float2*phi[[buffer(1)]],
    device float4*energyStress[[buffer(2)]],constant MultipoleMeshExtra&e[[buffer(3)]],
    constant PMECommand&c[[buffer(4)]],uint gid[[thread_position_in_grid]]){
    if(gid>=c.gridPointCount)return;uint3 q=decodeGrid(gid,c);
    int3 m=int3(signedMode(q.x,c.gridX),signedMode(q.y,c.gridY),signedMode(q.z,c.gridZ));
    // Omit Nyquist planes; their signed representatives are not a conjugate pair
    // on a skew cell. The declared active mode widths must lie below Nyquist.
    if(all(m==int3(0))||any(abs(m)>int3(e.halfWidths.xyz))||q.x==c.gridX/2||q.y==c.gridY/2||q.z==c.gridZ/2){
        phi[gid]=0;energyStress[2ul*gid]=0;energyStress[2ul*gid+1]=0;return;
    }
    float3 k=6.283185307179586f*(float(m.x)*c.reciprocalA.xyz+float(m.y)*c.reciprocalB.xyz+float(m.z)*c.reciprocalC.xyz);
    float k2=dot(k,k),a2=c.betaPerNM*c.betaPerNM;
    float b=sincPi(float(m.x)/float(c.gridX))*sincPi(float(m.y)/float(c.gridY))*sincPi(float(m.z)/float(c.gridZ));
    float b2=b*b,b4=b2*b2,b12=b4*b4*b4;
    float green=(12.566370614359172f/c.volumeNM3)*exp(-k2/(4*a2))/(k2*b12);
    phi[gid]=rho[gid]*(float(c.gridPointCount)*green);
    float energy=0.5f*green*dot(rho[gid],rho[gid]);float scale=2/k2+1/(2*a2);
    energyStress[2ul*gid]=float4(energy,energy*(-1+scale*k.x*k.x),energy*scale*k.x*k.y,energy*scale*k.x*k.z);
    energyStress[2ul*gid+1]=float4(energy*(-1+scale*k.y*k.y),energy*scale*k.y*k.z,energy*(-1+scale*k.z*k.z),0);
}
[[host_name("nvivo_pme_multipole_gather")]] kernel void multipoleGather(
    device const MultipoleSource*sources[[buffer(0)]],device const float2*potential[[buffer(1)]],
    device float4*out[[buffer(2)]],constant PMECommand&c[[buffer(3)]],uint gid[[thread_position_in_grid]]){
    if(gid>=c.particleCount)return;MultipoleSource s=sources[gid];float3 u=s.positionCharge.xyz;
    int3 base=int3(floor(u))-2;float3 t=u-floor(u);float4 wx[6],wy[6],wz[6];
    bspline6(t.x,wx);bspline6(t.y,wy);bspline6(t.z,wz);float m[10],lambda[10];momentsOf(s,m);
    for(uint a=0;a<10;++a)lambda[a]=0;float3 gradient=0;
    for(uint z=0;z<6;++z)for(uint y=0;y<6;++y)for(uint x=0;x<6;++x){
        uint index=gridIndex(wrapIndex(base.x+int(x),c.gridX),wrapIndex(base.y+int(y),c.gridY),wrapIndex(base.z+int(z),c.gridZ),c);
        float value=potential[index].x;
        for(uint a=0;a<10;++a){uint3 n=momentPowers[a];float f=momentWeights[a]*value;
            lambda[a]+=f*weightDerivative(n,wx[x],wy[y],wz[z]);
            for(uint b=0;b<3;++b){uint3 dn=n;dn[b]++;gradient[b]+=m[a]*f*weightDerivative(dn,wx[x],wy[y],wz[z]);}
        }
    }
    out[4ul*gid]=float4(-gradient,lambda[0]);out[4ul*gid+1]=float4(lambda[1],lambda[2],lambda[3],lambda[4]);
    out[4ul*gid+2]=float4(lambda[5],lambda[6],lambda[7],lambda[8]);out[4ul*gid+3]=float4(lambda[9],0,0,0);
}
} // namespace nvivo_pme
