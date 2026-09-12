#include "NumiVivoCore/NumiVivoOmicsMNN.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#pragma clang fp contract(off) reassociate(off)

extern "C" int32_t nvivo_omics_mnn_exact(
 const double *x,const intptr_t *levels,uint32_t n,uint32_t d,uint32_t k,uint32_t b,
 uint64_t xc,uint64_t lc,uint64_t budget,intptr_t *lower,double *ld,
 intptr_t *upper,double *ud,uint64_t nc,intptr_t *ln,intptr_t *un,uint64_t cc,
 uint64_t *evaluated,int32_t (*cancel)(void *),void *context) {
 if(!evaluated)return 1;
 *evaluated=0;
 if(!x||!levels||!lower||!ld||!upper||!ud||!ln||!un||n<2||n>2000000||d<1||d>64||k<1||k>100||b<2||b>128||xc<uint64_t(n)*d||lc<n||nc<uint64_t(n)*k||cc<n)return 1;
 if(cancel&&cancel(context))return 5;
 uint64_t counts[128]={},pairs=0,preceding=0;
 for(uint32_t i=0;i<n;++i){
  if((i&16383)==0&&cancel&&cancel(context))return 5;
  if(levels[i]<0||levels[i]>=b)return 1;
  ++counts[levels[i]];
 }
 for(uint32_t v=0;v<b;++v){if(counts[v]<k)return 1;pairs+=preceding*counts[v];preceding+=counts[v];}
 if(pairs*d>budget)return 2;
 for(uint64_t i=0;i<uint64_t(n)*d;++i){
  if((i&16383)==0&&cancel&&cancel(context))return 5;
  if(!std::isfinite(x[i]))return 4;
 }
 std::fill_n(lower,uint64_t(n)*k,intptr_t(-1));std::fill_n(upper,uint64_t(n)*k,intptr_t(-1));
 std::fill_n(ld,uint64_t(n)*k,std::numeric_limits<double>::infinity());
 std::fill_n(ud,uint64_t(n)*k,std::numeric_limits<double>::infinity());
 std::fill_n(ln,n,intptr_t(0));std::fill_n(un,n,intptr_t(0));
 auto worse=[&](intptr_t a,double ad,intptr_t z,double zd){
  return ad>zd||(ad==zd&&(levels[a]>levels[z]||(levels[a]==levels[z]&&a>z)));
 };
 auto offer=[&](uint32_t row,uint32_t candidate,double distance,intptr_t *ids,double *ds,intptr_t *sizes){
  size_t base=size_t(row)*k;
  if(sizes[row]<k){
   size_t slot=size_t(sizes[row]++);ids[base+slot]=candidate;ds[base+slot]=distance;
   while(slot){size_t parent=(slot-1)/2;if(!worse(ids[base+slot],ds[base+slot],ids[base+parent],ds[base+parent]))break;
    std::swap(ids[base+slot],ids[base+parent]);std::swap(ds[base+slot],ds[base+parent]);slot=parent;}
  }else{
   if(!worse(ids[base],ds[base],candidate,distance))return;
   ids[base]=candidate;ds[base]=distance;size_t slot=0;
   while(slot*2+1<k){size_t child=slot*2+1;
    if(child+1<k&&worse(ids[base+child+1],ds[base+child+1],ids[base+child],ds[base+child]))++child;
    if(!worse(ids[base+child],ds[base+child],ids[base+slot],ds[base+slot]))break;
    std::swap(ids[base+slot],ids[base+child]);std::swap(ds[base+slot],ds[base+child]);slot=child;}
  }
 };
 for(uint32_t i=0;i<n;++i){
  if(cancel&&cancel(context))return 5;
  for(uint32_t j=i+1;j<n;++j){
   if((j&4095)==0&&cancel&&cancel(context))return 5;
   if(levels[i]==levels[j])continue;
   double distance=0;
   for(uint32_t c=0;c<d;++c)distance+=std::abs(x[size_t(i)*d+c]-x[size_t(j)*d+c]);
   *evaluated+=d;
   if(!std::isfinite(distance))return 4;
   if(levels[i]<levels[j]){offer(i,j,distance,upper,ud,un);offer(j,i,distance,lower,ld,ln);}
   else{offer(j,i,distance,upper,ud,un);offer(i,j,distance,lower,ld,ln);}
  }
 }
 return cancel&&cancel(context)?5:0;
}
