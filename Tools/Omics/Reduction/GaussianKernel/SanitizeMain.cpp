#include "NumiVivoCore/NumiVivoOmicsGaussian.h"
#include <vector>
#include <cmath>
#include <cassert>
#include <cstdio>
#pragma clang fp contract(off) reassociate(off)
static int32_t stop(void *ctx) { return ++*static_cast<int *>(ctx)>=3; }
int main() {
 for (unsigned d: {1u,2u,20u,64u}) for (unsigned a: {1u,255u,256u,257u,513u}) {
  const unsigned q=32;std::vector<double>s(a*d),b(a*d),x(q*d),y(q*d),t(q),w(256*(d+q));
  for(unsigned i=0;i<a*d;++i){s[i]=std::sin(double(i))*.05;b[i]=std::cos(double(i));}
  for(unsigned i=0;i<q*d;++i)x[i]=std::cos(double(i+7))*.05;
  uint64_t evaluated=0;
  auto run=[&](uint64_t cap,uint64_t work, int32_t(*cancel)(void*)=nullptr,void *ctx=nullptr){return nvivo_omics_gaussian_bias(s.data(),b.data(),a,x.data(),q,d,15,s.size(),x.size(),work,y.data(),y.size(),t.data(),t.size(),w.data(),cap,&evaluated,cancel,ctx);};
  assert(run(w.size(),uint64_t(a)*q*d)==0);
  for(unsigned i=0;i<q;++i){std::vector<double>ref(d);double total=0;
   for(unsigned k=0;k<a;++k){double dist=0;for(unsigned j=0;j<d;++j){double diff=x[i*d+j]-s[k*d+j];dist+=diff*diff;}double wt=std::exp(-7.5*dist);total+=wt;for(unsigned j=0;j<d;++j)ref[j]+=wt*b[k*d+j];}
   std::vector<double>legacy(d);double legacyTotal=0;
   assert(nvivo_omics_gaussian_scalar(s.data(),b.data(),a,x.data()+i*d,d,15,s.size(),d,legacy.data(),d,&legacyTotal,nullptr,nullptr)==0);
   assert(legacyTotal==total);for(unsigned j=0;j<d;++j)assert(legacy[j]==ref[j]);
   assert(std::abs(total-t[i])<=1e-10*(1+total));for(unsigned j=0;j<d;++j)assert(std::abs(ref[j]/total-y[i*d+j])<1e-10);
  }
  std::vector<double>changed=x;std::vector<intptr_t>selected(q);for(unsigned i=0;i<q;++i)selected[i]=i;
  NVivoGaussianScalarReport phase{};
  auto apply=[&](uint64_t work, int32_t(*cancel)(void*)=nullptr,void *ctx=nullptr){return nvivo_omics_gaussian_scalar_apply(s.data(),b.data(),a,d,15,s.size(),changed.data(),q,changed.size(),selected.data(),q,work,&phase,cancel,ctx);};
  assert(apply(uint64_t(a)*q*d)==0);
  for(unsigned i=0;i<q*d;++i)assert(std::abs(changed[i]-(x[i]+y[i]))<1e-10);
  assert(phase.zero_weight_rows==0);
  assert(apply(uint64_t(a)*q*d-1)==2);
  selected.back()=q;assert(apply(uint64_t(a)*q*d)==1);selected.back()=q-1;
  int phaseCalls=0;assert(apply(uint64_t(a)*q*d,stop,&phaseCalls)==5);
  assert(run(w.size()-1,uint64_t(a)*q*d)==1);assert(run(w.size(),uint64_t(a)*q*d-1)==2);
  int calls=0;assert(run(w.size(),uint64_t(a)*q*d,stop,&calls)==5);
  s.back()=NAN;assert(run(w.size(),uint64_t(a)*q*d)==4);
 }
 puts("Native Gaussian ASan/UBSan and exhaustive scalar checks passed: 20 boundary cases.");
}
