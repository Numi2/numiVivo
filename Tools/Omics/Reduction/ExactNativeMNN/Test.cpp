#include "NumiVivoCore/NumiVivoOmicsMNN.h"
#include <vector>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <iostream>
struct Stop {int calls=0,limit=0;};
int32_t stop(void *p){auto&s=*static_cast<Stop*>(p);return ++s.calls>=s.limit;}
int main(){
 const uint32_t n=6,d=2,k=2,b=3;std::vector<double>x={0,0,0,0,1,0,2,0,1,0,2,0};std::vector<intptr_t>levels={2,0,1,2,0,1};
 std::vector<intptr_t>lo(n*k),up(n*k),lc(n),uc(n);std::vector<double>ld(n*k),ud(n*k);uint64_t ev=0;
 auto run=[&](uint64_t budget=24,uint64_t capacity=12,Stop*s=nullptr){return nvivo_omics_mnn_exact(x.data(),levels.data(),n,d,k,b,x.size(),levels.size(),budget,lo.data(),ld.data(),up.data(),ud.data(),capacity,lc.data(),uc.data(),lc.size(),&ev,s?stop:nullptr,s);};
 assert(run()==0&&ev==24);
 for(int dir=0;dir<2;++dir)for(uint32_t i=0;i<n;++i){
  std::vector<std::pair<double,intptr_t>>want,got;
  for(uint32_t j=0;j<n;++j)if(dir?levels[j]>levels[i]:levels[j]<levels[i])want.push_back({std::abs(x[2*i]-x[2*j])+std::abs(x[2*i+1]-x[2*j+1]),j});
  auto less=[&](auto a,auto z){return a.first<z.first||(a.first==z.first&&(levels[a.second]<levels[z.second]||(levels[a.second]==levels[z.second]&&a.second<z.second)));};
  std::sort(want.begin(),want.end(),less);if(want.size()>k)want.resize(k);
  auto&ids=dir?up:lo;auto&ds=dir?ud:ld;auto&counts=dir?uc:lc;
  assert(counts[i]==intptr_t(want.size()));for(intptr_t j=0;j<counts[i];++j)got.push_back({ds[i*k+j],ids[i*k+j]});std::sort(got.begin(),got.end(),less);assert(got==want);
 }
 assert(run(23)==2&&ev==0);assert(run(24,11)==1);
 Stop immediate{0,1};assert(run(24,12,&immediate)==5&&ev==0);
 Stop during{0,5};assert(run(24,12,&during)==5&&ev>0&&ev<24);
 auto saved=x;x[0]=std::numeric_limits<double>::quiet_NaN();assert(run()==4);x=saved;
 x[0]=std::numeric_limits<double>::max();x[2]=-std::numeric_limits<double>::max();assert(run()==4);x=saved;
 levels[0]=-1;assert(run()==1);
 std::cout<<"PASS interleaved-level exact ties, budget/capacity, immediate/mid-matching cancellation, nonfinite/overflow and invalid levels\n";
}
