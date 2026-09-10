// Standalone check of the exact development implementation under sanitizers.
#include "LocalNeighbors.cpp"
#include <cstdio>
int main(){
 for(uint32_t n:{2u,17u,120u}){
  uint32_t d=20,k=1;std::vector<double>x(size_t(n)*d),ds(size_t(n)*k*2);std::vector<uint32_t>levels(n);std::vector<int64_t>ids(size_t(n)*k*2);LocalReport r{};
  for(uint32_t i=0;i<n;++i){levels[i]=i%2;for(uint32_t j=0;j<d;++j)x[size_t(i)*d+j]=std::sin(double(i*j));}
  for(int kind=0;kind<2;++kind){
   int code=local_neighbors(x.data(),x.data(),levels.data(),n,n,d,k,2,kind,500000000,ids.data(),ds.data(),ids.size(),&r);if(code)return 1;
  }
  if(n>2 && local_neighbors(x.data(),x.data(),levels.data(),n,n,d,k,2,0,1,ids.data(),ds.data(),ids.size(),&r)!=2)return 2;
  x[0]=NAN;if(local_neighbors(x.data(),x.data(),levels.data(),n,n,d,k,2,0,500000000,ids.data(),ds.data(),ids.size(),&r)!=4)return 3;
 }
 std::puts("passed: ASan/UBSan prefix/suffix and external-query paths, nonfinite and budget failures");
}
