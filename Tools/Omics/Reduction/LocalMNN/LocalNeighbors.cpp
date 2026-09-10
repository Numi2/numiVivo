// Development bridge over the same vendored HNSW owner; no production API.
#include "ThirdParty/hnswlib/hnswlib.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <vector>
struct LocalReport { uint64_t construction,query,indexBytes; };
namespace {
class Space final: public hnswlib::SpaceInterface<double> {
 public:
    const double *source,*queries;uint32_t n,q,d;bool manhattan;uint64_t evaluations=0,maximum;int error=0;
    Space(const double*s,const double*x,uint32_t n,uint32_t q,uint32_t d,bool l1,uint64_t max):source(s),queries(x),n(n),q(q),d(d),manhattan(l1),maximum(max){}
    static double distance(const void *a,const void *b,const void *ctx) {
#pragma clang fp contract(off)
        auto &s=*static_cast<Space*>(const_cast<void*>(ctx));
        if(s.error)return INFINITY;
        if(s.evaluations>=s.maximum){s.error=2;return INFINITY;}++s.evaluations;
        uint32_t i,j;std::memcpy(&i,a,4);std::memcpy(&j,b,4);
        if(i>=s.n+s.q || j>=s.n+s.q){s.error=1;return INFINITY;}
        const double *x=i<s.n?s.source+size_t(i)*s.d:s.queries+size_t(i-s.n)*s.d;
        const double *y=j<s.n?s.source+size_t(j)*s.d:s.queries+size_t(j-s.n)*s.d;
        double value=0;
        for(uint32_t c=0;c<s.d;++c){double v=x[c]-y[c];value+=s.manhattan?std::abs(v):v*v;}
        if(!std::isfinite(value)){s.error=4;return INFINITY;}return value;
    }
    size_t get_data_size() override{return 4;}
    hnswlib::DISTFUNC<double> get_dist_func() override{return distance;}
    void *get_dist_func_param() override{return this;}
};

}
// kind=0: two level-filtered L1 outputs, lower then upper, each n*k entries.
// kind=1: external squared-Euclidean query to all source rows, q*k entries.
// All outputs invalid on nonzero status:1 args,2 budget,4 nonfinite,6 index,7 alloc.
extern "C" int local_neighbors(const double *source,const double *queries,const uint32_t *levels,
 uint32_t n,uint32_t q,uint32_t d,uint32_t k,uint32_t b,int kind,uint64_t budget,
 int64_t *ids,double *distances,uint64_t capacity,LocalReport *report){
    if(!source || !queries || !ids || !distances || !report || n<1 || n>2000000 || q<1 || q>2000000 ||
       d<1 || d>64 || k<1 || k>128 || k>n || budget<1 || budget>500000000 ||
       (kind!=0 && kind!=1) || capacity<uint64_t(q)*k*(kind==0?2:1) || (kind==0 && (!levels || n!=q || b<2 || b>128)))return 1;
    *report={};
    try {
        for(size_t i=0;i<size_t(n)*d;++i)if(!std::isfinite(source[i]))return 4;
        for(size_t i=0;i<size_t(q)*d;++i)if(!std::isfinite(queries[i]))return 4;
        std::vector<uint32_t> rank(n);
        if(kind==0){
            std::vector<uint32_t> counts(b),offset(b);
            for(uint32_t i=0;i<n;++i){if(levels[i]>=b)return 1;++counts[levels[i]];}
            for(uint32_t i=0;i<b;++i){if(counts[i]<k)return 1;if(i)offset[i]=offset[i-1]+counts[i-1];}
            for(uint32_t i=0;i<n;++i)rank[i]=offset[levels[i]]++;
        }else for(uint32_t i=0;i<n;++i)rank[i]=i;
        auto query=[&](Space &space,hnswlib::HierarchicalNSW<double> &index,uint32_t row,size_t offset)->int {
            uint32_t queryID=n+row;uint64_t before=space.evaluations;
            auto found=index.searchKnn(&queryID,k);report->query+=space.evaluations-before;if(space.error)return space.error;
            std::vector<std::pair<double,uint32_t>> nearest;
            while(!found.empty()){auto h=found.top();found.pop();if(h.second>=n || !std::isfinite(h.first) || h.first<0)return 6;nearest.emplace_back(h.first,uint32_t(h.second));}
            if(nearest.size()!=k)return 6;
            std::sort(nearest.begin(),nearest.end(),[&](auto a,auto b){return a.first<b.first || (a.first==b.first && rank[a.second]<rank[b.second]);});
            for(uint32_t j=0;j<k;++j){ids[offset+j]=nearest[j].second;distances[offset+j]=nearest[j].first;}
            return 0;
        };
        if(kind==0){
            // Only eligible levels enter each serial prefix/suffix index. This
            // avoids traversing a global graph whose results are mostly filtered.
            std::vector<std::vector<uint32_t>> rows(b);
            for(uint32_t row=0;row<n;++row)rows[levels[row]].push_back(row);
            for(int direction=0;direction<2;++direction){
                Space space(source,queries,n,q,d,true,budget);
                hnswlib::HierarchicalNSW<double> index(&space,n,16,200,7);index.setEf(128);
                for(uint32_t position=0;position<b;++position){
                    uint32_t level=direction?b-1-position:position;
                    for(uint32_t row:rows[level]){
                        const size_t offset=(size_t(direction)*q+row)*k;
                        if(position==0){for(uint32_t j=0;j<k;++j){ids[offset+j]=-1;distances[offset+j]=0;}}
                        else {int code=query(space,index,row,offset);if(code)return code;}
                    }
                    if(position+1<b)for(uint32_t row:rows[level]){
                        uint64_t before=space.evaluations;index.addPoint(&row,row);report->construction+=space.evaluations-before;if(space.error)return space.error;
                    }
                }
                report->indexBytes=std::max(report->indexBytes,uint64_t(index.indexFileSize()));
            }
        }else{
            Space space(source,queries,n,q,d,false,budget);
            hnswlib::HierarchicalNSW<double> index(&space,n,16,200,7);index.setEf(128);
            for(uint32_t row=0;row<n;++row){uint64_t before=space.evaluations;index.addPoint(&row,row);report->construction+=space.evaluations-before;if(space.error)return space.error;}
            report->indexBytes=index.indexFileSize();
            for(uint32_t row=0;row<q;++row){int code=query(space,index,row,size_t(row)*k);if(code)return code;}
        }
        return 0;
    }catch(const std::bad_alloc&){return 7;}catch(...){return 6;}
}
