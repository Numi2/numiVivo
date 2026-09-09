#include "NumiVivoCore/NumiVivoOmicsHNSW.h"
#include "ThirdParty/hnswlib/hnswlib.h"
#include <algorithm>
#include <bit>
#include <cmath>
#include <cerrno>
#include <cstdint>
#include <cstring>
#include <fcntl.h>
#include <list>
#include <stdexcept>
#include <sys/stat.h>
#include <unistd.h>
#include <unordered_map>

namespace {
constexpr uint32_t tileRows = 256;
uint32_t u32(const unsigned char *p) {
    return uint32_t(p[0]) | uint32_t(p[1])<<8 | uint32_t(p[2])<<16 | uint32_t(p[3])<<24;
}
uint64_t u64(const unsigned char *p) { return uint64_t(u32(p)) | uint64_t(u32(p+4))<<32; }
struct File {
    int fd;
    explicit File(const char *path): fd(open(path, O_RDONLY|O_NOFOLLOW|O_NONBLOCK)) {
        if (fd<0) throw std::runtime_error("score file");
    }
    ~File() { close(fd); }
};
struct Tile { uint32_t first; std::vector<double> values; };
class Scores final : public hnswlib::SpaceInterface<double> {
    File file;
    const NVivoHNSWOptions &options;
    size_t capacity;
    std::list<Tile> tiles;
    std::unordered_map<uint32_t, std::list<Tile>::iterator> positions;
    uint64_t cachedBytes = 0;
 public:
    NVivoHNSWReport &report;
    uint64_t evaluations = 0;
    int32_t error = 0;
    Scores(const char *path, const NVivoHNSWOptions &o, NVivoHNSWReport &r): file(path), options(o),
        capacity(o.score_cache_bytes/(uint64_t(tileRows)*o.dimensions*8)), report(r) {
        struct stat s{};
        if (fstat(file.fd,&s)!=0 || !S_ISREG(s.st_mode) || uint64_t(s.st_size)!=uint64_t(o.rows)*o.dimensions*16)
            throw std::runtime_error("score dimensions");
    }
    const double *row(uint32_t id) {
        if (id>=options.rows) throw std::runtime_error("score row");
        const uint32_t first=id/tileRows*tileRows;
        auto found=positions.find(first);
        if (found!=positions.end()) {
            ++report.score_cache_hits;
            tiles.splice(tiles.begin(),tiles,found->second);
            return tiles.front().values.data()+size_t(id-first)*options.dimensions;
        }
        if (tiles.size()==capacity) {
            cachedBytes-=tiles.back().values.size()*8; positions.erase(tiles.back().first); tiles.pop_back();
        }
        const size_t count=size_t(std::min(tileRows,options.rows-first))*options.dimensions;
        std::vector<unsigned char> bytes(count*16);
        size_t done=0;
        while (done<bytes.size()) {
            ssize_t got=pread(file.fd,bytes.data()+done,bytes.size()-done,off_t(uint64_t(first)*options.dimensions*16+done));
            if (got<0 && errno==EINTR) continue;
            if (got<=0) throw std::runtime_error("score read");
            done+=size_t(got); ++report.score_read_calls; report.score_read_bytes+=uint64_t(got);
        }
        Tile tile{first,std::vector<double>(count)};
        for (size_t i=0;i<count;++i) {
            const unsigned char *p=bytes.data()+i*16;
            const double value=std::bit_cast<double>(u64(p+8));
            if (u32(p)!=first+i/options.dimensions || u32(p+4)!=i%options.dimensions || !std::isfinite(value))
                throw std::runtime_error("score record");
            tile.values[i]=value;
        }
        cachedBytes+=count*8; report.peak_cached_score_bytes=std::max(report.peak_cached_score_bytes,cachedBytes);
        tiles.push_front(std::move(tile)); positions.emplace(first,tiles.begin());
        return tiles.front().values.data()+size_t(id-first)*options.dimensions;
    }
    static double distance(const void *a, const void *b, const void *context) {
#pragma clang fp contract(off)
        auto &self=*static_cast<Scores *>(const_cast<void *>(context));
        if (self.error) return INFINITY;
        if (self.evaluations>=self.options.maximum_distance_evaluations) { self.error=2; return INFINITY; }
        ++self.evaluations;
        uint32_t i,j; std::memcpy(&i,a,4); std::memcpy(&j,b,4);
        // At least two cache tiles are admitted. Touching the first makes it MRU,
        // so loading the second cannot invalidate its pointer.
        try {
            const double *x=self.row(i), *y=self.row(j);
            double value=0;
            for (uint32_t c=0;c<self.options.dimensions;++c) { const double delta=x[c]-y[c]; value+=delta*delta; }
            if (!std::isfinite(value)) { self.error=4; return INFINITY; }
            return value;
        } catch (const std::bad_alloc &) { self.error=7; return INFINITY; }
          catch (...) { self.error=3; return INFINITY; }
        // Never throw through hnswlib's borrowed visited-list lifetime. The
        // caller checks error after every insertion/query; no failed result is used.
    }
    size_t get_data_size() override { return sizeof(uint32_t); }
    hnswlib::DISTFUNC<double> get_dist_func() override { return distance; }
    void *get_dist_func_param() override { return this; }
};
}

int32_t nvivo_omics_hnsw_neighbors(const char *path, const NVivoHNSWOptions *o,
    uint32_t *indices, double *distances, uint64_t capacity, NVivoHNSWReport *report,
    NVivoHNSWCancel cancel, void *context) {
    if (!path || !o || !indices || !distances || !report || o->struct_size!=sizeof(*o) || o->abi_version!=1 ||
        o->rows<2 || o->rows>1'000'000 || o->dimensions<1 || o->dimensions>64 || o->neighbors<2 || o->neighbors>128 ||
        o->neighbors>o->rows || o->connections<8 || o->connections>64 || o->ef_construction<o->connections || o->ef_construction>512 ||
        o->ef_search<o->neighbors || o->ef_search>1024 || o->maximum_distance_evaluations<1 || o->maximum_distance_evaluations>2'000'000'000 ||
        o->score_cache_bytes<uint64_t(2)*tileRows*o->dimensions*8 || o->score_cache_bytes>67'108'864 ||
        uint64_t(o->rows)*o->neighbors>4'000'000 || capacity<uint64_t(o->rows)*o->neighbors) return 1;
    *report={};
    try {
        std::unique_ptr<Scores> scores;
        try { scores=std::make_unique<Scores>(path,*o,*report); }
        catch (const std::bad_alloc &) { return 7; }
        catch (...) { return 3; }
        hnswlib::HierarchicalNSW<double> index(scores.get(),o->rows,o->connections,o->ef_construction,size_t(o->seed));
        for (uint32_t row=0;row<o->rows;++row) {
            if (cancel && cancel(context)) return 5;
            index.addPoint(&row,row);
            if (scores->error) return scores->error;
        }
        report->construction_distances=scores->evaluations;
        report->index_storage_bytes=index.indexFileSize();
        index.setEf(o->ef_search);
        for (uint32_t row=0;row<o->rows;++row) {
            if (cancel && cancel(context)) return 5;
            auto found=index.searchKnn(&row,o->neighbors);
            if (scores->error) return scores->error;
            std::vector<std::pair<double,uint32_t>> nearest;
            while (!found.empty()) {
                const auto p=found.top(); found.pop();
                if (p.second>=o->rows || !std::isfinite(p.first) || p.first<0) return 6;
                if (p.second!=row) nearest.emplace_back(p.first,uint32_t(p.second));
            }
            std::sort(nearest.begin(),nearest.end());
            if (nearest.size()<o->neighbors-1) return 6;
            const size_t offset=size_t(row)*o->neighbors;
            indices[offset]=row; distances[offset]=0;
            for (uint32_t k=1;k<o->neighbors;++k) {
                if (k>1 && nearest[k-1].second==nearest[k-2].second) return 6;
                indices[offset+k]=nearest[k-1].second; distances[offset+k]=std::sqrt(nearest[k-1].first);
            }
        }
        report->query_distances=scores->evaluations-report->construction_distances;
        return 0;
    } catch (const std::bad_alloc &) { return 7; }
      catch (...) { return 6; }
}
