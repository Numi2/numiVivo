#include "NumiVivoCore/NumiVivoOmicsGaussian.h"
#define ACCELERATE_NEW_LAPACK
#include <Accelerate/Accelerate.h>
#include <algorithm>
#include <cmath>
#include <cstddef>
#pragma clang fp contract(off) reassociate(off)

extern "C" int32_t nvivo_omics_gaussian_bias(
    const double *source, const double *bias, uint32_t anchors,
    const double *queries, uint32_t query_rows, uint32_t dimensions, double sigma,
    uint64_t source_capacity, uint64_t query_capacity, uint64_t maximum_scalar_terms,
    double *delta, uint64_t delta_capacity, double *totals, uint64_t totals_capacity,
    double *workspace, uint64_t workspace_capacity, uint64_t *evaluated_scalar_terms,
    int32_t (*cancel)(void *), void *context) {
    if (!evaluated_scalar_terms) return 1;
    *evaluated_scalar_terms = 0;
    constexpr uint32_t tile = 256;
    const uint64_t ad = uint64_t(anchors)*dimensions, qd = uint64_t(query_rows)*dimensions;
    if (!source || !bias || !queries || !delta || !totals || !workspace ||
        !anchors || anchors > 200000000 || !query_rows || query_rows > 32 ||
        !dimensions || dimensions > 64 || !std::isfinite(sigma) || sigma < .001 || sigma > 1000 ||
        source_capacity < ad || query_capacity < qd || delta_capacity < qd || totals_capacity < query_rows ||
        workspace_capacity < uint64_t(tile)*(dimensions+query_rows)) return 1;
    if (ad*query_rows > maximum_scalar_terms) return 2;
    if (cancel && cancel(context)) return 5;
    for (uint64_t i=0;i<qd;++i) if (!std::isfinite(queries[i])) return 4;
    std::fill_n(delta,qd,0.); std::fill_n(totals,query_rows,0.);
    double *packed = workspace, *weights = workspace + tile*dimensions;
    for (uint32_t start=0;start<anchors;start+=tile) {
        if (cancel && cancel(context)) return 5;
        const uint32_t count = std::min(tile,anchors-start);
        for (uint32_t a=0;a<count;++a) for (uint32_t j=0;j<dimensions;++j) {
            const auto pos = uint64_t(start+a)*dimensions+j;
            if (!std::isfinite(source[pos]) || !std::isfinite(bias[pos])) return 4;
            packed[j*tile+a] = source[pos];
        }
        std::fill_n(weights,query_rows*count,0.);
        for (uint32_t q=0;q<query_rows;++q) {
            double *dist = weights+q*count;
            for (uint32_t j=0;j<dimensions;++j) {
                const double v = queries[q*dimensions+j];
                const double *column = packed+j*tile;
                // Vectorize independent anchors; never reassociate the dimension sum.
                #pragma clang loop vectorize(enable)
                for (uint32_t a=0;a<count;++a) {
                    const double difference = v-column[a];
                    dist[a] += difference*difference;
                }
            }
            for (uint32_t a=0;a<count;++a) dist[a] *= -.5*sigma;
        }
        const int length = int(query_rows*count);
        vvexp(weights,weights,&length);
        for (uint32_t q=0;q<query_rows;++q)
            for (uint32_t a=0;a<count;++a) totals[q] += weights[q*count+a];
        cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,int(query_rows),int(dimensions),int(count),
            1.,weights,int(count),bias+uint64_t(start)*dimensions,int(dimensions),1.,delta,int(dimensions));
    }
    *evaluated_scalar_terms = ad*query_rows;
    if (cancel && cancel(context)) return 5;
    for (uint32_t q=0;q<query_rows;++q) {
        // BLAS accumulation of subnormal products can amplify rounding after
        // normalization. Preserve the original scalar underflow semantics here.
        if (totals[q] < 1e-280) {
            if (ad > maximum_scalar_terms-*evaluated_scalar_terms) return 2;
            *evaluated_scalar_terms += ad;
            double *row = delta+q*dimensions;
            std::fill_n(row,dimensions,0.); totals[q]=0.;
            for (uint32_t a=0;a<anchors;++a) {
                if (a%tile==0 && cancel && cancel(context)) return 5;
                double squared=0.;
                for (uint32_t j=0;j<dimensions;++j) {
                    const double difference=queries[q*dimensions+j]-source[uint64_t(a)*dimensions+j];
                    squared += difference*difference;
                }
                const double weight=std::exp(-.5*sigma*squared);totals[q]+=weight;
                for (uint32_t j=0;j<dimensions;++j) row[j]+=weight*bias[uint64_t(a)*dimensions+j];
            }
        }
        if (!std::isfinite(totals[q])) return 4;
        for (uint32_t j=0;j<dimensions;++j) {
            double &v = delta[q*dimensions+j];
            if (!std::isfinite(v)) return 4;
            v = totals[q] > 0 ? v/totals[q] : 0.;
            if (!std::isfinite(v)) return 4;
        }
    }
    return 0;
}

// Legacy scalar Gaussian sum. FP contraction/reassociation remain disabled.
static int32_t scalarSum(
    const double *source, const double *bias, uint32_t anchors, const double *query,
    uint32_t dimensions, double sigma, uint64_t source_capacity, uint64_t query_capacity,
    double *numerator, uint64_t output_capacity, double *total,
    int32_t (*cancel)(void *), void *context, bool validate) {
    const uint64_t ad=uint64_t(anchors)*dimensions;
    if (!source || !bias || !query || !numerator || !total || !anchors || anchors>200000000 ||
        !dimensions || dimensions>64 || !std::isfinite(sigma) || sigma<.001 || sigma>1000 ||
        source_capacity<ad || query_capacity<dimensions || output_capacity<dimensions) return 1;
    if (cancel && cancel(context)) return 5;
    for (uint32_t j=0;j<dimensions;++j) if (!std::isfinite(query[j])) return 4;
    std::fill_n(numerator,dimensions,0.);*total=0.;
    for (uint32_t a=0;a<anchors;++a) {
        if (a%256==0 && cancel && cancel(context)) return 5;
        double squared=0.;
        for (uint32_t j=0;j<dimensions;++j) {
            const double value=source[uint64_t(a)*dimensions+j];
            if (validate && !std::isfinite(value)) return 4;
            const double difference=query[j]-value;squared+=difference*difference;
        }
        const double weight=std::exp(-.5*sigma*squared);*total+=weight;
        for (uint32_t j=0;j<dimensions;++j) {
            const double value=bias[uint64_t(a)*dimensions+j];
            if (validate && !std::isfinite(value)) return 4;
            numerator[j]+=weight*value;
        }
    }
    if (!std::isfinite(*total)) return 4;
    for (uint32_t j=0;j<dimensions;++j) if (!std::isfinite(numerator[j])) return 4;
    return 0;
}

extern "C" int32_t nvivo_omics_gaussian_scalar(
    const double *source, const double *bias, uint32_t anchors, const double *query,
    uint32_t dimensions, double sigma, uint64_t source_capacity, uint64_t query_capacity,
    double *numerator, uint64_t output_capacity, double *total,
    int32_t (*cancel)(void *), void *context) {
    return scalarSum(source,bias,anchors,query,dimensions,sigma,source_capacity,query_capacity,
        numerator,output_capacity,total,cancel,context,true);
}
extern "C" int32_t nvivo_omics_gaussian_scalar_apply(
    const double *source, const double *bias, uint32_t anchors, uint32_t dimensions,
    double sigma, uint64_t source_capacity, double *scores, uint32_t rows, uint64_t score_capacity,
    const intptr_t *selected, uint32_t selected_rows, uint64_t maximum_scalar_terms,
    NVivoGaussianScalarReport *report, int32_t (*cancel)(void *), void *context) {
    if (!source || !bias || !scores || !selected || !report || !anchors || anchors>200000000 ||
        !dimensions || dimensions>64 || !rows || rows>2000000 || !selected_rows || selected_rows>rows ||
        !std::isfinite(sigma) || sigma<.001 || sigma>1000 ||
        source_capacity<uint64_t(anchors)*dimensions || score_capacity<uint64_t(rows)*dimensions) return 1;
    if (uint64_t(anchors)*dimensions*selected_rows>maximum_scalar_terms) return 2;
    if (cancel && cancel(context)) return 5;
    // Validate immutable snapshots once per assembly phase, not once per query.
    for (uint64_t i=0;i<uint64_t(anchors)*dimensions;++i)
        if (!std::isfinite(source[i]) || !std::isfinite(bias[i])) return 4;
    for (uint32_t q=0;q<selected_rows;++q) if (selected[q]<0 || uint64_t(selected[q])>=rows) return 1;
    *report={0,INFINITY,0.,0.};
    double delta[64];
    for (uint32_t q=0;q<selected_rows;++q) {
        double *row=scores+uint64_t(selected[q])*dimensions,total=0.;
        const int32_t status=scalarSum(source,bias,anchors,row,dimensions,sigma,source_capacity,dimensions,
            delta,64,&total,cancel,context,false);
        if (status) return status;
        report->minimum_weight=std::min(report->minimum_weight,total);
        report->maximum_weight=std::max(report->maximum_weight,total);
        if (total==0.) ++report->zero_weight_rows;
        for (uint32_t j=0;j<dimensions;++j) {
            const double value=total>0. ? delta[j]/total : 0.;
            row[j]+=value;report->sum_squared+=value*value;
            if (!std::isfinite(row[j])) return 4;
        }
    }
    return 0;
}
