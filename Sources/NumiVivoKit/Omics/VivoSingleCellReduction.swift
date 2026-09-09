import Foundation

public struct VivoSingleCellReductionOptions: Codable, Sendable, Equatable {
    public var highlyVariableFeatures: Int = 2_000
    public var meanBins: Int = 20
    public var components: Int = 20
    public var maximumBasis: Int = 128
    public var relativeResidualTolerance: Double = 1e-6
    public var seed: UInt64 = 7
    public var retainProjectionCenters: Bool? = nil
    public init() {}
    private enum CodingKeys: String,CodingKey {
        case highlyVariableFeatures,meanBins,components,maximumBasis,relativeResidualTolerance,seed,retainProjectionCenters
    }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder,allowed: ["highlyVariableFeatures","meanBins","components","maximumBasis","relativeResidualTolerance","seed","retainProjectionCenters"])
        let c=try decoder.container(keyedBy: CodingKeys.self)
        highlyVariableFeatures=try c.decodeIfPresent(Int.self,forKey: .highlyVariableFeatures) ?? 2_000
        meanBins=try c.decodeIfPresent(Int.self,forKey: .meanBins) ?? 20
        components=try c.decodeIfPresent(Int.self,forKey: .components) ?? 20
        maximumBasis=try c.decodeIfPresent(Int.self,forKey: .maximumBasis) ?? 128
        relativeResidualTolerance=try c.decodeIfPresent(Double.self,forKey: .relativeResidualTolerance) ?? 1e-6
        seed=try c.decodeIfPresent(UInt64.self,forKey: .seed) ?? 7
        retainProjectionCenters=try c.decodeIfPresent(Bool.self,forKey: .retainProjectionCenters)
    }
    public func validate() throws {
        guard (1...10_000).contains(highlyVariableFeatures),(1...128).contains(meanBins),
              (1...64).contains(components),(components...256).contains(maximumBasis),
              relativeResidualTolerance.isFinite,(1e-10...1e-3).contains(relativeResidualTolerance) else {
            throw VivoOmicsError.invalid("sparse reduction options")
        }
    }
}
public struct VivoSingleCellVariableFeature: Codable, Sendable, Equatable {
    public let featureIndex: Int
    public let featureID: String
    public let meanNormalized: Double
    public let varianceNormalized: Double
    public let logMeanForBinning: Double
    public let meanBin: Int
    public let logDispersion: Double?
    public let normalizedDispersion: Double?
    public let selected: Bool
}
public struct VivoSingleCellReductionResult: Codable, Sendable, Equatable {
    public let method: String
    public let options: VivoSingleCellReductionOptions
    public let features: [VivoSingleCellVariableFeature]
    public let selectedFeatureIndices: [Int]
    public let cells: [VivoOmicsCellIdentity]
    /// Cell-major scores and selected-feature-major loadings, in source order.
    public let scores: [[Double]]
    public let loadings: [[Double]]
    public let explainedVariance: [Double]
    public let explainedVarianceRatio: [Double]
    public let relativeResiduals: [Double]
    public let maximumLoadingOrthogonalityError: Double
    public let basisSize: Int
    public let qualification: String
    public var projectionCenters: [Double]? = nil
}

/// Called only with the freshly validated common processing result. Neither
/// feature selection nor centering changes the source normalization denominator.
enum VivoSingleCellReduction {
    static func dot(_ a: [Double],_ b: [Double]) -> Double {
        var sum=0.0;for i in a.indices { sum+=a[i]*b[i] };return sum
    }
    static func norm(_ a: [Double]) -> Double { sqrt(dot(a,a)) }
    static func run(_ processed: VivoSingleCellProcessed,options: VivoSingleCellReductionOptions) throws -> VivoSingleCellReductionResult {
        try options.validate()
        let data=processed.dataset,normal=processed.normalized,n=data.cells.count,m=data.features.count
        guard n>options.components else { throw VivoOmicsError.invalid("PCA needs more cells than requested components") }
        var seen=[Int](repeating: 0,count: m),means=[Double](repeating: 0,count: m),m2=means
        for k in normal.values.indices {
            if k%65_536==0 { try Task.checkCancellation() }
            let j=normal.featureIndices[k],value=expm1(normal.values[k])
            seen[j]+=1;let delta=value-means[j];means[j]+=delta/Double(seen[j]);m2[j]+=delta*(value-means[j])
        }
        let selection = try selectFeatures(featureIDs: data.features.map(\.id), cells: n, seen: seen, nonzeroMeans: means, m2: m2, options: options)
        let selected = selection.selected, statistics = selection.statistics
        var local=[Int](repeating: -1,count: m)
        for (j,source) in selected.enumerated() { local[source]=j }
        var rowOffsets=[0],columns: [Int]=[],values: [Double]=[]
        for row in 0..<n {
            for k in normal.rowOffsets[row]..<normal.rowOffsets[row+1] {
                let j=local[normal.featureIndices[k]]
                if j>=0 { columns.append(j);values.append(normal.values[k]) }
            }
            rowOffsets.append(values.count)
        }
        let centers=selected.map { processed.features[$0].meanLogNormalized! }
        let totalVariance=selected.reduce(0) { $0+processed.features[$1].varianceLogNormalized! }
        guard totalVariance.isFinite,totalVariance>0 else { throw VivoOmicsError.invalid("PCA has no positive variance") }
        return try fit(cells: data.cells.map { .init(sampleID: $0.sampleID,barcode: $0.barcode) }, statistics: statistics,
            selected: selected, centers: centers, totalVariance: totalVariance, options: options,
            project: { v,shift in
                var result=[Double](repeating: -shift,count: n)
                for row in 0..<n { for k in rowOffsets[row]..<rowOffsets[row+1] { result[row]+=values[k]*v[columns[k]] } }
                return result
            }, transpose: { projected,initial in
                var result=initial
                for row in 0..<n { for k in rowOffsets[row]..<rowOffsets[row+1] { result[columns[k]]+=values[k]*projected[row] } }
                return result
            })
    }

    static func selectFeatures(featureIDs: [String], cells n: Int, seen: [Int], nonzeroMeans: [Double], m2: [Double],
                               options: VivoSingleCellReductionOptions) throws -> (statistics: [VivoSingleCellVariableFeature], selected: [Int]) {
        try options.validate()
        let m=featureIDs.count
        guard n>options.components, m>0, seen.count==m, nonzeroMeans.count==m, m2.count==m,
              seen.allSatisfy({ $0>=0 && $0<=n }) else { throw VivoOmicsError.invalid("HVG moments or axes") }
        var means=nonzeroMeans
        var variance=means,logMean=means,dispersion=[Double?](repeating: nil,count: m)
        for j in 0..<m {
            variance[j]=(m2[j]+means[j]*means[j]*Double(seen[j])*Double(n-seen[j])/Double(n))/Double(n-1)
            means[j]*=Double(seen[j])/Double(n)
            guard means[j].isFinite,variance[j].isFinite,variance[j]>=0 else { throw VivoOmicsError.invalid("nonfinite feature moments") }
            logMean[j]=log1p(means[j]==0 ? 1e-12:means[j])
            if variance[j]>0,means[j]>0 { dispersion[j]=log(variance[j]/means[j]) }
        }
        let minimum=logMean.min()!,maximum=logMean.max()!
        var edges: [Double]
        if minimum==maximum {
            let delta=minimum==0 ? 0.001:abs(minimum)*0.001
            edges=(0...options.meanBins).map { minimum-delta+2*delta*Double($0)/Double(options.meanBins) }
        } else {
            edges=(0...options.meanBins).map { minimum+(maximum-minimum)*Double($0)/Double(options.meanBins) }
            edges[0]-=(maximum-minimum)*0.001
        }
        var bins=[Int](repeating: 0,count: m),members=[[Double]](repeating: [],count: options.meanBins)
        for j in 0..<m {
            var bin=0;while bin+1<options.meanBins && logMean[j]>edges[bin+1] { bin+=1 }
            bins[j]=bin;if let d=dispersion[j] { members[bin].append(d) }
        }
        var averages=[Double](repeating: 0,count: options.meanBins),deviations=averages
        for b in members.indices where !members[b].isEmpty {
            let values=members[b],average=values.reduce(0,+)/Double(values.count)
            if values.count==1 { deviations[b]=average }
            else { averages[b]=average;deviations[b]=sqrt(values.reduce(0) { $0+pow($1-average,2) }/Double(values.count-1)) }
        }
        var standardized=[Double?](repeating: nil,count: m)
        for j in 0..<m {
            if let d=dispersion[j],deviations[bins[j]] != 0 {
                let value=(d-averages[bins[j]])/deviations[bins[j]]
                if value.isFinite { standardized[j]=value }
            }
        }
        let ranking=standardized.compactMap { $0 }.sorted(by: >)
        guard !ranking.isEmpty else { throw VivoOmicsError.invalid("no features with defined normalized dispersion") }
        let cutoff=ranking[min(options.highlyVariableFeatures,ranking.count)-1]
        let selected=(0..<m).filter { standardized[$0].map { $0>=cutoff } ?? false }
        guard selected.count>=options.components,selected.count<=10_000 else { throw VivoOmicsError.limit("PCA selected feature count, including cutoff ties") }
        let statistics=(0..<m).map { j in VivoSingleCellVariableFeature(featureIndex: j,featureID: featureIDs[j],
            meanNormalized: means[j],varianceNormalized: variance[j],logMeanForBinning: logMean[j],meanBin: bins[j],
            logDispersion: dispersion[j],normalizedDispersion: standardized[j],selected: standardized[j].map { $0>=cutoff } ?? false) }
        return (statistics,selected)
    }

    static func fit(cells: [VivoOmicsCellIdentity], statistics: [VivoSingleCellVariableFeature], selected: [Int],
                    centers: [Double], totalVariance: Double, options: VivoSingleCellReductionOptions,
                    project: ([Double],Double) throws -> [Double],
                    transpose: ([Double],[Double]) throws -> [Double]) throws -> VivoSingleCellReductionResult {
        try options.validate()
        let n=cells.count,width=selected.count
        guard n>options.components, width>=options.components, centers.count==width, centers.allSatisfy(\.isFinite),
              totalVariance.isFinite,totalVariance>0 else { throw VivoOmicsError.invalid("PCA moments or axes") }
        func scores(_ v: [Double]) throws -> [Double] { try project(v,dot(centers,v)) }
        func covariance(_ v: [Double]) throws -> [Double] {
            let projected=try scores(v),sum=projected.reduce(0,+)
            var result=try transpose(projected,centers.map { -$0*sum })
            for j in result.indices { result[j]/=Double(n-1) };return result
        }
        let capacity=min(options.maximumBasis,width)
        var state=options.seed
        func randomVector() -> [Double] {
            (0..<width).map { _ in
                state &+= 0x9E3779B97F4A7C15;var z=state
                z=(z^(z>>30)) &* 0xBF58476D1CE4E5B9;z=(z^(z>>27)) &* 0x94D049BB133111EB
                z ^= z>>31;return Double(z>>11)/9_007_199_254_740_992-0.5
            }
        }
        var basis: [[Double]]=[],q=randomVector(),h=[Double](repeating: 0,count: capacity*capacity)
        let initial=norm(q);q=q.map { $0/initial }
        for i in 0..<capacity {
            try Task.checkCancellation();basis.append(q)
            var z=try covariance(q);let scale=max(Double.leastNormalMagnitude,norm(z))
            var coefficients=[Double](repeating: 0,count: basis.count)
            for _ in 0..<2 { for j in basis.indices {
                let coefficient=dot(basis[j],z);coefficients[j]+=coefficient
                for f in z.indices { z[f]-=coefficient*basis[j][f] }
            } }
            for j in coefficients.indices { h[i*capacity+j]=coefficients[j];h[j*capacity+i]=coefficients[j] }
            if i+1==capacity { break }
            var length=norm(z)
            if length<=scale*1e-12 {
                // Restart in the orthogonal complement, including repeated eigenvalues.
                z=randomVector()
                for _ in 0..<2 { for vector in basis {
                    let coefficient=dot(vector,z);for f in z.indices { z[f]-=coefficient*vector[f] }
                } }
                length=norm(z)
            }
            guard length.isFinite,length>Double.leastNormalMagnitude else { throw VivoOmicsError.invalid("PCA basis breakdown") }
            q=z.map { $0/length }
        }
        let eigen=try symmetricEigen(h,n: capacity)
        var embedding=[[Double]](repeating: Array(repeating: 0,count: options.components),count: n)
        var loadings=[[Double]](repeating: Array(repeating: 0,count: options.components),count: width)
        var variances: [Double]=[],residuals: [Double]=[],vectors: [[Double]]=[]
        for component in 0..<options.components {
            try Task.checkCancellation()
            let mode=eigen.order[component],lambda=eigen.values[mode]
            guard lambda>totalVariance*1e-12 else { throw VivoOmicsError.invalid("requested PCA components exceed numerical rank") }
            var v=[Double](repeating: 0,count: width)
            for j in basis.indices { for f in v.indices { v[f]+=basis[j][f]*eigen.vectors[j*capacity+mode] } }
            let length=norm(v);v=v.map { $0/length }
            let largest=v.indices.max { abs(v[$0])<abs(v[$1]) }!
            if v[largest]<0 { v=v.map { -$0 } }
            let product=try covariance(v),residual=norm(v.indices.map { product[$0]-lambda*v[$0] })/lambda
            guard residual.isFinite,residual<=options.relativeResidualTolerance else {
                throw VivoOmicsError.invalid("PCA component \(component) residual \(residual) exceeds tolerance; increase maximumBasis")
            }
            let projected=try scores(v)
            for row in 0..<n { embedding[row][component]=projected[row] }
            for f in v.indices { loadings[f][component]=v[f] }
            variances.append(lambda);residuals.append(residual);vectors.append(v)
        }
        var orthogonality=0.0
        for i in vectors.indices { for j in vectors.indices { orthogonality=max(orthogonality,abs(dot(vectors[i],vectors[j])-(i==j ? 1:0))) } }
        guard orthogonality<1e-8 else { throw VivoOmicsError.invalid("PCA loadings lost orthogonality") }
        return .init(method: "sparse-seurat-dispersion-centered-krylov-PCA-v1",options: options,features: statistics,
            selectedFeatureIndices: selected,cells: cells,scores: embedding,
            loadings: loadings,explainedVariance: variances,explainedVarianceRatio: variances.map { $0/totalVariance },
            relativeResiduals: residuals,maximumLoadingOrthogonalityError: orthogonality,basisSize: capacity,
            qualification: "Descriptive unscaled log-normalized PCA; residual-qualified components, not donor integration or biological validation",
            projectionCenters: options.retainProjectionCenters == true ? centers : nil)
    }
    /// Only the bounded Krylov projection is dense, never cells by genes or the full covariance.
    static func symmetricEigen(_ input: [Double],n: Int) throws -> (values: [Double],vectors: [Double],order: [Int]) {
        var a=input,v=[Double](repeating: 0,count: n*n)
        for i in 0..<n { v[i*n+i]=1 }
        let scale=max(Double.leastNormalMagnitude,input.map(abs).max() ?? 0)
        var converged=false
        for _ in 0..<100 {
            var largest=0.0
            for p in 0..<n { for q in (p+1)..<n {
                let off=a[p*n+q];largest=max(largest,abs(off));if abs(off)<=scale*1e-13 { continue }
                let tau=(a[q*n+q]-a[p*n+p])/(2*off),t=(tau>=0 ? 1.0:-1.0)/(abs(tau)+hypot(1,tau))
                let c=1/sqrt(1+t*t),s=t*c
                a[p*n+p]-=t*off;a[q*n+q]+=t*off;a[p*n+q]=0;a[q*n+p]=0
                for k in 0..<n where k != p && k != q {
                    let x=a[k*n+p],y=a[k*n+q];a[k*n+p]=c*x-s*y;a[p*n+k]=a[k*n+p];a[k*n+q]=s*x+c*y;a[q*n+k]=a[k*n+q]
                }
                for k in 0..<n { let x=v[k*n+p],y=v[k*n+q];v[k*n+p]=c*x-s*y;v[k*n+q]=s*x+c*y }
            } }
            if largest<=scale*1e-13 { converged=true;break }
        }
        guard converged else { throw VivoOmicsError.invalid("PCA projected eigensolver did not converge") }
        let values=(0..<n).map { a[$0*n+$0] }
        return (values,v,(0..<n).sorted { values[$0]>values[$1] })
    }
}
