// Frozen from 36ccd858504ed643b21f2b622d81a3b261fc0950. Oracle only, never linked into production.
import Foundation
import NumiVivoKit

public struct VivoSingleCellEmbeddingResult: Codable, Sendable, Equatable {
    public let method: String
    public let options: VivoSingleCellEmbeddingOptions
    public let cells: [VivoOmicsCellIdentity]
    public let initialCoordinates: [[Double]]
    public let coordinates: [[Double]]
    public let curveA: Double
    public let curveB: Double
    public let curveSquaredError: Double
    public let retainedDirectedEdges: Int
    public let discardedDirectedEdges: Int
    public let edgeVisits: Int
    public let attractiveUpdates: Int
    public let negativeSamples: Int
    public let maximumCoordinateMagnitude: Double
    public let qualification: String
}

enum VivoLegacyEmbeddingReference {
    /// Fit UMAP's offset-exponential target using log-positive curve parameters.
    static func fitCurve(spread: Double,minimumDistance: Double) throws -> (a: Double,b: Double,error: Double) {
        let x = (0..<300).map { 3*spread*Double($0)/299 }
        let y = x.map { $0<minimumDistance ? 1 : exp(-($0-minimumDistance)/spread) }
        func evaluate(_ u: Double,_ v: Double) -> (error: Double,g0: Double,g1: Double,h00: Double,h01: Double,h11: Double) {
            let a = exp(u),b = exp(v)
            var error=0.0,g0=0.0,g1=0.0,h00=0.0,h01=0.0,h11=0.0
            for i in 1..<x.count {
                let q = 1/(1+a*pow(x[i],2*b)),r = q-y[i]
                let j0 = -q*(1-q),j1 = j0*2*b*log(x[i])
                error += r*r;g0 += j0*r;g1 += j1*r;h00 += j0*j0;h01 += j0*j1;h11 += j1*j1
            }
            return (error,g0,g1,h00,h01,h11)
        }
        var u=0.0,v=0.0,damping=1e-3
        for _ in 0..<500 {
            let e = evaluate(u,v)
            if hypot(e.g0,e.g1)<1e-8 { return (exp(u),exp(v),e.error) }
            let h00=e.h00+damping,h11=e.h11+damping,det=h00*h11-e.h01*e.h01
            guard det.isFinite,det>0 else { throw VivoOmicsError.invalid("UMAP curve fit system") }
            let du=(-h11*e.g0+e.h01*e.g1)/det,dv=(e.h01*e.g0-h00*e.g1)/det
            let nextU=max(-20,min(20,u+du)),nextV=max(-6,min(4,v+dv)),next=evaluate(nextU,nextV)
            if next.error<e.error { u=nextU;v=nextV;damping=max(1e-12,damping/3) }
            else { damping *= 10 }
            if !damping.isFinite { break }
        }
        throw VivoOmicsError.invalid("UMAP curve fit did not reach stationarity")
    }
    static func attractiveCoefficient(squaredDistance: Double,a: Double,b: Double) -> Double {
        squaredDistance>0 ? -2*a*b*pow(squaredDistance,b-1)/(1+a*pow(squaredDistance,b)) : 0
    }
    static func repulsiveCoefficient(squaredDistance: Double,a: Double,b: Double,strength: Double) -> Double {
        squaredDistance>0 ? 2*strength*b/((0.001+squaredDistance)*(1+a*pow(squaredDistance,b))) : 0
    }
    static func run(_ graph: VivoSingleCellNeighborGraph,scores: [[Double]],options: VivoSingleCellEmbeddingOptions) throws -> VivoSingleCellEmbeddingResult {
        try options.validate()
        let n=graph.cells.count,d=options.dimensions
        guard scores.count==n,scores.allSatisfy({ $0.count>=d }),graph.weights.max() != nil else {
            throw VivoOmicsError.invalid("embedding needs a nonempty graph and sufficient PCA dimensions")
        }
        let curve=try fitCurve(spread: options.spread,minimumDistance: options.minimumDistance)
        var state=options.seed
        func random() -> UInt64 {
            state &+= 0x9E3779B97F4A7C15;var z=state
            z=(z^(z>>30)) &* 0xBF58476D1CE4E5B9;z=(z^(z>>27)) &* 0x94D049BB133111EB
            return z^(z>>31)
        }
        var coordinates=[[Double]](repeating: Array(repeating: 0,count: d),count: n)
        for f in 0..<d {
            let column=scores.map { $0[f] },low=column.min()!,high=column.max()!
            guard high>low,low.isFinite,high.isFinite else { throw VivoOmicsError.invalid("embedding initialization has zero or nonfinite range") }
            for i in 0..<n { coordinates[i][f]=10*(scores[i][f]-low)/(high-low)+(Double(random()>>11)/9_007_199_254_740_992-0.5)*2e-4 }
        }
        let initial=coordinates,maximum=graph.weights.max()!
        var heads: [Int]=[],tails: [Int]=[],intervals: [Double]=[],bound=0.0
        for i in 0..<n { for edge in graph.rowOffsets[i]..<graph.rowOffsets[i+1] where graph.weights[edge]>=maximum/Double(options.epochs) {
            let interval=maximum/graph.weights[edge]
            heads.append(i);tails.append(graph.columnIndices[edge]);intervals.append(interval)
            bound += ceil(Double(options.epochs)/interval)*Double(1+options.negativeSampleRate)
        } }
        guard bound<=Double(options.maximumUpdates),Double(heads.count)*Double(options.epochs)<=Double(options.maximumUpdates) else { throw VivoOmicsError.limit("UMAP update budget before optimization") }
        let negativeIntervals=intervals.map { $0/Double(options.negativeSampleRate) }
        var nextPositive=intervals,nextNegative=negativeIntervals,attractive=0,negative=0,alpha=options.learningRate
        func distance(_ i: Int,_ j: Int) -> Double {
            var value=0.0;for f in 0..<d { let delta=coordinates[i][f]-coordinates[j][f];value += delta*delta };return value
        }
        for epoch in 0..<options.epochs {
            try Task.checkCancellation()
            for edge in heads.indices where nextPositive[edge]<=Double(epoch) {
                let i=heads[edge],j=tails[edge],coefficient=attractiveCoefficient(squaredDistance: distance(i,j),a: curve.a,b: curve.b)
                for f in 0..<d {
                    let gradient=max(-4,min(4,coefficient*(coordinates[i][f]-coordinates[j][f])))
                    coordinates[i][f] += alpha*gradient;coordinates[j][f] -= alpha*gradient
                }
                attractive += 1;nextPositive[edge] += intervals[edge]
                let samples=max(0,Int((Double(epoch)-nextNegative[edge])/negativeIntervals[edge]))
                for _ in 0..<samples {
                    let other=Int(random()%UInt64(n));negative += 1
                    if other==i { continue }
                    let coefficient=repulsiveCoefficient(squaredDistance: distance(i,other),a: curve.a,b: curve.b,strength: options.repulsionStrength)
                    for f in 0..<d { coordinates[i][f] += alpha*max(-4,min(4,coefficient*(coordinates[i][f]-coordinates[other][f]))) }
                }
                nextNegative[edge] += Double(samples)*negativeIntervals[edge]
                guard attractive+negative<=options.maximumUpdates else { throw VivoOmicsError.limit("UMAP update budget during optimization") }
            }
            guard coordinates.allSatisfy({ $0.allSatisfy(\.isFinite) }) else { throw VivoOmicsError.invalid("nonfinite embedding coordinates") }
            alpha=options.learningRate*(1-Double(epoch)/Double(options.epochs))
        }
        return .init(method: "sparse-UMAP-negative-sampling-Double-SplitMix64-v1",options: options,cells: graph.cells,
            initialCoordinates: initial,coordinates: coordinates,curveA: curve.a,curveB: curve.b,curveSquaredError: curve.error,
            retainedDirectedEdges: heads.count,discardedDirectedEdges: graph.weights.count-heads.count,edgeVisits: heads.count*options.epochs,attractiveUpdates: attractive,
            negativeSamples: negative,maximumCoordinateMagnitude: coordinates.flatMap { $0 }.map(abs).max() ?? 0,
            qualification: "Fixed-epoch descriptive UMAP-compatible coordinates; not optimizer convergence, clustering, donor integration or biological validation")
    }
}
