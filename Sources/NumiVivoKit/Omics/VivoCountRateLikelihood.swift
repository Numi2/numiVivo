import Foundation

/// Exact sufficient grouping for RATE inference with fixed NB2 dispersion.
/// All cells at the same full RNA depth share a probability parameter. Their
/// counts may be summed here; dispersion inference would need more information.
public struct VivoCountDepthStratum: Codable, Sendable, Equatable {
    public let libraryCounts: [UInt64]
    public let cellsPerLibrary: [Int]
    public let geneCountsPerLibrary: [UInt64]
    public init(libraryCounts: [UInt64],cellsPerLibrary: [Int],geneCountsPerLibrary: [UInt64]) {
        self.libraryCounts=libraryCounts;self.cellsPerLibrary=cellsPerLibrary;self.geneCountsPerLibrary=geneCountsPerLibrary
    }
}

public struct VivoCountRateLikelihood: Sendable {
    public let cells: Int
    public let geneCounts: UInt64
    public let libraryCounts: UInt64
    public let cellDispersion: Double
    public let maximumLikelihoodRateCPM: Double
    private let exposures: [Double]
    private let counts: [Double]
    private let ns: [Double]
    private let weights: [Double]
    private let logOffsets: [Double]
    private let mode: Double
    private let modeScore: Double
    private let modeProbabilities: [Double]

    public init(_ stratum: VivoCountDepthStratum,cellDispersion phi: Double) throws {
        let m=stratum.libraryCounts.count
        guard m>0,m<=1_000_000,stratum.cellsPerLibrary.count==m,stratum.geneCountsPerLibrary.count==m,
              phi.isFinite,phi==0 || (1e-8...100).contains(phi) else { throw VivoOmicsError.invalid("count rate likelihood dimensions or dispersion") }
        var previous: UInt64=0,total: UInt64=0,library: UInt64=0,n=0
        for i in 0..<m {
            let depth=stratum.libraryCounts[i],cells=stratum.cellsPerLibrary[i],y=stratum.geneCountsPerLibrary[i]
            guard depth>previous,depth<=1_000_000_000,cells>0,cells<=1_000_000 else { throw VivoOmicsError.invalid("count depth bins require sorted unique positive libraries and cells") }
            let capacity=depth*UInt64(cells)
            guard y<=capacity else { throw VivoOmicsError.invalid("gene counts exceed full RNA counts in depth bin") }
            total=try vivoOmicsSum(total,y);library=try vivoOmicsSum(library,capacity);n+=cells;previous=depth
        }
        guard n<=1_000_000,library<=1_000_000_000_000 else { throw VivoOmicsError.limit("count rate likelihood observation budget") }
        // Spell out numeric conversion: Double.init on UInt64 can select
        // the bitPattern initializer when passed as a function value.
        let e=stratum.libraryCounts.map { Double($0)/1e6 },y=stratum.geneCountsPerLibrary.map { Double($0) },numbers=stratum.cellsPerLibrary.map { Double($0) }
        let offsets=phi==0 ? [] : e.map { log(phi*$0) },w=phi==0 ? [] : y.indices.map { y[$0]+numbers[$0]/phi }
        func score(_ x: Double) -> Double {
            if phi==0 { return Double(total)-Double(library)/1e6*exp(x) }
            var value=0.0
            for i in y.indices { let q=VivoCountObservation.sigmoid(offsets[i]+x);value+=y[i]*(1-q)-numbers[i]*q/phi }
            return value
        }
        var mle=0.0,x=0.0,s=0.0
        if total>0 {
            if phi==0 { mle=Double(total)/Double(library)*1e6;x=log(mle);s=score(x) }
            else {
                var lo=log(Double(total)/Double(library)*1e6),hi=lo
                for _ in 0..<256 { if score(lo)>=0 { break };lo-=2 }
                for _ in 0..<256 { if score(hi)<=0 { break };hi+=2 }
                guard score(lo)>=0,score(hi)<=0 else { throw VivoOmicsError.invalid("count rate likelihood mode bracket") }
                for _ in 0..<80 { let mid=(lo+hi)/2;if score(mid)>0 { lo=mid } else { hi=mid } }
                x=(lo+hi)/2;mle=exp(x);s=score(x)
            }
        }
        guard mle.isFinite,mle>=0,mle<=1_000_000*(1+1e-12) else { throw VivoOmicsError.invalid("count rate mode outside RNA rate domain") }
        cells=n;geneCounts=total;libraryCounts=library;cellDispersion=phi;maximumLikelihoodRateCPM=mle
        exposures=e;counts=y;ns=numbers;weights=w;logOffsets=offsets;mode=x;modeScore=s
        modeProbabilities=phi==0 ? [] : offsets.map { VivoCountObservation.sigmoid($0+x) }
    }

    /// Log likelihood relative to its own maximum, omitting rate-independent
    /// constants. Negative infinity means impossible zero rate after counts.
    public func logRelativeLikelihood(rateCPM r: Double) throws -> Double {
        guard r.isFinite,r>=0,r<=1_000_000*(1+1e-12) else { throw VivoOmicsError.invalid("count likelihood rate") }
        if r==0 { return geneCounts==0 ? 0 : -.infinity }
        let phi=cellDispersion
        if geneCounts==0 {
            if phi==0 { return -Double(libraryCounts)/1e6*r }
            var value=0.0
            for i in exposures.indices { value-=ns[i]/phi*log1p(phi*exposures[i]*r) }
            return value
        }
        let d=log(r)-mode
        var value=modeScore*d
        if phi==0 { value-=Double(libraryCounts)/1e6*maximumLikelihoodRateCPM*VivoCountObservation.expm1Remainder(d) }
        else {
            for i in exposures.indices { value-=weights[i]*VivoCountObservation.softplusRemainder(logOffsets[i]+mode,modeProbabilities[i],d) }
        }
        guard value.isFinite else { throw VivoOmicsError.invalid("count rate likelihood arithmetic") }
        return value
    }
}
