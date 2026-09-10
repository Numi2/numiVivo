import Foundation

public enum VivoOmicsNBMomentMethod: String, Codable, Sendable {
    case direct, adaptive
}

public struct VivoOmicsNBMomentSummation: Codable, Sendable, Equatable {
    public let method: VivoOmicsNBMomentMethod
    public let coveredCounts: UInt64
    public let acceptedBlocks: Int
    public let probabilityErrorBound: Double
    public let firstRawMomentErrorBound: Double
    public let secondRawMomentErrorBound: Double
    public let meanErrorBound: Double
    public let varianceErrorBound: Double
}

public struct VivoOmicsNBDevianceMoments: Codable, Sendable, Equatable {
    public let mean: Double
    public let variance: Double
    public let devianceScale: Double
    public let degreesOfFreedom: Double
    public let firstCount: UInt64
    public let lastCount: UInt64
    public let evaluatedCounts: Int
    public let omittedProbabilityBound: Double
    public let meanTruncationBound: Double
    public let varianceTruncationBound: Double
    /// Absent for the original direct evaluator. Together with truncation
    /// bounds, these bound analytic approximation, not floating rounding.
    public internal(set) var summation: VivoOmicsNBMomentSummation? = nil
}

public struct VivoOmicsNBAdjustedResiduals: Codable, Sendable, Equatable {
    public let moments: [VivoOmicsNBDevianceMoments]
    public let leverage: [Double]
    public let unitDeviance: [Double]
    public let unitDegreesOfFreedom: [Double]
    public let deviance: Double
    public let degreesOfFreedom: Double
    public let quasiDispersion: Double?
}

/// Direct conditional moments for QL development. Tail bounds concern omitted
/// support, not interval bounds for floating-point rounding. No asymptotic or
/// table fallback is used if the explicit support-work bound is exhausted.
public enum VivoOmicsNBResidualAdjustment {
    private struct Sum {
        var value = 0.0, compensation = 0.0
        mutating func add(_ x: Double) {
            let y = x-compensation, next = value+y
            compensation = (next-value)-y; value = next
        }
    }
    private struct Tail {
        var mass: Double, first: Double, second: Double
        static let zero = Tail(mass: 0,first: 0,second: 0)
        static let unbounded = Tail(mass: .infinity,first: .infinity,second: .infinity)
    }
    public static func moments(mean mu: Double, dispersion a: Double,
        relativeTolerance: Double = 1e-10, maximumTerms: Int = 1_000_000,
        method: VivoOmicsNBMomentMethod = .direct) throws -> VivoOmicsNBDevianceMoments {
        guard mu.isFinite, mu >= 0, a.isFinite, a == 0 || (1e-8...100).contains(a),
              relativeTolerance.isFinite, (1e-13...1e-3).contains(relativeTolerance),
              (1...10_000_000).contains(maximumTerms), mu <= 9_007_199_254_740_992-Double(maximumTerms) else {
            throw VivoOmicsStatisticsError.invalid("deviance moment domain or work bound")
        }
        if mu == 0 {
            return .init(mean: 0,variance: 0,devianceScale: 0,degreesOfFreedom: 0,
                firstCount: 0,lastCount: 0,evaluatedCounts: 1,omittedProbabilityBound: 0,
                meanTruncationBound: 0,varianceTruncationBound: 0)
        }
        if method == .adaptive, mu >= 10_000, a > 0 {
            return try VivoOmicsNBAdaptiveMoments.evaluate(mean: mu,dispersion: a,
                relativeTolerance: relativeTolerance,maximumTerms: maximumTerms)
        }
        let r = a == 0 ? 0 : 1/a, q = a == 0 ? 0 : mu/(mu+r)
        let mode = UInt64(floor(a == 0 ? mu : max(0,mu*(1-a))))
        var low = mode, high = mode, lowWeight = 1.0, highWeight = 1.0, terms = 0
        var mass = Sum(), first = Sum(), second = Sum()
        func upward(_ k: UInt64) -> Double {
            a == 0 ? mu/Double(k+1) : q*(Double(k)+r)/Double(k+1)
        }
        func downward(_ k: UInt64) -> Double {
            a == 0 ? Double(k)/mu : Double(k)/(q*(Double(k)+r-1))
        }
        func add(_ k: UInt64, _ weight: Double) throws {
            let d = try VivoOmicsNegativeBinomial.unitDeviance(count: k,mean: mu,dispersion: a)
            mass.add(weight); first.add(weight*d); second.add(weight*d*d); terms += 1
        }
        let zeroDeviance = try VivoOmicsNegativeBinomial.unitDeviance(count: 0,mean: mu,dispersion: a)
        func leftTail() -> Tail {
            if low == 0 || lowWeight == 0 { return .zero }
            let rho = downward(low)
            if rho >= 1 { return .unbounded }
            let bound = lowWeight*rho/(1-rho)
            return .init(mass: bound,first: bound*zeroDeviance,second: bound*zeroDeviance*zeroDeviance)
        }
        func rightTail() -> Tail {
            if highWeight == 0 { return .zero }
            let rho = max(q,upward(high))
            if rho >= 1 { return .unbounded }
            let t = 1/(1-rho), s0 = rho*t, s1 = s0*t, s2 = s1*(1+rho)*t
            let k = Double(high), b1: Double, b2: Double
            if a > 0 {
                // D(k+j) <= A*(k+j)+B for NB. Combine this polynomial
                // envelope with the geometric upper bound on future masses.
                let slope = 2*log1p(r/mu), intercept = zeroDeviance
                let b = slope*k+intercept
                b1 = b*s0+slope*s1
                b2 = b*b*s0+2*b*slope*s1+slope*slope*s2
            } else {
                // log(x) <= x bounds the Poisson deviance by 2*k^2/mu+2*mu.
                let s3 = s1*(1+4*rho+rho*rho)*t*t
                let s4 = s1*(1+11*rho+11*rho*rho+rho*rho*rho)*t*t*t
                let k2 = k*k*s0+2*k*s1+s2
                let k4 = k*k*k*k*s0+4*k*k*k*s1+6*k*k*s2+4*k*s3+s4
                b1 = 2*k2/mu+2*mu*s0
                b2 = 4*k4/(mu*mu)+8*k2+4*mu*mu*s0
            }
            return .init(mass: highWeight*s0,first: highWeight*b1,second: highWeight*b2)
        }
        try add(mode,1)
        var extendLeft = low > 0
        while true {
            if terms == 1 || terms % 32 < 2 || terms >= maximumTerms {
                let l = leftTail(), u = rightTail(), total = l.mass+u.mass
                let m = first.value/mass.value, rawSecond = second.value/mass.value
                let variance = rawSecond-m*m
                let meanError = (l.first+u.first+abs(m)*total)/mass.value
                let secondError = (l.second+u.second+abs(rawSecond)*total)/mass.value
                let varianceError = secondError+(2*abs(m)+meanError)*meanError
                if m > 0, variance > 0, total/mass.value <= relativeTolerance,
                   meanError <= relativeTolerance*m, varianceError <= relativeTolerance*variance {
                    let scale = 2*m/variance, df = scale*m
                    guard scale.isFinite, df.isFinite else { throw VivoOmicsStatisticsError.invalid("nonfinite deviance moment scale") }
                    return .init(mean: m,variance: variance,devianceScale: scale,degreesOfFreedom: df,
                        firstCount: low,lastCount: high,evaluatedCounts: terms,
                        omittedProbabilityBound: total/mass.value,meanTruncationBound: meanError,
                        varianceTruncationBound: varianceError)
                }
                extendLeft = low > 0 && (l.mass > relativeTolerance*mass.value/8 ||
                    l.first > relativeTolerance*first.value/8 || l.second > relativeTolerance*second.value/8)
            }
            guard terms < maximumTerms else {
                throw VivoOmicsStatisticsError.invalid("deviance moments exhausted \(maximumTerms) support evaluations (\(low)...\(high))")
            }
            if extendLeft, low > 0 {
                lowWeight *= downward(low); low -= 1; try add(low,lowWeight)
            }
            if terms < maximumTerms {
                highWeight *= upward(high); high += 1; try add(high,highWeight)
            }
        }
    }

    /// Conditional residual stage at supplied fitted means and average QL scale.
    /// No new mean, trend, average scale, prior or hypothesis test is estimated.
    public static func adjustedResiduals(counts: [UInt64], means: [Double], design: [[Double]],
        dispersion: Double, averageQuasiDispersion: Double, relativeTolerance: Double = 1e-10,
        maximumTerms: Int = 1_000_000, method: VivoOmicsNBMomentMethod = .direct) throws -> VivoOmicsNBAdjustedResiduals {
        guard counts.count == means.count, counts.count == design.count,
              counts.allSatisfy({ $0 <= 9_007_199_254_740_992 }),
              means.allSatisfy({ $0.isFinite && $0 > 0 }), dispersion.isFinite, dispersion >= 0,
              averageQuasiDispersion.isFinite, averageQuasiDispersion > 0 else {
            throw VivoOmicsStatisticsError.invalid("adjusted residual shapes, means or QL scale")
        }
        let fittedDispersion = dispersion/averageQuasiDispersion
        let scaledMeans = means.map { $0/averageQuasiDispersion }
        guard scaledMeans.allSatisfy({ $0.isFinite && $0 > 0 }),
              fittedDispersion.isFinite, fittedDispersion == 0 || (1e-8...100).contains(fittedDispersion),
              dispersion == 0 || fittedDispersion > 0 else {
            throw VivoOmicsStatisticsError.invalid("QL scaling exceeds the representable mean or dispersion domain")
        }
        let weights = means.map { sqrt($0/(1+fittedDispersion*$0)) }
        let weighted = zip(design,weights).map { row,w in row.map { $0*w } }
        let leverage = try VivoOmicsQR(design: weighted).leverage
        var ms: [VivoOmicsNBDevianceMoments] = [], ds: [Double] = [], dfs: [Double] = []
        var deviance = Sum(), degrees = Sum()
        for i in counts.indices {
            let m = try moments(mean: scaledMeans[i],dispersion: dispersion,
                relativeTolerance: relativeTolerance,maximumTerms: maximumTerms,method: method)
            let complement = 1-leverage[i]
            let d: Double, df: Double
            if complement < 1e-4 { d = 0; df = 0 }
            else {
                d = try VivoOmicsNegativeBinomial.unitDeviance(count: counts[i],mean: means[i],dispersion: fittedDispersion)*m.devianceScale
                df = complement*m.degreesOfFreedom
            }
            guard d.isFinite, df.isFinite, d >= 0, df >= 0 else { throw VivoOmicsStatisticsError.invalid("nonfinite adjusted residual contribution") }
            ms.append(m);ds.append(d);dfs.append(df);deviance.add(d);degrees.add(df)
        }
        let quasiDispersion = degrees.value >= 1e-4 ? deviance.value/degrees.value : nil
        guard deviance.value.isFinite, degrees.value.isFinite,
              quasiDispersion == nil || quasiDispersion!.isFinite else {
            throw VivoOmicsStatisticsError.invalid("adjusted residual aggregate overflow")
        }
        return .init(moments: ms,leverage: leverage,unitDeviance: ds,unitDegreesOfFreedom: dfs,
            deviance: deviance.value,degreesOfFreedom: degrees.value,
            quasiDispersion: quasiDispersion)
    }
}
