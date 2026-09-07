import Foundation

/// Diagnostics for equal-length independent Markov chains. These are classical
/// split-R-hat and Geyer initial-positive-sequence estimates, not rank-normalized
/// tail diagnostics or proof that initialization bias has disappeared.
public struct VivoCorrelatedSamplingDiagnostics: Codable, Sendable, Equatable {
    public let chainCount: Int
    public let retainedSamplesPerChain: Int
    public let mean: Double
    public let sampleVariance: Double
    public let autocorrelationEffectiveSampleSize: Double
    public let monteCarloStandardError: Double
    public let splitRHat: Double
    public let autocorrelationSequenceTruncated: Bool
    public let interpretation: String
}

public enum VivoCorrelatedSamplingAnalysis {
    static let maximumAutocorrelationProducts = 10_000_000
    static let maximumRetainedSamples = 10_000_000
    public static let interpretation = "Equal-length independent chains analyzed with classical split-R-hat and a conservative Geyer initial-positive-sequence autocorrelation ESS. A bounded autocorrelation calculation that does not reach a nonpositive pair is marked truncated and cannot pass the sampling gate. Diagnostics apply only to the supplied scalar observable after the declared discard. They are not rank-normalized tail ESS, evidence of independent retained frames, or proof of equilibrium."

    public static func calculate(chains: [[Double]]) throws -> VivoCorrelatedSamplingDiagnostics {
        try calculate(chains:chains,maximumAutocorrelationProducts:maximumAutocorrelationProducts)
    }

    static func calculate(chains: [[Double]], maximumAutocorrelationProducts: Int) throws -> VivoCorrelatedSamplingDiagnostics {
        guard (2...128).contains(chains.count), let count = chains.first?.count,
              (32...1_000_000).contains(count), count.isMultiple(of: 2),
              chains.count*count <= maximumRetainedSamples,
              chains.allSatisfy({ $0.count == count && $0.allSatisfy(\.isFinite) }),
              maximumAutocorrelationProducts > 0 else {
            throw VivoChemistryError.invalid("correlated-series chain count, equal even length or finite samples")
        }
        let splitLength = count/2
        let splitChains = chains.flatMap { chain in
            [Array(chain[..<splitLength]),Array(chain[splitLength...])]
        }
        let splitCount = splitChains.count
        func mean(_ values: [Double]) -> Double {
            let anchor = values[0]
            return anchor+values.reduce(0) { $0+($1-anchor) }/Double(values.count)
        }
        func variance(_ values: [Double], around center: Double) -> Double {
            values.reduce(0) { $0+($1-center)*($1-center) }/Double(values.count-1)
        }
        let splitMeans = splitChains.map(mean)
        let splitVariances = zip(splitChains,splitMeans).map { entry in
            variance(entry.0,around:entry.1)
        }
        let within = splitVariances.reduce(0,+)/Double(splitCount)
        let meanOfMeans = mean(splitMeans)
        let between = Double(splitLength)*splitMeans.reduce(0) { $0+pow($1-meanOfMeans,2) }/Double(splitCount-1)
        let variancePlus = (Double(splitLength-1)/Double(splitLength))*within+between/Double(splitLength)
        let flattened = chains.flatMap { $0 }, overallMean = mean(flattened)
        let overallVariance = variance(flattened,around:overallMean)

        let total = Double(splitCount*splitLength)
        let rHat: Double
        let effective: Double
        var truncated = false
        if variancePlus == 0 {
            rHat = 1; effective = total
        } else if within == 0 {
            // Constant but mutually displaced chains are maximally non-mixed.
            rHat = Double.greatestFiniteMagnitude; effective = Double(splitCount)
        } else {
            rHat = max(1,sqrt(max(0,variancePlus/within)))
            var pairSums: [Double] = [], lag = 1, closed = false
            let maximumPairs=maximumAutocorrelationProducts/max(1,2*splitCount*splitLength)
            while lag+1 < splitLength && pairSums.count < maximumPairs {
                func rho(_ offset: Int) -> Double {
                    var covariance = 0.0
                    for (chain,center) in zip(splitChains,splitMeans) {
                        for i in 0..<(splitLength-offset) {
                            covariance += (chain[i]-center)*(chain[i+offset]-center)
                        }
                    }
                    covariance /= Double(splitCount*splitLength)
                    return min(1,1-(within-covariance)/variancePlus)
                }
                let pair = rho(lag)+rho(lag+1)
                guard pair.isFinite else { throw VivoChemistryError.convergence("correlated-series autocovariance overflow") }
                if pair <= 0 { closed=true;break }
                pairSums.append(pairSums.last.map { min($0,pair) } ?? pair)
                lag += 2
            }
            truncated = !closed && lag+1 < splitLength
            let autocorrelationTime = max(1,1+2*pairSums.reduce(0,+))
            effective = truncated ? Double(splitCount)
                : min(total,max(Double(splitCount),total/autocorrelationTime))
        }
        let mcse = sqrt(max(0,variancePlus)/effective)
        guard [overallMean,overallVariance,effective,mcse,rHat].allSatisfy(\.isFinite),
              effective > 0, effective <= total+1e-9, rHat >= 0 else {
            throw VivoChemistryError.convergence("correlated-series diagnostics overflow")
        }
        return .init(chainCount: chains.count,retainedSamplesPerChain: count,mean: overallMean,
            sampleVariance: overallVariance,autocorrelationEffectiveSampleSize: effective,
            monteCarloStandardError: mcse,splitRHat: rHat,
            autocorrelationSequenceTruncated:truncated,interpretation: interpretation)
    }
}
