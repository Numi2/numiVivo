import Foundation

public struct VivoOmicsLowessFit: Codable, Sendable, Equatable {
    public let fitted: [Double]
    public let robustnessWeights: [Double]
    public let completedRobustnessUpdates: Int
    public let anchorCount: Int
    public let neighborhoodVisits: Int
    public let constantFitAnchors: Int
    public let emptyWeightAnchors: Int
}

/// Robust local linear regression with tricube neighborhoods and Tukey weights.
/// Input coordinates are sorted stably, with fitted values restored to input order.
public enum VivoOmicsRobustLowess {
    public static func fit(x: [Double], y: [Double], span: Double = 0.5,
        robustnessIterations: Int = 3, deltaFraction: Double = 0.01,
        maximumNeighborhoodVisits: Int = 100_000_000) throws -> VivoOmicsLowessFit {
        guard x.count == y.count, x.count >= 2, x.count <= 100_000,
              x.allSatisfy(\.isFinite), y.allSatisfy(\.isFinite), span.isFinite, span > 0, span <= 1,
              (0...20).contains(robustnessIterations), deltaFraction.isFinite,
              (0...1).contains(deltaFraction), maximumNeighborhoodVisits > 0 else {
            throw VivoOmicsStatisticsError.invalid("LOWESS coordinates, span, iterations or work bound")
        }
        let order = x.indices.sorted { x[$0] == x[$1] ? $0 < $1 : x[$0] < x[$1] }
        let xs = order.map { x[$0] }, ys = order.map { y[$0] }, n = x.count
        let range = xs[n-1]-xs[0]
        guard range.isFinite else { throw VivoOmicsStatisticsError.invalid("LOWESS coordinate range overflow") }
        let delta = deltaFraction*range, neighbors = max(2,min(n,Int(floor(span*Double(n)))))
        var anchors = [0], at = 0
        while at < n-1 {
            var next = at+1
            while next < n, xs[next] <= xs[at]+delta { next += 1 }
            next = max(at+1,next-1)
            while next < n, xs[next] == xs[at] { next += 1 }
            if next >= n { break }
            anchors.append(next); at = next
        }
        var fitted = [Double](repeating: 0,count: n), robustness = [Double](repeating: 1,count: n)
        var visits = 0, constant = 0, empty = 0, updates = 0
        for pass in 0...robustnessIterations {
            try Task.checkCancellation()
            var left = 0
            for anchor in anchors {
                while left+neighbors < n, xs[anchor]-xs[left] > xs[left+neighbors]-xs[anchor] { left += 1 }
                let radius = max(xs[anchor]-xs[left],xs[left+neighbors-1]-xs[anchor])
                var lower = left, upper = left+neighbors
                while lower > 0, xs[anchor]-xs[lower-1] <= radius { lower -= 1 }
                while upper < n, xs[upper]-xs[anchor] <= radius { upper += 1 }
                guard upper-lower <= (maximumNeighborhoodVisits-visits)/2 else {
                    throw VivoOmicsStatisticsError.invalid("LOWESS exhausted \(maximumNeighborhoodVisits) neighborhood visits")
                }
                visits += 2*(upper-lower)
                func weight(_ j: Int) -> Double {
                    let distance = abs(xs[j]-xs[anchor])
                    if radius == 0 { return distance == 0 ? robustness[j] : 0 }
                    let t = min(1,distance/radius), kernel = 1-t*t*t
                    return kernel*kernel*kernel*robustness[j]
                }
                let scale = radius > 0 ? radius : 1
                var total = 0.0, meanX = 0.0, meanY = 0.0
                for j in lower..<upper {
                    let w = weight(j)
                    if w > 0 {
                        let next = total+w, fraction = w/next
                        meanX += fraction*((xs[j]-xs[anchor])/scale-meanX)
                        meanY += fraction*(ys[j]-meanY); total = next
                    }
                }
                if total == 0 { fitted[anchor] = ys[anchor]; empty += 1; continue }
                var variance = 0.0, covariance = 0.0
                for j in lower..<upper {
                    let w = weight(j)/total, dx = (xs[j]-xs[anchor])/scale-meanX
                    variance += w*dx*dx; covariance += w*dx*(ys[j]-meanY)
                }
                if variance > 0, sqrt(variance)*scale > 0.001*range {
                    fitted[anchor] = meanY-meanX*covariance/variance
                } else { fitted[anchor] = meanY; constant += 1 }
            }
            if anchors.count > 1 {
                for k in 1..<anchors.count {
                    let a = anchors[k-1], b = anchors[k], width = xs[b]-xs[a]
                    for j in (a+1)..<b {
                        let t = (xs[j]-xs[a])/width
                        fitted[j] = (1-t)*fitted[a]+t*fitted[b]
                    }
                }
            }
            if let last = anchors.last, last+1 < n { for j in (last+1)..<n { fitted[j] = fitted[last] } }
            guard fitted.allSatisfy(\.isFinite) else { throw VivoOmicsStatisticsError.invalid("nonfinite LOWESS prediction") }
            if pass == robustnessIterations { break }
            let absolute = ys.indices.map { abs(ys[$0]-fitted[$0]) }
            let cutoff = 6*(try VivoOmicsLinearStatistics.median(absolute))
            let meanAbsolute = absolute.reduce(0) { $0+$1/Double(n) }
            guard cutoff.isFinite, meanAbsolute.isFinite else { throw VivoOmicsStatisticsError.invalid("LOWESS residual scale overflow") }
            if cutoff <= 1e-7*meanAbsolute { break }
            robustness = absolute.map { value in
                if value >= cutoff { return 0 }
                let t = value/cutoff, weight = 1-t*t
                return weight*weight
            }
            updates += 1
        }
        var restored = fitted, restoredWeights = robustness
        for i in 0..<n { restored[order[i]] = fitted[i]; restoredWeights[order[i]] = robustness[i] }
        return .init(fitted: restored,robustnessWeights: restoredWeights,completedRobustnessUpdates: updates,
            anchorCount: anchors.count,neighborhoodVisits: visits,constantFitAnchors: constant,emptyWeightAnchors: empty)
    }
}
