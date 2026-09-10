import Foundation

public struct VivoOmicsPrecisionLowessFit: Codable, Sendable, Equatable {
    public let fitted: [Double]
    public let mode: String
    public let delta: Double
    public let anchorCount: Int
    public let neighborhoodVisits: Int
    public let constantFitAnchors: Int
    public let clippedLowIndices: [Int]
    public let clippedHighIndices: [Int]
}

/// A single precision-weighted local-linear pass; robust prior reweighting is external.
public enum VivoOmicsPrecisionLowess {
    public static func fit(x: [Double], y: [Double], weights: [Double], span: Double,
        maximumNeighborhoodVisits: Int = 100_000_000) throws -> VivoOmicsPrecisionLowessFit {
        let n = x.count
        guard (3...100_000).contains(n), y.count == n, weights.count == n,
              x.allSatisfy(\.isFinite), y.allSatisfy(\.isFinite),
              weights.allSatisfy({ $0.isFinite && $0 >= 0 }), weights.contains(where: { $0 > 0 }),
              span.isFinite, span > 0, span <= 1, maximumNeighborhoodVisits > 0 else {
            throw VivoOmicsStatisticsError.invalid("precision LOWESS coordinates, weights, span or work bound")
        }
        let low = weights.indices.filter { weights[$0] < 1e-8 }, high = weights.indices.filter { weights[$0] > 100 }
        let clipped = weights.map { min(100,max(1e-8,$0)) }
        if clipped.max()!-clipped.min()! < 1e-15 {
            let fit = try VivoOmicsRobustLowess.fit(x: x,y: y,span: span,robustnessIterations: 0,
                maximumNeighborhoodVisits: maximumNeighborhoodVisits)
            return .init(fitted: fit.fitted,mode: "equal-weight-lowess",delta: 0.01*(x.max()!-x.min()!),
                anchorCount: fit.anchorCount,neighborhoodVisits: fit.neighborhoodVisits,
                constantFitAnchors: fit.constantFitAnchors,clippedLowIndices: low,clippedHighIndices: high)
        }
        let order = x.indices.sorted { x[$0] == x[$1] ? $0 < $1 : x[$0] < x[$1] }
        let xs = order.map { x[$0] }, ys = order.map { y[$0] }, ws = order.map { clipped[$0] }
        let range = xs[n-1]-xs[0]
        guard range.isFinite else { throw VivoOmicsStatisticsError.invalid("precision LOWESS coordinate range overflow") }
        var visits = 0, constants = 0
        func regression(anchor: Int, lower: Int, upper: Int, radius: Double, local: Bool) throws -> Double {
            guard upper-lower <= (maximumNeighborhoodVisits-visits)/2 else {
                throw VivoOmicsStatisticsError.invalid("precision LOWESS exhausted neighborhood visits")
            }
            visits += 2*(upper-lower)
            let scale = radius > 0 ? radius : 1
            func weight(_ j: Int) -> Double {
                if !local || radius == 0 { return ws[j] }
                let distance = min(1,abs(xs[j]-xs[anchor])/radius)
                let kernel = 1-distance*distance*distance
                return ws[j]*kernel*kernel*kernel
            }
            var total = 0.0, mx = 0.0, my = 0.0
            for j in lower..<upper {
                let w = weight(j)
                if w > 0 {
                    let next = total+w, fraction = w/next
                    mx += fraction*((xs[j]-xs[anchor])/scale-mx)
                    my += fraction*(ys[j]-my); total = next
                }
            }
            guard total > 0 else { throw VivoOmicsStatisticsError.invalid("empty precision-weighted neighborhood") }
            var variance = 0.0, covariance = 0.0
            for j in lower..<upper {
                let w = weight(j)/total, dx = (xs[j]-xs[anchor])/scale-mx
                variance += w*dx*dx; covariance += w*dx*(ys[j]-my)
            }
            if variance > 0 { return my-mx*covariance/variance }
            constants += 1; return my
        }
        if Double(n) < 4+1/span {
            // Predict one global weighted line; two endpoints determine all values.
            let first = try regression(anchor: 0,lower: 0,upper: n,radius: range,local: false)
            let last = range > 0 ? try regression(anchor: n-1,lower: 0,upper: n,radius: range,local: false) : first
            var restored = [Double](repeating: first,count: n)
            if range > 0 { for i in 0..<n { restored[order[i]] = first+(last-first)*(xs[i]-xs[0])/range } }
            guard restored.allSatisfy(\.isFinite) else {
                throw VivoOmicsStatisticsError.invalid("nonfinite global precision LOWESS fit")
            }
            return .init(fitted: restored,mode: "global-weighted-line",delta: range,anchorCount: range > 0 ? 2 : 1,
                neighborhoodVisits: visits,constantFitAnchors: constants,clippedLowIndices: low,clippedHighIndices: high)
        }
        var delta = 0.0
        if n > 200 {
            let gaps = (1..<n).map { xs[$0]-xs[$0-1] }.sorted()
            var cumulative = gaps
            for i in 1..<gaps.count { cumulative[i] += cumulative[i-1] }
            delta = (0..<200).map { cumulative[n-2-$0]/Double(200-$0) }.min()!
        }
        var anchors = [0], at = 0
        while at < n-1 {
            var next = at+1
            while next < n, xs[next] <= xs[at]+delta { next += 1 }
            if next >= n { break }
            anchors.append(next); at = next
        }
        if xs[anchors.last!] < xs[n-1] { anchors.append(n-1) }
        let target = span*ws.reduce(0,+)
        var fitted = [Double](repeating: 0,count: n)
        for anchor in anchors {
            try Task.checkCancellation()
            var lower = anchor, upper = anchor+1, mass = ws[anchor]
            while mass < target, lower > 0 || upper < n {
                guard visits < maximumNeighborhoodVisits else { throw VivoOmicsStatisticsError.invalid("precision LOWESS exhausted neighborhood visits") }
                visits += 1
                if upper == n || (lower > 0 && xs[anchor]-xs[lower-1] <= xs[upper]-xs[anchor]) {
                    lower -= 1; mass += ws[lower]
                } else { mass += ws[upper]; upper += 1 }
            }
            let radius = max(xs[anchor]-xs[lower],xs[upper-1]-xs[anchor])
            while lower > 0, xs[anchor]-xs[lower-1] <= radius { lower -= 1 }
            while upper < n, xs[upper]-xs[anchor] <= radius { upper += 1 }
            fitted[anchor] = try regression(anchor: anchor,lower: lower,upper: upper,radius: radius,local: true)
        }
        if anchors.count > 1 {
            for k in 1..<anchors.count {
                let a = anchors[k-1], b = anchors[k], width = xs[b]-xs[a]
                for j in (a+1)..<b { let t = (xs[j]-xs[a])/width; fitted[j] = (1-t)*fitted[a]+t*fitted[b] }
            }
        }
        if let last = anchors.last, last+1 < n { for j in (last+1)..<n { fitted[j] = fitted[last] } }
        guard fitted.allSatisfy(\.isFinite) else { throw VivoOmicsStatisticsError.invalid("nonfinite precision LOWESS fit") }
        var restored = fitted
        for i in 0..<n { restored[order[i]] = fitted[i] }
        return .init(fitted: restored,mode: "weighted-neighborhood-lowess",delta: delta,anchorCount: anchors.count,
            neighborhoodVisits: visits,constantFitAnchors: constants,clippedLowIndices: low,clippedHighIndices: high)
    }
}
