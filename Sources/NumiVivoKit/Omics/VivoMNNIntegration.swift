import Foundation
import NumiVivoCore

/// Scale-aware latent integration. Algorithmic assembly follows Scanorama 1.7.4;
/// see Documentation/ThirdParty/Scanorama-LICENSE.txt. Native matching has explicit
/// grouped-source-index ties; legacy scalar and opt-in tiled reductions are explicit.
public struct VivoMNNIntegrationOptions: Codable, Sendable, Equatable {
    public enum Kernel: String, Codable, Sendable { case tiledGaussian }
    /// Nil retains the historical scalar reduction and canonical plan encoding.
    public var kernel: Kernel? = nil
    static let maximumNeighbors = 100
    /// Optional, complete observation-sample map for an explicitly protected
    /// categorical group. Nil preserves historical encoding and behavior.
    public var protectedSampleGroups: [String: String]? = nil
    public var covariate: VivoSingleCellIntegrationOptions.Covariate = .donor
    public var neighbors: Int = 20
    public var sigma: Double = 15
    public var minimumAlignment: Double = 0.1
    public var maximumWork: Int = 100_000_000_000
    public var maximumResidentBytes: Int = 536_870_912
    public init() {}
    private enum CodingKeys: String, CodingKey { case protectedSampleGroups, covariate, neighbors, sigma, minimumAlignment, maximumWork, maximumResidentBytes, kernel }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["protectedSampleGroups","covariate","neighbors","sigma","minimumAlignment","maximumWork","maximumResidentBytes","kernel"])
        let c = try decoder.container(keyedBy: CodingKeys.self)
        protectedSampleGroups = try c.decodeIfPresent([String: String].self, forKey: .protectedSampleGroups)
        kernel = try c.decodeIfPresent(Kernel.self, forKey: .kernel)
        covariate = try c.decodeIfPresent(VivoSingleCellIntegrationOptions.Covariate.self, forKey: .covariate) ?? .donor
        neighbors = try c.decodeIfPresent(Int.self, forKey: .neighbors) ?? 20
        sigma = try c.decodeIfPresent(Double.self, forKey: .sigma) ?? 15
        minimumAlignment = try c.decodeIfPresent(Double.self, forKey: .minimumAlignment) ?? 0.1
        maximumWork = try c.decodeIfPresent(Int.self, forKey: .maximumWork) ?? 100_000_000_000
        maximumResidentBytes = try c.decodeIfPresent(Int.self, forKey: .maximumResidentBytes) ?? 536_870_912
    }
    public func validate() throws {
        try VivoIntegrationProtection.validate(protectedSampleGroups)
        guard (1...Self.maximumNeighbors).contains(neighbors), sigma.isFinite, (0.001...1000).contains(sigma),
              minimumAlignment.isFinite, (0...1).contains(minimumAlignment),
              (1...10_000_000_000_000).contains(maximumWork),
              (1...4_294_967_296).contains(maximumResidentBytes) else { throw VivoOmicsError.invalid("MNN options") }
    }
}
public struct VivoMNNAlignment: Codable, Sendable, Equatable {
    public let firstLevel: Int
    public let secondLevel: Int
    public let anchorOffset: Int
    public let anchors: Int
    public let firstMatchedCells: Int
    public let secondMatchedCells: Int
    public let score: Double
}
public struct VivoMNNAssemblyStep: Codable, Sendable, Equatable {
    public let alignment: Int
    public let targetLevels: [Int]
    public let referenceLevels: [Int]
    public let anchors: Int
    public let correctedCells: Int
    public let zeroWeightCells: Int
    public let minimumWeight: Double
    public let maximumWeight: Double
    public let correctionRMS: Double
    public let status: String
}
public struct VivoMNNIntegrationReport: Codable, Sendable, Equatable {
    public let method: String
    public let cells: Int
    public let components: Int
    public let levels: [String]
    public let cellLevels: [Int]
    public let medianRowNorm: Double
    public let zeroDirectionCells: Int
    public let alignments: [VivoMNNAlignment]
    public let assemblyOrder: [Int]
    public let steps: [VivoMNNAssemblyStep]
    public let panoramas: [[Int]]
    public let distanceScalarTerms: Int
    public let kernelScalarTerms: Int
    public let estimatedMaximumLatentResidentBytes: Int
    public let qualification: String
}
struct VivoMNNAnchor: Sendable { let first: Int; let second: Int }
struct VivoMNNIntegrationSolution {
    let scores: [Double]
    let anchors: [VivoMNNAnchor]
    let anchorDistances: [Double]
    let report: VivoMNNIntegrationReport
}

enum VivoMNNIntegration {
    static func product(_ factors: [Int], limit: Int, reason: String) throws -> Int {
        var value = 1
        for factor in factors {
            let next = value.multipliedReportingOverflow(by: factor)
            guard factor >= 0, !next.overflow, next.partialValue <= limit else { throw VivoOmicsError.limit(reason) }
            value = next.partialValue
        }
        return value
    }
    /// Upper bound for the algorithm's latent arrays, including two anchor-coordinate
    /// buffers and heap/identity storage. Source metadata and JSON encoding are separate.
    static func admittedBytes(rows n: Int, dimensions d: Int, options o: VivoMNNIntegrationOptions) throws -> Int {
        try o.validate()
        guard (2...VivoPCAStorageLimits.maximumRows).contains(n), (1...64).contains(d) else { throw VivoOmicsError.invalid("MNN PCA axes") }
        return try product([n, 3*d*8 + 2*o.neighbors*16 + 80 + o.neighbors*40], limit: o.maximumResidentBytes, reason: "MNN latent memory budget")
    }
    // Keep pointer/callback bridging outside the scalar solver's optimization unit.
    @inline(never) private static func tiledBias(source: [Double], bias: [Double], queries: [Double],
        rows: Int, dimensions d: Int, sigma: Double, maximumWork: Int,
        deltas: inout [Double], totals: inout [Double], workspace: inout [Double]) throws -> Int {
        var evaluatedTerms: UInt64 = 0
        let status = source.withUnsafeBufferPointer { src in bias.withUnsafeBufferPointer { b in
            queries.withUnsafeBufferPointer { q in deltas.withUnsafeMutableBufferPointer { delta in
                totals.withUnsafeMutableBufferPointer { total in workspace.withUnsafeMutableBufferPointer { work in
                    nvivo_omics_gaussian_bias(src.baseAddress,b.baseAddress,UInt32(source.count/d),
                        q.baseAddress,UInt32(rows),UInt32(d),sigma,UInt64(src.count),UInt64(q.count),UInt64(maximumWork),
                        delta.baseAddress,UInt64(delta.count),total.baseAddress,UInt64(total.count),
                        work.baseAddress,UInt64(work.count),&evaluatedTerms,{ _ in Task.isCancelled ? 1 : 0 },nil)
                } }
            } }
        } }
        if status == 5 { throw CancellationError() }
        if status == 2 { throw VivoOmicsError.limit("MNN scalar underflow fallback work budget") }
        guard status == 0 else { throw VivoOmicsError.invalid("MNN tiled Gaussian status: \(status)") }
        return Int(evaluatedTerms)
    }
    static func run(cells: [VivoOmicsCellIdentity], scores x: [Double], dimensions d: Int,
                    samples: [VivoOmicsSample], options o: VivoMNNIntegrationOptions) throws -> VivoMNNIntegrationSolution {
        let n = cells.count, k = o.neighbors
        let admitted = try admittedBytes(rows: n, dimensions: d, options: o)
        guard x.count == n*d, x.allSatisfy(\.isFinite), Set(samples.map(\.id)).count == samples.count,
              Set(cells).count == n else { throw VivoOmicsError.invalid("MNN PCA identities or values") }
        let sampleMap = Dictionary(uniqueKeysWithValues: samples.map { ($0.id,$0) })
        let names = try cells.map { cell -> String in
            guard let sample = sampleMap[cell.sampleID], let name = o.covariate == .donor ? sample.donorID : sample.batchID,
                  !name.isEmpty, name != "unreported" else { throw VivoOmicsError.invalid("MNN requires known selected covariate") }
            return name
        }
        let levels = Set(names).sorted(), bCount = levels.count
        guard (2...128).contains(bCount) else { throw VivoOmicsError.invalid("MNN requires 2...128 covariate levels") }
        let levelMap = Dictionary(uniqueKeysWithValues: levels.enumerated().map { ($0.element,$0.offset) })
        let batch = names.map { levelMap[$0]! }
        var rows = [[Int]](repeating: [], count: bCount), conditionLevels: [String:Set<Int>] = [:]
        for i in 0..<n {
            rows[batch[i]].append(i)
            conditionLevels[sampleMap[cells[i].sampleID]!.condition, default: []].insert(batch[i])
        }
        guard rows.allSatisfy({ $0.count >= k }) else { throw VivoOmicsError.invalid("MNN covariate level has fewer cells than neighbors") }
        var connected: Set<Int> = [0]
        for _ in 0..<bCount {
            let before = connected.count
            for group in conditionLevels.values where !connected.isDisjoint(with: group) { connected.formUnion(group) }
            if connected.count == before { break }
        }
        guard connected.count == bCount else { throw VivoOmicsError.invalid("MNN covariate is confounded with condition") }
        try VivoIntegrationProtection.check(o.protectedSampleGroups, cells: cells,
            samples: sampleMap, covariate: o.covariate, levelIDs: levelMap)
        var pairCount = 0, preceding = 0
        for group in rows { pairCount += preceding*group.count; preceding += group.count }
        let distanceTerms = try product([pairCount,d], limit: o.maximumWork, reason: "MNN matching work budget")
        var norms = [Double](repeating: 0, count: n)
        for i in 0..<n {
            for j in 0..<d { norms[i] += x[i*d+j]*x[i*d+j] }
            norms[i] = sqrt(norms[i])
        }
        guard norms.allSatisfy(\.isFinite) else { throw VivoOmicsError.invalid("MNN PCA norm overflow") }
        let orderedNorms = norms.sorted()
        let scale = n%2 == 1 ? orderedNorms[n/2] : orderedNorms[n/2-1]/2 + orderedNorms[n/2]/2
        guard scale.isFinite, scale > 0 else { throw VivoOmicsError.invalid("MNN needs positive median PCA norm") }
        var unit = x, corrected = x
        for i in 0..<n { for j in 0..<d {
            unit[i*d+j] = norms[i] > 0 ? x[i*d+j]/norms[i] : 0
            corrected[i*d+j] /= scale
        } }
        guard corrected.allSatisfy(\.isFinite) else { throw VivoOmicsError.invalid("MNN scaled PCA overflow") }
        var rank = [Int](repeating: 0, count: n)
        for (position,i) in rows.flatMap({ $0 }).enumerated() { rank[i] = position }
        var lower = [Int](repeating: -1, count: n*k), upper = lower
        var lowerD = [Double](repeating: .infinity, count: n*k), upperD = lowerD
        var lowerN = [Int](repeating: 0, count: n), upperN = lowerN
        func worse(_ a: Int, _ ad: Double, _ b: Int, _ bd: Double) -> Bool { ad > bd || (ad == bd && rank[a] > rank[b]) }
        func offer(row: Int, candidate: Int, distance: Double, ids: inout [Int], ds: inout [Double], counts: inout [Int]) {
            let base = row*k
            if counts[row] < k {
                var slot = counts[row]; counts[row] += 1
                ids[base+slot] = candidate; ds[base+slot] = distance
                while slot > 0 {
                    let parent = (slot-1)/2
                    if !worse(ids[base+slot],ds[base+slot],ids[base+parent],ds[base+parent]) { break }
                    ids.swapAt(base+slot,base+parent); ds.swapAt(base+slot,base+parent); slot = parent
                }
            } else {
                if !worse(ids[base],ds[base],candidate,distance) { return }
                ids[base] = candidate; ds[base] = distance; var slot = 0
                while slot*2+1 < k {
                    var child = slot*2+1
                    if child+1 < k && worse(ids[base+child+1],ds[base+child+1],ids[base+child],ds[base+child]) { child += 1 }
                    if !worse(ids[base+child],ds[base+child],ids[base+slot],ds[base+slot]) { break }
                    ids.swapAt(base+slot,base+child); ds.swapAt(base+slot,base+child); slot = child
                }
            }
        }
        for first in 0..<bCount { for second in (first+1)..<bCount {
            for i in rows[first] {
                try Task.checkCancellation()
                for j in rows[second] {
                    var distance = 0.0
                    for column in 0..<d { distance += abs(unit[i*d+column]-unit[j*d+column]) }
                    offer(row: i, candidate: j, distance: distance, ids: &upper, ds: &upperD, counts: &upperN)
                    offer(row: j, candidate: i, distance: distance, ids: &lower, ds: &lowerD, counts: &lowerN)
                }
            }
        } }
        var pairs = [[VivoMNNAnchor]](repeating: [], count: bCount*bCount)
        for i in 0..<n { for slot in 0..<upperN[i] {
            let j = upper[i*k+slot]
            if (0..<lowerN[j]).contains(where: { lower[j*k+$0] == i }) {
                pairs[batch[i]*bCount+batch[j]].append(.init(first: i, second: j))
            }
        } }
        var anchors: [VivoMNNAnchor] = [], anchorDistances: [Double] = [], alignments: [VivoMNNAlignment] = []
        for first in 0..<bCount { for second in (first+1)..<bCount {
            let key = first*bCount+second
            pairs[key].sort { $0.first == $1.first ? $0.second < $1.second : $0.first < $1.first }
            let matched = pairs[key]; if matched.isEmpty { continue }
            let firstCount = Set(matched.map(\.first)).count, secondCount = Set(matched.map(\.second)).count
            alignments.append(.init(firstLevel: first, secondLevel: second, anchorOffset: anchors.count, anchors: matched.count,
                firstMatchedCells: firstCount, secondMatchedCells: secondCount,
                score: max(Double(firstCount)/Double(rows[first].count), Double(secondCount)/Double(rows[second].count))))
            anchors += matched
            for p in matched {
                var distance = 0.0
                for column in 0..<d { distance += abs(unit[p.first*d+column]-unit[p.second*d+column]) }
                anchorDistances.append(distance)
            }
        } }
        // Release nearest-neighbor heaps before assembly-coordinate buffers.
        lower.removeAll(); upper.removeAll(); lowerD.removeAll(); upperD.removeAll(); unit.removeAll()
        let order = alignments.indices.filter { alignments[$0].score > o.minimumAlignment }.sorted {
            let a = alignments[$0], b = alignments[$1]
            if a.score != b.score { return a.score > b.score }
            if a.firstLevel != b.firstLevel { return a.firstLevel > b.firstLevel }
            return a.secondLevel > b.secondLevel
        }
        var maximumLatentBytes = admitted
        var panoramas: [[Int]] = [], visits = [Int](repeating: 0, count: bCount), steps: [VivoMNNAssemblyStep] = [], kernelTerms = 0
        func oriented(_ source: Int, _ reference: Int) -> [VivoMNNAnchor] {
            if source < reference { return pairs[source*bCount+reference] }
            return pairs[reference*bCount+source].map { .init(first: $0.second, second: $0.first) }
        }
        func apply(targetLevels: [Int], referenceLevels: [Int], matched: [VivoMNNAnchor], alignment: Int) throws -> VivoMNNAssemblyStep {
            guard !matched.isEmpty else { throw VivoOmicsError.invalid("MNN assembly has no anchors") }
            let selected = targetLevels.flatMap { rows[$0] }
            let cost = try product([selected.count,matched.count,d], limit: o.maximumWork, reason: "MNN kernel work budget")
            guard cost <= o.maximumWork-distanceTerms-kernelTerms else { throw VivoOmicsError.limit("MNN total work budget") }
            kernelTerms += cost
            let bufferBytes = try product([matched.count,d,16], limit: o.maximumResidentBytes, reason: "MNN anchor memory budget")
            let tileBytes = o.kernel == nil ? 0 : (256*(d+32) + 2*32*d + 32)*8
            guard tileBytes <= o.maximumResidentBytes-admitted, bufferBytes <= o.maximumResidentBytes-admitted-tileBytes else { throw VivoOmicsError.limit("MNN assembly memory budget") }
            maximumLatentBytes = max(maximumLatentBytes,admitted+bufferBytes+tileBytes)
            // Duplicate anchor records are retained (including panorama merge pairs).
            // Snapshot both coordinates and differences before mutating any query row.
            var source = [Double](repeating: 0, count: matched.count*d), bias = source
            for (a,p) in matched.enumerated() { for j in 0..<d {
                source[a*d+j] = corrected[p.first*d+j]
                bias[a*d+j] = corrected[p.second*d+j]-corrected[p.first*d+j]
            } }
            var zero = 0, minWeight = Double.infinity, maxWeight = 0.0, sumSquared = 0.0
            if o.kernel == .tiledGaussian {
                var queries = [Double](repeating: 0, count: 32*d)
                var deltas = queries, totals = [Double](repeating: 0, count: 32)
                var workspace = [Double](repeating: 0, count: 256*(d+32))
                for start in stride(from: 0, to: selected.count, by: 32) {
                    try Task.checkCancellation()
                    let count = min(32,selected.count-start)
                    for q in 0..<count { for j in 0..<d { queries[q*d+j] = corrected[selected[start+q]*d+j] } }
                    let baseTerms = count*matched.count*d
                    let allowedTerms = baseTerms + o.maximumWork-distanceTerms-kernelTerms
                    let evaluatedTerms = try tiledBias(source: source,bias: bias,queries: queries,
                        rows: count,dimensions: d,sigma: o.sigma,maximumWork: allowedTerms,
                        deltas: &deltas,totals: &totals,workspace: &workspace)
                    kernelTerms += evaluatedTerms-baseTerms
                    for q in 0..<count {
                        let total = totals[q], i = selected[start+q]
                        minWeight = min(minWeight,total); maxWeight = max(maxWeight,total)
                        if total == 0 { zero += 1 }
                        for j in 0..<d {
                            let value = deltas[q*d+j]
                            corrected[i*d+j] += value; sumSquared += value*value
                        }
                    }
                }
            } else {
                var scalar = NVivoGaussianScalarReport()
                let status = source.withUnsafeBufferPointer { src in bias.withUnsafeBufferPointer { b in
                    selected.withUnsafeBufferPointer { indices in corrected.withUnsafeMutableBufferPointer { values in
                        nvivo_omics_gaussian_scalar_apply(src.baseAddress,b.baseAddress,UInt32(matched.count),UInt32(d),
                            o.sigma,UInt64(src.count),values.baseAddress,UInt32(n),UInt64(values.count),indices.baseAddress,
                            UInt32(indices.count),UInt64(cost),&scalar,{ _ in Task.isCancelled ? 1 : 0 },nil)
                    } }
                } }
                if status == 5 { throw CancellationError() }
                guard status == 0 else { throw VivoOmicsError.invalid("MNN scalar phase status: \(status)") }
                zero = Int(scalar.zero_weight_rows); minWeight = scalar.minimum_weight
                maxWeight = scalar.maximum_weight; sumSquared = scalar.sum_squared
            }
            return .init(alignment: alignment, targetLevels: targetLevels, referenceLevels: referenceLevels,
                anchors: matched.count, correctedCells: selected.count, zeroWeightCells: zero, minimumWeight: minWeight,
                maximumWeight: maxWeight, correctionRMS: sqrt(sumSquared/Double(selected.count))*scale, status: "corrected")
        }
        for alignment in order {
            var i = alignments[alignment].firstLevel, j = alignments[alignment].secondLevel
            visits[i] += 1; visits[j] += 1
            if visits[i] > 3 && visits[j] > 3 {
                steps.append(.init(alignment: alignment, targetLevels: [], referenceLevels: [], anchors: 0, correctedCells: 0,
                    zeroWeightCells: 0, minimumWeight: 0, maximumWeight: 0, correctionRMS: 0, status: "both-level-visit-counts-exceed-three")); continue
            }
            var pi = panoramas.firstIndex { $0.contains(i) }, pj = panoramas.firstIndex { $0.contains(j) }
            if pi == nil && pj == nil {
                if rows[i].count < rows[j].count { swap(&i,&j) }
                panoramas.append([i]); pi = panoramas.count-1
            }
            if pi == nil {
                let reference = panoramas[pj!]
                let match = reference.flatMap { oriented(i,$0) }
                steps.append(try apply(targetLevels: [i], referenceLevels: reference, matched: match, alignment: alignment)); panoramas[pj!].append(i)
            } else if pj == nil {
                let reference = panoramas[pi!]
                let match = reference.flatMap { oriented(j,$0) }
                steps.append(try apply(targetLevels: [j], referenceLevels: reference, matched: match, alignment: alignment)); panoramas[pi!].append(j)
            } else {
                let target = panoramas[pi!], reference = panoramas[pj!], one = oriented(i,j)
                steps.append(try apply(targetLevels: target, referenceLevels: reference, matched: one+one, alignment: alignment))
                if pi! != pj! { panoramas[pi!] += reference; panoramas.remove(at: pj!) }
            }
        }
        for b in 0..<bCount where !panoramas.contains(where: { $0.contains(b) }) { panoramas.append([b]) }
        for i in corrected.indices { corrected[i] *= scale }
        guard corrected.allSatisfy(\.isFinite) else { throw VivoOmicsError.invalid("MNN corrected PCA nonfinite") }
        let report = VivoMNNIntegrationReport(method: o.kernel == nil ? "mutual-nearest-neighbor-panorama-median-norm-Double-v1" : "mutual-nearest-neighbor-panorama-median-norm-tiled-Gaussian-Double-v1", cells: n, components: d,
            levels: levels, cellLevels: batch, medianRowNorm: scale, zeroDirectionCells: norms.filter { $0 == 0 }.count,
            alignments: alignments, assemblyOrder: order, steps: steps, panoramas: panoramas,
            distanceScalarTerms: distanceTerms, kernelScalarTerms: kernelTerms, estimatedMaximumLatentResidentBytes: maximumLatentBytes,
            qualification: (o.kernel == nil ? "" : "All-anchor 32-query/256-anchor direct-distance tiles with Accelerate vector exp and BLAS FP64 reduction, with budgeted scalar fallback below 1e-280 total weight; framework internal workspace is outside the latent estimate. ")+"Native exact Manhattan mutual neighbors on row-normalized PCA; Gaussian panorama correction in global median-row-norm units; explicit grouped-source-index ties. Original counts/PCA retained. Labels and conditions do not enter correction; condition is an eligibility check only. Latent arrays are resident and budgeted, not a million-cell/out-of-core or Metal qualification. Source metadata and artifact encoding memory are outside the latent estimate. Zero kernel weights preserve the affected row and remain reported. Replay/numerical agreement is separate from biological or prospective validation.")
        return .init(scores: corrected, anchors: anchors, anchorDistances: anchorDistances, report: report)
    }
}
