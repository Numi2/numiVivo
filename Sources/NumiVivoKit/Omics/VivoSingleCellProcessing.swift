import Foundation

public struct VivoOmicsCellIdentity: Codable, Sendable, Equatable, Hashable {
    public let sampleID: String
    public let barcode: String
    public init(sampleID: String, barcode: String) { self.sampleID = sampleID; self.barcode = barcode }
    private enum CodingKeys: String, CodingKey { case sampleID, barcode }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["sampleID", "barcode"])
        let c = try decoder.container(keyedBy: CodingKeys.self)
        sampleID = try c.decode(String.self, forKey: .sampleID); barcode = try c.decode(String.self, forKey: .barcode)
    }
}
public enum VivoOmicsMissingMitochondrialPolicy: String, Codable, Sendable { case rejectCell, retainCell }
public struct VivoSingleCellFilterPolicy: Codable, Sendable, Equatable {
    public var minimumCounts: UInt64 = 1
    public var maximumCounts: UInt64?
    public var minimumDetectedFeatures: Int = 1
    public var maximumDetectedFeatures: Int?
    public var maximumMitochondrialFraction: Double?
    public var missingMitochondrial: VivoOmicsMissingMitochondrialPolicy = .rejectCell
    public var includedSampleIDs: [String]?
    public var excludedCells: [VivoOmicsCellIdentity] = []
    public init() {}
    private enum CodingKeys: String, CodingKey {
        case minimumCounts, maximumCounts, minimumDetectedFeatures, maximumDetectedFeatures
        case maximumMitochondrialFraction, missingMitochondrial, includedSampleIDs, excludedCells
    }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["minimumCounts", "maximumCounts", "minimumDetectedFeatures",
            "maximumDetectedFeatures", "maximumMitochondrialFraction", "missingMitochondrial", "includedSampleIDs", "excludedCells"])
        let c = try decoder.container(keyedBy: CodingKeys.self)
        minimumCounts = try c.decodeIfPresent(UInt64.self, forKey: .minimumCounts) ?? 1
        maximumCounts = try c.decodeIfPresent(UInt64.self, forKey: .maximumCounts)
        minimumDetectedFeatures = try c.decodeIfPresent(Int.self, forKey: .minimumDetectedFeatures) ?? 1
        maximumDetectedFeatures = try c.decodeIfPresent(Int.self, forKey: .maximumDetectedFeatures)
        maximumMitochondrialFraction = try c.decodeIfPresent(Double.self, forKey: .maximumMitochondrialFraction)
        missingMitochondrial = try c.decodeIfPresent(VivoOmicsMissingMitochondrialPolicy.self, forKey: .missingMitochondrial) ?? .rejectCell
        includedSampleIDs = try c.decodeIfPresent([String].self, forKey: .includedSampleIDs)
        excludedCells = try c.decodeIfPresent([VivoOmicsCellIdentity].self, forKey: .excludedCells) ?? []
    }
    public func validate() throws {
        guard minimumDetectedFeatures >= 0, maximumCounts.map({ $0 >= minimumCounts }) ?? true,
              maximumDetectedFeatures.map({ $0 >= minimumDetectedFeatures }) ?? true,
              maximumMitochondrialFraction.map({ $0.isFinite && (0...1).contains($0) }) ?? true,
              excludedCells.count <= 100_000, Set(excludedCells).count == excludedCells.count else {
            throw VivoOmicsError.invalid("quality-filter bounds or duplicate excluded cells")
        }
        if let ids = includedSampleIDs {
            guard !ids.isEmpty, ids.count <= 100_000, Set(ids).count == ids.count, ids.allSatisfy(vivoOmicsID) else {
                throw VivoOmicsError.invalid("included sample IDs must be bounded, nonempty and unique")
            }
        }
    }
}
public enum VivoSingleCellRejectionReason: String, Codable, Sendable {
    case belowMinimumCounts, aboveMaximumCounts, belowMinimumDetectedFeatures, aboveMaximumDetectedFeatures
    case aboveMaximumMitochondrialFraction, missingMitochondrialFraction, sampleNotSelected, explicitlyExcluded
}
public struct VivoSingleCellDecision: Codable, Sendable, Equatable {
    public let sourceCellIndex: Int
    public let identity: VivoOmicsCellIdentity
    public let accepted: Bool
    public let reasons: [VivoSingleCellRejectionReason]
    public let quality: VivoCellQuality
    public let mitochondrialThresholdNotEvaluated: Bool
}
public struct VivoSingleCellFeatureStatistics: Codable, Sendable, Equatable {
    public let featureIndex: Int
    public let featureID: String
    public let totalCounts: UInt64
    public let detectedCells: Int
    public let meanLogNormalized: Double?
    public let varianceLogNormalized: Double?
}
public struct VivoSingleCellProcessed: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let policy: VivoSingleCellFilterPolicy
    public let decisions: [VivoSingleCellDecision]
    /// Selected row -> original row. All features are retained, preventing a
    /// feature-selection-dependent change to library-size denominators.
    public let sourceCellIndices: [Int]
    public let dataset: VivoSingleCellDataset
    public let normalized: VivoLogNormalizedCounts
    public let features: [VivoSingleCellFeatureStatistics]
    public let pseudobulk: VivoPseudobulkCounts
}
public enum VivoSingleCellProcessing {
    public static func run(_ source: VivoSingleCellDataset, policy: VivoSingleCellFilterPolicy = .init(),
                           normalizationTarget: Double = 10_000, limits: VivoOmicsLimits = .init()) throws -> VivoSingleCellProcessed {
        try source.validate(limits: limits); try policy.validate()
        let quality = try VivoSingleCellAnalysis.quality(source, limits: limits)
        let sampleIDs = Set(source.samples.map(\.id)), selectedSamples = policy.includedSampleIDs.map(Set.init)
        guard selectedSamples.map({ $0.isSubset(of: sampleIDs) }) ?? true else { throw VivoOmicsError.invalid("filter names an unknown sample") }
        let identities = source.cells.map { VivoOmicsCellIdentity(sampleID: $0.sampleID, barcode: $0.barcode) }
        let excluded = Set(policy.excludedCells)
        guard excluded.isSubset(of: Set(identities)) else { throw VivoOmicsError.invalid("filter excludes an unknown cell identity") }
        var decisions: [VivoSingleCellDecision] = [], selected: [Int] = []
        var offsets = [0], columns: [Int] = [], counts: [UInt64] = []
        for row in source.cells.indices {
            if row % 4096 == 0 { try Task.checkCancellation() }
            let qc = quality[row]
            var reasons: [VivoSingleCellRejectionReason] = []
            if qc.totalCounts < policy.minimumCounts { reasons.append(.belowMinimumCounts) }
            if let maximum = policy.maximumCounts, qc.totalCounts > maximum { reasons.append(.aboveMaximumCounts) }
            if qc.detectedFeatures < policy.minimumDetectedFeatures { reasons.append(.belowMinimumDetectedFeatures) }
            if let maximum = policy.maximumDetectedFeatures, qc.detectedFeatures > maximum { reasons.append(.aboveMaximumDetectedFeatures) }
            var missing = false
            if let maximum = policy.maximumMitochondrialFraction {
                if let fraction = qc.mitochondrialFraction {
                    if fraction > maximum { reasons.append(.aboveMaximumMitochondrialFraction) }
                } else {
                    missing = true
                    if policy.missingMitochondrial == .rejectCell { reasons.append(.missingMitochondrialFraction) }
                }
            }
            if let ids = selectedSamples, !ids.contains(source.cells[row].sampleID) { reasons.append(.sampleNotSelected) }
            if excluded.contains(identities[row]) { reasons.append(.explicitlyExcluded) }
            decisions.append(.init(sourceCellIndex: row, identity: identities[row], accepted: reasons.isEmpty,
                                   reasons: reasons, quality: qc, mitochondrialThresholdNotEvaluated: missing))
            if reasons.isEmpty {
                selected.append(row)
                let range = source.matrix.rowOffsets[row]..<source.matrix.rowOffsets[row + 1]
                columns.append(contentsOf: source.matrix.featureIndices[range]); counts.append(contentsOf: source.matrix.counts[range])
                offsets.append(counts.count)
            }
        }
        let dataset = VivoSingleCellDataset(id: source.id, evidence: source.evidence, sourceDescription: source.sourceDescription,
            countUnit: source.countUnit, samples: source.samples, features: source.features, cells: selected.map { source.cells[$0] },
            matrix: .init(cellCount: selected.count, featureCount: source.features.count, rowOffsets: offsets, featureIndices: columns, counts: counts))
        try dataset.validate(limits: limits)
        let normalized = try VivoSingleCellAnalysis.logNormalize(dataset, targetSum: normalizationTarget, limits: limits)
        let number = source.features.count, cells = selected.count
        var totals = [UInt64](repeating: 0, count: number), detected = [Int](repeating: 0, count: number)
        var means = [Double](repeating: 0, count: number), m2 = [Double](repeating: 0, count: number)
        for k in counts.indices {
            let j = columns[k], value = normalized.values[k]
            totals[j] = try vivoOmicsSum(totals[j], counts[k]); detected[j] += 1
            let delta = value - means[j]; means[j] += delta / Double(detected[j]); m2[j] += delta * (value - means[j])
        }
        let features = source.features.indices.map { j -> VivoSingleCellFeatureStatistics in
            // Combine the nonzero Welford accumulator with the implicit zeros.
            let mean = cells == 0 ? nil : means[j] * Double(detected[j]) / Double(cells)
            let variance: Double? = cells < 2 ? nil :
                (m2[j] + means[j] * means[j] * Double(detected[j]) * Double(cells - detected[j]) / Double(cells)) / Double(cells - 1)
            return .init(featureIndex: j, featureID: source.features[j].id, totalCounts: totals[j], detectedCells: detected[j],
                         meanLogNormalized: mean, varianceLogNormalized: variance)
        }
        return .init(schemaVersion: 1, policy: policy, decisions: decisions, sourceCellIndices: selected,
                     dataset: dataset, normalized: normalized, features: features,
                     pseudobulk: try VivoSingleCellAnalysis.pseudobulk(dataset, limits: limits))
    }
}
