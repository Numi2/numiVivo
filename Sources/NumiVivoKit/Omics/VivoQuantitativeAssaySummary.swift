import Foundation

/// Per-feature coverage and finite-value moments for a quantitative assay.
/// Missing values are absent sparse entries; an explicit zero contributes to
/// `measuredZeroCount` and to the moments.
public struct VivoQuantitativeFeatureSummary: Codable, Sendable, Equatable {
    public let assayID: String
    public let featureID: String
    public let assayRowCount: Int
    public let measuredValueCount: Int
    public let missingValueCount: Int
    public let measuredZeroCount: Int
    public let minimum: Double?
    public let maximum: Double?
    public let mean: Double?
    /// Sample variance when at least two values are measured. It is nil for
    /// zero or one measured value rather than treating missing uncertainty as
    /// zero.
    public let variance: Double?

    public init(assayID: String, featureID: String, assayRowCount: Int,
                measuredValueCount: Int, missingValueCount: Int,
                measuredZeroCount: Int, minimum: Double?, maximum: Double?,
                mean: Double?, variance: Double?) {
        self.assayID = assayID
        self.featureID = featureID
        self.assayRowCount = assayRowCount
        self.measuredValueCount = measuredValueCount
        self.missingValueCount = missingValueCount
        self.measuredZeroCount = measuredZeroCount
        self.minimum = minimum
        self.maximum = maximum
        self.mean = mean
        self.variance = variance
    }
}

/// Assay-level coverage totals. `observationCount` is the common dataset
/// axis, while `assayRowCount` records how many rows the modality supplied.
public struct VivoQuantitativeAssayCoverage: Codable, Sendable, Equatable {
    public let assayID: String
    public let kind: VivoQuantitativeAssayKind
    public let observationCount: Int
    public let assayRowCount: Int
    public let featureCount: Int
    public let measuredValueCount: Int
    public let missingValueCount: Int
    public let measuredZeroCount: Int

    public init(assayID: String, kind: VivoQuantitativeAssayKind,
                observationCount: Int, assayRowCount: Int, featureCount: Int,
                measuredValueCount: Int, missingValueCount: Int,
                measuredZeroCount: Int) {
        self.assayID = assayID
        self.kind = kind
        self.observationCount = observationCount
        self.assayRowCount = assayRowCount
        self.featureCount = featureCount
        self.measuredValueCount = measuredValueCount
        self.missingValueCount = missingValueCount
        self.measuredZeroCount = measuredZeroCount
    }
}

/// Deterministic, provenance-bound descriptive output for quantitative assays.
/// This report describes what was measured and how much is missing; it does
/// not normalize, impute, test, fit or predict a biological outcome.
public struct VivoQuantitativeAssaySummary: Codable, Sendable, Equatable {
    public let method: String
    public let datasetFingerprint: VivoFingerprint
    public let observationCount: Int
    public let assayCount: Int
    public let coverage: [VivoQuantitativeAssayCoverage]
    public let features: [VivoQuantitativeFeatureSummary]
    public let qualification: String

    public init(method: String = VivoQuantitativeAssaySummaries.method,
                datasetFingerprint: VivoFingerprint, observationCount: Int,
                assayCount: Int, coverage: [VivoQuantitativeAssayCoverage],
                features: [VivoQuantitativeFeatureSummary],
                qualification: String) {
        self.method = method
        self.datasetFingerprint = datasetFingerprint
        self.observationCount = observationCount
        self.assayCount = assayCount
        self.coverage = coverage
        self.features = features
        self.qualification = qualification
    }

    public func validate() throws {
        guard method == VivoQuantitativeAssaySummaries.method,
              observationCount >= 1, assayCount >= 1,
              coverage.count == assayCount, !features.isEmpty,
              !qualification.isEmpty, qualification.utf8.count <= 16_384,
              Set(coverage.map(\.assayID)).count == coverage.count,
              Set(features.map { "\($0.assayID)\u{1f}\($0.featureID)" }).count == features.count else {
            throw VivoOmicsError.invalid("quantitative assay summary schema or identity")
        }
        for item in coverage {
            guard vivoOmicsID(item.assayID), item.observationCount == observationCount,
                  item.assayRowCount >= 0, item.assayRowCount <= observationCount,
                  item.featureCount > 0, item.measuredValueCount >= 0,
                  item.measuredZeroCount >= 0,
                  item.measuredZeroCount <= item.measuredValueCount,
                  item.missingValueCount >= 0 else {
                throw VivoOmicsError.invalid("quantitative assay coverage")
            }
            let possible = item.assayRowCount.multipliedReportingOverflow(by: item.featureCount)
            guard !possible.overflow,
                  item.measuredValueCount <= possible.partialValue else {
                throw VivoOmicsError.invalid("quantitative assay coverage bounds")
            }
            let globalPossible = observationCount.multipliedReportingOverflow(by: item.featureCount)
            let total = item.measuredValueCount.addingReportingOverflow(item.missingValueCount)
            guard !globalPossible.overflow, !total.overflow,
                  total.partialValue == globalPossible.partialValue else {
                throw VivoOmicsError.invalid("quantitative assay coverage total")
            }
        }
        let coverageByAssay = Dictionary(uniqueKeysWithValues: coverage.map { ($0.assayID, $0) })
        var featuresByAssay: [String: Int] = [:]
        var featureIDsByAssay: [String: Set<String>] = [:]
        var measuredByAssay: [String: Int] = [:]
        var zerosByAssay: [String: Int] = [:]
        for item in features {
            let total = item.measuredValueCount.addingReportingOverflow(item.missingValueCount)
            guard vivoOmicsID(item.assayID), vivoOmicsID(item.featureID),
                  let assayCoverage = coverageByAssay[item.assayID],
                  item.assayRowCount == assayCoverage.assayRowCount,
                  item.assayRowCount >= 0, item.assayRowCount <= observationCount,
                  item.measuredValueCount >= 0, item.missingValueCount >= 0,
                  item.measuredZeroCount >= 0,
                  item.measuredZeroCount <= item.measuredValueCount,
                  !total.overflow, total.partialValue == observationCount,
                  item.measuredValueCount <= item.assayRowCount else {
                throw VivoOmicsError.invalid("quantitative feature summary counts")
            }
            if item.measuredValueCount == 0 {
                guard item.minimum == nil, item.maximum == nil, item.mean == nil, item.variance == nil else {
                    throw VivoOmicsError.invalid("empty quantitative feature summary moments")
                }
            } else {
                guard let minimum = item.minimum, let maximum = item.maximum,
                      let mean = item.mean, minimum.isFinite, maximum.isFinite,
                      mean.isFinite, minimum <= maximum,
                      item.variance.map({ $0.isFinite && $0 >= 0 }) ?? true,
                      item.measuredValueCount >= 2 || item.variance == nil else {
                    throw VivoOmicsError.invalid("quantitative feature summary moments")
                }
            }
            featuresByAssay[item.assayID, default: 0] += 1
            featureIDsByAssay[item.assayID, default: []].insert(item.featureID)
            let measuredTotal = measuredByAssay[item.assayID, default: 0].addingReportingOverflow(item.measuredValueCount)
            let zeroTotal = zerosByAssay[item.assayID, default: 0].addingReportingOverflow(item.measuredZeroCount)
            guard !measuredTotal.overflow, !zeroTotal.overflow else {
                throw VivoOmicsError.limit("quantitative summary feature totals")
            }
            measuredByAssay[item.assayID] = measuredTotal.partialValue
            zerosByAssay[item.assayID] = zeroTotal.partialValue
        }
        guard featuresByAssay.values.allSatisfy({ $0 > 0 }),
              Set(featuresByAssay.keys) == Set(coverageByAssay.keys),
              featuresByAssay.allSatisfy({ pair in pair.value == coverageByAssay[pair.key]?.featureCount }),
              featureIDsByAssay.values.allSatisfy({ $0.count > 0 }),
              coverage.allSatisfy({ item in
                  measuredByAssay[item.assayID] == item.measuredValueCount &&
                  zerosByAssay[item.assayID] == item.measuredZeroCount
              }) else {
            throw VivoOmicsError.invalid("quantitative summary feature coverage")
        }
    }
}

public enum VivoQuantitativeAssaySummaries {
    public static let method = "quantitative-descriptive-summary-v1"

    /// Computes bounded descriptive coverage in the dataset's declared assay
    /// spaces. Sparse row absence and explicit measured zeros remain distinct.
    public static func summarize(_ dataset: VivoQuantitativeAssayDataset) throws -> VivoQuantitativeAssaySummary {
        try dataset.validate()
        let datasetFingerprint = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(dataset))
        var coverage: [VivoQuantitativeAssayCoverage] = []
        var features: [VivoQuantitativeFeatureSummary] = []
        coverage.reserveCapacity(dataset.assays.count)
        for assay in dataset.assays {
            try Task.checkCancellation()
            var measuredTotal = 0
            var missingTotal = 0
            var zeroTotal = 0
            var measured = [Int](repeating: 0, count: assay.matrix.featureCount)
            var zeros = [Int](repeating: 0, count: assay.matrix.featureCount)
            var minima = [Double](repeating: .infinity, count: assay.matrix.featureCount)
            var maxima = [Double](repeating: -.infinity, count: assay.matrix.featureCount)
            var means = [Double](repeating: 0, count: assay.matrix.featureCount)
            var m2 = [Double](repeating: 0, count: assay.matrix.featureCount)
            for row in 0..<assay.matrix.observationCount {
                if row % 4_096 == 0 { try Task.checkCancellation() }
                let start = assay.matrix.rowOffsets[row], end = assay.matrix.rowOffsets[row + 1]
                measuredTotal += end - start
                for cursor in start..<end {
                    let feature = assay.matrix.featureIndices[cursor]
                    let value = assay.matrix.values[cursor]
                    measured[feature] += 1
                    if value == 0 {
                        zeros[feature] += 1
                        zeroTotal += 1
                    }
                    minima[feature] = min(minima[feature], value)
                    maxima[feature] = max(maxima[feature], value)
                    let count = measured[feature]
                    let delta = value - means[feature]
                    guard delta.isFinite else {
                        throw VivoOmicsError.invalid("quantitative summary moment overflow")
                    }
                    means[feature] += delta / Double(count)
                    let correction = value - means[feature]
                    m2[feature] += delta * correction
                    guard means[feature].isFinite, m2[feature].isFinite else {
                        throw VivoOmicsError.invalid("quantitative summary moment overflow")
                    }
                }
            }
            let possible = dataset.observations.count.multipliedReportingOverflow(by: assay.matrix.featureCount)
            guard !possible.overflow, measuredTotal <= possible.partialValue else {
                throw VivoOmicsError.limit("quantitative summary observation-feature bounds")
            }
            missingTotal = possible.partialValue - measuredTotal
            coverage.append(.init(assayID: assay.id, kind: assay.kind,
                                  observationCount: dataset.observations.count,
                                  assayRowCount: assay.matrix.observationCount,
                                  featureCount: assay.matrix.featureCount,
                                  measuredValueCount: measuredTotal,
                                  missingValueCount: missingTotal,
                                  measuredZeroCount: zeroTotal))
            for featureIndex in assay.features.indices {
                let count = measured[featureIndex]
                let variance: Double? = count >= 2 ? max(0, m2[featureIndex] / Double(count - 1)) : nil
                features.append(.init(assayID: assay.id, featureID: assay.features[featureIndex].id,
                                      assayRowCount: assay.matrix.observationCount,
                                      measuredValueCount: count,
                                      missingValueCount: dataset.observations.count - count,
                                      measuredZeroCount: zeros[featureIndex],
                                      minimum: count > 0 ? minima[featureIndex] : nil,
                                      maximum: count > 0 ? maxima[featureIndex] : nil,
                                      mean: count > 0 ? means[featureIndex] : nil,
                                      variance: variance))
            }
        }
        let result = VivoQuantitativeAssaySummary(
            datasetFingerprint: datasetFingerprint,
            observationCount: dataset.observations.count,
            assayCount: dataset.assays.count,
            coverage: coverage,
            features: features,
            qualification: "Descriptive quantitative-assay coverage and finite moments only. Missing values remain distinct from measured zeros; no normalization, imputation, differential test, joint model or biological outcome prediction is performed.")
        try result.validate()
        return result
    }
}

public extension VivoQuantitativeAssayDataset {
    func descriptiveSummary() throws -> VivoQuantitativeAssaySummary {
        try VivoQuantitativeAssaySummaries.summarize(self)
    }
}
