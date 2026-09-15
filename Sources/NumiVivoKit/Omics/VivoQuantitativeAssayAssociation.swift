import Foundation

/// Identifies two declared quantitative features for a pairwise complete-case
/// association.  The observation axis is inherited from the dataset, so an
/// absent assay row or sparse entry remains missing rather than becoming zero.
public struct VivoQuantitativeAssayFeaturePair: Codable, Sendable, Equatable, Hashable {
    public let leftAssayID: String
    public let leftFeatureID: String
    public let rightAssayID: String
    public let rightFeatureID: String

    public init(leftAssayID: String, leftFeatureID: String,
                rightAssayID: String, rightFeatureID: String) {
        self.leftAssayID = leftAssayID
        self.leftFeatureID = leftFeatureID
        self.rightAssayID = rightAssayID
        self.rightFeatureID = rightFeatureID
    }

    private enum CodingKeys: String, CodingKey {
        case leftAssayID, leftFeatureID, rightAssayID, rightFeatureID
    }

    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: [
            "leftAssayID", "leftFeatureID", "rightAssayID", "rightFeatureID"
        ])
        let values = try decoder.container(keyedBy: CodingKeys.self)
        leftAssayID = try values.decode(String.self, forKey: .leftAssayID)
        leftFeatureID = try values.decode(String.self, forKey: .leftFeatureID)
        rightAssayID = try values.decode(String.self, forKey: .rightAssayID)
        rightFeatureID = try values.decode(String.self, forKey: .rightFeatureID)
    }

    public func validate() throws {
        guard [leftAssayID, leftFeatureID, rightAssayID, rightFeatureID].allSatisfy(vivoOmicsID),
              leftAssayID != rightAssayID || leftFeatureID != rightFeatureID else {
            throw VivoOmicsError.invalid("quantitative association feature pair")
        }
    }
}

public enum VivoQuantitativeAssayAssociationStatus: String, Codable, Sendable {
    case computed
    case insufficientOverlap
    case zeroVariance
}

/// Pairwise complete-case moments for two declared quantitative features.
/// Variances and covariance require at least two complete observations;
/// Pearson correlation is reported only when both features vary.  This is a
/// descriptive association, not a causal or outcome model.
public struct VivoQuantitativeAssayAssociation: Codable, Sendable, Equatable {
    public let leftAssayID: String
    public let leftFeatureID: String
    public let rightAssayID: String
    public let rightFeatureID: String
    public let overlapCount: Int
    public let missingPairCount: Int
    public let leftMean: Double?
    public let rightMean: Double?
    public let leftVariance: Double?
    public let rightVariance: Double?
    public let covariance: Double?
    public let pearsonCorrelation: Double?
    public let status: VivoQuantitativeAssayAssociationStatus

    public init(pair: VivoQuantitativeAssayFeaturePair, overlapCount: Int,
                missingPairCount: Int, leftMean: Double?, rightMean: Double?,
                leftVariance: Double?, rightVariance: Double?, covariance: Double?,
                pearsonCorrelation: Double?, status: VivoQuantitativeAssayAssociationStatus) {
        self.leftAssayID = pair.leftAssayID
        self.leftFeatureID = pair.leftFeatureID
        self.rightAssayID = pair.rightAssayID
        self.rightFeatureID = pair.rightFeatureID
        self.overlapCount = overlapCount
        self.missingPairCount = missingPairCount
        self.leftMean = leftMean
        self.rightMean = rightMean
        self.leftVariance = leftVariance
        self.rightVariance = rightVariance
        self.covariance = covariance
        self.pearsonCorrelation = pearsonCorrelation
        self.status = status
    }

    private enum CodingKeys: String, CodingKey {
        case leftAssayID, leftFeatureID, rightAssayID, rightFeatureID
        case overlapCount, missingPairCount
        case leftMean, rightMean, leftVariance, rightVariance, covariance
        case pearsonCorrelation, status
    }

    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: [
            "leftAssayID", "leftFeatureID", "rightAssayID", "rightFeatureID",
            "overlapCount", "missingPairCount", "leftMean", "rightMean",
            "leftVariance", "rightVariance", "covariance", "pearsonCorrelation", "status"
        ])
        let values = try decoder.container(keyedBy: CodingKeys.self)
        leftAssayID = try values.decode(String.self, forKey: .leftAssayID)
        leftFeatureID = try values.decode(String.self, forKey: .leftFeatureID)
        rightAssayID = try values.decode(String.self, forKey: .rightAssayID)
        rightFeatureID = try values.decode(String.self, forKey: .rightFeatureID)
        overlapCount = try values.decode(Int.self, forKey: .overlapCount)
        missingPairCount = try values.decode(Int.self, forKey: .missingPairCount)
        leftMean = try values.decodeIfPresent(Double.self, forKey: .leftMean)
        rightMean = try values.decodeIfPresent(Double.self, forKey: .rightMean)
        leftVariance = try values.decodeIfPresent(Double.self, forKey: .leftVariance)
        rightVariance = try values.decodeIfPresent(Double.self, forKey: .rightVariance)
        covariance = try values.decodeIfPresent(Double.self, forKey: .covariance)
        pearsonCorrelation = try values.decodeIfPresent(Double.self, forKey: .pearsonCorrelation)
        status = try values.decode(VivoQuantitativeAssayAssociationStatus.self, forKey: .status)
    }

    private var pair: VivoQuantitativeAssayFeaturePair {
        .init(leftAssayID: leftAssayID, leftFeatureID: leftFeatureID,
              rightAssayID: rightAssayID, rightFeatureID: rightFeatureID)
    }

    fileprivate func validate(observationCount: Int) throws {
        try pair.validate()
        guard observationCount >= 1, overlapCount >= 0, missingPairCount >= 0,
              overlapCount <= observationCount, missingPairCount <= observationCount else {
            throw VivoOmicsError.invalid("quantitative association counts")
        }
        let total = overlapCount.addingReportingOverflow(missingPairCount)
        guard !total.overflow, total.partialValue == observationCount else {
            throw VivoOmicsError.invalid("quantitative association coverage")
        }

        func finite(_ value: Double?) -> Bool { value.map(\.isFinite) ?? true }
        guard finite(leftMean), finite(rightMean), finite(leftVariance), finite(rightVariance),
              finite(covariance), finite(pearsonCorrelation) else {
            throw VivoOmicsError.invalid("nonfinite quantitative association moment")
        }

        switch status {
        case .insufficientOverlap:
            guard overlapCount < 2,
                  leftVariance == nil, rightVariance == nil,
                  covariance == nil, pearsonCorrelation == nil else {
                throw VivoOmicsError.invalid("insufficient quantitative association moments")
            }
            if overlapCount == 0 {
                guard leftMean == nil, rightMean == nil else {
                    throw VivoOmicsError.invalid("empty quantitative association moments")
                }
            } else {
                guard leftMean != nil, rightMean != nil else {
                    throw VivoOmicsError.invalid("single-observation association means")
                }
            }
        case .zeroVariance:
            guard overlapCount >= 2, let leftVariance, let rightVariance, covariance != nil,
                  leftVariance >= 0, rightVariance >= 0, (leftVariance == 0 || rightVariance == 0),
                  leftMean != nil, rightMean != nil, pearsonCorrelation == nil else {
                throw VivoOmicsError.invalid("zero-variance quantitative association")
            }
        case .computed:
            guard overlapCount >= 2, let leftVariance, let rightVariance, let covariance,
                  let correlation = pearsonCorrelation, leftVariance > 0, rightVariance > 0,
                  correlation >= -1, correlation <= 1,
                  leftMean != nil, rightMean != nil else {
                throw VivoOmicsError.invalid("computed quantitative association moments")
            }
            let scale = max(leftVariance, rightVariance)
            let bound = scale * sqrt((leftVariance / scale) * (rightVariance / scale))
            guard bound.isFinite, abs(covariance) <= bound * (1 + 1e-9) else {
                throw VivoOmicsError.invalid("quantitative association covariance bound")
            }
        }
    }
}

/// Fingerprint-bound pairwise association report over quantitative assays.
/// Complete observations are selected independently for each requested pair;
/// no missing-value imputation or assay-specific normalization is performed.
public struct VivoQuantitativeAssayAssociationSummary: Codable, Sendable, Equatable {
    public let method: String
    public let datasetFingerprint: VivoFingerprint
    public let observationCount: Int
    public let associations: [VivoQuantitativeAssayAssociation]
    public let qualification: String

    public init(method: String = VivoQuantitativeAssayAssociations.method,
                datasetFingerprint: VivoFingerprint, observationCount: Int,
                associations: [VivoQuantitativeAssayAssociation], qualification: String) {
        self.method = method
        self.datasetFingerprint = datasetFingerprint
        self.observationCount = observationCount
        self.associations = associations
        self.qualification = qualification
    }

    private enum CodingKeys: String, CodingKey {
        case method, datasetFingerprint, observationCount, associations, qualification
    }

    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: [
            "method", "datasetFingerprint", "observationCount", "associations", "qualification"
        ])
        let values = try decoder.container(keyedBy: CodingKeys.self)
        method = try values.decode(String.self, forKey: .method)
        datasetFingerprint = try values.decode(VivoFingerprint.self, forKey: .datasetFingerprint)
        observationCount = try values.decode(Int.self, forKey: .observationCount)
        associations = try values.decode([VivoQuantitativeAssayAssociation].self, forKey: .associations)
        qualification = try values.decode(String.self, forKey: .qualification)
    }

    public func validate() throws {
        guard method == VivoQuantitativeAssayAssociations.method,
              observationCount >= 1, !associations.isEmpty,
              associations.count <= VivoQuantitativeAssayAssociations.maximumPairs,
              !qualification.isEmpty, qualification.utf8.count <= 16_384 else {
            throw VivoOmicsError.invalid("quantitative association summary schema")
        }
        var seen = Set<String>()
        for association in associations {
            try association.validate(observationCount: observationCount)
            let left = association.leftAssayID + "\u{1f}" + association.leftFeatureID
            let right = association.rightAssayID + "\u{1f}" + association.rightFeatureID
            let key = left < right ? left + "\u{1e}" + right : right + "\u{1e}" + left
            guard seen.insert(key).inserted else {
                throw VivoOmicsError.invalid("duplicate quantitative association pair")
            }
        }
    }
}

public enum VivoQuantitativeAssayAssociations {
    public static let method = "quantitative-pairwise-association-v1"
    public static let maximumPairs = 4_096
    private static let maximumPairOperations = 50_000_000

    private struct IndexedAssay {
        let assay: VivoQuantitativeAssaySpace
        let rowsByObservation: [Int: Int]
        let featuresByID: [String: Int]
    }

    private static func sparseValue(_ matrix: VivoSparseValues, row: Int, feature: Int) -> Double? {
        var low = matrix.rowOffsets[row]
        var high = matrix.rowOffsets[row + 1]
        while low < high {
            let middle = low + (high - low) / 2
            let candidate = matrix.featureIndices[middle]
            if candidate < feature {
                low = middle + 1
            } else {
                high = middle
            }
        }
        guard low < matrix.rowOffsets[row + 1], matrix.featureIndices[low] == feature else { return nil }
        return matrix.values[low]
    }

    public static func summarize(_ dataset: VivoQuantitativeAssayDataset,
                                 pairs: [VivoQuantitativeAssayFeaturePair]) throws -> VivoQuantitativeAssayAssociationSummary {
        try dataset.validate()
        guard !pairs.isEmpty, pairs.count <= maximumPairs else {
            throw VivoOmicsError.limit("quantitative association pair count")
        }
        let operations = pairs.count.multipliedReportingOverflow(by: dataset.observations.count)
        guard !operations.overflow, operations.partialValue <= maximumPairOperations else {
            throw VivoOmicsError.limit("quantitative association observation-pair operations")
        }

        var indexed: [String: IndexedAssay] = [:]
        indexed.reserveCapacity(dataset.assays.count)
        for assay in dataset.assays {
            var rows: [Int: Int] = [:]
            rows.reserveCapacity(assay.observationIndices.count)
            for (row, observation) in assay.observationIndices.enumerated() {
                rows[observation] = row
            }
            let features = Dictionary(uniqueKeysWithValues: assay.features.enumerated().map { ($0.element.id, $0.offset) })
            indexed[assay.id] = IndexedAssay(assay: assay, rowsByObservation: rows, featuresByID: features)
        }

        var associations: [VivoQuantitativeAssayAssociation] = []
        associations.reserveCapacity(pairs.count)
        for pair in pairs {
            try Task.checkCancellation()
            try pair.validate()
            guard let left = indexed[pair.leftAssayID], let right = indexed[pair.rightAssayID],
                  let leftFeature = left.featuresByID[pair.leftFeatureID],
                  let rightFeature = right.featuresByID[pair.rightFeatureID] else {
                throw VivoOmicsError.invalid("quantitative association references unknown assay or feature")
            }

            var overlap = 0
            var leftMean = 0.0, rightMean = 0.0
            var leftM2 = 0.0, rightM2 = 0.0, coM2 = 0.0
            for observation in dataset.observations.indices {
                if observation % 4_096 == 0 { try Task.checkCancellation() }
                guard let leftRow = left.rowsByObservation[observation],
                      let rightRow = right.rowsByObservation[observation],
                      let x = sparseValue(left.assay.matrix, row: leftRow, feature: leftFeature),
                      let y = sparseValue(right.assay.matrix, row: rightRow, feature: rightFeature) else {
                    continue
                }
                overlap += 1
                let deltaX = x - leftMean
                let deltaY = y - rightMean
                guard deltaX.isFinite, deltaY.isFinite else {
                    throw VivoOmicsError.invalid("quantitative association moment overflow")
                }
                leftMean += deltaX / Double(overlap)
                rightMean += deltaY / Double(overlap)
                leftM2 += deltaX * (x - leftMean)
                rightM2 += deltaY * (y - rightMean)
                coM2 += deltaX * (y - rightMean)
                guard leftMean.isFinite, rightMean.isFinite, leftM2.isFinite,
                      rightM2.isFinite, coM2.isFinite else {
                    throw VivoOmicsError.invalid("quantitative association moment overflow")
                }
            }

            let missing = dataset.observations.count - overlap
            let leftMeanValue = overlap > 0 ? leftMean : nil
            let rightMeanValue = overlap > 0 ? rightMean : nil
            var leftVariance: Double?
            var rightVariance: Double?
            var covariance: Double?
            var correlation: Double?
            var status: VivoQuantitativeAssayAssociationStatus
            if overlap < 2 {
                leftVariance = nil
                rightVariance = nil
                covariance = nil
                correlation = nil
                status = .insufficientOverlap
            } else {
                leftVariance = max(0, leftM2 / Double(overlap - 1))
                rightVariance = max(0, rightM2 / Double(overlap - 1))
                covariance = coM2 / Double(overlap - 1)
                guard leftVariance!.isFinite, rightVariance!.isFinite, covariance!.isFinite else {
                    throw VivoOmicsError.invalid("quantitative association covariance overflow")
                }
                if leftVariance == 0 || rightVariance == 0 {
                    correlation = nil
                    status = .zeroVariance
                } else {
                    let scale = max(leftVariance!, rightVariance!)
                    let denominator = scale * sqrt((leftVariance! / scale) * (rightVariance! / scale))
                    guard denominator.isFinite, denominator > 0 else {
                        throw VivoOmicsError.invalid("quantitative association correlation overflow")
                    }
                    let rawCorrelation = covariance! / denominator
                    guard rawCorrelation.isFinite, abs(rawCorrelation) <= 1 + 1e-9 else {
                        throw VivoOmicsError.invalid("quantitative association correlation bound")
                    }
                    correlation = min(1, max(-1, rawCorrelation))
                    status = .computed
                }
            }
            associations.append(.init(pair: pair, overlapCount: overlap, missingPairCount: missing,
                                      leftMean: leftMeanValue, rightMean: rightMeanValue,
                                      leftVariance: leftVariance, rightVariance: rightVariance,
                                      covariance: covariance, pearsonCorrelation: correlation,
                                      status: status))
        }

        let result = VivoQuantitativeAssayAssociationSummary(
            datasetFingerprint: try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(dataset)),
            observationCount: dataset.observations.count, associations: associations,
            qualification: "Pairwise complete-observation moments and Pearson association over declared quantitative assays only. Missing rows and sparse entries are excluded without imputation; measured zeros are included. No normalization, assay-specific transform, causal model, differential test or biological outcome prediction is performed.")
        try result.validate()
        return result
    }
}

public extension VivoQuantitativeAssayDataset {
    func pairwiseAssociationSummary(_ pairs: [VivoQuantitativeAssayFeaturePair]) throws -> VivoQuantitativeAssayAssociationSummary {
        try VivoQuantitativeAssayAssociations.summarize(self, pairs: pairs)
    }
}
