import Foundation

/// Continuous assay families that share the observation/provenance axis with
/// count assays.  They are intentionally separate from `VivoAssayKind`: raw
/// intensities and concentrations must not be relabeled as molecule counts.
public enum VivoQuantitativeAssayKind: String, Codable, Sendable {
    case proteomics
    case metabolomics
    case spatialImaging
}

/// Canonical sparse storage for finite real-valued measurements.  An absent
/// entry means unmeasured/missing; an explicit zero is retained as a measured
/// zero.  Values may be negative when the declared unit is a transformed scale.
public struct VivoSparseValues: Codable, Sendable, Equatable {
    public let observationCount: Int
    public let featureCount: Int
    public let rowOffsets: [Int]
    public let featureIndices: [Int]
    public let values: [Double]

    public init(observationCount: Int, featureCount: Int, rowOffsets: [Int],
                featureIndices: [Int], values: [Double]) {
        self.observationCount = observationCount
        self.featureCount = featureCount
        self.rowOffsets = rowOffsets
        self.featureIndices = featureIndices
        self.values = values
    }

    public func validate(limits: VivoOmicsLimits = .init()) throws {
        try limits.validate()
        guard observationCount >= 0, observationCount <= limits.maximumCells,
              featureCount > 0, featureCount <= limits.maximumFeatures,
              values.count <= limits.maximumNonzeros,
              rowOffsets.count == observationCount + 1,
              rowOffsets.first == 0, rowOffsets.last == values.count,
              featureIndices.count == values.count else {
            throw VivoOmicsError.limit("sparse real matrix shape or nonzero count")
        }
        for row in 0..<observationCount {
            let start = rowOffsets[row], end = rowOffsets[row + 1]
            guard start >= 0, end >= start, end <= values.count else {
                throw VivoOmicsError.invalid("sparse real row offsets")
            }
            var previous = -1
            for index in start..<end {
                let feature = featureIndices[index]
                guard feature > previous, feature < featureCount,
                      values[index].isFinite else {
                    throw VivoOmicsError.invalid("sparse real indices must be sorted, unique and finite")
                }
                previous = feature
            }
        }
    }
}

/// A feature space for proteomics, metabolomics or spatial-image measurements.
/// The row map preserves partial assay coverage without padding missing rows.
public struct VivoQuantitativeAssaySpace: Codable, Sendable, Equatable {
    public let id: String
    public let kind: VivoQuantitativeAssayKind
    public let featureNamespace: String
    public let unit: String
    public let sourceDescription: String
    public let features: [VivoAssayFeature]
    public let observationIndices: [Int]
    public let matrix: VivoSparseValues

    public init(id: String, kind: VivoQuantitativeAssayKind, featureNamespace: String,
                unit: String, sourceDescription: String, features: [VivoAssayFeature],
                observationIndices: [Int], matrix: VivoSparseValues) {
        self.id = id; self.kind = kind; self.featureNamespace = featureNamespace
        self.unit = unit; self.sourceDescription = sourceDescription; self.features = features
        self.observationIndices = observationIndices; self.matrix = matrix
    }
}

/// Shared identity/provenance container for non-count omics measurements.
/// This is an interchange foundation only: no imputation, normalization,
/// differential testing or biological outcome inference is performed here.
public struct VivoQuantitativeAssayDataset: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let id: String
    public let evidence: VivoOmicsEvidence
    public let sourceDescription: String
    public let samples: [VivoOmicsSample]
    public let observations: [VivoAssayObservation]
    public let spatialFrames: [VivoSpatialFrame]
    public let assays: [VivoQuantitativeAssaySpace]

    public init(schemaVersion: Int = 1, id: String, evidence: VivoOmicsEvidence,
                sourceDescription: String, samples: [VivoOmicsSample],
                observations: [VivoAssayObservation], spatialFrames: [VivoSpatialFrame],
                assays: [VivoQuantitativeAssaySpace]) {
        self.schemaVersion = schemaVersion; self.id = id; self.evidence = evidence
        self.sourceDescription = sourceDescription; self.samples = samples
        self.observations = observations; self.spatialFrames = spatialFrames; self.assays = assays
    }

    public static var limits: VivoOmicsLimits {
        var result = VivoOmicsLimits()
        result.maximumFeatures = 200_000
        result.maximumNonzeros = 32_000_000
        result.maximumInputBytes = 1_073_741_824
        return result
    }

    public func validate() throws {
        let limits = Self.limits
        guard schemaVersion == 1, vivoOmicsID(id), !sourceDescription.isEmpty,
              sourceDescription.utf8.count <= 16_384, !samples.isEmpty,
              samples.count <= limits.maximumCells, !observations.isEmpty,
              observations.count <= limits.maximumCells, !assays.isEmpty,
              assays.count <= 16, spatialFrames.count <= 64 else {
            throw VivoOmicsError.invalid("quantitative assay schema, identity or bounds")
        }
        try vivoOmicsValidateIdentities(samples: samples, cells: observations.map(\.identity))

        var frames: [String: VivoSpatialFrame] = [:]
        for frame in spatialFrames {
            guard vivoOmicsID(frame.id), frames.updateValue(frame, forKey: frame.id) == nil,
                  (2...3).contains(frame.axes.count), Set(frame.axes).count == frame.axes.count,
                  frame.axes.allSatisfy(vivoOmicsID), !frame.sourceDescription.isEmpty,
                  frame.sourceDescription.utf8.count <= 16_384 else {
                throw VivoOmicsError.invalid("quantitative spatial frame")
            }
        }
        for observation in observations {
            if let position = observation.position {
                guard let frame = frames[position.frameID],
                      position.coordinates.count == frame.axes.count,
                      position.coordinates.allSatisfy(\.isFinite) else {
                    throw VivoOmicsError.invalid("quantitative spatial position/frame mismatch")
                }
            }
        }

        var assayIDs = Set<String>()
        var remainingFeatures = limits.maximumFeatures
        var remainingValues = limits.maximumNonzeros
        for assay in assays {
            guard vivoOmicsID(assay.id), assayIDs.insert(assay.id).inserted,
                  vivoOmicsID(assay.featureNamespace), vivoOmicsID(assay.unit),
                  !assay.sourceDescription.isEmpty,
                  assay.sourceDescription.utf8.count <= 16_384,
                  assay.observationIndices.count == assay.matrix.observationCount,
                  assay.features.count == assay.matrix.featureCount,
                  Set(assay.observationIndices).count == assay.observationIndices.count,
                  assay.observationIndices.allSatisfy({ observations.indices.contains($0) }),
                  assay.features.count <= remainingFeatures,
                  assay.matrix.values.count <= remainingValues else {
                throw VivoOmicsError.invalid("quantitative assay identity, row map or bounds")
            }
            remainingFeatures -= assay.features.count
            remainingValues -= assay.matrix.values.count
            try assay.matrix.validate(limits: limits)

            var featureIDs = Set<String>()
            for feature in assay.features {
                guard vivoOmicsID(feature.id), vivoOmicsID(feature.name),
                      featureIDs.insert(feature.id).inserted else {
                    throw VivoOmicsError.invalid("quantitative feature identity")
                }
                if let interval = feature.interval {
                    guard vivoOmicsID(interval.contig), interval.start < interval.end else {
                        throw VivoOmicsError.invalid("quantitative genomic interval")
                    }
                }
            }
        }
    }

    /// Returns nil for both an unmeasured row and a missing sparse value.  A
    /// measured zero is returned as `0`, so callers never need to infer missing
    /// values from a numeric sentinel.
    public func value(assayID: String, observationIndex: Int, featureID: String) throws -> Double? {
        try validate()
        guard observations.indices.contains(observationIndex),
              let assay = assays.first(where: { $0.id == assayID }),
              let feature = assay.features.firstIndex(where: { $0.id == featureID }) else {
            throw VivoOmicsError.invalid("unknown quantitative assay, observation or feature")
        }
        guard let row = assay.observationIndices.firstIndex(of: observationIndex) else { return nil }
        for index in assay.matrix.rowOffsets[row]..<assay.matrix.rowOffsets[row + 1]
            where assay.matrix.featureIndices[index] == feature {
            return assay.matrix.values[index]
        }
        return nil
    }
}
