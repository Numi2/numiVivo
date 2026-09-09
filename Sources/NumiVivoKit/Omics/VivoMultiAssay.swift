import Foundation

/// Independent feature spaces share observation identities without padding missing
/// measurements with biological zeros. Count matrices retain exact UInt64 values.
public enum VivoAssayKind: String, Codable, Sendable { case rna, antibodyCapture, chromatinAccessibility, guideCapture }
public enum VivoAssayCountUnit: String, Codable, Sendable { case umiCount, readCount, fragmentCount }
public enum VivoObservationKind: String, Codable, Sendable { case cell, spot }
public enum VivoSpatialUnit: String, Codable, Sendable { case micrometer, pixel }

public struct VivoGenomicInterval: Codable, Sendable, Equatable {
    public let contig: String
    public let start: UInt64
    public let end: UInt64
    public init(contig: String, start: UInt64, end: UInt64) {
        self.contig = contig
        self.start = start
        self.end = end
    }
    private enum CodingKeys: String, CodingKey { case contig, start, end }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["contig", "start", "end"])
        let v = try decoder.container(keyedBy: CodingKeys.self)
        contig = try v.decode(String.self, forKey: .contig)
        start = try v.decode(UInt64.self, forKey: .start)
        end = try v.decode(UInt64.self, forKey: .end)
    }
}

public struct VivoAssayFeature: Codable, Sendable, Equatable {
    public let id: String
    public let name: String
    public let interval: VivoGenomicInterval?
    public init(id: String, name: String, interval: VivoGenomicInterval?) {
        self.id = id
        self.name = name
        self.interval = interval
    }
    private enum CodingKeys: String, CodingKey { case id, name, interval }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["id", "name", "interval"])
        let v = try decoder.container(keyedBy: CodingKeys.self)
        id = try v.decode(String.self, forKey: .id)
        name = try v.decode(String.self, forKey: .name)
        interval = try v.decodeIfPresent(VivoGenomicInterval.self, forKey: .interval)
    }
}

public struct VivoSpatialFrame: Codable, Sendable, Equatable {
    public let id: String
    public let unit: VivoSpatialUnit
    public let axes: [String]
    public let sourceDescription: String
    public init(id: String, unit: VivoSpatialUnit, axes: [String], sourceDescription: String) {
        self.id = id
        self.unit = unit
        self.axes = axes
        self.sourceDescription = sourceDescription
    }
    private enum CodingKeys: String, CodingKey { case id, unit, axes, sourceDescription }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["id", "unit", "axes", "sourceDescription"])
        let v = try decoder.container(keyedBy: CodingKeys.self)
        id = try v.decode(String.self, forKey: .id)
        unit = try v.decode(VivoSpatialUnit.self, forKey: .unit)
        axes = try v.decode([String].self, forKey: .axes)
        sourceDescription = try v.decode(String.self, forKey: .sourceDescription)
    }
}

public struct VivoSpatialPosition: Codable, Sendable, Equatable {
    public let frameID: String
    public let coordinates: [Double]
    public init(frameID: String, coordinates: [Double]) {
        self.frameID = frameID
        self.coordinates = coordinates
    }
    private enum CodingKeys: String, CodingKey { case frameID, coordinates }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["frameID", "coordinates"])
        let v = try decoder.container(keyedBy: CodingKeys.self)
        frameID = try v.decode(String.self, forKey: .frameID)
        coordinates = try v.decode([Double].self, forKey: .coordinates)
    }
}

public struct VivoAssayObservation: Codable, Sendable, Equatable {
    public let identity: VivoOmicsCell
    public let kind: VivoObservationKind
    public let position: VivoSpatialPosition?
    public init(identity: VivoOmicsCell, kind: VivoObservationKind, position: VivoSpatialPosition?) {
        self.identity = identity
        self.kind = kind
        self.position = position
    }
    private enum CodingKeys: String, CodingKey { case identity, kind, position }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["identity", "kind", "position"])
        let v = try decoder.container(keyedBy: CodingKeys.self)
        identity = try v.decode(VivoOmicsCell.self, forKey: .identity)
        kind = try v.decode(VivoObservationKind.self, forKey: .kind)
        position = try v.decodeIfPresent(VivoSpatialPosition.self, forKey: .position)
    }
}

public struct VivoAssaySpace: Codable, Sendable, Equatable {
    public let id: String
    public let kind: VivoAssayKind
    public let featureNamespace: String
    public let countUnit: VivoAssayCountUnit
    public let genomeAssembly: String?
    public let sourceDescription: String
    public let features: [VivoAssayFeature]
    public let observationIndices: [Int]
    public let matrix: VivoSparseCounts
    public init(id: String, kind: VivoAssayKind, featureNamespace: String, countUnit: VivoAssayCountUnit, genomeAssembly: String?, sourceDescription: String, features: [VivoAssayFeature], observationIndices: [Int], matrix: VivoSparseCounts) {
        self.id = id
        self.kind = kind
        self.featureNamespace = featureNamespace
        self.countUnit = countUnit
        self.genomeAssembly = genomeAssembly
        self.sourceDescription = sourceDescription
        self.features = features
        self.observationIndices = observationIndices
        self.matrix = matrix
    }
    private enum CodingKeys: String, CodingKey { case id, kind, featureNamespace, countUnit, genomeAssembly, sourceDescription, features, observationIndices, matrix }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["id", "kind", "featureNamespace", "countUnit", "genomeAssembly", "sourceDescription", "features", "observationIndices", "matrix"])
        let v = try decoder.container(keyedBy: CodingKeys.self)
        id = try v.decode(String.self, forKey: .id)
        kind = try v.decode(VivoAssayKind.self, forKey: .kind)
        featureNamespace = try v.decode(String.self, forKey: .featureNamespace)
        countUnit = try v.decode(VivoAssayCountUnit.self, forKey: .countUnit)
        genomeAssembly = try v.decodeIfPresent(String.self, forKey: .genomeAssembly)
        sourceDescription = try v.decode(String.self, forKey: .sourceDescription)
        features = try v.decode([VivoAssayFeature].self, forKey: .features)
        observationIndices = try v.decode([Int].self, forKey: .observationIndices)
        matrix = try v.decode(VivoSparseCounts.self, forKey: .matrix)
    }
}

public struct VivoMultiAssayDataset: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let id: String
    public let evidence: VivoOmicsEvidence
    public let sourceDescription: String
    public let samples: [VivoOmicsSample]
    public let observations: [VivoAssayObservation]
    public let spatialFrames: [VivoSpatialFrame]
    public let assays: [VivoAssaySpace]
    public init(schemaVersion: Int, id: String, evidence: VivoOmicsEvidence, sourceDescription: String, samples: [VivoOmicsSample], observations: [VivoAssayObservation], spatialFrames: [VivoSpatialFrame], assays: [VivoAssaySpace]) {
        self.schemaVersion = schemaVersion
        self.id = id
        self.evidence = evidence
        self.sourceDescription = sourceDescription
        self.samples = samples
        self.observations = observations
        self.spatialFrames = spatialFrames
        self.assays = assays
    }
    private enum CodingKeys: String, CodingKey { case schemaVersion, id, evidence, sourceDescription, samples, observations, spatialFrames, assays }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["schemaVersion", "id", "evidence", "sourceDescription", "samples", "observations", "spatialFrames", "assays"])
        let v = try decoder.container(keyedBy: CodingKeys.self)
        schemaVersion = try v.decode(Int.self, forKey: .schemaVersion)
        id = try v.decode(String.self, forKey: .id)
        evidence = try v.decode(VivoOmicsEvidence.self, forKey: .evidence)
        sourceDescription = try v.decode(String.self, forKey: .sourceDescription)
        samples = try v.decode([VivoOmicsSample].self, forKey: .samples)
        observations = try v.decode([VivoAssayObservation].self, forKey: .observations)
        spatialFrames = try v.decode([VivoSpatialFrame].self, forKey: .spatialFrames)
        assays = try v.decode([VivoAssaySpace].self, forKey: .assays)
    }
}

public struct VivoTenXAssayMapping: Codable, Sendable, Equatable {
    public let featureType: String
    public let id: String
    public let kind: VivoAssayKind
    public let featureNamespace: String
    public let countUnit: VivoAssayCountUnit
    public let genomeAssembly: String?
    public init(featureType: String, id: String, kind: VivoAssayKind, featureNamespace: String, countUnit: VivoAssayCountUnit, genomeAssembly: String?) {
        self.featureType = featureType
        self.id = id
        self.kind = kind
        self.featureNamespace = featureNamespace
        self.countUnit = countUnit
        self.genomeAssembly = genomeAssembly
    }
    private enum CodingKeys: String, CodingKey { case featureType, id, kind, featureNamespace, countUnit, genomeAssembly }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["featureType", "id", "kind", "featureNamespace", "countUnit", "genomeAssembly"])
        let v = try decoder.container(keyedBy: CodingKeys.self)
        featureType = try v.decode(String.self, forKey: .featureType)
        id = try v.decode(String.self, forKey: .id)
        kind = try v.decode(VivoAssayKind.self, forKey: .kind)
        featureNamespace = try v.decode(String.self, forKey: .featureNamespace)
        countUnit = try v.decode(VivoAssayCountUnit.self, forKey: .countUnit)
        genomeAssembly = try v.decodeIfPresent(String.self, forKey: .genomeAssembly)
    }
}

public struct VivoTenXMultiAssayPlan: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let id: String
    public let evidence: VivoOmicsEvidence
    public let sourceDescription: String
    public let sample: VivoOmicsSample
    public let assays: [VivoTenXAssayMapping]
    public init(schemaVersion: Int, id: String, evidence: VivoOmicsEvidence, sourceDescription: String, sample: VivoOmicsSample, assays: [VivoTenXAssayMapping]) {
        self.schemaVersion = schemaVersion
        self.id = id
        self.evidence = evidence
        self.sourceDescription = sourceDescription
        self.sample = sample
        self.assays = assays
    }
    private enum CodingKeys: String, CodingKey { case schemaVersion, id, evidence, sourceDescription, sample, assays }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["schemaVersion", "id", "evidence", "sourceDescription", "sample", "assays"])
        let v = try decoder.container(keyedBy: CodingKeys.self)
        schemaVersion = try v.decode(Int.self, forKey: .schemaVersion)
        id = try v.decode(String.self, forKey: .id)
        evidence = try v.decode(VivoOmicsEvidence.self, forKey: .evidence)
        sourceDescription = try v.decode(String.self, forKey: .sourceDescription)
        sample = try v.decode(VivoOmicsSample.self, forKey: .sample)
        assays = try v.decode([VivoTenXAssayMapping].self, forKey: .assays)
    }
}

extension VivoMultiAssayDataset {
    public static var limits: VivoOmicsLimits {
        var l = VivoOmicsLimits()
        l.maximumFeatures = 200_000; l.maximumNonzeros = 20_000_000
        l.maximumInputBytes = 1_073_741_824
        return l
    }
    public func validate() throws {
        let limits = Self.limits
        guard schemaVersion == 1, vivoOmicsID(id), !sourceDescription.isEmpty,
              sourceDescription.utf8.count <= 16_384, !assays.isEmpty, assays.count <= 16,
              observations.count <= limits.maximumCells, spatialFrames.count <= 64 else {
            throw VivoOmicsError.invalid("multi-assay schema, identity or bounds")
        }
        guard !samples.isEmpty, samples.count <= limits.maximumCells else { throw VivoOmicsError.invalid("multi-assay sample bounds") }
        try vivoOmicsValidateIdentities(samples: samples, cells: observations.map(\.identity))
        var frameByID: [String: VivoSpatialFrame] = [:]
        for f in spatialFrames {
            guard vivoOmicsID(f.id), frameByID.updateValue(f, forKey: f.id) == nil,
                  (2...3).contains(f.axes.count), Set(f.axes).count == f.axes.count,
                  f.axes.allSatisfy(vivoOmicsID), !f.sourceDescription.isEmpty,
                  f.sourceDescription.utf8.count <= 16_384 else { throw VivoOmicsError.invalid("spatial frame") }
        }
        for o in observations {
            if let p = o.position {
                guard let f = frameByID[p.frameID], p.coordinates.count == f.axes.count,
                      p.coordinates.allSatisfy(\.isFinite) else { throw VivoOmicsError.invalid("spatial position/frame mismatch") }
            }
        }
        var ids = Set<String>(), remaining = limits.maximumNonzeros, remainingFeatures = limits.maximumFeatures
        for a in assays {
            guard vivoOmicsID(a.id), ids.insert(a.id).inserted, vivoOmicsID(a.featureNamespace),
                  !a.sourceDescription.isEmpty, a.sourceDescription.utf8.count <= 16_384,
                  a.genomeAssembly.map(vivoOmicsID) ?? true,
                  a.features.count == a.matrix.featureCount, a.observationIndices.count == a.matrix.cellCount,
                  Set(a.observationIndices).count == a.observationIndices.count,
                  a.observationIndices.allSatisfy({ observations.indices.contains($0) }),
                  a.matrix.counts.count <= remaining, a.features.count <= remainingFeatures else { throw VivoOmicsError.invalid("assay identities, row map or aggregate bounds") }
            remaining -= a.matrix.counts.count; remainingFeatures -= a.features.count
            try a.matrix.validate(limits: limits)
            if a.kind == .chromatinAccessibility {
                guard a.genomeAssembly != nil, a.countUnit != .umiCount else { throw VivoOmicsError.invalid("accessibility requires assembly and fragment/read units") }
            } else if a.countUnit == .fragmentCount {
                throw VivoOmicsError.invalid("fragment units require an accessibility assay")
            }
            var featureIDs = Set<String>()
            for f in a.features {
                guard vivoOmicsID(f.id), vivoOmicsID(f.name), featureIDs.insert(f.id).inserted else { throw VivoOmicsError.invalid("assay feature identity") }
                if let interval = f.interval {
                    guard a.genomeAssembly != nil, vivoOmicsID(interval.contig), interval.start < interval.end else { throw VivoOmicsError.invalid("zero-based half-open genomic interval") }
                }
                if a.kind == .chromatinAccessibility, f.interval == nil { throw VivoOmicsError.invalid("accessibility feature lacks genomic interval") }
            }
        }
    }
    /// Returns nil for an unmeasured assay/observation, zero only for measured
    /// sparse absence. Feature IDs live within the selected feature space.
    public func count(assayID: String, observationIndex: Int, featureID: String) throws -> UInt64? {
        try validate()
        guard observations.indices.contains(observationIndex), let a = assays.first(where: { $0.id == assayID }),
              let column = a.features.firstIndex(where: { $0.id == featureID }) else { throw VivoOmicsError.invalid("unknown assay, observation or feature") }
        guard let row = a.observationIndices.firstIndex(of: observationIndex) else { return nil }
        for k in a.matrix.rowOffsets[row]..<a.matrix.rowOffsets[row + 1] where a.matrix.featureIndices[k] == column { return a.matrix.counts[k] }
        return 0
    }
}
