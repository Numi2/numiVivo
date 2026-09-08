import Foundation

public struct VivoSingleCellLibraryPaths: Codable, Sendable, Equatable {
    public let matrix: String
    public let features: String
    public let barcodes: String
    public let metadata: VivoSingleCellImport
    private enum CodingKeys: String, CodingKey { case matrix, features, barcodes, metadata }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["matrix", "features", "barcodes", "metadata"])
        let values = try decoder.container(keyedBy: CodingKeys.self)
        matrix = try values.decode(String.self, forKey: .matrix)
        features = try values.decode(String.self, forKey: .features)
        barcodes = try values.decode(String.self, forKey: .barcodes)
        metadata = try values.decode(VivoSingleCellImport.self, forKey: .metadata)
    }
    public init(matrix: String, features: String, barcodes: String, metadata: VivoSingleCellImport) {
        self.matrix = matrix; self.features = features; self.barcodes = barcodes; self.metadata = metadata
    }
}
public struct VivoSingleCellManifest: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let id: String
    public let sourceDescription: String
    public let libraries: [VivoSingleCellLibraryPaths]
    public let normalizationTarget: Double?
    public let limits: VivoOmicsLimits?
    private enum CodingKeys: String, CodingKey { case schemaVersion, id, sourceDescription, libraries, normalizationTarget, limits }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["schemaVersion", "id", "sourceDescription", "libraries", "normalizationTarget", "limits"])
        let values = try decoder.container(keyedBy: CodingKeys.self)
        schemaVersion = try values.decode(Int.self, forKey: .schemaVersion)
        id = try values.decode(String.self, forKey: .id)
        sourceDescription = try values.decode(String.self, forKey: .sourceDescription)
        libraries = try values.decode([VivoSingleCellLibraryPaths].self, forKey: .libraries)
        normalizationTarget = try values.decodeIfPresent(Double.self, forKey: .normalizationTarget)
        limits = try values.decodeIfPresent(VivoOmicsLimits.self, forKey: .limits)
    }
    public init(id: String, sourceDescription: String, libraries: [VivoSingleCellLibraryPaths],
                normalizationTarget: Double? = nil, limits: VivoOmicsLimits? = nil) {
        self.schemaVersion = 1; self.id = id; self.sourceDescription = sourceDescription
        self.libraries = libraries; self.normalizationTarget = normalizationTarget; self.limits = limits
    }
    public func admittedLimits(ceiling: VivoOmicsLimits = .init()) throws -> VivoOmicsLimits {
        let requested = limits ?? .init()
        try requested.validate(); try ceiling.validate()
        guard schemaVersion == 1, vivoOmicsID(id), !sourceDescription.isEmpty,
              sourceDescription.utf8.count <= 16_384, !libraries.isEmpty, libraries.count <= 64 else {
            throw VivoOmicsError.invalid("manifest schema, identity or library count")
        }
        guard requested.maximumCells <= ceiling.maximumCells, requested.maximumFeatures <= ceiling.maximumFeatures,
              requested.maximumNonzeros <= ceiling.maximumNonzeros, requested.maximumInputBytes <= ceiling.maximumInputBytes,
              requested.maximumLineBytes <= ceiling.maximumLineBytes else {
            throw VivoOmicsError.limit("manifest exceeds caller admission ceiling")
        }
        if let target = normalizationTarget, !target.isFinite || target <= 0 {
            throw VivoOmicsError.invalid("normalization target")
        }
        for library in libraries {
            for path in [library.matrix, library.features, library.barcodes] {
                let parts = path.split(separator: "/", omittingEmptySubsequences: false)
                guard path.utf8.count <= 4096, !parts.isEmpty,
                      parts.allSatisfy({ !$0.isEmpty && $0 != "." && $0 != ".." &&
                          !$0.contains("\\") && !$0.contains("\0") && $0.utf8.count <= 240 }) else {
                    throw VivoOmicsError.invalid("source paths must be safe relative paths within the manifest directory")
                }
            }
        }
        return requested
    }
}
public struct VivoSingleCellLibraryBytes: Codable, Sendable, Equatable {
    public let matrix: Data
    public let features: Data
    public let barcodes: Data
    public init(matrix: Data, features: Data, barcodes: Data) {
        self.matrix = matrix; self.features = features; self.barcodes = barcodes
    }
}

/// Exact original manifest and source bytes, including any gzip encoding.
/// The existing artifact store retains sources independently of mutable paths.
public struct VivoSingleCellInputBundle: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let manifest: Data
    public let libraries: [VivoSingleCellLibraryBytes]
    public init(manifest: Data, libraries: [VivoSingleCellLibraryBytes]) {
        self.schemaVersion = 1; self.manifest = manifest; self.libraries = libraries
    }
}
public struct VivoSingleCellReport: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let numericalProfile: String
    public let dataset: VivoSingleCellDataset
    public let quality: [VivoCellQuality]
    public let normalized: VivoLogNormalizedCounts?
    public let pseudobulk: VivoPseudobulkCounts
}
public enum VivoSingleCellCampaign {
    public static let numericalProfile = "singlecell-exact-counts-fp64-log1p-v1"
    public static let maximumManifestBytes = 2 * 1_024 * 1_024
    public static func manifest(from bytes: Data, ceiling: VivoOmicsLimits = .init()) throws -> VivoSingleCellManifest {
        guard bytes.count <= maximumManifestBytes else { throw VivoOmicsError.limit("manifest bytes") }
        let manifest = try JSONDecoder().decode(VivoSingleCellManifest.self, from: bytes)
        _ = try manifest.admittedLimits(ceiling: ceiling)
        return manifest
    }
    public static func evaluate(_ input: VivoSingleCellInputBundle,
                                ceiling: VivoOmicsLimits = .init()) throws -> VivoSingleCellReport {
        let plan = try manifest(from: input.manifest, ceiling: ceiling)
        let limits = try plan.admittedLimits(ceiling: ceiling)
        guard input.schemaVersion == 1, input.libraries.count == plan.libraries.count else {
            throw VivoOmicsError.invalid("input bundle schema or library count")
        }
        var remaining = limits.maximumInputBytes
        for data in [input.manifest] + input.libraries.flatMap({ [$0.matrix, $0.features, $0.barcodes] }) {
            guard data.count <= remaining else { throw VivoOmicsError.limit("aggregate campaign source bytes") }
            remaining -= data.count
        }
        // The raw-source and expanded-source allowances are independent and
        // apply globally across all libraries, not once per compressed file.
        var remainingExpanded = limits.maximumInputBytes - input.manifest.count
        func expand(_ data: Data) throws -> Data {
            let result = try VivoOmicsSourceDecoder.decode(data, maximumExpandedBytes: remainingExpanded)
            remainingExpanded -= result.count; return result
        }
        var datasets: [VivoSingleCellDataset] = []
        var remainingCells = limits.maximumCells, remainingNNZ = limits.maximumNonzeros
        for (index, raw) in input.libraries.enumerated() {
            try Task.checkCancellation()
            var local = limits; local.maximumCells = max(1, remainingCells); local.maximumNonzeros = remainingNNZ
            let matrix = try expand(raw.matrix), features = try expand(raw.features), barcodes = try expand(raw.barcodes)
            let dataset = try VivoMatrixMarketCounts.decode(matrix: matrix, features: features, barcodes: barcodes,
                metadata: plan.libraries[index].metadata, limits: local)
            guard dataset.cells.count <= remainingCells else { throw VivoOmicsError.limit("aggregate campaign cells") }
            remainingCells -= dataset.cells.count; remainingNNZ -= dataset.matrix.counts.count
            datasets.append(dataset)
        }
        let dataset = try VivoSingleCellAnalysis.concatenate(datasets, id: plan.id, sourceDescription: plan.sourceDescription, limits: limits)
        let quality = try VivoSingleCellAnalysis.quality(dataset, limits: limits)
        let normalized = try plan.normalizationTarget.map { try VivoSingleCellAnalysis.logNormalize(dataset, targetSum: $0, limits: limits) }
        let bulk = try VivoSingleCellAnalysis.pseudobulk(dataset, limits: limits)
        try Task.checkCancellation()
        return .init(schemaVersion: 1, numericalProfile: numericalProfile, dataset: dataset,
                     quality: quality, normalized: normalized, pseudobulk: bulk)
    }
    /// Checks reconstruction from original input. It does not certify the assay,
    /// the supplied sample metadata, independent replication, or biology.
    public static func verify(_ report: VivoSingleCellReport, input: VivoSingleCellInputBundle,
                              ceiling: VivoOmicsLimits = .init()) throws {
        guard report.schemaVersion == 1, report.numericalProfile == numericalProfile,
              try evaluate(input, ceiling: ceiling) == report else {
            throw VivoOmicsError.invalid("report differs from native reconstruction under the declared profile")
        }
    }
}
