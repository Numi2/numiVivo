import Foundation

public struct VivoCellQuality: Codable, Sendable, Equatable {
    public let sampleID: String
    public let barcode: String
    public let totalCounts: UInt64
    public let detectedFeatures: Int
    public let mitochondrialCounts: UInt64
    public let mitochondrialFeatureCount: Int
    /// Nil for an empty cell or absent feature annotations, not a fabricated zero percent.
    public let mitochondrialFraction: Double?
}
/// Cell-only QC shared by streamed PCA and pseudobulk; no group aggregation.
final class VivoSingleCellQualityAccumulator {
    let metadata: VivoSingleCellCountMetadata
    private var totals: [UInt64], mitochondrial: [UInt64], detected: [Int]
    private(set) var nonzeros = 0
    init(_ metadata: VivoSingleCellCountMetadata) {
        self.metadata = metadata
        totals = .init(repeating: 0, count: metadata.cells.count)
        mitochondrial = totals; detected = .init(repeating: 0, count: metadata.cells.count)
    }
    func add(row: Int, feature: Int, count: UInt64) throws {
        totals[row] = try vivoOmicsSum(totals[row], count)
        detected[row] += 1; nonzeros += 1
        if metadata.features[feature].mitochondrial {
            mitochondrial[row] = try vivoOmicsSum(mitochondrial[row], count)
        }
    }
    func finish() -> [VivoCellQuality] {
        let mitoFeatures = metadata.features.filter(\.mitochondrial).count
        return metadata.cells.indices.map { i in
            .init(sampleID: metadata.cells[i].sampleID, barcode: metadata.cells[i].barcode,
                totalCounts: totals[i], detectedFeatures: detected[i], mitochondrialCounts: mitochondrial[i],
                mitochondrialFeatureCount: mitoFeatures,
                mitochondrialFraction: totals[i] == 0 || mitoFeatures == 0 ? nil : Double(mitochondrial[i]) / Double(totals[i]))
        }
    }
}
public struct VivoLogNormalizedCounts: Codable, Sendable, Equatable {
    public let method: String
    public let targetSum: Double
    public let values: [Double]
    /// Same sparsity/row ordering as the source; never usable as raw counts.
    public let rowOffsets: [Int]
    public let featureIndices: [Int]
    public let featureCount: Int
}
public struct VivoPseudobulkGroup: Codable, Sendable, Equatable {
    public let biologicalReplicateID: String
    public let donorID: String?
    public let condition: String
    public let organism: String
    public let cellGroup: String?
    public let sampleIDs: [String]
    public let batchIDs: [String]
    public let sourceCellIndices: [Int]
}
public struct VivoPseudobulkCounts: Codable, Sendable, Equatable {
    public let method: String
    public let countUnit: VivoOmicsCountUnit
    public let groups: [VivoPseudobulkGroup]
    public let featureIDs: [String]
    /// Rows are pseudobulk groups, not cells. The shared sparse representation
    /// uses cellCount as its row-count field; this wrapper declares row semantics.
    public let matrix: VivoSparseCounts
}

public enum VivoSingleCellAnalysis {
    /// Concatenate libraries with an identical ordered feature dictionary. There
    /// is no implicit gene-ID remapping, unit conversion or evidence promotion.
    public static func concatenate(_ inputs: [VivoSingleCellDataset], id: String,
                                   sourceDescription: String, limits: VivoOmicsLimits = .init()) throws -> VivoSingleCellDataset {
        try limits.validate()
        guard let first = inputs.first, inputs.count <= limits.maximumCells else {
            throw VivoOmicsError.invalid("at least one bounded input library is required")
        }
        var remainingCells = limits.maximumCells, remainingNNZ = limits.maximumNonzeros
        var samples: [VivoOmicsSample] = [], cells: [VivoOmicsCell] = []
        var offsets = [0], columns: [Int] = [], counts: [UInt64] = []
        for input in inputs {
            try input.validate(limits: limits)
            guard input.features == first.features, input.countUnit == first.countUnit,
                  input.evidence == first.evidence else {
                throw VivoOmicsError.invalid("libraries differ in ordered features, count units or evidence class")
            }
            guard input.cells.count <= remainingCells, input.matrix.counts.count <= remainingNNZ,
                  input.samples.count <= limits.maximumCells - samples.count else {
                throw VivoOmicsError.limit("combined cells, samples or nonzeros")
            }
            remainingCells -= input.cells.count; remainingNNZ -= input.matrix.counts.count
            let base = counts.count
            offsets += input.matrix.rowOffsets.dropFirst().map { $0 + base }
            columns += input.matrix.featureIndices; counts += input.matrix.counts
            samples += input.samples; cells += input.cells
        }
        let output = VivoSingleCellDataset(id: id, evidence: first.evidence, sourceDescription: sourceDescription,
            countUnit: first.countUnit, samples: samples, features: first.features, cells: cells,
            matrix: .init(cellCount: cells.count, featureCount: first.features.count,
                          rowOffsets: offsets, featureIndices: columns, counts: counts))
        try output.validate(limits: limits); return output
    }

    public static func quality(_ dataset: VivoSingleCellDataset, limits: VivoOmicsLimits = .init()) throws -> [VivoCellQuality] {
        try dataset.validate(limits: limits)
        let mitochondrialFeatures = dataset.features.filter(\.mitochondrial).count
        return try dataset.cells.indices.map { row in
            let cell = dataset.cells[row], start = dataset.matrix.rowOffsets[row], end = dataset.matrix.rowOffsets[row + 1]
            var total: UInt64 = 0, mitochondrial: UInt64 = 0
            for k in start..<end {
                let count = dataset.matrix.counts[k]
                total = try vivoOmicsSum(total, count)
                if dataset.features[dataset.matrix.featureIndices[k]].mitochondrial { mitochondrial = try vivoOmicsSum(mitochondrial, count) }
            }
            return .init(sampleID: cell.sampleID, barcode: cell.barcode, totalCounts: total,
                         detectedFeatures: end - start, mitochondrialCounts: mitochondrial, mitochondrialFeatureCount: mitochondrialFeatures,
                         mitochondrialFraction: total == 0 || mitochondrialFeatures == 0 ? nil : Double(mitochondrial) / Double(total))
        }
    }
    public static func logNormalize(_ dataset: VivoSingleCellDataset, targetSum: Double = 10_000,
                                    limits: VivoOmicsLimits = .init()) throws -> VivoLogNormalizedCounts {
        guard targetSum.isFinite, targetSum > 0 else { throw VivoOmicsError.invalid("normalization target must be finite and positive") }
        let qc = try quality(dataset, limits: limits)
        var values: [Double] = []; values.reserveCapacity(dataset.matrix.counts.count)
        for row in dataset.cells.indices {
            for k in dataset.matrix.rowOffsets[row]..<dataset.matrix.rowOffsets[row + 1] {
                // Divide first: an individual count/total is bounded by one.
                let value = log1p((Double(dataset.matrix.counts[k]) / Double(qc[row].totalCounts)) * targetSum)
                guard value.isFinite else { throw VivoOmicsError.invalid("nonfinite normalized value") }
                values.append(value)
            }
        }
        return .init(method: "library-size-log1p-v1", targetSum: targetSum, values: values,
                     rowOffsets: dataset.matrix.rowOffsets, featureIndices: dataset.matrix.featureIndices,
                     featureCount: dataset.matrix.featureCount)
    }
    /// Technical libraries sharing the explicitly supplied replicate, condition
    /// and cell group are pooled. Donor IDs remain visible for repeated measures.
    /// This is aggregation, not a test of differential expression or independence.
    public static func pseudobulk(_ dataset: VivoSingleCellDataset, limits: VivoOmicsLimits = .init()) throws -> VivoPseudobulkCounts {
        try dataset.validate(limits: limits)
        struct Key: Hashable { let replicate: String; let condition: String; let group: String? }
        let samples = Dictionary(uniqueKeysWithValues: dataset.samples.map { ($0.id, $0) })
        var members: [Key: [Int]] = [:]
        for (index, cell) in dataset.cells.enumerated() {
            guard let sample = samples[cell.sampleID] else { throw VivoOmicsError.invalid("missing sample") }
            members[.init(replicate: sample.biologicalReplicateID, condition: sample.condition, group: cell.group), default: []].append(index)
        }
        let keys = members.keys.sorted {
            if $0.replicate != $1.replicate { return $0.replicate < $1.replicate }
            if $0.condition != $1.condition { return $0.condition < $1.condition }
            if $0.group == nil { return $1.group != nil }
            if $1.group == nil { return false }
            return $0.group! < $1.group!
        }
        var groups: [VivoPseudobulkGroup] = [], offsets = [0], columns: [Int] = [], values: [UInt64] = []
        for key in keys {
            let rows = members[key]!
            var sums: [Int: UInt64] = [:], sampleIDs = Set<String>(), batches = Set<String>()
            let representative = samples[dataset.cells[rows[0]].sampleID]!
            for row in rows {
                let sample = samples[dataset.cells[row].sampleID]!
                sampleIDs.insert(sample.id); batches.insert(sample.batchID)
                for k in dataset.matrix.rowOffsets[row]..<dataset.matrix.rowOffsets[row + 1] {
                    let column = dataset.matrix.featureIndices[k]
                    sums[column] = try vivoOmicsSum(sums[column, default: 0], dataset.matrix.counts[k])
                }
            }
            for column in sums.keys.sorted() { columns.append(column); values.append(sums[column]!) }
            offsets.append(values.count)
            groups.append(.init(biologicalReplicateID: key.replicate, donorID: representative.donorID,
                condition: key.condition, organism: representative.organism, cellGroup: key.group,
                sampleIDs: sampleIDs.sorted(), batchIDs: batches.sorted(), sourceCellIndices: rows))
        }
        let matrix = VivoSparseCounts(cellCount: groups.count, featureCount: dataset.features.count,
                                      rowOffsets: offsets, featureIndices: columns, counts: values)
        try matrix.validate(limits: limits)
        return .init(method: "raw-sum-by-replicate-condition-group-v1", countUnit: dataset.countUnit,
                     groups: groups, featureIDs: dataset.features.map(\.id), matrix: matrix)
    }
}
