import Foundation

/// New-directory publication only. Files are rendered before filesystem writes,
/// written through the shared rooted no-follow reader/writer, then moved from a
/// sibling staging directory. An existing destination is never overwritten.
public enum VivoOmicsDirectoryExport {
    public static func write(_ contents: [String: Data], to destination: URL) throws {
        guard destination.isFileURL, !contents.isEmpty else { throw VivoOmicsError.invalid("export destination or empty file set") }
        let url = destination.standardizedFileURL
        guard !FileManager.default.fileExists(atPath: url.path),
              (try? FileManager.default.attributesOfItem(atPath: url.path)) == nil else {
            throw VivoOmicsError.invalid("export destination already exists")
        }
        var bytes = 0
        for data in contents.values {
            guard data.count <= 128 * 1_024 * 1_024 - bytes else { throw VivoOmicsError.limit("export exceeds 128 MiB") }
            bytes += data.count
        }
        let staging = url.deletingLastPathComponent().appendingPathComponent(".numivivo-export-" + UUID().uuidString)
        try FileManager.default.createDirectory(at: staging, withIntermediateDirectories: false)
        defer { try? FileManager.default.removeItem(at: staging) }
        let files = try VivoRootedFileStore(rootURL: staging, createIfNeeded: false)
        for path in contents.keys.sorted() {
            try Task.checkCancellation()
            guard try files.writeFile(contents[path]!, relative: path, immutable: true) else { throw VivoOmicsError.invalid("duplicate staged export path") }
        }
        try Task.checkCancellation()
        try FileManager.default.moveItem(at: staging, to: url)
    }
}

public enum VivoSingleCellMEXExchange {
    /// Preserves exact UInt64 counts, sample/barcode identities, annotations and
    /// ordered features. The manifest groups rows by sample; the sidecar gives
    /// exported row -> source row so interleaved source order is not lost.
    public static func files(_ dataset: VivoSingleCellDataset) throws -> [String: Data] {
        try dataset.validate()
        guard dataset.samples.count <= 64 else { throw VivoOmicsError.limit("MEX manifest exports at most 64 sample libraries") }
        var output: [String: Data] = [:], remaining = 64 * 1_024 * 1_024
        func append(_ text: String, to data: inout Data) throws {
            let count = text.utf8.count
            guard count <= remaining else { throw VivoOmicsError.limit("MEX export exceeds source-byte admission limit") }
            remaining -= count; data.append(contentsOf: text.utf8)
        }
        var features = Data()
        for feature in dataset.features { try append("\(feature.id)\t\(feature.name)\tGene Expression\n", to: &features) }
        output["features.tsv"] = features
        var rowsBySample: [String: [Int]] = [:]
        for (row, cell) in dataset.cells.enumerated() { rowsBySample[cell.sampleID, default: []].append(row) }
        var libraries: [VivoSingleCellLibraryPaths] = [], lineage: [Int] = []
        let mitochondrial = dataset.features.filter(\.mitochondrial).map(\.id)
        for (sampleIndex, sample) in dataset.samples.enumerated() {
            try Task.checkCancellation()
            let rows = rowsBySample[sample.id] ?? [], directory = String(format: "library-%04d", sampleIndex)
            // Every library references the same feature file, but native source
            // admission counts each read, including repeated feature dictionaries.
            if sampleIndex > 0 {
                guard features.count <= remaining else { throw VivoOmicsError.limit("repeated feature dictionaries exceed import allowance") }
                remaining -= features.count
            }
            let nnz = rows.reduce(0) { $0 + dataset.matrix.rowOffsets[$1 + 1] - dataset.matrix.rowOffsets[$1] }
            var matrix = Data(), barcodes = Data(), groups: [String: String] = [:]
            try append("%%MatrixMarket matrix coordinate integer general\n\(dataset.features.count) \(rows.count) \(nnz)\n", to: &matrix)
            for (column, row) in rows.enumerated() {
                let cell = dataset.cells[row]
                try append(cell.barcode + "\n", to: &barcodes)
                if let group = cell.group { groups[cell.barcode] = group }
                for k in dataset.matrix.rowOffsets[row]..<dataset.matrix.rowOffsets[row + 1] {
                    try append("\(dataset.matrix.featureIndices[k] + 1) \(column + 1) \(dataset.matrix.counts[k])\n", to: &matrix)
                }
            }
            let matrixPath = directory + "/matrix.mtx", barcodePath = directory + "/barcodes.tsv"
            output[matrixPath] = matrix; output[barcodePath] = barcodes
            libraries.append(.init(matrix: matrixPath, features: "features.tsv", barcodes: barcodePath,
                metadata: .init(datasetID: sample.id, evidence: dataset.evidence, sourceDescription: dataset.sourceDescription,
                    countUnit: dataset.countUnit, sample: sample, mitochondrialFeatureIDs: mitochondrial, cellGroups: groups)))
            lineage.append(contentsOf: rows)
        }
        let manifest = try VivoCanonicalJSON.encode(VivoSingleCellManifest(id: dataset.id, sourceDescription: dataset.sourceDescription, libraries: libraries))
        guard manifest.count <= remaining, manifest.count <= VivoSingleCellCampaign.maximumManifestBytes else {
            throw VivoOmicsError.limit("export manifest exceeds source or manifest allowance")
        }
        output["manifest.json"] = manifest
        struct Lineage: Codable {
            let schema: String
            let sourceDataset: VivoFingerprint
            let sourceCellIndices: [Int]
            let featureOrder: String
        }
        output["source-cell-indices.json"] = try VivoCanonicalJSON.encode(Lineage(schema: "numivivo.org/mex-row-lineage/v1",
            sourceDataset: VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(dataset)), sourceCellIndices: lineage,
            featureOrder: "unchanged from source dataset"))
        return output
    }
}

public enum VivoSingleCellAnalysisTables {
    public static func files(_ report: VivoSingleCellCohortReport, receipt: VivoSingleCellAnalysisReceipt) throws -> [String: Data] {
        var result: [String: Data] = [:], remaining = 128 * 1_024 * 1_024
        func line(_ values: [String], to data: inout Data) throws {
            guard values.allSatisfy({ !$0.contains("\t") && !$0.contains("\n") && !$0.contains("\r") }) else { throw VivoOmicsError.invalid("unescaped tabular value") }
            let text = values.joined(separator: "\t") + "\n"
            guard text.utf8.count <= remaining else { throw VivoOmicsError.limit("tabular export exceeds 128 MiB") }
            remaining -= text.utf8.count; data.append(contentsOf: text.utf8)
        }
        func number(_ value: Double?) -> String { value.map { String($0) } ?? "" }
        var quality = Data()
        try line(["source_cell", "sample", "barcode", "accepted", "counts", "detected_features", "mitochondrial_fraction", "reasons", "mitochondrial_threshold_not_evaluated"], to: &quality)
        for row in report.processed.decisions {
            try line([String(row.sourceCellIndex), row.identity.sampleID, row.identity.barcode, String(row.accepted),
                String(row.quality.totalCounts), String(row.quality.detectedFeatures), number(row.quality.mitochondrialFraction),
                row.reasons.map(\.rawValue).joined(separator: ";"), String(row.mitochondrialThresholdNotEvaluated)], to: &quality)
        }
        result["cell-quality.tsv"] = quality
        var features = Data()
        try line(["feature_index", "feature_id", "counts", "detected_cells", "mean_log_normalized", "variance_log_normalized"], to: &features)
        for row in report.processed.features {
            try line([String(row.featureIndex), row.featureID, String(row.totalCounts), String(row.detectedCells),
                      number(row.meanLogNormalized), number(row.varianceLogNormalized)], to: &features)
        }
        result["feature-statistics.tsv"] = features
        struct ContrastFile: Codable { let id: String; let method: String; let statistics: String; let design: String; let diagnostics: String? }
        var contrastFiles: [ContrastFile] = []
        for (index, contrast) in report.contrasts.enumerated() {
            try Task.checkCancellation()
            let prefix = String(format: "contrast-%04d", index)
            var statistics = Data(), design = Data()
            try line(["feature_index", "feature_id", "status", "counts", "expressing_pseudobulks", "mean_normalized_count", "log2_effect", "standard_error", "t", "df", "interval_lower", "interval_upper", "p_value", "BH_adjusted_p"] + (contrast.negativeBinomial == nil ? [] : ["z"]), to: &statistics)
            for row in contrast.features {
                try line([String(row.featureIndex), row.featureID, row.status.rawValue, String(row.totalCounts), String(row.expressingPseudobulks),
                    String(row.meanNormalizedCount), number(row.log2FoldChange), number(row.standardError), number(row.tStatistic),
                    number(row.degreesOfFreedom), number(row.intervalLower), number(row.intervalUpper), number(row.pValue), number(row.adjustedPValue)] + (contrast.negativeBinomial == nil ? [] : [number(row.zStatistic)]), to: &statistics)
            }
            try line(["source_pseudobulk", "replicate", "donor", "condition", "size_factor", "library_counts"] + contrast.design.columnNames, to: &design)
            for (row, observation) in contrast.design.observations.enumerated() {
                try line([String(contrast.design.sourcePseudobulkIndices[row]), observation.biologicalReplicateID, observation.donorID ?? "",
                    observation.condition, String(contrast.design.sizeFactorValues[row]), String(contrast.design.libraryCounts[row])] +
                    contrast.design.rows[row].map { String($0) }, to: &design)
            }
            result[prefix + ".tsv"] = statistics; result[prefix + "-design.tsv"] = design
            let diagnosticsPath = contrast.negativeBinomial.map { _ in prefix + "-nb-diagnostics.json" }
            if let diagnostics = contrast.negativeBinomial, let path = diagnosticsPath {
                result[path] = try VivoCanonicalJSON.encode(diagnostics)
            }
            contrastFiles.append(.init(id: contrast.request.id, method: contrast.method, statistics: prefix + ".tsv", design: prefix + "-design.tsv", diagnostics: diagnosticsPath))
        }
        struct Index: Codable {
            let schema: String
            let receipt: VivoSingleCellAnalysisReceipt
            let contrasts: [ContrastFile]
            let missingValues: String
        }
        result["index.json"] = try VivoCanonicalJSON.encode(Index(schema: "numivivo.org/singlecell-analysis-tables/v1", receipt: receipt,
            contrasts: contrastFiles, missingValues: "empty TSV fields mean not available or not tested, never zero"))
        return result
    }
}

extension VivoSingleCellCampaignIO {
    public static func readDocument(_ url: URL, maximumBytes: Int) throws -> Data {
        guard url.isFileURL else { throw VivoOmicsError.invalid("document must be a local file") }
        let path = url.standardizedFileURL
        let files = try VivoRootedFileStore(rootURL: path.deletingLastPathComponent(), createIfNeeded: false)
        return try files.readFile(path.lastPathComponent, maximumBytes: maximumBytes)
    }
}
