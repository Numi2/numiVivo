import Foundation
import Testing
@testable import NumiVivoKit

@Suite struct SingleCellInterchangeTests {
    @Test func h5adImportMappingRejectsAmbiguousPathsAndColumns() throws {
        let source = try VivoSingleCellExamples.pairedCounts()
        func mapping(matrixPath: String = "X", sampleColumn: String = "sample",
                     samples: [VivoOmicsSample] = source.samples,
                     mitochondrial: [String] = []) -> VivoH5ADImportPlan {
            .init(id: "mapping", evidence: .synthetic, sourceDescription: "mapping fixture",
                  countUnit: .umiCount, matrixPath: matrixPath, samples: samples,
                  sampleColumn: sampleColumn, mitochondrialFeatureIDs: mitochondrial)
        }
        try mapping(matrixPath: "layers/counts").validate()
        for invalid in [
            mapping(matrixPath: "layers/counts/extra"),
            mapping(sampleColumn: "obs/sample"),
            mapping(samples: []),
            mapping(mitochondrial: [source.features[0].id, source.features[0].id])
        ] {
            #expect(throws: (any Error).self) { try invalid.validate() }
        }
    }

    @Test func h5adDerivedPlansRejectInvalidMappingsBeforeSourceWork() throws {
        let source = try VivoSingleCellExamples.pairedCounts()
        let invalid = VivoH5ADImportPlan(id: source.id, evidence: source.evidence,
            sourceDescription: source.sourceDescription, countUnit: source.countUnit,
            matrixPath: "layers/counts/extra", samples: source.samples, sampleColumn: "sample")
        let program = VivoSingleCellProgramDefinition(id: "fixture", organism: "human",
            featureNamespace: "fixture-id", sourceURI: "urn:numivivo:fixture", sourceVersion: "1",
            sourceDescription: "Numerical fixture only", members: [.init(featureID: source.features[0].id)])
        let options = VivoSingleCellProgramOptions(definitions: [program])
        #expect(throws: (any Error).self) { try VivoH5ADPseudobulkPlan(mapping: invalid).validate() }
        #expect(throws: (any Error).self) { try VivoH5ADPCAPlan(mapping: invalid).validate() }
        #expect(throws: (any Error).self) { try VivoH5ADProgramPlan(mapping: invalid, programs: options).validate() }
        #expect(throws: (any Error).self) { try VivoH5ADPCAQueryPlan(mapping: invalid, featureNamespace: "fixture-id").validate() }
    }

    @Test func gzipChecksumsMembersAndExpansionLimits() throws {
        let plain = Data("native compressed count fixture\n".utf8)
        let gzip = Data(base64Encoded: "H4sIAAAAAAAC/8tLLMksS1VIzs8tKEotLk5NATJL80oU0jIrSkqLUrkA65t8niAAAAA=")!
        #expect(try VivoOmicsSourceDecoder.decode(gzip, maximumExpandedBytes: plain.count) == plain)
        #expect(try VivoOmicsSourceDecoder.decode(gzip + gzip, maximumExpandedBytes: 2 * plain.count) == plain + plain)
        #expect(throws: (any Error).self) { try VivoOmicsSourceDecoder.decode(gzip, maximumExpandedBytes: plain.count - 1) }
        #expect(throws: (any Error).self) { try VivoOmicsSourceDecoder.decode(Data(gzip.dropLast()), maximumExpandedBytes: 1000) }
        #expect(throws: (any Error).self) { try VivoOmicsSourceDecoder.decode(gzip + Data([0]), maximumExpandedBytes: 1000) }
        var corrupt = gzip; corrupt[corrupt.count - 8] ^= 1
        #expect(throws: (any Error).self) { try VivoOmicsSourceDecoder.decode(corrupt, maximumExpandedBytes: 1000) }
    }
    @Test func gzipReadableSnapshotIsPrivateAndPlainSourcesStayInPlace() throws {
        let root = FileManager.default.temporaryDirectory.appendingPathComponent("numivivo-h5ad-readable-" + UUID().uuidString)
        try FileManager.default.createDirectory(at: root, withIntermediateDirectories: false)
        defer { try? FileManager.default.removeItem(at: root) }
        let plain = Data("native compressed count fixture\n".utf8)
        let gzip = Data(base64Encoded: "H4sIAAAAAAAC/8tLLMksS1VIzs8tKEotLk5NATJL80oU0jIrSkqLUrkA65t8niAAAAA=")!
        let plainURL = root.appendingPathComponent("plain.h5ad"), gzipURL = root.appendingPathComponent("wrapped.h5ad")
        try plain.write(to: plainURL); try gzip.write(to: gzipURL)
        let plainResult = try VivoSingleCellH5AD.withReadableSnapshot(plainURL, limits: .init()) { readable in
            #expect(readable == plainURL)
            return try Data(contentsOf: readable)
        }
        #expect(plainResult == plain)
        var decodedURL: URL?
        let gzipResult = try VivoSingleCellH5AD.withReadableSnapshot(gzipURL, limits: .init()) { readable in
            decodedURL = readable
            #expect(readable != gzipURL)
            return try Data(contentsOf: readable)
        }
        #expect(gzipResult == plain)
        if let decodedURL { #expect(!FileManager.default.fileExists(atPath: decodedURL.path)) }
    }
    @Test func nativeMEXRoundTripPreservesCountsAnnotationsAndIdentities() throws {
        let root = FileManager.default.temporaryDirectory.appendingPathComponent("numivivo-mex-" + UUID().uuidString)
        try FileManager.default.createDirectory(at: root, withIntermediateDirectories: true)
        defer { try? FileManager.default.removeItem(at: root) }
        let source = try VivoSingleCellExamples.pairedCounts()
        let destination = root.appendingPathComponent("export")
        let files = try VivoSingleCellMEXExchange.files(source)
        try VivoOmicsDirectoryExport.write(files, to: destination)
        let input = try VivoSingleCellCampaignIO.snapshot(manifestURL: destination.appendingPathComponent("manifest.json"))
        let imported = try VivoSingleCellCampaign.evaluate(input)
        #expect(imported.dataset == source)
        #expect(files["source-cell-indices.json"] != nil)
        #expect(throws: (any Error).self) { try VivoOmicsDirectoryExport.write(files, to: destination) }
        #expect(try VivoSingleCellCampaignIO.snapshot(manifestURL: destination.appendingPathComponent("manifest.json")) == input)
    }
    private func legacySamples() -> [VivoOmicsSample] {
        [.init(id: "guide-a", biologicalReplicateID: "replicate-a", condition: "treated", batchID: "batch-a", organism: "human"),
         .init(id: "guide-b", biologicalReplicateID: "replicate-b", condition: "control", batchID: "batch-b", organism: "human")]
    }
    private func legacyMapping() -> VivoH5ADImportPlan {
        .init(id: "legacy-h5ad-fixture", evidence: .synthetic, sourceDescription: "Legacy H5AD interoperability fixture",
              countUnit: .umiCount, matrixPath: "X", samples: legacySamples(), sampleColumn: "sample")
    }
    /// Write the AnnData 0.7 dataframe dialect directly: an untagged root,
    /// direct string indexes, an object-reference categorical column, and an
    /// untagged Float32 dense X.  The fixture is deliberately unlike the
    /// project's current writer so it exercises the compatibility reader.
    private func writeLegacyH5AD(_ url: URL, values: [Float], codes: [Int16] = [0, 1],
                                 labels: [String] = ["guide-a", "guide-b"], rootEncoding: Bool = false,
                                 invalidCategoryReferenceType: Bool = false) throws {
        guard values.count == 4, codes.count == 2 else { throw VivoOmicsError.invalid("legacy H5AD fixture shape") }
        try VivoHDF5.lock.withLock {
            let h = try VivoHDF5(), file = try h.file(url.path, create: true)
            defer { h.close(file, "H5Fclose") }
            if rootEncoding { try h.encoding(file, "anndata", "0.1.0") }
            let obs = try h.group(file, "obs"); defer { h.close(obs, "H5Gclose") }
            try h.encoding(obs, "dataframe", "0.1.0")
            try h.writeStrings(obs, "_index", ["cell_barcode"], attribute: true, scalar: true)
            try h.writeStrings(obs, "column-order", ["sample"], attribute: true)
            try h.writeStrings(obs, "cell_barcode", ["cell-a", "cell-b"])
            try codes.withUnsafeBytes {
                try h.write(obs, "sample", type: h.native("NATIVE_SHORT"), dimensions: [2], attribute: false, buffer: $0.baseAddress)
            }
            let categories = try h.group(file, "uns"); defer { h.close(categories, "H5Gclose") }
            try h.writeStrings(categories, "sample_categories", labels)
            let sample = try h.dataset(obs, "sample"); defer { h.close(sample, "H5Dclose") }
            if invalidCategoryReferenceType {
                var invalid: UInt64 = 0
                try withUnsafeBytes(of: &invalid) {
                    try h.write(sample, "categories", type: h.native("NATIVE_ULLONG"), dimensions: [], attribute: true, buffer: $0.baseAddress)
                }
            } else {
                try h.legacyWriteObjectReference(sample, "categories", file: file, path: "uns/sample_categories")
            }
            let variable = try h.group(file, "var"); defer { h.close(variable, "H5Gclose") }
            try h.encoding(variable, "dataframe", "0.1.0")
            try h.writeStrings(variable, "_index", ["gene_id"], attribute: true, scalar: true)
            try h.writeStrings(variable, "column-order", [], attribute: true)
            try h.writeStrings(variable, "gene_id", ["gene-a", "gene-b"])
            try values.withUnsafeBytes {
                try h.write(file, "X", type: h.native("NATIVE_FLOAT"), dimensions: [2, 2], attribute: false, buffer: $0.baseAddress)
            }
        }
    }
    @Test func legacyDataframeH5ADStreamsExactCountsAndPreflightsWithoutAStore() throws {
        let root = FileManager.default.temporaryDirectory.appendingPathComponent("numivivo-h5ad-legacy-" + UUID().uuidString)
        try FileManager.default.createDirectory(at: root, withIntermediateDirectories: false)
        defer { try? FileManager.default.removeItem(at: root) }
        let source = root.appendingPathComponent("legacy.h5ad")
        try writeLegacyH5AD(source, values: [1, 0, 0, 2])
        let preflight = try VivoSingleCellH5AD.preflightCountMatrix(source)
        #expect(preflight.format == "numivivo.org/h5ad-count-matrix-preflight/v1")
        #expect(preflight.matrixPath == "X" && preflight.rows == 2 && preflight.features == 2)
        #expect(preflight.nonzeros == 2 && preflight.countRecordBytes == 32)
        #expect(preflight.retainedSourceAndCountBytes == preflight.sourceBytes + preflight.countRecordBytes)
        let document = try VivoSingleCellH5AD.read(source, plan: legacyMapping())
        #expect(document.dataset.matrix.counts == [1, 2])
        let store = root.appendingPathComponent("store")
        let implementation = try VivoFingerprint(bytes: Array(repeating: 9, count: 32))
        let receipt = try VivoH5ADCountStore.publish(source: source, plan: legacyMapping(), implementation: implementation, to: store)
        #expect(receipt.entries == 2)
        #expect(try VivoH5ADCountStore.verify(store, implementation: implementation) == receipt)
    }
    @Test func legacyDataframeH5ADRejectsInvalidFloatCountsAndCategoryMappings() throws {
        let root = FileManager.default.temporaryDirectory.appendingPathComponent("numivivo-h5ad-legacy-invalid-" + UUID().uuidString)
        try FileManager.default.createDirectory(at: root, withIntermediateDirectories: false)
        defer { try? FileManager.default.removeItem(at: root) }
        for (name, values) in [("fractional", [Float(1.5), 0, 0, 2]), ("negative", [Float(-1), 0, 0, 2]), ("nan", [Float.nan, 0, 0, 2]), ("above-exact-count-limit", [Float(18_014_398_509_481_984), 0, 0, 2])] {
            let source = root.appendingPathComponent(name + ".h5ad")
            try writeLegacyH5AD(source, values: values)
            #expect(throws: (any Error).self) { try VivoSingleCellH5AD.preflightCountMatrix(source) }
        }
        let badCode = root.appendingPathComponent("bad-code.h5ad")
        try writeLegacyH5AD(badCode, values: [1, 0, 0, 2], codes: [0, 2])
        #expect(throws: (any Error).self) { try VivoSingleCellH5AD.read(badCode, plan: legacyMapping()) }
        let badReference = root.appendingPathComponent("bad-reference.h5ad")
        try writeLegacyH5AD(badReference, values: [1, 0, 0, 2], invalidCategoryReferenceType: true)
        #expect(throws: (any Error).self) { try VivoSingleCellH5AD.read(badReference, plan: legacyMapping()) }
        let mixed = root.appendingPathComponent("mixed-root.h5ad")
        try writeLegacyH5AD(mixed, values: [1, 0, 0, 2], rootEncoding: true)
        #expect(throws: (any Error).self) { try VivoSingleCellH5AD.preflightCountMatrix(mixed) }
    }
    @Test func nativeH5ADWriterRoundTripsThroughCountReaderAndPreservesSourceBytes() throws {
        let root = FileManager.default.temporaryDirectory.appendingPathComponent("numivivo-h5ad-roundtrip-" + UUID().uuidString)
        try FileManager.default.createDirectory(at: root, withIntermediateDirectories: false)
        defer { try? FileManager.default.removeItem(at: root) }
        let source = try VivoSingleCellExamples.pairedCounts()
        let native = root.appendingPathComponent("native.h5ad")
        try VivoSingleCellH5AD.write(source, to: native)
        let mapping = VivoH5ADImportPlan(id: source.id, evidence: source.evidence,
            sourceDescription: source.sourceDescription, countUnit: source.countUnit,
            matrixPath: "X", samples: source.samples, sampleColumn: "sample",
            barcodeColumn: "barcode", groupColumn: "group", featureNameColumn: "name",
            mitochondrialFeatureIDs: ["g31"])
        let document = try VivoSingleCellH5AD.read(native, plan: mapping)
        #expect(document.dataset == source)
        let exported = root.appendingPathComponent("exported.h5ad")
        try document.exportOriginal(to: exported)
        #expect(try Data(contentsOf: exported) == Data(contentsOf: native))
    }
    @Test func rawExchangeDoesNotRoundCountsAboveDoubleIntegerPrecision() throws {
        let fixture = try VivoSingleCellExamples.pairedCounts()
        var counts = fixture.matrix.counts; counts[0] = 9_007_199_254_740_993
        let source = VivoSingleCellDataset(id: fixture.id, evidence: fixture.evidence, sourceDescription: fixture.sourceDescription,
            countUnit: fixture.countUnit, samples: fixture.samples, features: fixture.features, cells: fixture.cells,
            matrix: .init(cellCount: fixture.matrix.cellCount, featureCount: fixture.matrix.featureCount,
                rowOffsets: fixture.matrix.rowOffsets, featureIndices: fixture.matrix.featureIndices, counts: counts))
        let files = try VivoSingleCellMEXExchange.files(source)
        let manifest = try VivoSingleCellCampaign.manifest(from: files["manifest.json"]!)
        let input = VivoSingleCellInputBundle(manifest: files["manifest.json"]!, libraries: manifest.libraries.map {
            .init(matrix: files[$0.matrix]!, features: files[$0.features]!, barcodes: files[$0.barcodes]!)
        })
        #expect(try VivoSingleCellCampaign.evaluate(input).dataset == source)
        #expect(throws: (any Error).self) {
            try VivoPseudobulkDifferentialExpression.run(source, contrast: VivoSingleCellExamples.pairedPlan().contrasts[0])
        }
    }
    @Test func completeNativeAnalysisReconstructsAndTablesRetainDesign() throws {
        let source = try VivoSingleCellExamples.pairedCounts(), plan = VivoSingleCellExamples.pairedPlan()
        let report = try VivoSingleCellCohortAnalysis.run(source, plan: plan)
        #expect(report.processed.dataset.cells.count == 24)
        #expect(report.contrasts.count == 1 && report.contrasts[0].design.residualDegreesOfFreedom == 5)
        #expect(try VivoSingleCellCohortAnalysis.run(source, plan: plan) == report)
        let fingerprint = try VivoFingerprint(bytes: [UInt8](repeating: 1, count: 32))
        let tables = try VivoSingleCellAnalysisTables.files(report, receipt: .init(request: fingerprint, result: fingerprint, implementation: fingerprint))
        #expect(tables["cell-quality.tsv"] != nil && tables["contrast-0000-design.tsv"] != nil && tables["index.json"] != nil)
        #expect(report.contrasts[0].features.prefix(3).allSatisfy { ($0.log2FoldChange ?? 0) > 1.5 })
    }
}
