import Foundation
import Testing
@testable import NumiVivoKit

@Suite struct SingleCellInterchangeTests {
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
