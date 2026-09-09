import Foundation
import Testing
@testable import NumiVivoKit

@Suite struct SingleCellPCABundleTests {
    @Test func sharedStreamedQualityMatchesResidentQuality() throws {
        let data = try VivoSingleCellExamples.pairedCounts()
        let metadata = VivoSingleCellCountMetadata(id: data.id, evidence: data.evidence, sourceDescription: data.sourceDescription,
            countUnit: data.countUnit, samples: data.samples, features: data.features, cells: data.cells)
        let qc = VivoSingleCellQualityAccumulator(metadata)
        #expect(qc.finish().allSatisfy { $0.totalCounts == 0 && $0.mitochondrialFraction == nil })
        for row in data.cells.indices {
            for k in data.matrix.rowOffsets[row]..<data.matrix.rowOffsets[row + 1] {
                try qc.add(row: row, feature: data.matrix.featureIndices[k], count: data.matrix.counts[k])
            }
        }
        #expect(qc.finish() == (try VivoSingleCellAnalysis.quality(data)))
        #expect(qc.nonzeros == data.matrix.counts.count)
    }
    @Test func densePCARecordsRetainSignedAndExtremeFiniteValues() throws {
        let root = FileManager.default.temporaryDirectory.appendingPathComponent(UUID().uuidString)
        try FileManager.default.createDirectory(at: root, withIntermediateDirectories: false)
        defer { try? FileManager.default.removeItem(at: root) }
        let values: [[Double]] = [[-0.0, -2, .leastNonzeroMagnitude], [1, .greatestFiniteMagnitude, 0]]
        let url = root.appendingPathComponent("scores.bin")
        let hash = try VivoH5ADPCA.writeMatrix(values, columns: 3, to: url)
        #expect(try VivoH5ADCountStore.fingerprint(url) == hash)
        let records = try VivoWindowedCountRecords(url, entries: 6)
        for i in 0..<6 {
            let v = try records.record(i)
            #expect(v.row == i / 3 && v.feature == i % 3 && v.bits == values[i / 3][i % 3].bitPattern)
        }
        #expect(throws: (any Error).self) { try VivoH5ADPCA.writeMatrix([[.nan]], columns: 1, to: root.appendingPathComponent("nan")) }
        #expect(throws: (any Error).self) { try VivoH5ADPCA.writeMatrix([[1]], columns: 2, to: root.appendingPathComponent("width")) }
    }
    @Test func standalonePlanRequiresProjectionCenters() throws {
        let data = try VivoSingleCellExamples.pairedCounts()
        let mapping = VivoH5ADImportPlan(id: data.id, evidence: data.evidence, sourceDescription: data.sourceDescription,
            countUnit: data.countUnit, matrixPath: "X", samples: data.samples, sampleColumn: "sample")
        var plan = VivoH5ADPCAPlan(mapping: mapping)
        try plan.validate()
        #expect(plan.reduction.pca.retainProjectionCenters == true)
        #expect(try VivoCanonicalJSON.decode(VivoH5ADPCAPlan.self, from: VivoCanonicalJSON.encode(plan)) == plan)
        plan.reduction.pca.retainProjectionCenters = false
        #expect(throws: (any Error).self) { try plan.validate() }
    }
}
