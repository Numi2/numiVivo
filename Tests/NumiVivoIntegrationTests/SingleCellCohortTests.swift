import Foundation
import Testing
@testable import NumiVivoKit

@Suite struct SingleCellCohortTests {
    static func fixture(confoundedBatch: Bool = false) throws -> VivoSingleCellDataset {
        var samples: [VivoOmicsSample] = [], cells: [VivoOmicsCell] = []
        var offsets = [0], columns: [Int] = [], counts: [UInt64] = []
        let features = (0..<32).map { VivoOmicsFeature(id: "g\($0)", name: "gene\($0)", mitochondrial: $0 == 31) }
        for donor in 0..<6 { for condition in 0..<2 {
            let name = condition == 0 ? "control" : "treated", sample = "s\(donor)-\(condition)"
            samples.append(.init(id: sample, biologicalReplicateID: "r\(donor)", donorID: "d\(donor)",
                condition: name, batchID: confoundedBatch ? name : "shared", organism: "synthetic-organism"))
            for cell in 0..<3 {
                cells.append(.init(barcode: "c\(cell)", sampleID: sample, group: "declared-group"))
                if cell < 2 { for gene in 0..<32 {
                    let raw = 50 + gene * 3 + donor * 7 + cell * 2 + (donor * 11 + gene * 5 + condition * gene * 3 + cell * 7) % 17
                    columns.append(gene); counts.append(UInt64(raw * (condition == 1 && gene < 3 ? 4 : 1)))
                } }
                offsets.append(counts.count)
            }
        } }
        let value = VivoSingleCellDataset(id: "synthetic-paired", evidence: .synthetic, sourceDescription: "Synthetic test fixture, not measurements",
            countUnit: .umiCount, samples: samples, features: features, cells: cells,
            matrix: .init(cellCount: cells.count, featureCount: features.count, rowOffsets: offsets, featureIndices: columns, counts: counts))
        try value.validate(); return value
    }
    static func contrast() -> VivoOmicsExpressionContrast {
        var value = VivoOmicsExpressionContrast(id: "paired", controlCondition: "control", treatmentCondition: "treated",
            cellGroup: "declared-group", design: .pairedDonors)
        value.minimumCellsPerPseudobulk = 2
        return value
    }
    @Test func qualitySelectionPreservesCountsAndSourceIdentity() throws {
        let source = try Self.fixture(), processed = try VivoSingleCellProcessing.run(source)
        #expect(source.cells.count == 36 && processed.dataset.cells.count == 24)
        #expect(processed.decisions.filter { !$0.accepted }.count == 12)
        #expect(processed.sourceCellIndices == source.cells.indices.filter { $0 % 3 != 2 })
        #expect(processed.dataset.features == source.features)
        #expect(processed.dataset.matrix.counts == source.matrix.counts)
        #expect(processed.features.allSatisfy { $0.detectedCells == 24 && ($0.varianceLogNormalized ?? -1) >= 0 })
    }
    @Test func pairedModeratedExpressionUsesDonorsAndRecoversDeclaredSignal() throws {
        let processed = try VivoSingleCellProcessing.run(Self.fixture())
        let result = try VivoPseudobulkDifferentialExpression.run(processed.dataset, contrast: Self.contrast())
        #expect(result.design.controlReplicates == 6 && result.design.treatmentReplicates == 6)
        #expect(result.design.residualDegreesOfFreedom == 5)
        #expect(result.testedFeatures == 32 && result.variancePrior != nil)
        #expect(result.features.prefix(3).allSatisfy { ($0.log2FoldChange ?? 0) > 1.5 })
        #expect(result.features.allSatisfy { (0...1).contains($0.adjustedPValue ?? -1) })
        var wrong = Self.contrast(); wrong.design = .independentReplicates
        #expect(throws: (any Error).self) { try VivoPseudobulkDifferentialExpression.run(processed.dataset, contrast: wrong) }
    }
    @Test func incompletePairsAndConfoundedBatchAreRejected() throws {
        let source = try Self.fixture()
        var filter = VivoSingleCellFilterPolicy(); filter.includedSampleIDs = source.samples.map(\.id).filter { $0 != "s0-1" }
        let incomplete = try VivoSingleCellProcessing.run(source, policy: filter)
        #expect(throws: (any Error).self) { try VivoPseudobulkDifferentialExpression.run(incomplete.dataset, contrast: Self.contrast()) }
        #expect(throws: (any Error).self) { try VivoPseudobulkDifferentialExpression.run(Self.fixture(confoundedBatch: true), contrast: Self.contrast()) }
    }
    @Test func qrStudentAndMultipleTestingMatchAnalyticalFixtures() throws {
        let design = [[1.0, 0], [1, 0], [1, 0], [1, 1], [1, 1], [1, 1]]
        let fit = try VivoOmicsQR(design: design).fit([1, 2, 3, 4, 5, 6], contrast: [0, 1])
        #expect(abs(fit.effect - 3) < 1e-12)
        #expect(abs(fit.residualVariance - 1) < 1e-12)
        #expect(abs(fit.contrastVarianceScale - 2.0 / 3) < 1e-12)
        #expect(abs(try VivoOmicsLinearStatistics.studentTwoSidedP(t: 2.7764451051977987, degreesOfFreedom: 4) - 0.05) < 1e-10)
        #expect(abs(try VivoOmicsLinearStatistics.studentCriticalValue(degreesOfFreedom: 4) - 2.7764451051977987) < 1e-9)
        let q = try VivoOmicsLinearStatistics.benjaminiHochberg([0.01, 0.04, 0.03, 0.8])
        #expect(zip(q, [0.04, 0.05333333333333334, 0.05333333333333334, 0.8]).allSatisfy { abs($0 - $1) < 1e-14 })
    }
    @Test func unknownScientificOptionsCannotBeIgnored() throws {
        let bytes = Data("{\"maximumMitochondrialFractoin\":0.1}".utf8)
        #expect(throws: (any Error).self) { try JSONDecoder().decode(VivoSingleCellFilterPolicy.self, from: bytes) }
    }
}
