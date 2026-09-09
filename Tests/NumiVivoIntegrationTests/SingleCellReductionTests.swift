import Foundation
import Testing
@testable import NumiVivoKit

@Suite struct SingleCellReductionTests {
    @Test func centeredScoresRecoverVarianceAndRepeatExactly() throws {
        let data = try SingleCellCohortTests.fixture()
        var options = VivoSingleCellReductionOptions()
        options.highlyVariableFeatures = 32; options.meanBins = 2
        options.components = 8; options.maximumBasis = 32
        let plan = VivoSingleCellAnalysisPlan(id: "reduction-test", reduction: options)
        let report = try VivoSingleCellCohortAnalysis.run(data, plan: plan)
        let result = try #require(report.reduction)
        #expect(result == (try VivoSingleCellCohortAnalysis.run(data, plan: plan)).reduction)
        #expect(result.cells.count == 24 && result.selectedFeatureIndices.count == 32)
        #expect(result.maximumLoadingOrthogonalityError < 1e-10)
        for j in 0..<options.components {
            let scores = result.scores.map { $0[j] }
            #expect(abs(scores.reduce(0,+)) < 1e-9)
            let variance = scores.reduce(0) { $0 + $1*$1 } / Double(scores.count-1)
            #expect(abs(variance-result.explainedVariance[j]) < 1e-9)
            #expect(result.relativeResiduals[j] < options.relativeResidualTolerance)
        }
        // This fixture has more selected genes than retained cells. A random
        // feature-space start contains nullspace; the full basis must handle it.
        #expect(result.basisSize > result.cells.count)
        var insufficient = options; insufficient.maximumBasis = options.components
        #expect(throws: (any Error).self) {
            try VivoSingleCellCohortAnalysis.run(data, plan: .init(id: "short-basis", reduction: insufficient))
        }
    }
    @Test func undefinedDispersionAndRepeatedProjectedEigenvalues() throws {
        let cells = (0..<4).map { VivoOmicsCell(barcode: "c\($0)", sampleID: "s", group: "test") }
        let data = VivoSingleCellDataset(id: "repeated-spectrum", evidence: .synthetic,
            sourceDescription: "Analytical repeated-eigenvalue fixture, not experimental data", countUnit: .umiCount,
            samples: [.init(id: "s", biologicalReplicateID: "r", donorID: "d", condition: "test", batchID: "b", organism: "synthetic")],
            features: (0..<4).map { .init(id: "g\($0)", name: "g\($0)", mitochondrial: false) }, cells: cells,
            matrix: .init(cellCount: 4, featureCount: 4, rowOffsets: [0,1,2,3,4], featureIndices: [0,1,2,3], counts: [10,10,10,10]))
        // Equal gene dispersions are deliberately undefined under Seurat bin
        // normalization, so this must reject rather than invent an HVG ranking.
        var options = VivoSingleCellReductionOptions(); options.components = 2; options.highlyVariableFeatures = 4
        #expect(throws: (any Error).self) {
            try VivoSingleCellCohortAnalysis.run(data, plan: .init(id: "undefined-hvg", reduction: options))
        }
        let eigen = try VivoSingleCellReduction.symmetricEigen([2,0,0,0,2,0,0,0,0], n: 3)
        #expect(eigen.order.map { eigen.values[$0] } == [2,2,0])
    }
    @Test func fullSparseBasisRecoversRepeatedSpectrum() throws {
        var cells: [VivoOmicsCell] = [], counts: [UInt64] = [], columns: [Int] = [], offsets = [0]
        for a in 0..<4 { for b in 0..<4 {
            cells.append(.init(barcode: "c\(a)-\(b)", sampleID: "s", group: "test"))
            for j in 0..<8 {
                columns.append(j)
                counts.append(j < 4 ? (j == a ? 100 : 10) : (j-4 == b ? 200 : 20))
            }
            offsets.append(counts.count)
        } }
        let data = VivoSingleCellDataset(id: "repeated-covariance", evidence: .synthetic,
            sourceDescription: "Analytical independent categorical blocks, not experimental data", countUnit: .umiCount,
            samples: [.init(id: "s", biologicalReplicateID: "r", donorID: "d", condition: "test", batchID: "b", organism: "synthetic")],
            features: (0..<8).map { .init(id: "g\($0)", name: "g\($0)", mitochondrial: false) }, cells: cells,
            matrix: .init(cellCount: 16, featureCount: 8, rowOffsets: offsets, featureIndices: columns, counts: counts))
        var options = VivoSingleCellReductionOptions()
        options.highlyVariableFeatures = 8; options.meanBins = 1; options.components = 6; options.maximumBasis = 8
        let result = try #require(VivoSingleCellCohortAnalysis.run(data, plan: .init(id: "repeated", reduction: options)).reduction)
        let first = pow(log1p(100.0/390*10000)-log1p(10.0/390*10000),2)*4/15
        let second = pow(log1p(200.0/390*10000)-log1p(20.0/390*10000),2)*4/15
        let expected = ([Double](repeating: first,count: 3)+[Double](repeating: second,count: 3)).sorted(by: >)
        for j in expected.indices { #expect(abs(result.explainedVariance[j]-expected[j]) < 1e-10) }
        options.components = 7
        #expect(throws: (any Error).self) {
            try VivoSingleCellCohortAnalysis.run(data, plan: .init(id: "excess-rank", reduction: options))
        }
    }
    @Test func invalidOptionsAndLegacyEncoding() throws {
        #expect(throws: (any Error).self) {
            try JSONDecoder().decode(VivoSingleCellReductionOptions.self, from: Data("{\"componentz\":2}".utf8))
        }
        var options = VivoSingleCellReductionOptions(); options.components = 300
        #expect(throws: (any Error).self) { try options.validate() }
        let report = try VivoSingleCellCohortAnalysis.run(SingleCellCohortTests.fixture(), plan: .init(id: "legacy"))
        let object = try #require(JSONSerialization.jsonObject(with: VivoCanonicalJSON.encode(report)) as? [String: Any])
        #expect(object["reduction"] == nil)
        #expect((object["plan"] as? [String: Any])?["reduction"] == nil)
    }
}
