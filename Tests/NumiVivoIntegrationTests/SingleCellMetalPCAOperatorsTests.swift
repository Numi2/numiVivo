import Foundation
import Testing
@testable import NumiVivoKit

@Suite struct SingleCellMetalPCAOperatorsTests {
    @Test(.enabled(if: ProcessInfo.processInfo.environment["NUMIVIVO_TEST_METAL"] == "1"))
    func sparseProjectAndTransposeMatchCPUOracle() throws {
        let rowOffsets = [0, 2, 3, 5, 7]
        let columns = [0, 2, 1, 0, 2, 0, 1]
        let values = [0.1, 0.4, 0.2, 0.3, 0.5, 0.7, 0.9]
        let operators = try VivoMetalSparsePCAOperators(rowOffsets: rowOffsets, columns: columns,
            values: values, rows: 4, columnCount: 3)
        let vector = [0.25, -0.5, 0.75], shift = 0.2
        let project = try operators.project(vector, shift: shift)
        let expectedProject = (0..<4).map { row in
            values[rowOffsets[row]..<rowOffsets[row + 1]].enumerated().reduce(-shift) { total, item in
                total + item.element * vector[columns[rowOffsets[row] + item.offset]]
            }
        }
        #expect(zip(project, expectedProject).allSatisfy { abs($0 - $1) <= 1e-5 })

        let projected = [0.5, -0.25, 0.75, 0.125], initial = [0.4, -0.2, 0.8]
        let transpose = try operators.transpose(projected, initial: initial)
        var expectedTranspose = initial
        for row in 0..<4 {
            for index in rowOffsets[row]..<rowOffsets[row + 1] {
                expectedTranspose[columns[index]] += values[index] * projected[row]
            }
        }
        #expect(zip(transpose, expectedTranspose).allSatisfy { abs($0 - $1) <= 1e-5 })
    }

    @Test(.enabled(if: ProcessInfo.processInfo.environment["NUMIVIVO_TEST_METAL"] == "1"))
    func residentReductionRecordsOptInMetalOperators() throws {
        let data = try SingleCellCohortTests.fixture()
        var options = VivoSingleCellReductionOptions()
        options.highlyVariableFeatures = 32
        options.meanBins = 2
        options.components = 8
        options.maximumBasis = 32
        options.relativeResidualTolerance = 2e-4
        options.featureStatisticsBackend = .metalFP32
        options.pcaOperatorsBackend = .metalFP32
        let report = try VivoSingleCellCohortAnalysis.run(data,
            plan: .init(id: "metal-pca-operators", reduction: options))
        let reduction = try #require(report.reduction)
        #expect(reduction.method == "sparse-seurat-dispersion-centered-krylov-PCA-v1+metal-feature-statistics+metal-pca-operators")
        #expect(reduction.relativeResiduals.allSatisfy { $0 <= options.relativeResidualTolerance })
        #expect(reduction.scores.flatMap { $0 }.allSatisfy { $0.isFinite })
        #expect(reduction.qualification.contains("opt-in Metal PCA-operators profile"))
    }

    @Test func pcaOperatorBackendIsOptionalAndUnknownKeysRemainRejected() throws {
        let options = VivoSingleCellReductionOptions()
        #expect(options.pcaOperatorsBackend == nil)
        #expect(throws: (any Error).self) {
            try JSONDecoder().decode(VivoSingleCellReductionOptions.self,
                from: Data(#"{"pcaOperatorsBacken":"metalFP32"}"#.utf8))
        }
    }
}
