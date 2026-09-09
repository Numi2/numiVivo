import Foundation
import Testing
@testable import NumiVivoKit

@Suite struct SingleCellAdaptiveIntegrationTests {
    @Test func heterogeneousRidgeMatchesIndependentDenseSolve() throws {
        let fit = try VivoSingleCellIntegration.ridgeFit(masses: [2,3,5], sums: [[4,2],[9,-3],[0,10]], penalties: [0.4,1.2,2.5])
        // NumPy solve of the full intercept + three indicator normal system.
        let expected = [[1.1333333333333335,0.9833333333333336],[0.7222222222222221,0.013888888888888689],[1.333333333333333,-1.4166666666666667],[-0.7555555555555556,0.6777777777777776]]
        for (row, oracle) in zip([fit.intercept] + fit.effects, expected) {
            for (a,b) in zip(row,oracle) { #expect(abs(a-b) < 1e-12) }
        }
        #expect(fit.residual < 1e-12)
        let constant = try VivoSingleCellIntegration.ridgeFit(masses: [2,3,5], sums: [[14],[21],[35]], penalties: [0.4,1.2,2.5])
        #expect(abs(constant.intercept[0] - 7) < 1e-12)
        #expect(constant.effects.allSatisfy { abs($0[0]) < 1e-12 })
        for penalties in [[0.4,0,2.5],[0.4,-1,2.5],[0.4,.infinity,2.5],[0.4,1.2]] {
            #expect(throws: (any Error).self) { try VivoSingleCellIntegration.ridgeFit(masses: [2,3,5], sums: [[4],[9],[0]], penalties: penalties) }
        }
    }
    @Test func expectedMassPenaltiesAndMappedTrajectoryAreExplicit() throws {
        let input = SingleCellIntegrationTests.fixture()
        let samples = (0..<2).map { VivoOmicsSample(id: "s\($0)", biologicalReplicateID: "d\($0)", donorID: "d\($0)", condition: "shared", batchID: "unreported", organism: "human") }
        let root = FileManager.default.temporaryDirectory.appendingPathComponent(UUID().uuidString)
        try FileManager.default.createDirectory(at: root, withIntermediateDirectories: false)
        defer { try? FileManager.default.removeItem(at: root) }
        for seed: UInt64 in [7,19,41] {
            var options = VivoSingleCellIntegrationOptions(); options.clusters = 3; options.seed = seed
            options.ridge = 0.2; options.ridgeScaling = .expectedClusterBatchMass
            let resident = try VivoSingleCellIntegration.run(input, samples: samples, options: options)
            let source = try VivoIntegrationMatrix(rows: 40, columns: 3, scratch: root)
            for i in 0..<40 { try source.setRow(i, input.scores[i]) }
            let file = try VivoSingleCellIntegration.run(cells: input.cells, x: source, samples: samples, options: options) {
                try VivoIntegrationMatrix(rows: $0, columns: $1, scratch: root)
            }
            #expect(try file.materialize(cells: input.cells, options: options) == resident)
            let penalties = try #require(resident.ridgePenalties)
            #expect(penalties.count == 3 && penalties.allSatisfy { $0.count == 2 })
            for c in 0..<3 {
                let mass = resident.memberships.reduce(0) { $0 + $1[c] }
                for b in 0..<2 { #expect(abs(penalties[c][b] - 0.2 * mass * 20 / 40) < 1e-12) }
            }
            #expect(resident.method == "diversity-soft-clustering-expected-mass-ridge-Double-v1")
        }
    }
    @Test func legacyDefaultsAndUnknownModesRemainUnambiguous() throws {
        let options = VivoSingleCellIntegrationOptions()
        let raw = try VivoCanonicalJSON.encode(options)
        #expect(!String(decoding: raw, as: UTF8.self).contains("ridgeScaling"))
        #expect(try VivoCanonicalJSON.decode(VivoSingleCellIntegrationOptions.self, from: raw) == options)
        #expect(throws: (any Error).self) { try VivoCanonicalJSON.decode(VivoSingleCellIntegrationOptions.self, from: Data("{\"ridgeScaling\":\"automatic\"}".utf8)) }
        var adaptive = options; adaptive.ridgeScaling = .expectedClusterBatchMass; adaptive.ridge = 0
        #expect(throws: (any Error).self) { try adaptive.validate() }
    }
}
