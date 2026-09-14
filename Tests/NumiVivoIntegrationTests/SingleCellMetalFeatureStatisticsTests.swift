import Foundation
import Testing
@testable import NumiVivoKit

@Suite struct SingleCellMetalFeatureStatisticsTests {
    @Test(.enabled(if: ProcessInfo.processInfo.environment["NUMIVIVO_TEST_METAL"] == "1"))
    func sparseMomentsMatchCPUAndReductionUsesExplicitMetalProfile() throws {
        let data = try SingleCellCohortTests.fixture()
        let processed = try VivoSingleCellProcessing.run(data)
        var seen = [Int](repeating: 0, count: processed.dataset.features.count)
        var means = [Double](repeating: 0, count: seen.count)
        var m2 = means
        for k in processed.normalized.values.indices {
            let feature = processed.normalized.featureIndices[k]
            let value = expm1(processed.normalized.values[k])
            seen[feature] += 1
            let delta = value - means[feature]
            means[feature] += delta / Double(seen[feature])
            m2[feature] += delta * (value - means[feature])
        }
        let metal = try VivoMetalFeatureStatistics.run(values: processed.normalized.values,
            featureIndices: processed.normalized.featureIndices, featureCount: processed.dataset.features.count)
        #expect(metal.seen == seen)
        for feature in seen.indices where seen[feature] > 0 {
            #expect(abs(metal.means[feature] - means[feature]) <= max(2e-4, abs(means[feature]) * 2e-5))
            #expect(abs(metal.m2[feature] - m2[feature]) <= max(5e-3, abs(m2[feature]) * 5e-5))
        }

        var options = VivoSingleCellReductionOptions()
        options.highlyVariableFeatures = 32
        options.meanBins = 2
        options.components = 8
        options.maximumBasis = 32
        options.featureStatisticsBackend = .metalFP32
        let report = try VivoSingleCellCohortAnalysis.run(data,
            plan: .init(id: "metal-feature-statistics", reduction: options))
        let reduction = try #require(report.reduction)
        #expect(reduction.method == "sparse-seurat-dispersion-centered-krylov-PCA-v1+metal-feature-statistics")
        #expect(reduction.qualification.contains("opt-in Metal FP32 sparse feature-statistics profile"))
        #expect(reduction.selectedFeatureIndices.count == options.highlyVariableFeatures)
    }

    @Test func backendIsOptionalAndUnknownReductionKeysRemainRejected() throws {
        let options = VivoSingleCellReductionOptions()
        #expect(options.featureStatisticsBackend == nil)
        #expect(throws: (any Error).self) {
            try JSONDecoder().decode(VivoSingleCellReductionOptions.self,
                from: Data(#"{"featureStatisticsBacken":"metalFP32"}"#.utf8))
        }
        let encoded = try VivoCanonicalJSON.encode(options)
        #expect(!String(decoding: encoded, as: UTF8.self).contains("featureStatisticsBackend"))
    }
}
