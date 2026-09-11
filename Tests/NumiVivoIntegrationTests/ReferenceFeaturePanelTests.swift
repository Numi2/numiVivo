import Foundation
import Testing
@testable import NumiVivoKit

@Suite struct ReferenceFeaturePanelTests {
    @Test func fullPanelPreservesLegacyNumericsAndEncodingIsOptional() throws {
        let processed = try VivoSingleCellProcessing.run(VivoSingleCellExamples.pairedCounts())
        var options = VivoSingleCellReductionOptions()
        options.components = 2; options.highlyVariableFeatures = 4; options.maximumBasis = 32
        let original = try VivoSingleCellReduction.run(processed, options: options)
        #expect(!String(decoding: try JSONEncoder().encode(options), as: UTF8.self).contains("featurePanel"))
        options.featurePanel = processed.dataset.features.map(\.id).reversed()
        let explicit = try VivoSingleCellReduction.run(processed, options: options)
        #expect(explicit.features == original.features)
        #expect(explicit.scores == original.scores && explicit.loadings == original.loadings)
        #expect(explicit.explainedVariance == original.explainedVariance)
    }

    @Test func subsetSelectionKeepsFullSourceNormalization() throws {
        let processed = try VivoSingleCellProcessing.run(VivoSingleCellExamples.pairedCounts())
        let data = processed.dataset
        var options = VivoSingleCellReductionOptions()
        options.components = 2; options.highlyVariableFeatures = 12; options.maximumBasis = 16
        options.featurePanel = (0..<16).map { "g\($0)" }; options.retainProjectionCenters = true
        let result = try VivoSingleCellReduction.run(processed, options: options)
        #expect(result.selectedFeatureIndices.allSatisfy { $0 < 16 })
        for item in result.features where item.featureIndex >= 16 {
            #expect(item.meanBin == -1 && item.normalizedDispersion == nil && !item.selected)
        }
        let centers = try #require(result.projectionCenters)
        for (j, source) in result.selectedFeatureIndices.enumerated() {
            var fullMean = 0.0, restrictedMean = 0.0
            for row in data.cells.indices {
                let range = data.matrix.rowOffsets[row]..<data.matrix.rowOffsets[row + 1]
                let total = range.reduce(UInt64(0)) { $0 + data.matrix.counts[$1] }
                let subset = range.filter { data.matrix.featureIndices[$0] < 16 }.reduce(UInt64(0)) { $0 + data.matrix.counts[$1] }
                let count = range.filter { data.matrix.featureIndices[$0] == source }.reduce(UInt64(0)) { $0 + data.matrix.counts[$1] }
                fullMean += log1p(Double(count) / Double(total) * 10_000) / Double(data.cells.count)
                restrictedMean += log1p(Double(count) / Double(subset) * 10_000) / Double(data.cells.count)
            }
            #expect(abs(centers[j] - fullMean) < 1e-12)
            #expect(abs(centers[j] - restrictedMean) > 0.1)
        }
    }

    @Test func invalidAndUnmeasuredPanelsReject() throws {
        let processed = try VivoSingleCellProcessing.run(VivoSingleCellExamples.pairedCounts())
        var options = VivoSingleCellReductionOptions(); options.components = 2
        for panel in [[], ["g1"], ["g1", "g1"], ["g1", ""], ["g1", "unmeasured"]] {
            options.featurePanel = panel
            #expect(throws: (any Error).self) { try VivoSingleCellReduction.run(processed, options: options) }
        }
        #expect(throws: (any Error).self) {
            try JSONDecoder().decode(VivoSingleCellReductionOptions.self, from: Data(#"{"featurePannel":["g1","g2"]}"#.utf8))
        }
    }
}
