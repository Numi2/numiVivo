import Foundation
import Testing
@testable import NumiVivoKit

@Suite struct SingleCellPipelineTests {
    @Test func cohortPlanRunsNBReductionGraphAndIntegrationAsOneReport() throws {
        let data = try SingleCellNBCohortTests.fixture()
        var contrast = VivoOmicsExpressionContrast(id: "pipeline-nb", controlCondition: "ctrl",
            treatmentCondition: "stim", design: .pairedDonors)
        contrast.model = .negativeBinomial
        contrast.minimumCellsPerPseudobulk = 1
        var nb = VivoOmicsNBCohortOptions()
        nb.trend = .mean
        contrast.negativeBinomialOptions = nb

        var reduction = VivoSingleCellReductionOptions()
        reduction.highlyVariableFeatures = 24
        reduction.meanBins = 4
        reduction.components = 4
        reduction.maximumBasis = 24
        reduction.relativeResidualTolerance = 1e-4

        var integration = VivoSingleCellIntegrationOptions()
        integration.clusters = 3
        integration.maximumIterations = 2

        var neighbors = VivoSingleCellNeighborOptions()
        neighbors.neighbors = 5
        neighbors.representation = .integrated

        var clustering = VivoSingleCellClusteringOptions()
        clustering.maximumSweeps = 20
        var embedding = VivoSingleCellEmbeddingOptions()
        embedding.epochs = 20
        embedding.negativeSampleRate = 1

        let plan = VivoSingleCellAnalysisPlan(id: "complete-native-pipeline", contrasts: [contrast],
            reduction: reduction, neighbors: neighbors, clustering: clustering,
            embedding: embedding, integration: integration)
        let report = try VivoSingleCellCohortAnalysis.run(data, plan: plan)
        let result = try #require(report.contrasts.first)
        #expect(result.negativeBinomial != nil)
        #expect(report.reduction?.scores.count == report.processed.dataset.cells.count)
        #expect(report.integration?.scores.count == report.processed.dataset.cells.count)
        #expect(report.neighbors?.neighborIndices.count == report.processed.dataset.cells.count * neighbors.neighbors)
        #expect(report.clustering?.labels.count == report.processed.dataset.cells.count)
        #expect(report.embedding?.coordinates.count == report.processed.dataset.cells.count)
        #expect(report.neighbors?.options.representation == .integrated)
        #expect(report.integration?.maximumRidgeResidual ?? .infinity < 1e-10)
        #expect(try VivoSingleCellCohortAnalysis.run(data, plan: plan) == report)
    }
}
