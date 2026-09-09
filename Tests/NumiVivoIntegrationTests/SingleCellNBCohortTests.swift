import Foundation
import Testing
@testable import NumiVivoKit

@Suite struct SingleCellNBCohortTests {
    @Test func parametricTrendRecoversDeclaredMeanRelationship() throws {
        let means = (0..<64).map { exp(Double($0)/10+1) }
        let dispersions = means.map { 0.04 + 2.0/$0 }
        let trend = try VivoOmicsNBCohort.fitTrend(means: means,dispersions: dispersions,
            featureIndices: Array(0..<64),residualDF: 7,options: .init())
        #expect(abs(trend.intercept-0.04) < 1e-6)
        #expect(abs(trend.inverseMeanCoefficient-2) < 1e-5)
        #expect(trend.priorVarianceFloorReached)
        #expect(trend.priorLogVariance == 0.25)
        #expect(throws: (any Error).self) { try VivoOmicsNBCohort.fitTrend(means: Array(repeating: 10,count: 64),dispersions: dispersions,featureIndices: Array(0..<64),residualDF: 7,options: .init()) }
    }
    static func fixture() throws -> VivoSingleCellDataset {
        var samples: [VivoOmicsSample] = [], cells: [VivoOmicsCell] = []
        var offsets = [0], indices: [Int] = [], counts: [UInt64] = []
        for donor in 0..<6 { for condition in 0..<2 {
            let row = donor*2+condition, id = "sample-\(row)"
            samples.append(.init(id: id,biologicalReplicateID: "d\(donor)",donorID: "d\(donor)",condition: condition == 0 ? "ctrl" : "stim",batchID: "shared",organism: "synthetic-organism"))
            cells.append(.init(barcode: "cell",sampleID: id))
            for gene in 0..<96 {
                let residual = exp(0.9*sin(Double((gene+1)*(row+3))*1.731))
                let base = Double(40+gene*8) * exp(Double(donor)/8) * residual
                let value = base * (gene < 3 && condition == 1 ? 8 : 1)
                indices.append(gene); counts.append(UInt64(value.rounded()))
            }
            offsets.append(counts.count)
        } }
        return .init(id: "nb-controls",evidence: .synthetic,sourceDescription: "Deterministic numerical regression, not experimental evidence",countUnit: .umiCount,
            samples: samples,features: (0..<96).map { .init(id: "g\($0)",name: "g\($0)") },cells: cells,
            matrix: .init(cellCount: 12,featureCount: 96,rowOffsets: offsets,featureIndices: indices,counts: counts))
    }
    @Test func cohortUsesNBAndRetainsInfluenceAndReplicationGates() throws {
        let data = try Self.fixture()
        var contrast = VivoOmicsExpressionContrast(id: "nb",controlCondition: "ctrl",treatmentCondition: "stim",design: .pairedDonors)
        contrast.model = .negativeBinomial; contrast.minimumCellsPerPseudobulk = 1
        var options = VivoOmicsNBCohortOptions(); options.trend = .mean
        contrast.negativeBinomialOptions = options
        let result = try VivoPseudobulkDifferentialExpression.run(data,contrast: contrast)
        #expect(result.negativeBinomial != nil && result.variancePrior == nil)
        #expect(result.testedFeatures > 80)
        #expect(result.features.prefix(3).allSatisfy { ($0.log2FoldChange ?? 0) > 2 })
        #expect(result.features.allSatisfy { $0.tStatistic == nil && $0.degreesOfFreedom == nil })
        #expect(result.features.filter { $0.status == .tested }.allSatisfy { $0.zStatistic != nil && $0.pValue != nil })
        options.maximumCooksDistance = 1e-15; contrast.negativeBinomialOptions = options
        let excluded = try VivoPseudobulkDifferentialExpression.run(data,contrast: contrast)
        #expect(excluded.testedFeatures == 0)
        #expect(excluded.features.allSatisfy { $0.pValue == nil && $0.adjustedPValue == nil })
        contrast.design = .independentReplicates
        #expect(throws: (any Error).self) { try VivoPseudobulkDifferentialExpression.run(data,contrast: contrast) }
    }
    @Test func optionsCannotBeSilentlyAppliedToWrongModel() throws {
        var contrast = SingleCellCohortTests.contrast()
        contrast.negativeBinomialOptions = .init()
        #expect(throws: (any Error).self) { try contrast.validate() }
        contrast.model = .negativeBinomial; contrast.priorCount = 1
        #expect(throws: (any Error).self) { try contrast.validate() }
        let bytes = Data("{\"trend\":\"mean\",\"fakeOption\":true}".utf8)
        #expect(throws: (any Error).self) { try JSONDecoder().decode(VivoOmicsNBCohortOptions.self,from: bytes) }
    }
}
