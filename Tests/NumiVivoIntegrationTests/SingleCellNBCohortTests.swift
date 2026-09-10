import Foundation
import Testing
@testable import NumiVivoKit

@Suite struct SingleCellNBCohortTests {
    @Test func singletonBatchDoesNotRequireUnrequestedInfluenceDiagnostics() throws {
        let source = try Self.fixture()
        let samples = source.samples.enumerated().map { i,s in
            VivoOmicsSample(id: s.id,biologicalReplicateID: s.id,donorID: s.id,
                condition: s.condition,batchID: i == 0 ? "singleton" : "shared",organism: s.organism)
        }
        let data = VivoSingleCellDataset(id: source.id,evidence: source.evidence,sourceDescription: source.sourceDescription,
            countUnit: source.countUnit,samples: samples,features: source.features,cells: source.cells,matrix: source.matrix)
        var request = VivoOmicsExpressionContrast(id: "singleton",controlCondition: "ctrl",treatmentCondition: "stim",design: .independentReplicates)
        request.model = .negativeBinomial;request.minimumCellsPerPseudobulk = 1;request.adjustForBatch = true
        let methods: [VivoOmicsNBTestMethod?] = [nil,.likelihoodRatio]
        for method in methods {
            var options = VivoOmicsNBCohortOptions();options.trend = .mean;options.testMethod = method
            request.negativeBinomialOptions = options
            let result = try VivoPseudobulkDifferentialExpression.run(data,contrast: request)
            let diagnostics = try #require(result.negativeBinomial?.features)
            let missing = result.features.filter { $0.status == .tested && diagnostics[$0.featureIndex].finalFit?.cooksDistances == nil }
            #expect(missing.count > 50)
            #expect(missing.allSatisfy { $0.pValue != nil && $0.adjustedPValue != nil })
            #expect(result.design.residualDegreesOfFreedom == 9)
            options.maximumCooksDistance = 1e9;request.negativeBinomialOptions = options
            let requested = try VivoPseudobulkDifferentialExpression.run(data,contrast: request)
            for f in missing {
                #expect(requested.features[f.featureIndex].status == .numericalFailure)
                #expect(requested.features[f.featureIndex].pValue == nil)
                #expect(requested.negativeBinomial?.features[f.featureIndex].finalFit == diagnostics[f.featureIndex].finalFit)
            }
        }
    }
    @Test func explicitLikelihoodRatioUsesSameDispersionAndOwnBHFamily() throws {
        let data = try Self.fixture()
        var request = VivoOmicsExpressionContrast(id: "lrt",controlCondition: "ctrl",treatmentCondition: "stim",design: .pairedDonors)
        request.model = .negativeBinomial; request.minimumCellsPerPseudobulk = 1
        var options = VivoOmicsNBCohortOptions(); options.trend = .mean
        request.negativeBinomialOptions = options
        let wald = try VivoPseudobulkDifferentialExpression.run(data,contrast: request)
        let encoded = try JSONEncoder().encode(options)
        #expect(!(String(data: encoded,encoding: .utf8) ?? "").contains("testMethod"))
        options.testMethod = .likelihoodRatio; request.negativeBinomialOptions = options
        #expect(try JSONDecoder().decode(VivoOmicsNBCohortOptions.self,from: JSONEncoder().encode(options)) == options)
        let lrt = try VivoPseudobulkDifferentialExpression.run(data,contrast: request)
        #expect(lrt.method.hasSuffix("LRT-v1") && lrt.testedFeatures == wald.testedFeatures)
        #expect(lrt.negativeBinomial?.trend == wald.negativeBinomial?.trend)
        let diagnostics = try #require(lrt.negativeBinomial?.features), old = try #require(wald.negativeBinomial?.features)
        for (a,b) in zip(diagnostics,old) { #expect(a.finalFit == b.finalFit) }
        let tested = lrt.features.filter { $0.status == .tested }
        let bh = try VivoOmicsLinearStatistics.benjaminiHochberg(tested.map { $0.pValue! })
        for (i,f) in tested.enumerated() {
            #expect(f.pValue == diagnostics[f.featureIndex].likelihoodRatioFit?.pValue)
            #expect(f.adjustedPValue == bh[i])
            #expect(f.zStatistic == wald.features[f.featureIndex].zStatistic)
            #expect(f.intervalLower == wald.features[f.featureIndex].intervalLower)
        }
        options.maximumCooksDistance = 1e-15; request.negativeBinomialOptions = options
        let excluded = try VivoPseudobulkDifferentialExpression.run(data,contrast: request)
        #expect(excluded.testedFeatures == 0)
        #expect(excluded.negativeBinomial!.features.allSatisfy { $0.likelihoodRatioFit == nil })
    }
    @Test func gammaTrendRecoversMeanAndRejectsUnidentifiedDesign() throws {
        let means=(0..<64).map { exp(Double($0)/10+1) }
        let values=means.map { 0.04+2/$0 }
        var options=VivoOmicsNBCohortOptions();options.trend = .gammaParametric
        let fit=try VivoOmicsNBCohort.fitTrend(means: means,dispersions: values,featureIndices: Array(0..<64),residualDF: 7,options: options)
        #expect(abs(fit.intercept-0.04)<1e-7)
        #expect(abs(fit.inverseMeanCoefficient-2)<1e-6)
        #expect(fit.trendFitFeatureIndices == Array(0..<64))
        #expect(throws: (any Error).self) { try VivoOmicsNBCohort.fitTrend(means: Array(repeating: 10,count: 64),dispersions: values,featureIndices: Array(0..<64),residualDF: 7,options: options) }
    }
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
        options.effectPriorStandardDeviationLog2 = 1; contrast.negativeBinomialOptions = options
        let shrunk = try VivoPseudobulkDifferentialExpression.run(data,contrast: contrast)
        #expect(shrunk.features == result.features)
        #expect(shrunk.negativeBinomial?.trend == result.negativeBinomial?.trend)
        #expect(shrunk.negativeBinomial?.effectShrinkage?.convergedFeatures == result.testedFeatures)
        #expect(shrunk.negativeBinomial?.effectShrinkage?.failedFeatures == 0)
        #expect(result.negativeBinomial?.effectShrinkage == nil)
        var empirical = options; empirical.effectPriorStandardDeviationLog2 = nil
        empirical.effectPriorEstimation = .weightedUpperQuantile; contrast.negativeBinomialOptions = empirical
        let learned = try VivoPseudobulkDifferentialExpression.run(data,contrast: contrast)
        #expect(learned.features == result.features)
        #expect(learned.negativeBinomial?.trend == result.negativeBinomial?.trend)
        #expect(learned.negativeBinomial?.effectPriorEstimationError == nil)
        let estimate = try #require(learned.negativeBinomial?.effectPriorEstimate)
        #expect(estimate.referenceFeatureIndices.count >= 20)
        #expect(learned.negativeBinomial?.effectShrinkage?.priorStandardDeviationLog2 == estimate.priorStandardDeviationLog2)
        #expect(learned.negativeBinomial?.effectShrinkage?.failedFeatures == 0)
        empirical.maximumCooksDistance = 1e-15; contrast.negativeBinomialOptions = empirical
        let unavailable = try VivoPseudobulkDifferentialExpression.run(data,contrast: contrast)
        #expect(unavailable.testedFeatures == 0)
        #expect(unavailable.negativeBinomial?.effectPriorEstimate == nil)
        #expect(unavailable.negativeBinomial?.effectPriorEstimationError != nil)
        #expect(unavailable.negativeBinomial?.effectShrinkage == nil)
        options.maximumCooksDistance = 1e-15; contrast.negativeBinomialOptions = options
        let excluded = try VivoPseudobulkDifferentialExpression.run(data,contrast: contrast)
        #expect(excluded.testedFeatures == 0)
        #expect(excluded.features.allSatisfy { $0.pValue == nil && $0.adjustedPValue == nil })
        #expect(excluded.negativeBinomial?.effectShrinkage?.eligibleFeatures == 0)
        #expect(excluded.negativeBinomial!.features.allSatisfy { $0.effectShrinkageFit == nil })
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
        let decoded = try JSONDecoder().decode(VivoOmicsNBCohortOptions.self,from: Data("{\"effectPriorStandardDeviationLog2\":1}".utf8))
        #expect(decoded.effectPriorStandardDeviationLog2 == 1)
        #expect(try JSONDecoder().decode(VivoOmicsNBCohortOptions.self,from: JSONEncoder().encode(decoded)) == decoded)
        for sd in [0.0,-1,0.001,101,.infinity,.nan] {
            var invalid = decoded; invalid.effectPriorStandardDeviationLog2 = sd
            #expect(throws: (any Error).self) { try invalid.validate() }
        }
        var both = decoded; both.effectPriorEstimation = .weightedUpperQuantile
        #expect(throws: (any Error).self) { try both.validate() }
        let learned = try JSONDecoder().decode(VivoOmicsNBCohortOptions.self,from: Data("{\"effectPriorEstimation\":\"weightedUpperQuantile\"}".utf8))
        #expect(try JSONDecoder().decode(VivoOmicsNBCohortOptions.self,from: JSONEncoder().encode(learned)) == learned)
    }
    @Test func empiricalPriorUsesWeightedQuantileAndRetainsDomainBoundaries() throws {
        let x = (0..<20).map { Double($0)/10 }
        func estimate(_ x: [Double], _ means: [Double]) throws -> VivoOmicsNBEffectPriorEstimate {
            try VivoOmicsNBCohort.estimateEffectPrior(effectsLog2: x,means: means,
                trendDispersions: Array(repeating: 0.1,count: x.count),featureIndices: Array(x.indices))
        }
        let equal = try estimate(x,Array(repeating: 100,count: 20))
        #expect(abs(equal.absoluteEffectQuantileLog2-1.805) < 1e-12)
        #expect(abs(equal.effectiveReferenceFeatures-20) < 1e-12)
        let reversed = try estimate(x.reversed().map { -$0 },Array(repeating: 100,count: 20))
        #expect(abs(equal.priorStandardDeviationLog2-reversed.priorStandardDeviationLog2) < 1e-12)
        let lowCountOutlier = try estimate(x+[9],Array(repeating: 100,count: 20)+[1e-5])
        #expect(lowCountOutlier.absoluteEffectQuantileLog2 < 2)
        let excluded = try estimate(x+[10],Array(repeating: 100,count: 21))
        #expect(excluded.excludedHighEffectFeatureIndices == [20])
        #expect(excluded.priorStandardDeviationLog2 == equal.priorStandardDeviationLog2)
        let zero = try estimate(Array(repeating: 0,count: 20),Array(repeating: 100,count: 20))
        #expect(zero.standardDeviationFloorReached && zero.priorStandardDeviationLog2 == 0.01)
        #expect(throws: (any Error).self) { try estimate(Array(x.prefix(19)),Array(repeating: 100,count: 19)) }
        #expect(throws: (any Error).self) { try estimate(x,Array(repeating: .nan,count: 20)) }
    }
}
