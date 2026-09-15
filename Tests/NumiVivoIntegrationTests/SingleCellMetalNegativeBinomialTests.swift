import Foundation
import Testing
@testable import NumiVivoKit

@Suite struct SingleCellMetalNegativeBinomialTests {
    @Test(.enabled(if: ProcessInfo.processInfo.environment["NUMIVIVO_TEST_METAL"] == "1"))
    func batchedVariableObjectiveMatchesFP64CPUOracle() throws {
        let counts: [UInt64] = [0, 3, 12, 7, 1, 19, 2, 0, 31, 5]
        let means = [0.25, 1.2, 4.5, 8.0, 2.1, 17.0, 0.75, 3.0, 29.0, 6.5]
        let dispersion = 0.2
        let metal = try VivoMetalNegativeBinomialLikelihood(counts: counts, dispersion: dispersion)
        let value = try metal.evaluate(means: means)
        let expected = zip(counts, means).reduce(0.0) { total, pair in
            let y = Double(pair.0), mu = pair.1
            return total + y * log(mu) - (y + 1 / dispersion) * log1p(dispersion * mu)
        }
        #expect(abs(value - expected) <= max(5e-5, abs(expected) * 5e-6))
    }

    @Test(.enabled(if: ProcessInfo.processInfo.environment["NUMIVIVO_TEST_METAL"] == "1"))
    func packedIndependentObjectivesMatchIndividualMetalAndFP64Oracle() throws {
        let inputs: [([UInt64], [Double], Double)] = [
            ([0, 3, 12, 7, 1, 19], [0.25, 1.2, 4.5, 8.0, 2.1, 17.0], 0.2),
            ([9, 2, 0, 31, 5, 15, 1, 8], [3.0, 0.75, 1.1, 29.0, 6.5, 11.0, 0.5, 4.2], 0.01),
            ([2, 11, 4, 18, 0, 6, 23, 9, 3], [0.4, 8.5, 2.0, 13.0, 1.0, 7.2, 19.0, 5.5, 2.5], 1.7)
        ]
        let requests = try inputs.map { counts, means, dispersion in
            try VivoMetalNegativeBinomialLikelihood.BatchRequest(counts: counts,
                means: means, dispersion: dispersion)
        }
        let packed = try VivoMetalNegativeBinomialLikelihood.evaluateBatch(requests)
        let replay = try VivoMetalNegativeBinomialLikelihood.evaluateBatch(requests)
        #expect(packed.count == inputs.count)
        #expect(packed == replay)
        for index in inputs.indices {
            let (counts, means, dispersion) = inputs[index]
            let individual = try VivoMetalNegativeBinomialLikelihood(counts: counts, dispersion: dispersion)
            let scalar = try individual.evaluate(means: means)
            let expected = zip(counts, means).reduce(0.0) { total, pair in
                let y = Double(pair.0), mu = pair.1
                return total + y * log(mu) - (y + 1 / dispersion) * log1p(dispersion * mu)
            }
            #expect(abs(packed[index] - scalar) <= max(5e-5, abs(scalar) * 5e-6))
            #expect(abs(packed[index] - expected) <= max(5e-5, abs(expected) * 5e-6))
        }
    }

    @Test(.enabled(if: ProcessInfo.processInfo.environment["NUMIVIVO_TEST_METAL"] == "1"))
    func optInFitUsesMetalObjectiveAndRetainsExactReport() throws {
        let counts: [UInt64] = [11, 20, 12, 55, 43, 61]
        let design = [[1.0, 0], [1, 0], [1, 0], [1, 1], [1, 1], [1, 1]]
        let offsets = [0.0, 0.2, 0.1, 0.3, 0.1, 0.2]
        let cpu = try VivoOmicsNegativeBinomial.fit(counts: counts, design: design, offsets: offsets,
            contrast: [0, 1], dispersion: 0.1)
        let metal = try VivoOmicsNegativeBinomial.fit(counts: counts, design: design, offsets: offsets,
            contrast: [0, 1], dispersion: 0.1, backend: .metalFP32)
        #expect(metal.backend == .metalFP32)
        #expect(metal.converged)
        #expect(abs((metal.effect ?? .nan) - (cpu.effect ?? .nan)) < 1e-4)
        #expect(abs(metal.logLikelihood - cpu.logLikelihood) < 1e-6)
        let encoded = try VivoCanonicalJSON.encode(metal)
        #expect(String(decoding: encoded, as: UTF8.self).contains(#""backend":"metalFP32""#))
    }

    @Test(.enabled(if: ProcessInfo.processInfo.environment["NUMIVIVO_TEST_METAL"] == "1"))
    func objectiveRejectsCountsOutsideFP32Contract() throws {
        #expect(throws: (any Error).self) {
            _ = try VivoMetalNegativeBinomialLikelihood(counts: [16_777_217], dispersion: 0.2)
        }
        #expect(throws: (any Error).self) {
            _ = try VivoMetalNegativeBinomialLikelihood.BatchRequest(counts: [1, 2],
                means: [1.0], dispersion: 0.2)
        }
    }

    @Test func backendIsAbsentFromDefaultFitEncoding() throws {
        let counts: [UInt64] = [2, 5, 10, 12, 20, 23]
        let design = counts.map { _ in [1.0] }
        let fit = try VivoOmicsNegativeBinomial.fit(counts: counts, design: design,
            offsets: Array(repeating: 0, count: counts.count), contrast: [1], dispersion: 0.2)
        #expect(fit.backend == nil)
        let encoded = try VivoCanonicalJSON.encode(fit)
        #expect(!String(decoding: encoded, as: UTF8.self).contains("backend"))
    }

    @Test func legacyExecutionProfileDecodesWithoutTheV2Audit() throws {
        let legacy = """
        {"schemaVersion":"numi.vivo.metal-nb-execution-profile.v1","cohortEvaluationWallClockNanoseconds":1,"fixedDispersionFitAttempts":2,"completedFixedDispersionFits":2,"fixedDispersionFitWallClockNanoseconds":3,"metalObjectiveAttempts":4,"completedMetalObjectiveEvaluations":4,"metalObjectiveWallClockNanoseconds":5,"metalObjectiveBatches":6,"metalObjectiveObservations":7,"fp64ToFP32ConversionNanoseconds":8,"sharedBufferWriteNanoseconds":9,"commandEncodingAndSubmissionNanoseconds":10,"commandCompletionNanoseconds":11,"cpuReductionNanoseconds":12,"cpuExactLikelihoodEvaluations":13,"cpuExactLikelihoodObservations":14,"cpuExactLikelihoodNanoseconds":15,"qualification":"legacy"}
        """
        let profile = try JSONDecoder().decode(VivoMetalNBExecutionProfile.self,
            from: Data(legacy.utf8))
        #expect(profile.schemaVersion == "numi.vivo.metal-nb-execution-profile.v1")
        #expect(profile.lineSearchCandidateEvaluations == nil)
        #expect(profile.metalAcceptedCPURejectedLineSearchCandidates == nil)
        #expect(profile.cpuAcceptedMetalRejectedLineSearchCandidates == nil)
        #expect(profile.maximumAbsoluteNormalizedLineSearchMarginDifference == nil)
    }

    @Test(.enabled(if: ProcessInfo.processInfo.environment["NUMIVIVO_TEST_METAL"] == "1"))
    func cohortAnalysisCarriesExplicitMetalNBProfile() throws {
        let data = try SingleCellNBCohortTests.fixture()
        var contrast = VivoOmicsExpressionContrast(id: "metal-nb", controlCondition: "ctrl",
            treatmentCondition: "stim", design: .pairedDonors)
        contrast.model = .negativeBinomial
        contrast.minimumCellsPerPseudobulk = 1
        var options = VivoOmicsNBCohortOptions()
        options.trend = .mean
        options.backend = .metalFP32
        contrast.negativeBinomialOptions = options
        let plan = VivoSingleCellAnalysisPlan(id: "metal-nb-plan", contrasts: [contrast])
        let report = try VivoSingleCellCohortAnalysis.run(data, plan: plan)
        let result = try #require(report.contrasts.first)
        let diagnostics = try #require(result.negativeBinomial?.features)
        let fits = diagnostics.compactMap(\.finalFit)
        let execution = try #require(result.negativeBinomial?.metalExecution)
        #expect(result.negativeBinomial != nil)
        #expect(fits.count >= 20)
        #expect(fits.allSatisfy { $0.backend == .metalFP32 })
        #expect(execution.schemaVersion == "numi.vivo.metal-nb-execution-profile.v2")
        #expect(execution.cohortEvaluationWallClockNanoseconds > 0)
        #expect(execution.fixedDispersionFitAttempts > 0)
        #expect(execution.completedFixedDispersionFits > 0)
        #expect(execution.metalObjectiveAttempts > 0)
        #expect(execution.completedMetalObjectiveEvaluations > 0)
        #expect(execution.metalObjectiveBatches > 0)
        #expect(execution.metalObjectiveObservations > 0)
        #expect(execution.cpuExactLikelihoodEvaluations > 0)
        #expect(execution.cpuExactLikelihoodObservations > 0)
        let candidateEvaluations = try #require(execution.lineSearchCandidateEvaluations)
        let metalAcceptedCPURejected = try #require(execution.metalAcceptedCPURejectedLineSearchCandidates)
        let cpuAcceptedMetalRejected = try #require(execution.cpuAcceptedMetalRejectedLineSearchCandidates)
        let maximumMarginDifference = try #require(execution.maximumAbsoluteNormalizedLineSearchMarginDifference)
        #expect(candidateEvaluations > 0)
        #expect(metalAcceptedCPURejected >= 0)
        #expect(cpuAcceptedMetalRejected >= 0)
        #expect(maximumMarginDifference.isFinite)
        let cohortDiagnostics = try #require(result.negativeBinomial)
        let encodedDiagnostics = try VivoCanonicalJSON.encode(cohortDiagnostics)
        #expect(!String(decoding: encodedDiagnostics, as: UTF8.self).contains("metalExecution"))
        let encodedExecution = try VivoCanonicalJSON.encode(execution)
        #expect(String(decoding: encodedExecution, as: UTF8.self).contains("cohortEvaluationWallClockNanoseconds"))
        #expect(String(decoding: encodedExecution, as: UTF8.self).contains("lineSearchCandidateEvaluations"))
        #expect(result.features.filter { $0.status == .tested }.count == result.testedFeatures)
        let encodedPlan = try VivoCanonicalJSON.encode(report.plan)
        #expect(String(decoding: encodedPlan, as: UTF8.self).contains(#""backend":"metalFP32""#))
        options.testMethod = .likelihoodRatio
        contrast.negativeBinomialOptions = options
        let lrt = try VivoSingleCellCohortAnalysis.run(data,
            plan: .init(id: "metal-nb-lrt-plan", contrasts: [contrast]))
        let lrtResult = try #require(lrt.contrasts.first)
        let lrtDiagnostics = try #require(lrtResult.negativeBinomial?.features)
        #expect(lrtResult.method.hasSuffix("LRT-v1"))
        #expect(lrtResult.features.filter { $0.status == .tested }.count == lrtResult.testedFeatures)
        #expect(lrtDiagnostics.compactMap(\.likelihoodRatioFit).count >= 20)
        var adjustedQL = options
        adjustedQL.testMethod = .quasiLikelihoodAdjusted
        contrast.negativeBinomialOptions = adjustedQL
        #expect(throws: (any Error).self) {
            _ = try VivoSingleCellCohortAnalysis.run(data,
                plan: .init(id: "metal-ql-rejected", contrasts: [contrast]))
        }
    }
}
