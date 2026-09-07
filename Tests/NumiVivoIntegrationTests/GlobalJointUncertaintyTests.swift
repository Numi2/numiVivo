import Foundation
import Testing
@testable import NumiVivoKit

@Suite(.serialized) struct GlobalJointUncertaintyTests {
    private func model(_ id: String = "a", probability: Double = 1) throws -> VivoJointProbabilityModel {
        .init(identifier: id, probability: probability, evidenceFingerprint: try PrecisionSamplingFixtures.id(id),
              interpretation: "explicit finite distribution fixture, not an inferred posterior")
    }
    private func request(_ draws: [VivoJointParameterDraw], parameters: [String] = ["x","y"],
                         models: [VivoJointProbabilityModel]? = nil) throws -> VivoGlobalJointUncertaintyRequest {
        .init(contextFingerprint: try PrecisionSamplingFixtures.id("joint-fixture"),parameterIdentifiers: parameters,
            models: try models ?? [model()],draws: draws,unresolvedComponents: ["electronic model discrepancy"],
            assumptions: ["finite support is supplied, not inferred"])
    }
    @Test func nonlinearJointPushForwardPreservesPerfectDependence() throws {
        let input = try request([
            .init(identifier: "low",modelIdentifier: "a",dependenceBlockIdentifier: "low",values: [-1,-1]),
            .init(identifier: "high",modelIdentifier: "a",dependenceBlockIdentifier: "high",values: [1,1])])
        let result = try VivoGlobalJointUncertainty.propagate(input,observableIdentifiers: ["nonlinear","difference"]) { x,_,_ in
            .init(values: [exp(x[0]+x[1]), x[0]-x[1]],chargedPrimitiveWork: 3)
        }
        let expected = cosh(2.0)
        #expect(abs(try #require(result.summaries[0].mean)-expected) < 1e-12)
        #expect(result.summaries[1].variance == 0)
        #expect(abs(expected-pow(cosh(1),2)) > 1) // independent marginal resampling would change the answer
        #expect(result.unresolvedComponents == input.unresolvedComponents)
        #expect(result.assumptions == input.assumptions)
        #expect(try #require(result.covariance)[0,1] == 0)
    }
    @Test func modelProbabilitiesDoNotDependOnDrawCounts() throws {
        let models = try [model("a",probability: 0.9),model("b",probability: 0.1)]
        var draws = [VivoJointParameterDraw(identifier: "a0",modelIdentifier: "a",dependenceBlockIdentifier: "a0",values: [0])]
        for i in 0..<9 { draws.append(.init(identifier: "b\(i)",modelIdentifier: "b",dependenceBlockIdentifier: "b\(i)",values: [10])) }
        let input = try request(draws,parameters: ["x"],models: models)
        let result = try VivoGlobalJointUncertainty.propagate(input,observableIdentifiers: ["x"]) { x,_,_ in .init(values: x.map(Optional.some),chargedPrimitiveWork: 1) }
        let s = result.summaries[0]
        #expect(abs(try #require(s.mean)-1) < 1e-12)
        #expect(abs(try #require(s.variance)-9) < 1e-12)
        #expect(abs(try #require(s.withinModelVariance)) < 1e-12)
        #expect(abs(try #require(s.betweenModelVariance)-9) < 1e-12)
        #expect(s.median == 0 && s.quantile975 == 10)
    }
    @Test func undefinedObservablesAreNotDroppedFromTheMeasure() throws {
        let input = try request([
            .init(identifier: "0",modelIdentifier: "a",dependenceBlockIdentifier: "same",values: [0]),
            .init(identifier: "1",modelIdentifier: "a",dependenceBlockIdentifier: "same",values: [1])],parameters: ["x"])
        let result = try VivoGlobalJointUncertainty.propagate(input,observableIdentifiers: ["log"]) { x,_,_ in
            .init(values: [x[0] == 0 ? nil : log(x[0])],chargedPrimitiveWork: 1)
        }
        #expect(result.summaries[0].definedProbability == 0.5)
        #expect(result.summaries[0].mean == 0 && result.covariance == nil)
        #expect(result.declaredBlockEffectiveCount == 1 && result.weightEffectiveDraws == 2)
    }
    @Test func gaussianSingularCovarianceIsPreservedWithoutJitter() throws {
        let covariance = try VivoQMMatrix(rows: 2,columns: 2,values: [1,1,1,1])
        let a = try VivoJointGaussianDraws.generate(mean: [2,2],covariance: covariance,count: 64,seed: 0,modelIdentifier: "joint")
        let b = try VivoJointGaussianDraws.generate(mean: [2,2],covariance: covariance,count: 64,seed: 0,modelIdentifier: "joint")
        #expect(a == b && a.allSatisfy { abs($0.values[0]-$0.values[1]) < 1e-12 })
        let tiny = try VivoQMMatrix(rows: 2,columns: 2,values: [1e-24,1e-24,1e-24,1e-24])
        let small = try VivoJointGaussianDraws.generate(mean: [0,0],covariance: tiny,count: 64,seed: 0,modelIdentifier: "tiny")
        #expect(small.allSatisfy { abs($0.values[0]-$0.values[1]) < 1e-24 })
        let indefinite = try VivoQMMatrix(rows: 2,columns: 2,values: [1,2,2,1])
        #expect(throws: (any Error).self) { try VivoJointGaussianDraws.generate(mean: [0,0],covariance: indefinite,count: 3,seed: 1,modelIdentifier: "bad") }
    }
    @Test func budgetInvalidDrawAndForwardFailureCannotBeSilentlySkipped() throws {
        var input = try request([
            .init(identifier: "0",modelIdentifier: "a",dependenceBlockIdentifier: "0",values: [0]),
            .init(identifier: "1",modelIdentifier: "a",dependenceBlockIdentifier: "1",values: [1])],parameters: ["x"])
        #expect(throws: (any Error).self) {
            try VivoGlobalJointUncertainty.propagate(input,observableIdentifiers: ["x"]) { x,_,_ in
                if x[0] > 0 { throw VivoChemistryError.convergence("forward model failed") }
                return .init(values: x.map(Optional.some),chargedPrimitiveWork: 1)
            }
        }
        input.maximumPrimitiveWork = 1
        #expect(throws: (any Error).self) { try VivoGlobalJointUncertainty.propagate(input,observableIdentifiers: ["x"]) { x,_,_ in .init(values: x.map(Optional.some),chargedPrimitiveWork: 1) } }
        input.maximumPrimitiveWork = 100_000
        input.draws[1] = .init(identifier: "1",modelIdentifier: "a",dependenceBlockIdentifier: "1",values: [.nan])
        #expect(throws: (any Error).self) { try input.validate() }
    }
}
