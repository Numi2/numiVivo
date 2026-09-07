import Foundation
import Testing
@testable import NumiVivoKit

@Suite(.serialized) struct GlobalKineticUncertaintyTests {
    @Test func fullNonlinearSurvivalUsesQualifiedBaselineAndJointDraws() throws {
        let pathway = try QMMMChemicalExchangeNetworkTests().replicatedPathway(state: "A",executionPrefix: "global-fixture",seedBase: 4000)
        let k = pathway.result.geometricMeanRatePerSecond
        let state = VivoQMMMChemicalExchangeState(identifier: "A",initialPopulation: 1,populationOrigin: .assumed,
            populationEvidence: .init(source: "synthetic fixture",locator: "initial A"),pathways: [pathway])
        let network = VivoQMMMChemicalExchangeNetworkRequest(identifier: "joint-nonlinear",states: [state],exchangeEdges: [],observationTimesSeconds: [0,1/k])
        let scenarios = [VivoJointKineticScenario(modelIdentifier: "m",network: network)]
        let parameters: [VivoJointKineticParameter] = [.pathwayLogPMFFactor(state: 0,pathway: 0),.pathwayLogNuclearFactor(state: 0,pathway: 0)]
        let context = try VivoGlobalKineticUncertainty.context(scenarios: scenarios,parameters: parameters)
        let ensemble = VivoGlobalJointUncertaintyRequest(contextFingerprint: context,parameterIdentifiers: parameters.map(\.identifier),
            models: [.init(identifier: "m",probability: 1,evidenceFingerprint: pathway.result.evidenceFingerprint,interpretation: "synthetic baseline model")],
            draws: [.init(identifier: "a",modelIdentifier: "m",dependenceBlockIdentifier: "a",values: [log(2),log(2)]),
                    .init(identifier: "b",modelIdentifier: "m",dependenceBlockIdentifier: "b",values: [log(0.5),log(0.5)])],
            unresolvedComponents: ["chemical model error"],assumptions: ["deliberate perfect dependence fixture"])
        let request = VivoGlobalKineticUncertaintyRequest(scenarios: scenarios,parameters: parameters,ensemble: ensemble)
        let result = try VivoGlobalKineticUncertainty.calculate(request)
        #expect(abs(try #require(result.summaries[5].mean)-(exp(-4)+exp(-0.25))/2) < 1e-10)
        #expect(abs(try #require(result.summaries[7].mean)/k-2.125) < 1e-10)
        #expect(result.summaries.count == 10 && result.summaries[9].mean == result.summaries[5].mean)
        #expect(result.unresolvedComponents == ["chemical model error"])
        #expect(result.assumptions.contains { $0.contains("baseline classical transmission is assumed") })
        try VivoGlobalKineticUncertainty.validate(result,request: request)
        var transplanted = ensemble; transplanted.contextFingerprint = try PrecisionSamplingFixtures.id("other-Hamiltonian")
        #expect(throws: (any Error).self) { try VivoGlobalKineticUncertainty.calculate(.init(scenarios: scenarios,parameters: parameters,ensemble: transplanted)) }
    }
}
