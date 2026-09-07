import Foundation
import Testing
@testable import NumiVivoKit

private enum NCMCMetalTestStop: Error { case requested }

@Suite(.serialized) struct ConstantPHNCMCMetalTests {
    private func fp(_ value: String) throws -> VivoFingerprint {
        try VivoCanonicalJSON.fingerprint(Data(value.utf8))
    }
    private func f32(_ value: Double) -> Double { Double(Float(value)) }
    private func system(_ id: String, length: Double) throws -> VivoClassicalSystem {
        .init(identifier: id, structureFingerprint: try fp("ncmc-two-atoms"), particles: [
            .init(index: 0, atomIndex: 0, typeIdentifier: "X", massDa: 12, chargeE: 0, sigmaNM: 0, epsilonKJPerMol: 0),
            .init(index: 1, atomIndex: 1, typeIdentifier: "X", massDa: 12, chargeE: 0, sigmaNM: 0, epsilonKJPerMol: 0)],
            bonds: [.init(a: 0, b: 1, lengthNM: length, forceConstant: 1000)])
    }
    private func configuration() -> VivoMDConfiguration {
        .init(timeStepPS: 0.0005, cutoffNM: 1, neighborSkinNM: 0, electrostatics: .cutoff,
              ensemble: .nve, thermostat: .none, targetTemperatureK: nil, frictionPerPS: nil,
              neighborListEnabled: false, randomSeed: 818)
    }
    private func initial() throws -> VivoConstantPHPhysicalState {
        try .init(stepIndex: 0, timePS: 0, positionsNM: [.init(-f32(0.12),0,0), .init(f32(0.12),0,0)],
                  velocitiesNMPerPS: [.init(f32(0.08),0,0), .init(-f32(0.08),0,0)], periodicCell: nil)
    }
    private func specs() throws -> [VivoConstantPHMetalStateSpecification] {
        let a = try system("endpoint-A", length: 0.2), b = try system("endpoint-B", length: 0.3)
        // This provider is a synthetic external quadratic potential, not a QM
        // validation claim. It verifies that the full provider energy and force,
        // not only retained classical terms, participate in alchemical switching.
        let id = try fp("B-external-quadratic"), systemID = try b.fingerprint()
        let extra = try VivoMDCandidateForceProvider(fingerprint: id, retainedSystemFingerprint: systemID,
            boundary: .finiteCluster, supportsCellMoves: false, maximumAcceptedResidual: 1e-8,
            molecularConnectivitySystem: b) { geometry in
                let x = geometry.particlePositionsNM[0].x
                return try .init(providerFingerprint: id, geometry: geometry,
                    additionalEnergyKJPerMol: -2 + 2.5*x*x,
                    physicalParticleForcesKJPerMolNM: [.init(-5*x,0,0), .zero],
                    derivativeMethod: "synthetic quadratic derivative", convergenceResidual: 0, requiredResidual: 1e-8)
            }
        return [.init(identifier: "A", system: a, configuration: configuration()),
                .init(identifier: "B", system: b, configuration: configuration(), forceProvider: extra)]
    }
    private func setup() throws -> VivoConstantPHNCMCMetalSetup {
        try VivoConstantPHNCMCMetalFactory.make(specifications: specs(), samplingTemperatureK: 300,
            configuration: .init(schedule: .init(lambdas: [0,0.125,0.5,1], propagationStepsPerLambda: 2), timeStepPS: 0.0005))
    }
    private func kinetic(_ state: VivoConstantPHPhysicalState) -> Double {
        state.velocitiesNMPerPS.reduce(0) { $0 + 6*$1.squaredNorm }
    }

    @Test func nativeSwitchMeasuresFullWorkAndFollowsTheReverseProtocol() async throws {
        let setup = try setup(), input = try initial()
        let a = setup.executableStates[0], b = setup.executableStates[1]
        let u0 = try await a.potentialEnergyKJPerMol(input)
        let forward = try await setup.switchEngine.switchCandidate(input, a, b)
        let u1 = try await b.potentialEnergyKJPerMol(forward.finalPhysicalState)
        try forward.validate(engine: setup.switchEngine, initial: input, from: a, to: b, maximumAbsoluteWorkKJPerMol: 1e8)
        let deltaH = (u1-u0) + kinetic(forward.finalPhysicalState)-kinetic(input)
        #expect(abs(forward.totalNonequilibriumWorkKJPerMol - deltaH) < 1e-7)
        #expect(forward.switchingSteps == 6 && forward.finalPhysicalState.stepIndex == 6)
        #expect(forward.logReverseOverForwardPathProbability == 0)
        #expect(abs(forward.totalNonequilibriumWorkKJPerMol - (u1-u0)) > 1e-7)
        var reverseInput = forward.finalPhysicalState
        reverseInput.velocitiesNMPerPS = reverseInput.velocitiesNMPerPS.map { $0 * -1 }
        let reverse = try await setup.switchEngine.switchCandidate(reverseInput, b, a)
        for i in input.positionsNM.indices {
            #expect((reverse.finalPhysicalState.positionsNM[i] - input.positionsNM[i]).norm < 1e-6)
            #expect((reverse.finalPhysicalState.velocitiesNMPerPS[i] + input.velocitiesNMPerPS[i]).norm < 1e-5)
        }
        #expect(abs(forward.totalNonequilibriumWorkKJPerMol + reverse.totalNonequilibriumWorkKJPerMol) < 1e-4)
    }

    @Test func samplerRestartsFromAnAcceptedTransactionWithoutResettingRNGOrHamiltonian() async throws {
        let setup = try setup(), input = try initial()
        let evidence = VivoKineticEvidence(source: "synthetic native NCMC fixture", locator: "ConstantPHNCMCMetalTests")
        let definitions = [
            VivoConstantPHStateDefinition(identifier: "A", boundProtonOffset: 0, referenceSemigrandBiasKJPerMol: 0,
                origin: .assumed, evidence: evidence, neighbors: ["B"]),
            .init(identifier: "B", boundProtonOffset: 1, referenceSemigrandBiasKJPerMol: 0,
                origin: .assumed, evidence: evidence, neighbors: ["A"])]
        let config = VivoConstantPHNCMCConfiguration(base: .init(identifier: "native-NCMC-restart", temperatureK: 300,
            referencePH: 7, targetPH: 7, states: definitions, initialStateIdentifier: "A",
            mdStepsPerAttempt: 1, attemptCount: 2, seed: 0))
        let direct = try VivoConstantPHNCMC(configuration: config, initialPhysicalState: input,
            executableStates: setup.executableStates, switchEngine: setup.switchEngine)
        let expected = try await direct.run()
        let interrupted = try VivoConstantPHNCMC(configuration: config, initialPhysicalState: input,
            executableStates: setup.executableStates, switchEngine: setup.switchEngine)
        do {
            _ = try await interrupted.run { checkpoint in
                if checkpoint.completedAttempts == 1 { throw NCMCMetalTestStop.requested }
            }
            Issue.record("expected explicit interruption")
        } catch NCMCMetalTestStop.requested {}
        let checkpoint = await interrupted.checkpoint()
        #expect(checkpoint.completedAttempts == 1)
        let resumed = try VivoConstantPHNCMC(configuration: config, initialPhysicalState: input,
            executableStates: setup.executableStates, switchEngine: setup.switchEngine, checkpoint: checkpoint)
        #expect(try await resumed.run() == expected)
    }

    @Test func factoryRejectsUnbudgetedSwitchesAndDoesNotSilentlyRoundIndependentState() async throws {
        #expect(throws: (any Error).self) {
            _ = try VivoConstantPHNCMCMetalFactory.make(specifications: specs(), samplingTemperatureK: 300,
                configuration: .init(schedule: .linear(intervals: 4, propagationStepsPerLambda: 2), maximumEndpointEvaluations: 1))
        }
        let setup = try setup()
        var unrounded = try initial(); unrounded.positionsNM[0].x = 0.1
        await #expect(throws: (any Error).self) {
            _ = try await setup.switchEngine.switchCandidate(unrounded, setup.executableStates[0], setup.executableStates[1])
        }
        var invalid = try specs(); var changed = invalid[1].system; changed.particles[0].massDa = 13
        invalid[1] = .init(identifier: "B", system: changed, configuration: configuration())
        #expect(throws: (any Error).self) {
            _ = try VivoConstantPHNCMCMetalFactory.make(specifications: invalid, samplingTemperatureK: 300,
                configuration: .init(schedule: .linear(intervals: 2, propagationStepsPerLambda: 1)))
        }
    }
}
