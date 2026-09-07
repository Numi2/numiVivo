import Foundation
import Testing
@testable import NumiVivoKit

@Suite(.serialized) struct NCMCProtocolKernelTests {
    private func fingerprint(_ text: String) throws -> VivoFingerprint {
        try VivoCanonicalJSON.fingerprint(Data(text.utf8))
    }
    private func physical() throws -> VivoConstantPHPhysicalState {
        try .init(stepIndex: 0, timePS: 0, positionsNM: [.init(0.25, 0, 0)],
                  velocitiesNMPerPS: [.zero], periodicCell: nil)
    }

    @Test func perturbationWorkTelescopesForCoordinateIndependentLinearHamiltonian() async throws {
        let manifold = try fingerprint("lambda-manifold")
        let hamiltonian = VivoNCMCLambdaHamiltonian(
            fingerprint: try fingerprint("linear-lambda-hamiltonian"),
            physicalManifoldFingerprint: manifold,
            fromStateIdentifier: "A", toStateIdentifier: "B",
            completePotentialEnergyKJPerMol: { _, lambda in 7.5 * lambda + 3 },
            propagate: { state, _, steps in
                let next = try VivoConstantPHPhysicalState(stepIndex: state.stepIndex + steps,
                    timePS: state.timePS + Double(steps) * 0.001,
                    positionsNM: state.positionsNM,
                    velocitiesNMPerPS: state.velocitiesNMPerPS,
                    periodicCell: state.periodicCell)
                return .init(finalPhysicalState: next, shadowWorkKJPerMol: 0,
                             steps: steps, method: "synthetic identity propagation")
            })
        let schedule = try VivoNCMCLambdaSchedule.linear(intervals: 5, propagationStepsPerLambda: 2)
        let result = try await VivoNCMCProtocolKernel.run(initial: physical(), hamiltonian: hamiltonian, schedule: schedule)
        #expect(abs(result.perturbationWorkKJPerMol - 7.5) < 1e-12)
        #expect(result.shadowWorkKJPerMol == 0)
        #expect(result.totalNonequilibriumWorkKJPerMol == result.perturbationWorkKJPerMol)
        #expect(result.steps.count == 5)
        #expect(result.finalPhysicalState.stepIndex == 10)
        #expect(abs(result.steps.map(\.perturbationWorkKJPerMol).reduce(0,+) - 7.5) < 1e-12)
    }

    @Test func shadowWorkAndPathProbabilityAccumulateSeparately() async throws {
        let manifold = try fingerprint("shadow-manifold")
        let hamiltonian = VivoNCMCLambdaHamiltonian(
            fingerprint: try fingerprint("shadow-hamiltonian"),
            physicalManifoldFingerprint: manifold,
            fromStateIdentifier: "A", toStateIdentifier: "B",
            completePotentialEnergyKJPerMol: { state, lambda in
                state.positionsNM[0].x * lambda
            },
            propagate: { state, _, steps in
                let next = try VivoConstantPHPhysicalState(stepIndex: state.stepIndex + steps,
                    timePS: state.timePS + Double(steps) * 0.001,
                    positionsNM: state.positionsNM.map { $0 + .init(0.1,0,0) },
                    velocitiesNMPerPS: state.velocitiesNMPerPS,
                    periodicCell: state.periodicCell)
                return .init(finalPhysicalState: next, shadowWorkKJPerMol: 0.2,
                             logReverseOverForwardPathProbability: -0.05,
                             steps: steps, method: "synthetic stochastic segment")
            })
        let schedule = try VivoNCMCLambdaSchedule.linear(intervals: 4, propagationStepsPerLambda: 1)
        let result = try await VivoNCMCProtocolKernel.run(initial: physical(), hamiltonian: hamiltonian, schedule: schedule)
        #expect(abs(result.shadowWorkKJPerMol - 0.8) < 1e-12)
        #expect(abs(result.logReverseOverForwardPathProbability + 0.2) < 1e-12)
        #expect(result.totalNonequilibriumWorkKJPerMol == result.perturbationWorkKJPerMol + result.shadowWorkKJPerMol)
        #expect(result.finalPhysicalState.positionsNM[0].x == 0.65)
    }

    @Test func scheduleRejectsNonmonotoneLambdaAndPropagationWorkOverflow() throws {
        let bad = VivoNCMCLambdaSchedule(lambdas: [0,0.5,0.4,1], propagationStepsPerLambda: 1)
        #expect(throws: (any Error).self) { try bad.validate() }
        let excessive = VivoNCMCLambdaSchedule(lambdas: (0...1000).map { Double($0)/1000 },
                                               propagationStepsPerLambda: 200_000)
        #expect(throws: (any Error).self) { try excessive.validate() }
    }
}
