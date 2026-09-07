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

    @Test func deterministicShadowWorkAccumulatesAcrossSymmetricIntervals() async throws {
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
                             logReverseOverForwardPathProbability: 0,
                             steps: steps, method: "synthetic deterministic bookkeeping segment")
            })
        let schedule = try VivoNCMCLambdaSchedule.linear(intervals: 4, propagationStepsPerLambda: 1)
        let result = try await VivoNCMCProtocolKernel.run(initial: physical(), hamiltonian: hamiltonian, schedule: schedule)
        #expect(abs(result.shadowWorkKJPerMol - 0.8) < 1e-12)
        #expect(result.logReverseOverForwardPathProbability == 0)
        #expect(result.totalNonequilibriumWorkKJPerMol == result.perturbationWorkKJPerMol + result.shadowWorkKJPerMol)
        #expect(abs(result.finalPhysicalState.positionsNM[0].x - 0.65) < 1e-12)
    }

    @Test func scheduleRejectsNonmonotoneLambdaAndPropagationWorkOverflow() throws {
        let bad = VivoNCMCLambdaSchedule(lambdas: [0,0.5,0.4,1], propagationStepsPerLambda: 1)
        #expect(throws: (any Error).self) { try bad.validate() }
        let excessive = VivoNCMCLambdaSchedule(lambdas: (0...1000).map { Double($0)/1000 },
                                               propagationStepsPerLambda: 200_000)
        #expect(throws: (any Error).self) { try excessive.validate() }
    }

    @Test func midpointProtocolHasAnActualReverseAndAntisymmetricWork() async throws {
        let manifold = try fingerprint("round-trip-manifold")
        func hamiltonian(reverse: Bool) throws -> VivoNCMCLambdaHamiltonian {
            let reverseDirection = reverse
            @Sendable func potential(_ state: VivoConstantPHPhysicalState, _ lambda: Double) -> Double {
                let l = reverseDirection ? 1 - lambda : lambda
                let x = state.positionsNM[0].x
                return 0.5 * (1 + 3*l) * x*x + 0.2*l*x
            }
            @Sendable func force(_ x: Double, _ lambda: Double) -> Double {
                let l = reverseDirection ? 1 - lambda : lambda
                return -(1 + 3*l)*x - 0.2*l
            }
            return .init(fingerprint: try fingerprint("round-trip-\(reverse)"), physicalManifoldFingerprint: manifold,
                fromStateIdentifier: reverse ? "B" : "A", toStateIdentifier: reverse ? "A" : "B",
                completePotentialEnergyKJPerMol: { state, l in potential(state, l) },
                propagate: { state, l, steps in
                    let dt = 0.03
                    let initialEnergy = potential(state, l) + 0.5 * state.velocitiesNMPerPS[0].squaredNorm
                    var next = state
                    for _ in 0..<steps {
                        next.velocitiesNMPerPS[0].x += 0.5 * dt * force(next.positionsNM[0].x, l)
                        next.positionsNM[0].x += dt * next.velocitiesNMPerPS[0].x
                        next.velocitiesNMPerPS[0].x += 0.5 * dt * force(next.positionsNM[0].x, l)
                    }
                    next.stepIndex += steps; next.timePS += Double(steps)*dt
                    let finalEnergy = potential(next, l) + 0.5 * next.velocitiesNMPerPS[0].squaredNorm
                    return .init(finalPhysicalState: next, shadowWorkKJPerMol: finalEnergy-initialEnergy,
                                 steps: steps, method: "analytic oscillator velocity-Verlet fixture")
                })
        }
        var initial = try physical(); initial.velocitiesNMPerPS[0].x = 0.6
        let schedule = VivoNCMCLambdaSchedule(lambdas: [0,0.125,0.5,0.875,1], propagationStepsPerLambda: 3)
        let forward = try await VivoNCMCProtocolKernel.run(initial: initial,
            hamiltonian: hamiltonian(reverse: false), schedule: schedule)
        var reversed = forward.finalPhysicalState
        reversed.velocitiesNMPerPS = reversed.velocitiesNMPerPS.map { $0 * -1 }
        let backward = try await VivoNCMCProtocolKernel.run(initial: reversed,
            hamiltonian: hamiltonian(reverse: true), schedule: schedule.reversed())
        #expect((backward.finalPhysicalState.positionsNM[0] - initial.positionsNM[0]).norm < 1e-13)
        #expect((backward.finalPhysicalState.velocitiesNMPerPS[0] + initial.velocitiesNMPerPS[0]).norm < 1e-13)
        #expect(abs(forward.perturbationWorkKJPerMol + backward.perturbationWorkKJPerMol) < 1e-13)
        #expect(abs(forward.shadowWorkKJPerMol + backward.shadowWorkKJPerMol) < 1e-13)
        #expect(abs(forward.shadowWorkKJPerMol) > 1e-8)
        #expect(forward.steps.allSatisfy {
            abs($0.perturbationWorkKJPerMol - $0.prePropagationPerturbationWorkKJPerMol - $0.postPropagationPerturbationWorkKJPerMol) < 1e-13
        })
    }

    @Test func rawStochasticPathProbabilityIsNotCountedTwiceAsWork() async throws {
        let h = VivoNCMCLambdaHamiltonian(fingerprint: try fingerprint("unqualified-heat"),
            physicalManifoldFingerprint: try fingerprint("manifold"), fromStateIdentifier: "A", toStateIdentifier: "B",
            completePotentialEnergyKJPerMol: { _, l in l }, propagate: { input, _, steps in
                var output = input; output.stepIndex += steps; output.timePS += Double(steps)*0.001
                return .init(finalPhysicalState: output, shadowWorkKJPerMol: 0.1,
                    logReverseOverForwardPathProbability: 0.2, steps: steps, method: "unqualified stochastic kernel")
            })
        await #expect(throws: (any Error).self) {
            _ = try await VivoNCMCProtocolKernel.run(initial: physical(), hamiltonian: h,
                schedule: .linear(intervals: 2, propagationStepsPerLambda: 1))
        }
    }

    @Test func wrongStepCountAndCellChangesCannotMasqueradeAsFixedManifoldPropagation() throws {
        let input = try physical()
        var wrong = input; wrong.stepIndex = 7; wrong.timePS = 0.001
        let result = VivoNCMCPropagationResult(finalPhysicalState: wrong, shadowWorkKJPerMol: 0, steps: 1, method: "bad clock")
        #expect(throws: (any Error).self) { try result.validate(initial: input, expectedSteps: 1, maximumAbsoluteWorkKJPerMol: 10) }
        wrong.stepIndex = 1
        wrong.periodicCell = .init(a: .init(2,0,0), b: .init(0,2,0), c: .init(0,0,2))
        let changedCell = VivoNCMCPropagationResult(finalPhysicalState: wrong, shadowWorkKJPerMol: 0, steps: 1, method: "unaccounted cell move")
        #expect(throws: (any Error).self) { try changedCell.validate(initial: input, expectedSteps: 1, maximumAbsoluteWorkKJPerMol: 10) }
    }
}
