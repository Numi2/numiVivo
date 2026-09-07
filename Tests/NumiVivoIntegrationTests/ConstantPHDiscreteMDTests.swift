import Foundation
import Testing
@testable import NumiVivoKit

@Suite(.serialized) struct ConstantPHDiscreteMDTests {
    private func fingerprint(_ text: String) throws -> VivoFingerprint {
        try VivoCanonicalJSON.fingerprint(Data(text.utf8))
    }

    private func evidence(_ label: String) -> VivoKineticEvidence {
        .init(source: "synthetic constant-pH regression", locator: label)
    }

    private func physical() throws -> VivoConstantPHPhysicalState {
        try .init(stepIndex: 0, timePS: 0,
                  positionsNM: [.init(0, 0, 0)],
                  velocitiesNMPerPS: [.zero], periodicCell: nil)
    }

    private func executable(_ identifier: String, manifold: VivoFingerprint,
                            energy: Double) throws -> VivoConstantPHExecutableState {
        VivoConstantPHExecutableState(identifier: identifier,
            physicalManifoldFingerprint: manifold,
            hamiltonianFingerprint: try fingerprint("hamiltonian-\(identifier)"),
            potentialEnergyKJPerMol: { _ in energy },
            propagate: { state, steps in
                try VivoConstantPHPhysicalState(stepIndex: state.stepIndex + steps,
                    timePS: state.timePS + Double(steps) * 0.001,
                    positionsNM: state.positionsNM,
                    velocitiesNMPerPS: state.velocitiesNMPerPS,
                    periodicCell: state.periodicCell)
            })
    }

    @Test func onePHUnitProducesExpectedTenfoldPopulationRatio() async throws {
        let states = [
            VivoConstantPHStateDefinition(identifier: "protonated", boundProtonOffset: 1,
                referenceSemigrandBiasKJPerMol: 0, origin: .assumed,
                evidence: evidence("equal reference free energy"), neighbors: ["deprotonated"]),
            VivoConstantPHStateDefinition(identifier: "deprotonated", boundProtonOffset: 0,
                referenceSemigrandBiasKJPerMol: 0, origin: .assumed,
                evidence: evidence("equal reference free energy"), neighbors: ["protonated"])
        ]
        let configuration = VivoConstantPHConfiguration(identifier: "two-state-pH-regression",
            temperatureK: 300, referencePH: 7, targetPH: 8, states: states,
            initialStateIdentifier: "protonated", mdStepsPerAttempt: 1,
            attemptCount: 50_000, seed: 20260907)
        let manifold = try fingerprint("one-particle-common-manifold")
        let sampler = try VivoConstantPHDiscreteMD(configuration: configuration,
            initialPhysicalState: physical(), executableStates: [
                executable("protonated", manifold: manifold, energy: 0),
                executable("deprotonated", manifold: manifold, energy: 0)
            ])
        let result = try await sampler.run()
        let protonated = try #require(result.populations["protonated"])
        let deprotonated = try #require(result.populations["deprotonated"])
        #expect(abs(deprotonated / protonated - 10) < 0.5)
        #expect(result.acceptanceFraction > 0.15 && result.acceptanceFraction < 0.25)
        #expect(result.finalCheckpoint.physicalState.stepIndex == 50_000)
        #expect(result.finalCheckpoint.completedAttempts == 50_000)
        #expect(result.finalCheckpoint.attempts.count == 50_000)
    }

    @Test func equalStatesAtReferencePHRemainEquiprobableAndReplayExactly() async throws {
        let states = [
            VivoConstantPHStateDefinition(identifier: "a", boundProtonOffset: 0,
                referenceSemigrandBiasKJPerMol: 0, origin: .assumed,
                evidence: evidence("a"), neighbors: ["b"]),
            VivoConstantPHStateDefinition(identifier: "b", boundProtonOffset: 1,
                referenceSemigrandBiasKJPerMol: 0, origin: .assumed,
                evidence: evidence("b"), neighbors: ["a"])
        ]
        let configuration = VivoConstantPHConfiguration(identifier: "deterministic-replay",
            temperatureK: 310, referencePH: 6.5, targetPH: 6.5, states: states,
            initialStateIdentifier: "a", mdStepsPerAttempt: 2,
            attemptCount: 1000, seed: 77)
        let manifold = try fingerprint("replay-manifold")
        let executableStates = try [executable("a", manifold: manifold, energy: 5),
                                    executable("b", manifold: manifold, energy: 5)]
        let first = try await VivoConstantPHDiscreteMD(configuration: configuration,
            initialPhysicalState: physical(), executableStates: executableStates).run()
        let second = try await VivoConstantPHDiscreteMD(configuration: configuration,
            initialPhysicalState: physical(), executableStates: executableStates).run()
        #expect(first == second)
        #expect(first.acceptanceFraction == 1)
        #expect(first.finalCheckpoint.physicalState.stepIndex == 2000)
        #expect(abs((first.populations["a"] ?? 0) - 0.5) < 0.001)
        #expect(abs((first.populations["b"] ?? 0) - 0.5) < 0.001)
    }

    @Test func asymmetricDegreeUsesHastingsProposalCorrection() async throws {
        let states = [
            VivoConstantPHStateDefinition(identifier: "hub", boundProtonOffset: 0,
                referenceSemigrandBiasKJPerMol: 0, origin: .assumed,
                evidence: evidence("hub"), neighbors: ["left", "right"]),
            VivoConstantPHStateDefinition(identifier: "left", boundProtonOffset: 0,
                referenceSemigrandBiasKJPerMol: 0, origin: .assumed,
                evidence: evidence("left"), neighbors: ["hub"]),
            VivoConstantPHStateDefinition(identifier: "right", boundProtonOffset: 0,
                referenceSemigrandBiasKJPerMol: 0, origin: .assumed,
                evidence: evidence("right"), neighbors: ["hub"])
        ]
        let configuration = VivoConstantPHConfiguration(identifier: "hastings-degree",
            temperatureK: 300, referencePH: 7, targetPH: 7, states: states,
            initialStateIdentifier: "hub", mdStepsPerAttempt: 1,
            attemptCount: 4000, seed: 9)
        let manifold = try fingerprint("hastings-manifold")
        let sampler = try VivoConstantPHDiscreteMD(configuration: configuration,
            initialPhysicalState: physical(), executableStates: try states.map {
                try executable($0.identifier, manifold: manifold, energy: 0)
            })
        let result = try await sampler.run()
        let hubOutgoing = try #require(result.finalCheckpoint.attempts.first { $0.fromStateIdentifier == "hub" })
        #expect(abs(hubOutgoing.logProposalRatio - log(2)) < 1e-14)
        let leafOutgoing = try #require(result.finalCheckpoint.attempts.first { $0.fromStateIdentifier != "hub" })
        #expect(abs(leafOutgoing.logProposalRatio + log(2)) < 1e-14)
        for value in result.populations.values { #expect(abs(value - 1.0 / 3.0) < 0.04) }
    }

    @Test func rejectsNonreciprocalGraphAndMixedPhysicalManifold() async throws {
        let invalid = [
            VivoConstantPHStateDefinition(identifier: "a", boundProtonOffset: 0,
                referenceSemigrandBiasKJPerMol: 0, origin: .assumed,
                evidence: evidence("a"), neighbors: ["b"]),
            VivoConstantPHStateDefinition(identifier: "b", boundProtonOffset: 0,
                referenceSemigrandBiasKJPerMol: 0, origin: .assumed,
                evidence: evidence("b"), neighbors: ["a", "c"]),
            VivoConstantPHStateDefinition(identifier: "c", boundProtonOffset: 0,
                referenceSemigrandBiasKJPerMol: 0, origin: .assumed,
                evidence: evidence("c"), neighbors: ["b"])
        ]
        var nonreciprocal = invalid
        nonreciprocal[2].neighbors = ["a"]
        let badConfiguration = VivoConstantPHConfiguration(identifier: "bad-graph",
            temperatureK: 300, referencePH: 7, targetPH: 7,
            states: nonreciprocal, initialStateIdentifier: "a",
            mdStepsPerAttempt: 1, attemptCount: 2, seed: 1)
        #expect(throws: (any Error).self) { try badConfiguration.validate() }

        let validStates = [
            VivoConstantPHStateDefinition(identifier: "a", boundProtonOffset: 0,
                referenceSemigrandBiasKJPerMol: 0, origin: .assumed,
                evidence: evidence("a"), neighbors: ["b"]),
            VivoConstantPHStateDefinition(identifier: "b", boundProtonOffset: 0,
                referenceSemigrandBiasKJPerMol: 0, origin: .assumed,
                evidence: evidence("b"), neighbors: ["a"])
        ]
        let validConfiguration = VivoConstantPHConfiguration(identifier: "mixed-manifold",
            temperatureK: 300, referencePH: 7, targetPH: 7,
            states: validStates, initialStateIdentifier: "a",
            mdStepsPerAttempt: 1, attemptCount: 2, seed: 1)
        let a = try executable("a", manifold: fingerprint("m1"), energy: 0)
        let b = try executable("b", manifold: fingerprint("m2"), energy: 0)
        #expect(throws: (any Error).self) {
            _ = try VivoConstantPHDiscreteMD(configuration: validConfiguration,
                initialPhysicalState: physical(), executableStates: [a, b])
        }
    }
}