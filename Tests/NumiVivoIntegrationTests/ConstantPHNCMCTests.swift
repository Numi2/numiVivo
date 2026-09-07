import Foundation
import Testing
@testable import NumiVivoKit

@Suite(.serialized) struct ConstantPHNCMCTests {
    private func fingerprint(_ text: String) throws -> VivoFingerprint {
        try VivoCanonicalJSON.fingerprint(Data(text.utf8))
    }
    private func evidence() -> VivoKineticEvidence {
        .init(source: "synthetic NCMC contract fixture", locator: "ConstantPHNCMCTests.swift")
    }
    private func definitions() -> [VivoConstantPHStateDefinition] {
        [
            .init(identifier: "A", boundProtonOffset: 0, referenceSemigrandBiasKJPerMol: 0,
                  origin: .assumed, evidence: evidence(), neighbors: ["B"]),
            .init(identifier: "B", boundProtonOffset: 0, referenceSemigrandBiasKJPerMol: 0,
                  origin: .assumed, evidence: evidence(), neighbors: ["A"])
        ]
    }
    private func physical() throws -> VivoConstantPHPhysicalState {
        try .init(stepIndex: 0, timePS: 0,
                  positionsNM: [.zero], velocitiesNMPerPS: [.zero], periodicCell: nil)
    }
    private func executable(_ id: String, manifold: VivoFingerprint) throws -> VivoConstantPHExecutableState {
        VivoConstantPHExecutableState(identifier: id,
            physicalManifoldFingerprint: manifold,
            hamiltonianFingerprint: try fingerprint("hamiltonian-\(id)"),
            potentialEnergyKJPerMol: { _ in 0 },
            propagate: { input, steps in
                try VivoConstantPHPhysicalState(stepIndex: input.stepIndex + steps,
                    timePS: input.timePS + Double(steps) * 0.001,
                    positionsNM: input.positionsNM.map { $0 + .init(Double(steps), 0, 0) },
                    velocitiesNMPerPS: input.velocitiesNMPerPS,
                    periodicCell: input.periodicCell)
            })
    }

    @Test func zeroWorkSymmetricSwitchIsAcceptedAndCommitsSwitchEndpoint() async throws {
        let manifold = try fingerprint("common-manifold")
        let states = try [executable("A", manifold: manifold), executable("B", manifold: manifold)]
        let base = VivoConstantPHConfiguration(identifier: "zero-work", temperatureK: 300,
            referencePH: 7, targetPH: 7, states: definitions(), initialStateIdentifier: "A",
            mdStepsPerAttempt: 1, attemptCount: 1, seed: 0)
        let cfg = VivoConstantPHNCMCConfiguration(base: base)
        let engineID = try fingerprint("zero-work-engine")
        let engine = VivoConstantPHNCMCSwitchEngine(fingerprint: engineID, physicalManifoldFingerprint: manifold) { input, from, to in
            let final = try VivoConstantPHPhysicalState(stepIndex: input.stepIndex + 4,
                timePS: input.timePS + 0.004,
                positionsNM: input.positionsNM.map { $0 + .init(0, 2, 0) },
                velocitiesNMPerPS: input.velocitiesNMPerPS,
                periodicCell: input.periodicCell)
            return try .init(switchEngineFingerprint: engineID,
                physicalManifoldFingerprint: manifold,
                fromStateIdentifier: from.identifier, toStateIdentifier: to.identifier,
                initialPhysicalState: input, finalPhysicalState: final,
                protocolWorkKJPerMol: 0, shadowWorkKJPerMol: 0,
                switchingSteps: 4, interpretation: "synthetic reversible zero-work switch")
        }
        let sampler = try VivoConstantPHNCMC(configuration: cfg, initialPhysicalState: physical(),
                                             executableStates: states, switchEngine: engine)
        let result = try await sampler.run()
        #expect(result.finalCheckpoint.currentStateIdentifier == "B")
        #expect(result.finalCheckpoint.acceptedMoves == 1)
        #expect(result.finalCheckpoint.physicalState.stepIndex == 5)
        #expect(result.finalCheckpoint.physicalState.positionsNM[0] == .init(1, 2, 0))
        #expect(result.finalCheckpoint.attempts[0].logAcceptanceProbability == 0)
        #expect(result.finalCheckpoint.attempts[0].endpointPotentialEnergyDifferenceKJPerMol == 0)
        let hamiltonians = Dictionary(uniqueKeysWithValues: states.map { ($0.identifier, $0.hamiltonianFingerprint) })
        try VivoConstantPHNCMC.validate(checkpoint: result.finalCheckpoint, configuration: cfg,
            manifold: manifold, hamiltonians: hamiltonians, switchEngineFingerprint: engineID)
        var encoded = try #require(JSONSerialization.jsonObject(with: VivoCanonicalJSON.encode(result.finalCheckpoint)) as? [String: Any])
        var attempts = try #require(encoded["attempts"] as? [[String: Any]])
        attempts[0]["accepted"] = false; encoded["attempts"] = attempts
        let tampered = try VivoCanonicalJSON.decode(VivoConstantPHNCMCCheckpoint.self,
            from: JSONSerialization.data(withJSONObject: encoded))
        #expect(throws: (any Error).self) {
            try VivoConstantPHNCMC.validate(checkpoint: tampered, configuration: cfg,
                manifold: manifold, hamiltonians: hamiltonians, switchEngineFingerprint: engineID)
        }
    }

    @Test func rejectedHighWorkCandidateRollsBackSwitchEndpointButKeepsPreSwitchMD() async throws {
        let manifold = try fingerprint("rollback-manifold")
        let states = try [executable("A", manifold: manifold), executable("B", manifold: manifold)]
        let base = VivoConstantPHConfiguration(identifier: "high-work", temperatureK: 300,
            referencePH: 7, targetPH: 7, states: definitions(), initialStateIdentifier: "A",
            mdStepsPerAttempt: 2, attemptCount: 1, seed: 29)
        let cfg = VivoConstantPHNCMCConfiguration(base: base)
        let engineID = try fingerprint("high-work-engine")
        let engine = VivoConstantPHNCMCSwitchEngine(fingerprint: engineID, physicalManifoldFingerprint: manifold) { input, from, to in
            let final = try VivoConstantPHPhysicalState(stepIndex: input.stepIndex + 10,
                timePS: input.timePS + 0.010,
                positionsNM: input.positionsNM.map { $0 + .init(100, 0, 0) },
                velocitiesNMPerPS: input.velocitiesNMPerPS,
                periodicCell: input.periodicCell)
            return try .init(switchEngineFingerprint: engineID,
                physicalManifoldFingerprint: manifold,
                fromStateIdentifier: from.identifier, toStateIdentifier: to.identifier,
                initialPhysicalState: input, finalPhysicalState: final,
                protocolWorkKJPerMol: 10_000, shadowWorkKJPerMol: 0,
                switchingSteps: 10, interpretation: "synthetic high-work rejected switch")
        }
        var initial = try physical(); initial.velocitiesNMPerPS[0] = .init(0.25,-0.5,1)
        let sampler = try VivoConstantPHNCMC(configuration: cfg, initialPhysicalState: initial,
                                             executableStates: states, switchEngine: engine)
        let result = try await sampler.run()
        #expect(result.finalCheckpoint.physicalState.velocitiesNMPerPS[0] == initial.velocitiesNMPerPS[0] * -1)
        #expect(result.finalCheckpoint.currentStateIdentifier == "A")
        #expect(result.finalCheckpoint.acceptedMoves == 0)
        #expect(result.finalCheckpoint.physicalState.stepIndex == 2)
        #expect(result.finalCheckpoint.physicalState.positionsNM[0] == .init(2, 0, 0))
        #expect(result.finalCheckpoint.attempts[0].proposedPhysicalStateFingerprint != (try result.finalCheckpoint.physicalState.fingerprint()))
        #expect(result.finalCheckpoint.attempts[0].totalNonequilibriumWorkKJPerMol == 10_000)
    }

    @Test func switchResultRejectsMissingOrNonfiniteWork() async throws {
        let manifold = try fingerprint("invalid-work-manifold")
        let a = try executable("A", manifold: manifold), b = try executable("B", manifold: manifold)
        let engineID = try fingerprint("invalid-work-engine")
        let engine = VivoConstantPHNCMCSwitchEngine(fingerprint: engineID, physicalManifoldFingerprint: manifold) { input, from, to in
            try .init(switchEngineFingerprint: engineID, physicalManifoldFingerprint: manifold,
                fromStateIdentifier: from.identifier, toStateIdentifier: to.identifier,
                initialPhysicalState: input, finalPhysicalState: input,
                protocolWorkKJPerMol: .nan, shadowWorkKJPerMol: 0,
                switchingSteps: 1, interpretation: "invalid synthetic work")
        }
        let initial = try physical()
        let value = try await engine.switchCandidate(initial, a, b)
        #expect(throws: (any Error).self) {
            try value.validate(engine: engine, initial: initial, from: a, to: b,
                               maximumAbsoluteWorkKJPerMol: 1e6)
        }
    }
}
