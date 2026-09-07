import Foundation
import Testing
@testable import NumiVivoKit

/// Identity/admission tests only. Factories construct closures, but no Metal
/// runtime, physical propagation, electronic evaluator or stochastic run starts.
@Suite(.serialized) struct MetalSamplingIdentityTests {
    private func system() throws -> VivoClassicalSystem {
        try .init(identifier: "metal-sampling-identity",structureFingerprint:
            VivoCanonicalJSON.fingerprint(Data("metal-sampling-identity-structure".utf8)),particles: [
                .init(index: 0,atomIndex: 0,typeIdentifier: "X",massDa: 12,chargeE: 0,sigmaNM: 0,epsilonKJPerMol: 0)])
    }
    private func dynamics() -> VivoMDConfiguration {
        .init(timeStepPS: 0.001,cutoffNM: 1,neighborSkinNM: 0,electrostatics: .cutoff,
            ensemble: .nve,thermostat: .none,targetTemperatureK: nil,frictionPerPS: nil,
            neighborListEnabled: false,randomSeed: 991)
    }
    private func physical() throws -> VivoConstantPHPhysicalState {
        try .init(stepIndex: 4,timePS: 0.004,positionsNM: [.zero],velocitiesNMPerPS: [.zero],periodicCell: nil)
    }
    private func specifications() throws -> [VivoConstantPHMetalStateSpecification] {
        let system = try system(), configuration = dynamics()
        return ["A","B"].map { .init(identifier: $0,system: system,configuration: configuration) }
    }
    private func sampling() -> VivoConstantPHConfiguration {
        let evidence = VivoKineticEvidence(source: "synthetic identity fixture",locator: "MetalSamplingIdentityTests")
        let states: [VivoConstantPHStateDefinition] = [
            .init(identifier: "A",boundProtonOffset: 0,referenceSemigrandBiasKJPerMol: 0,
                  origin: .assumed,evidence: evidence,neighbors: ["B"]),
            .init(identifier: "B",boundProtonOffset: 1,referenceSemigrandBiasKJPerMol: 0,
                  origin: .assumed,evidence: evidence,neighbors: ["A"])]
        return .init(identifier: "identity-resume",temperatureK: 300,referencePH: 7,targetPH: 7,
            states: states,initialStateIdentifier: "A",mdStepsPerAttempt: 1,attemptCount: 2,seed: 71)
    }
    // A nil field encodes the exact old identity: synthesized Encodable omits
    // absent optional keys. This is a reconstructed historical hash, not a
    // random wrong fingerprint that could miss omission of the contract.
    private func constantPHIdentity(_ specification: VivoConstantPHMetalStateSpecification,
                                    contract: String?) throws -> VivoFingerprint {
        struct Identity: Encodable {
            let schema: String; let numericalContract: String?
            let system: VivoFingerprint; let execution: VivoFingerprint; let provider: VivoFingerprint?
        }
        return try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(Identity(
            schema: "numivivo.org/constant-ph-metal-hamiltonian/v1",numericalContract: contract,
            system: specification.system.fingerprint(),execution: VivoMDCandidateForceProvider.executionFingerprint(
                configuration: specification.configuration,provider: specification.forceProvider),
            provider: specification.forceProvider?.fingerprint)))
    }
    private func oldStates(_ current: [VivoConstantPHExecutableState],
                           specifications: [VivoConstantPHMetalStateSpecification]) throws -> [VivoConstantPHExecutableState] {
        try zip(current,specifications).map { state,specification in
            .init(identifier: state.identifier,physicalManifoldFingerprint: state.physicalManifoldFingerprint,
                  hamiltonianFingerprint: try constantPHIdentity(specification,contract: nil),
                  potentialEnergyKJPerMol: { _ in throw VivoChemistryError.invalid("unexpected identity-test evaluation") },
                  propagate: { _,_ in throw VivoChemistryError.invalid("unexpected identity-test propagation") })
        }
    }

    @Test func constantPHRejectsHistoricalMetalCheckpointBeforeRestampingPhysicalState() async throws {
        let specs = try specifications(), config = sampling(), input = try physical()
        let current = try VivoConstantPHMetalStateFactory.make(specifications: specs,samplingTemperatureK: 300)
        let legacy = try oldStates(current,specifications: specs)
        for (state,specification) in zip(current,specs) {
            #expect(state.hamiltonianFingerprint == (try constantPHIdentity(specification,contract: VivoMDExecutionIdentity.current)))
            #expect(state.hamiltonianFingerprint != (try constantPHIdentity(specification,contract: nil)))
            #expect(state.hamiltonianFingerprint != (try constantPHIdentity(specification,contract: "retired-MD-contract")))
        }
        // Generic CPU/external states remain valid under their own supplied IDs.
        let oldSampler = try VivoConstantPHDiscreteMD(configuration: config,initialPhysicalState: input,executableStates: legacy)
        let oldCheckpoint = try VivoCanonicalJSON.decode(VivoConstantPHCheckpoint.self,
            from: VivoCanonicalJSON.encode(await oldSampler.checkpoint()))
        _ = try VivoConstantPHDiscreteMD(configuration: config,initialPhysicalState: input,
            executableStates: legacy,checkpoint: oldCheckpoint)
        #expect(throws: (any Error).self) {
            _ = try VivoConstantPHDiscreteMD(configuration: config,initialPhysicalState: input,
                executableStates: current,checkpoint: oldCheckpoint)
        }
        let newSampler = try VivoConstantPHDiscreteMD(configuration: config,initialPhysicalState: input,executableStates: current)
        let newCheckpoint = await newSampler.checkpoint()
        _ = try VivoConstantPHDiscreteMD(configuration: config,initialPhysicalState: input,
            executableStates: current,checkpoint: newCheckpoint)
    }

    @Test func ncmcInheritsChangedEndpointAndEngineIdentities() async throws {
        let specs = try specifications(), config = VivoConstantPHNCMCConfiguration(base: sampling()), input = try physical()
        let settings = VivoNCMCMetalConfiguration(schedule: .init(lambdas: [0,0.5,1],propagationStepsPerLambda: 1),timeStepPS: 0.0005)
        let setup = try VivoConstantPHNCMCMetalFactory.make(specifications: specs,samplingTemperatureK: 300,configuration: settings)
        let legacy = try oldStates(setup.executableStates,specifications: specs)
        let manifold = try #require(legacy.first?.physicalManifoldFingerprint)
        var switching = dynamics(); switching.timeStepPS = settings.timeStepPS
        struct EngineIdentity: Encodable {
            let schema: String; let endpoints: [String: VivoFingerprint]; let manifold: VivoFingerprint
            let settings: VivoNCMCMetalConfiguration; let switching: VivoMDConfiguration
        }
        func engineID(_ states: [VivoConstantPHExecutableState]) throws -> VivoFingerprint {
            try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(EngineIdentity(
                schema: "numivivo.org/ncmc-metal-engine/v1",
                endpoints: Dictionary(uniqueKeysWithValues: states.map { ($0.identifier,$0.hamiltonianFingerprint) }),
                manifold: manifold,settings: settings,switching: switching)))
        }
        let oldEngineID = try engineID(legacy)
        #expect(setup.switchEngine.fingerprint == (try engineID(setup.executableStates)))
        #expect(setup.switchEngine.fingerprint != oldEngineID)
        let oldEngine = VivoConstantPHNCMCSwitchEngine(fingerprint: oldEngineID,physicalManifoldFingerprint: manifold) { _,_,_ in
            throw VivoChemistryError.invalid("unexpected identity-test switching")
        }
        let oldSampler = try VivoConstantPHNCMC(configuration: config,initialPhysicalState: input,
            executableStates: legacy,switchEngine: oldEngine)
        let oldCheckpoint = try VivoCanonicalJSON.decode(VivoConstantPHNCMCCheckpoint.self,
            from: VivoCanonicalJSON.encode(await oldSampler.checkpoint()))
        _ = try VivoConstantPHNCMC(configuration: config,initialPhysicalState: input,
            executableStates: legacy,switchEngine: oldEngine,checkpoint: oldCheckpoint)
        // Isolate each guard, so neither identity can accidentally stop binding.
        #expect(throws: (any Error).self) {
            _ = try VivoConstantPHNCMC(configuration: config,initialPhysicalState: input,
                executableStates: setup.executableStates,switchEngine: oldEngine,checkpoint: oldCheckpoint)
        }
        #expect(throws: (any Error).self) {
            _ = try VivoConstantPHNCMC(configuration: config,initialPhysicalState: input,
                executableStates: legacy,switchEngine: setup.switchEngine,checkpoint: oldCheckpoint)
        }
        let newSampler = try VivoConstantPHNCMC(configuration: config,initialPhysicalState: input,
            executableStates: setup.executableStates,switchEngine: setup.switchEngine)
        let newCheckpoint = await newSampler.checkpoint()
        _ = try VivoConstantPHNCMC(configuration: config,initialPhysicalState: input,
            executableStates: setup.executableStates,switchEngine: setup.switchEngine,checkpoint: newCheckpoint)
    }

    private func nuclearSpecification() throws -> VivoNuclearMetalSpecification {
        let system = try system()
        return .init(system: system,initialState: .init(systemFingerprint: try system.fingerprint(),positionsNM: [.zero]),
            dynamics: dynamics())
    }
    private func nuclearDefinition(_ current: VivoNuclearPotentialDefinition,
                                   identity: VivoFingerprint) throws -> VivoNuclearPotentialDefinition {
        try .init(hamiltonianFingerprint: identity,atomIndices: current.atomIndices,particleIndices: current.particleIndices,
            massesDa: current.massesDa,periodicCell: current.periodicCell,periodicMoleculeGroups: current.periodicMoleculeGroups,
            coordinateEvaluation: current.coordinateEvaluation)
    }
    @Test func nuclearSpecificationRejectsOldRingCheckpointAndForceLabels() throws {
        let specification = try nuclearSpecification(), current = try specification.definition()
        struct Identity: Encodable {
            let schema: String; let numericalContract: String?; let specification: VivoNuclearMetalSpecification
        }
        func identity(_ contract: String?) throws -> VivoFingerprint {
            try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(Identity(
                schema: "numivivo.org/nuclear-metal-complete-potential/v1;explicit-FP32-coordinate-projection",
                numericalContract: contract,specification: specification)))
        }
        #expect(current.hamiltonianFingerprint == (try identity(VivoMDExecutionIdentity.current)))
        let legacy = try nuclearDefinition(current,identity: identity(nil))
        #expect(current != legacy)
        let config = VivoRingPolymerConfiguration(temperatureK: 300,beadCount: 2,integrationSteps: 1)
        let oldCheckpoint = try VivoRingPolymerCheckpoint(definition: legacy,configuration: config,seed: 13,
            beadPositionsNM: [[.zero],[.zero]])
        // The old artifact remains self-consistent; only current-Metal admission
        // rejects it. Generic numerical-potential schemas are not globally gated.
        try oldCheckpoint.validate()
        #expect(throws: (any Error).self) {
            try VivoRingPolymerWorkflowRequest(potential: specification,checkpoint: oldCheckpoint,sweeps: 1).admit(.init())
        }
        let currentCheckpoint = try VivoRingPolymerCheckpoint(definition: current,configuration: config,seed: 13,
            beadPositionsNM: oldCheckpoint.beadPositionsNM)
        try VivoRingPolymerWorkflowRequest(potential: specification,checkpoint: currentCheckpoint,sweeps: 1).admit(.init())
        let oldLabel = try VivoNuclearPotentialEvaluation(definition: legacy,positionsNM: [.zero],energyKJPerMol: 0,
            forcesKJPerMolNM: [.zero],normalizedConvergenceResidual: 0)
        try oldLabel.validate(definition: legacy,positionsNM: [.zero])
        #expect(throws: (any Error).self) { try oldLabel.validate(definition: current,positionsNM: [.zero]) }
        let currentLabel = try VivoNuclearPotentialEvaluation(definition: current,positionsNM: [.zero],energyKJPerMol: 0,
            forcesKJPerMolNM: [.zero],normalizedConvergenceResidual: 0)
        try currentLabel.validate(definition: current,positionsNM: [.zero])
    }

    @Test func directNuclearMetalIdentityBindsTheCurrentEvaluator() throws {
        let system = try system(), configuration = dynamics()
        struct Identity: Encodable {
            let system: VivoFingerprint; let execution: VivoFingerprint
            let projection: String; let numericalContract: String?
        }
        func identity(_ contract: String?) throws -> VivoFingerprint {
            try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(Identity(system: system.fingerprint(),
                execution: VivoMDCandidateForceProvider.executionFingerprint(configuration: configuration,provider: nil),
                projection: "explicit-nearest-FP32-coordinate-projection;complete-Metal-plus-BO-potential/v1",
                numericalContract: contract)))
        }
        let current = try VivoNuclearMetalPotentialIdentity.fingerprint(system: system,configuration: configuration,provider: nil)
        #expect(current == (try identity(VivoMDExecutionIdentity.current)))
        #expect(current != (try identity(nil)))
        #expect(current != (try identity("retired-MD-contract")))
    }
}
