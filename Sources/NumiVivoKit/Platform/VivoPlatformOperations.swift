import Foundation

public struct VivoWorkflowStructureImport: Codable, Sendable, Equatable {
    public var format: VivoStructureFormat
    public var policy: VivoStructureImportPolicy
    public var sourceIdentifier: String
    public init(format: VivoStructureFormat, policy: VivoStructureImportPolicy = .init(), sourceIdentifier: String) {
        self.format = format; self.policy = policy; self.sourceIdentifier = sourceIdentifier
    }
}
public struct VivoWorkflowContactAnalysis: Codable, Sendable, Equatable {
    public var selectedAtoms: [UInt32]
    public var configuration: VivoContactAnalysisConfiguration
    public init(selectedAtoms: [UInt32], configuration: VivoContactAnalysisConfiguration = .init()) {
        self.selectedAtoms = selectedAtoms; self.configuration = configuration
    }
}
public struct VivoWorkflowClassicalPreparation: Codable, Sendable, Equatable {
    public var assignment: VivoParameterAssignment
    public var options: VivoMDPreparationOptions
    public init(assignment: VivoParameterAssignment, options: VivoMDPreparationOptions = .init()) { self.assignment = assignment; self.options = options }
}
public struct VivoWorkflowHamiltonianPreparation: Codable, Sendable, Equatable {
    public var reference: VivoSCFConfiguration
    public init(reference: VivoSCFConfiguration = .init()) { self.reference = reference }
}

/// Thin adapters only. Existing production engines remain numerical authority;
/// no CLI-specific chemistry, restart format, or alternative force evaluator.
public enum VivoPlatformOperations {
    static func decode<T: Decodable>(_ type: T.Type, _ value: VivoJSONValue) throws -> T { try VivoCanonicalJSON.decode(type, from: VivoCanonicalJSON.encode(value)) }
    static func json<T: Encodable>(_ value: T) throws -> VivoJSONValue { try VivoCanonicalJSON.decode(VivoJSONValue.self, from: VivoCanonicalJSON.encode(value)) }
    static func input<T: Decodable>(_ type: T.Type, _ name: String, _ data: [String: Data]) throws -> T {
        guard let bytes = data[name] else { throw VivoChemistryError.invalid("missing operation input \(name)") }
        return try VivoCanonicalJSON.decode(type, from: bytes)
    }
    static func empty(_ value: VivoJSONValue) throws { guard value == .object([:]) else { throw VivoChemistryError.invalid("operation expects an empty configuration object") } }
    static func pure(identifier: String, id: VivoFingerprint, inputs: [String: String], outputs: [VivoChemistryTaskOutput],
                     summary: String, configure: @escaping @Sendable (VivoJSONValue) throws -> Void,
                     calculate: @escaping @Sendable (VivoJSONValue, [String: Data], VivoChemistryBudget) throws -> [String: Data]) -> VivoWorkflowDefinition {
        let operation = VivoChemistryOperation(identifier: identifier, version: "1", implementationFingerprint: id, outputs: outputs,
            execute: { cfg, data, budget in
                guard Set(data.keys) == Set(inputs.keys) else { throw VivoChemistryError.invalid("operation input ports") }
                try configure(cfg); let result = try calculate(cfg, data, budget)
                guard Set(result.keys) == Set(outputs.map(\.name)) else { throw VivoChemistryError.invalid("operation output ports") }
                return result
            }, validateOutputs: { cfg, data, actual, budget in
                guard Set(data.keys) == Set(inputs.keys), Set(actual.keys) == Set(outputs.map(\.name)) else {
                    throw VivoChemistryError.invalid("validated operation port mismatch")
                }
                try configure(cfg)
                guard try calculate(cfg, data, budget) == actual else { throw VivoChemistryError.invalid("platform output reconstruction differs") }
            })
        return .init(operation: operation, inputKinds: inputs, summary: summary,
            validationScope: "native input validation and deterministic output reconstruction; scientific scope remains the wrapped engine's scope",
            validateConfiguration: configure)
    }
    public static func registry(implementationFingerprint id: VivoFingerprint) throws -> VivoWorkflowRegistry {
        var definitions = foundations(id)
        let base = VivoElectronicWorkflowOperations.self, advanced = VivoAdvancedChemistryOperations.self
        func add(_ operation: VivoChemistryOperation, _ inputs: [String: String], _ summary: String,
                 _ validation: @escaping @Sendable (VivoJSONValue) throws -> Void) {
            definitions.append(.init(operation: operation, inputKinds: inputs, summary: summary,
                validationScope: "existing native operation output validator; interface availability is not enlarged-system or chemical validation",
                validateConfiguration: validation))
        }
        add(base.integrals(implementationFingerprint: id), ["system": "vivo.electronic-system", "basis": "vivo.gaussian-basis"],
            "All-electron normalized Cartesian Gaussian AO integrals.", empty)
        add(base.hartreeFock(implementationFingerprint: id), ["integrals": "vivo.ao-integrals"], "Restricted or unrestricted Hartree-Fock.") {
            try decode(VivoSCFConfiguration.self, $0).validate()
        }
        add(base.hamiltonian(implementationFingerprint: id), ["integrals": "vivo.ao-integrals", "scf": "vivo.hartree-fock"],
            "Common RHF-frame embedded Hamiltonian with scalar and source identities preserved.") {
            let settings = try decode(VivoSCFConfiguration.self, $0); try settings.validate()
            guard settings.reference == .restricted else { throw VivoChemistryError.unsupported("RHF Hamiltonian requires a restricted reference") }
        }
        add(base.manyBody(implementationFingerprint: id), ["hamiltonian": "vivo.embedded-hamiltonian"], "Explicit common classical many-body solver selection.") {
            _ = try decode(VivoManyBodySolverRequest.self, $0)
        }
        add(advanced.manyBody(implementationFingerprint: id), ["hamiltonian": "vivo.embedded-hamiltonian"], "Direct CI, tensor CCSD or state-averaged CASSCF.") {
            _ = try decode(VivoAdvancedManyBodyRequest.self, $0)
        }
        add(advanced.eccDMET(implementationFingerprint: id), ["hamiltonian": "vivo.embedded-hamiltonian"], "Existing bounded ECC-DMET implementation and validator.") {
            _ = try decode(VivoECCDMETConfiguration.self, $0)
        }
        add(advanced.smoothSolvent(implementationFingerprint: id), ["integrals": "vivo.ao-integrals"], "Smooth variational C-PCM with the existing restricted electronic reference.") {
            let cfg = try decode(VivoAdvancedChemistryOperations.SmoothSolventConfiguration.self, $0); try cfg.solvent.validate(); try cfg.scf.validate()
        }
        add(VivoECCPathWorkflowOperations.sharedPath(implementationFingerprint: id), ["request": "vivo.molecular-ecc-path-request"],
            "Shared-orbital native molecular ECC electronic profile; not thermal chemistry or a kinetic rate.", empty)
        add(VivoNuclearWorkflowOperations.qualifyStationaryPoint(implementationFingerprint: id), ["request": "vivo.nuclear-point-request"],
            "Native optimization, finite-difference nuclear derivatives, Hessian index and explicit harmonic thermochemistry.", empty)
        add(VivoNuclearWorkflowOperations.characterize(implementationFingerprint: id), ["request": "vivo.reaction-characterization-request"],
            "Stationary-point reaction characterization through the existing qualified geometry contract.", empty)
        add(VivoReactionConnectivityOperations.validate(implementationFingerprint: id), ["request": "vivo.mapped-reaction-connectivity-request"],
            "Mapped endpoint/descent connectivity, preserving atom and model identities.", empty)
        add(VivoReactionConvergenceOperations.barrier(implementationFingerprint: id), ["request": "vivo.barrier-convergence-request"],
            "Explicit basis/solvent/embedding/nuclear barrier-refinement evidence.", empty)
        add(VivoReactionConvergenceOperations.compareSolvationCycle(implementationFingerprint: id), ["request": "vivo.solvation-cycle-comparison-request"],
            "Declared reaction-cycle contribution comparison.", empty)
        definitions += VivoPlatformMechanismOperations.definitions(implementationFingerprint: id)
        definitions += VivoPlatformMDOperations.definitions(implementationFingerprint: id)
        definitions += VivoPlatformQMMMOperations.definitions(implementationFingerprint: id)
        definitions += VivoPlatformQMMMFreeEnergyOperations.definitions(implementationFingerprint: id)
        definitions += VivoPlatformRefinementOperations.definitions(implementationFingerprint: id)
        definitions += VivoPlatformAdaptiveOperations.definitions(implementationFingerprint: id)
        return try .init(implementationFingerprint: id, definitions: definitions)
    }
    private static func foundations(_ id: VivoFingerprint) -> [VivoWorkflowDefinition] {
        [
            pure(identifier: "vivo.platform.structure-import", id: id, inputs: ["source": "vivo.structure-text"],
                outputs: [.init(name: "structure", kind: "vivo.molecular-structure")], summary: "PDB/mmCIF/SDF/MOL2 through the shared importer, retaining source and diagnostics.",
                configure: { value in
                    let cfg = try decode(VivoWorkflowStructureImport.self, value)
                    guard !cfg.sourceIdentifier.isEmpty, cfg.sourceIdentifier.utf8.count <= 4096 else { throw VivoChemistryError.invalid("structure source identifier") }
                }, calculate: { value, data, _ in
                    let cfg = try decode(VivoWorkflowStructureImport.self, value)
                    guard let source = data["source"], let text = String(data: source, encoding: .utf8) else { throw VivoChemistryError.invalid("structure text must be UTF-8") }
                    let imported = try VivoStructureImporter(policy: cfg.policy).importText(text, format: cfg.format, sourceIdentifier: cfg.sourceIdentifier)
                    return ["structure": try VivoCanonicalJSON.encode(imported.structure)]
                }),
            pure(identifier: "vivo.platform.interaction-geometry", id: id, inputs: ["structure": "vivo.molecular-structure"],
                outputs: [.init(name: "analysis", kind: "vivo.interaction-analysis")], summary: "Existing selected-atom contact and distance analysis.",
                configure: { value in try decode(VivoWorkflowContactAnalysis.self, value).configuration.validate() }, calculate: { value, data, _ in
                    let cfg = try decode(VivoWorkflowContactAnalysis.self, value)
                    let structure = try input(VivoMolecularStructure.self, "structure", data)
                    return ["analysis": try VivoCanonicalJSON.encode(VivoStructureAnalyzer(configuration: cfg.configuration).analyze(structure, selectedAtomIndices: cfg.selectedAtoms))]
                }),
            pure(identifier: "vivo.platform.classical-prepare", id: id, inputs: ["structure": "vivo.molecular-structure", "parameters": "vivo.classical-parameter-set"],
                outputs: [.init(name: "system", kind: "vivo.classical-system")], summary: "Prepare the existing parameterized MD system without native-parameter substitution.",
                configure: { _ = try decode(VivoWorkflowClassicalPreparation.self, $0) }, calculate: { value, data, _ in
                    let cfg = try decode(VivoWorkflowClassicalPreparation.self, value)
                    let structure = try input(VivoMolecularStructure.self, "structure", data), parameters = try input(VivoClassicalParameterSet.self, "parameters", data)
                    return ["system": try VivoCanonicalJSON.encode(VivoMDPreparation.prepare(structure, parameters: parameters, assignment: cfg.assignment, options: cfg.options))]
                }),
            pure(identifier: "vivo.platform.chemical-state-free-energy", id: id, inputs: ["request": "vivo.chemical-state-free-energy-request"],
                outputs: [.init(name: "ensemble", kind: "vivo.chemical-state-free-energy-result")], summary: "Stable chemical-state log-sum-exp with contribution provenance and group aggregation.",
                configure: empty, calculate: { _, data, _ in
                    let request = try input(VivoChemicalStateFreeEnergyRequest.self, "request", data)
                    return ["ensemble": try VivoCanonicalJSON.encode(VivoChemicalStateFreeEnergy.calculate(request))]
                }),
            pure(identifier: "vivo.platform.event-qualification", id: id, inputs: ["request": "vivo.event-qualification-request"],
                outputs: [.init(name: "qualification", kind: "vivo.event-qualification-result")], summary: "Replicated escape/association event qualification with censoring and milestone controls.",
                configure: empty, calculate: { _, data, _ in
                    let request = try input(VivoEventQualificationRequest.self, "request", data)
                    return ["qualification": try VivoCanonicalJSON.encode(VivoEventQualification.qualify(request))]
                }),
            pure(identifier: "vivo.platform.free-energy-qualification", id: id, inputs: ["request": "vivo.free-energy-qualification-request"],
                outputs: [.init(name: "qualification", kind: "vivo.free-energy-qualification-result")], summary: "Existing replicated free-energy overlap and uncertainty qualification.",
                configure: empty, calculate: { _, data, _ in
                    let request = try input(VivoFreeEnergyQualificationRequest.self, "request", data)
                    return ["qualification": try VivoCanonicalJSON.encode(VivoFreeEnergyQualification.qualify(request))]
                }),
            pure(identifier: "vivo.platform.kinetic-pack", id: id, inputs: ["pack": "vivo.kinetic-pack"],
                outputs: [.init(name: "compiled", kind: "vivo.compiled-kinetics")], summary: "Validate and compile context/evidence-qualified kinetic packs, without fabricating rates.",
                configure: empty, calculate: { _, data, _ in
                    let pack = try input(VivoKineticPack.self, "pack", data)
                    return ["compiled": try VivoCanonicalJSON.encode(VivoKineticCompiler.compile(pack))]
                }),
            pure(identifier: "vivo.platform.target-occupancy", id: id, inputs: ["pack": "vivo.kinetic-pack", "experiment": "vivo.target-occupancy-experiment"],
                outputs: [.init(name: "simulation", kind: "vivo.target-occupancy-simulation")], summary: "Native coupled binding/covalent occupancy dynamics from the declared pack.",
                configure: empty, calculate: { _, data, _ in
                    let pack = try input(VivoKineticPack.self, "pack", data), experiment = try input(VivoTargetOccupancyExperiment.self, "experiment", data)
                    let compiled = try VivoKineticCompiler.compile(pack)
                    return ["simulation": try VivoCanonicalJSON.encode(VivoTargetOccupancySimulator().run(experiment, compiled: compiled))]
                })
        ]
    }
}
