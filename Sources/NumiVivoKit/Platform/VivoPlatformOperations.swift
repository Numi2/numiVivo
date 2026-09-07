import Foundation

public struct VivoWorkflowStructureImport: Codable, Sendable, Equatable {
    public var format: VivoStructureFormat
    public var identifier: String
    public init(format: VivoStructureFormat, identifier: String = "molecule") { self.format = format; self.identifier = identifier }
}
public struct VivoWorkflowStructureSlice: Codable, Sendable, Equatable {
    public var atomIndices: [UInt32]
    public init(atomIndices: [UInt32]) { self.atomIndices = atomIndices }
}
public struct VivoWorkflowElectronicSystem: Codable, Sendable, Equatable {
    public var conformerIndex: Int
    public var alphaElectrons: Int
    public var betaElectrons: Int
    public init(conformerIndex: Int = 0, alphaElectrons: Int, betaElectrons: Int) {
        self.conformerIndex = conformerIndex; self.alphaElectrons = alphaElectrons; self.betaElectrons = betaElectrons
    }
}

/// Explicit allowlist of adapters to existing scientific authorities. Catalog
/// membership means the interface is executable, never that chemistry is certified.
public enum VivoPlatformOperations {
    static func decode<T: Decodable>(_ type: T.Type, _ cfg: VivoJSONValue) throws -> T {
        try VivoCanonicalJSON.decode(type, from: VivoCanonicalJSON.encode(cfg))
    }
    static func input<T: Decodable>(_ type: T.Type, _ name: String, _ values: [String: Data]) throws -> T {
        guard let bytes = values[name] else { throw VivoChemistryError.invalid("missing platform input \(name)") }
        return try VivoCanonicalJSON.decode(type, from: bytes)
    }
    static func empty(_ cfg: VivoJSONValue) throws {
        guard cfg == .object([:]) else { throw VivoChemistryError.invalid("this operation has no configuration fields") }
    }
    static func json<T: Encodable>(_ value: T) throws -> VivoJSONValue {
        try VivoCanonicalJSON.decode(VivoJSONValue.self, from: VivoCanonicalJSON.encode(value))
    }
    static func pure(identifier: String, id: VivoFingerprint, inputs: [String: String], outputs: [VivoChemistryTaskOutput],
                     summary: String, configure: @escaping @Sendable (VivoJSONValue) throws -> Void,
                     calculate: @escaping @Sendable (VivoJSONValue, [String: Data], VivoChemistryBudget) throws -> [String: Data]) -> VivoWorkflowDefinition {
        let operation = VivoChemistryOperation(identifier: identifier, version: "1", implementationFingerprint: id, outputs: outputs,
            execute: { cfg, data, budget in
                guard Set(data.keys) == Set(inputs.keys) else { throw VivoChemistryError.invalid("platform input slots") }
                try configure(cfg); return try calculate(cfg, data, budget)
            }, validateOutputs: { cfg, data, actual, budget in
                guard Set(data.keys) == Set(inputs.keys), Set(actual.keys) == Set(outputs.map(\.name)) else {
                    throw VivoChemistryError.invalid("platform output slots")
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
        add(advanced.manyBody(implementationFingerprint: id), ["hamiltonian": "vivo.embedded-hamiltonian"],
            "Matrix-free CI, tensor CCSD or multistate CASSCF through their existing solver interface.") {
            _ = try decode(VivoAdvancedManyBodyRequest.self, $0)
        }
        add(advanced.eccDMET(implementationFingerprint: id), ["hamiltonian": "vivo.embedded-hamiltonian"],
            "Existing bounded ECC-DMET implementation and validator.") {
            _ = try decode(VivoECCDMETConfiguration.self, $0)
        }
        add(advanced.ri(implementationFingerprint: id), ["system": "vivo.electronic-system", "basis": "vivo.gaussian-basis"],
            "Native RI integrals with an explicit auxiliary basis.") { _ = try decode(VivoRIWorkflowConfiguration.self, $0) }
        add(advanced.riHF(implementationFingerprint: id), ["integrals": "vivo.ri-integrals"], "Factorized HF for the RI-MP2 path.") {
            try decode(VivoSCFConfiguration.self, $0).validate()
        }
        add(advanced.riMP2(implementationFingerprint: id), ["reference": "vivo.ri-hf"], "Factorized MP2 without a complete MO ERI tensor.", empty)
        add(advanced.smoothSCF(implementationFingerprint: id), ["integrals": "vivo.ao-integrals"], "Existing smooth-CPCM SCF operation.") {
            _ = try decode(VivoSmoothSCFWorkflowConfiguration.self, $0)
        }
        add(VivoReactionQualificationWorkflow.operation(implementationFingerprint: id), ["request": "vivo.reaction-calculation-request"],
            "General nuclear, solvent, embedding and reaction requests; no paper data required.", empty)

        definitions.append(pure(identifier: "vivo.platform.target-reference", id: id,
            inputs: ["experiment": "vivo.target-engagement-experiment"],
            outputs: [.init(name: "result", kind: "vivo.target-engagement-result")],
            summary: "Exposure-driven target fractions through the existing FP64 reference; assumptions retained.",
            configure: { try decode(VivoTargetEngagementNumerics.self, $0).validate() }, calculate: { cfg, inputs, _ in
                let experiment = try input(VivoTargetEngagementExperiment.self, "experiment", inputs)
                let result = try VivoTargetEngagementReference.run(experiment, numerics: decode(VivoTargetEngagementNumerics.self, cfg))
                return ["result": try VivoCanonicalJSON.encode(result)]
            }))
        definitions.append(pure(identifier: "vivo.platform.finite-drug", id: id,
            inputs: ["experiment": "vivo.finite-drug-experiment"],
            outputs: [.init(name: "result", kind: "vivo.finite-drug-result")],
            summary: "Existing mass-balanced finite-drug reaction operator; no inferred pharmacological rate.",
            configure: empty, calculate: { _, inputs, _ in
                ["result": try VivoCanonicalJSON.encode(VivoFiniteDrugRunRecord.run(input(VivoFiniteDrugExperiment.self, "experiment", inputs)))]
            }))
        definitions += VivoPlatformMechanismOperations.definitions(implementationFingerprint: id)
        definitions += VivoPlatformMDOperations.definitions(implementationFingerprint: id)
        definitions.append(VivoPlatformSnapshotOperations.definition(implementationFingerprint: id))
        definitions += VivoPlatformQMMMOperations.definitions(implementationFingerprint: id)
        // QM/MM module owns its free-energy adapters; register each authority once.
        definitions += VivoPlatformRefinementOperations.definitions(implementationFingerprint: id)
        definitions += VivoPlatformAdaptiveOperations.definitions(implementationFingerprint: id)
        return try .init(implementationFingerprint: id, definitions: definitions)
    }

    private static func foundations(_ id: VivoFingerprint) -> [VivoWorkflowDefinition] {
        var definitions: [VivoWorkflowDefinition] = []
        definitions.append(pure(identifier: "vivo.platform.structure-import", id: id, inputs: ["source": "vivo.structure-source-text"],
            outputs: [.init(name: "structure", kind: "vivo.molecular-structure-document")], summary: "Single-structure text import with explicit format and source mapping.",
            configure: { cfg in
                let c = try decode(VivoWorkflowStructureImport.self, cfg)
                guard !c.identifier.isEmpty, c.identifier.utf8.count <= 1024 else { throw VivoChemistryError.invalid("structure import identity") }
            }, calculate: { cfg, inputs, _ in
                let c = try decode(VivoWorkflowStructureImport.self, cfg), text = try input(String.self, "source", inputs)
                if c.format == .sdf, text.components(separatedBy: "$$$$").filter({ !$0.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty }).count > 1 {
                    throw VivoChemistryError.unsupported("multiple SDF records require separate workflow inputs; first-record truncation is not automatic")
                }
                if c.format == .smiles, text.split(whereSeparator: \.isNewline).filter({ !$0.trimmingCharacters(in: .whitespaces).isEmpty }).count != 1 {
                    throw VivoChemistryError.unsupported("multi-record SMILES requires separate workflow inputs")
                }
                if c.format == .mol2, text.split(whereSeparator: \.isNewline).filter({ $0.trimmingCharacters(in: .whitespaces) == "@<TRIPOS>MOLECULE" }).count > 1 {
                    throw VivoChemistryError.unsupported("multi-record MOL2 requires separate workflow inputs")
                }
                let imported = try VivoStructureCodec.decode(Data(text.utf8), format: c.format, identifier: c.identifier)
                let document = try VivoMolecularStructureDocument(structure: imported.structure,
                    sourceFingerprint: imported.sourceFingerprint, sourceFormat: imported.format)
                return ["structure": try document.canonicalData()]
            }))
        definitions.append(pure(identifier: "vivo.platform.structure-slice", id: id,
            inputs: ["structure": "vivo.molecular-structure-document"],
            outputs: [.init(name: "structure", kind: "vivo.molecular-structure-document"), .init(name: "mapping", kind: "vivo.structure-rewrite")],
            summary: "Mapped atom subset with explicit forward/reverse identity maps.", configure: { cfg in
                let c = try decode(VivoWorkflowStructureSlice.self, cfg)
                guard !c.atomIndices.isEmpty, Set(c.atomIndices).count == c.atomIndices.count else { throw VivoChemistryError.invalid("slice requires unique atom indices") }
            }, calculate: { cfg, inputs, _ in
                let source = try input(VivoMolecularStructureDocument.self, "structure", inputs); try source.validate()
                let rewrite = try VivoStructureSlicer.slice(source.structure, atomIndices: decode(VivoWorkflowStructureSlice.self, cfg).atomIndices)
                try rewrite.validate(originalAtomCount: source.structure.atoms.count)
                let document = try VivoMolecularStructureDocument(structure: rewrite.structure, sourceFingerprint: source.structureFingerprint)
                return ["structure": try document.canonicalData(), "mapping": try VivoCanonicalJSON.encode(rewrite)]
            }))
        definitions.append(pure(identifier: "vivo.platform.structure-electronic-system", id: id,
            inputs: ["structure": "vivo.molecular-structure-document"], outputs: [.init(name: "system", kind: "vivo.electronic-system")],
            summary: "Explicit conformer and electron sector to mapped nuclear coordinates; nanometres converted to Bohr.", configure: { cfg in
                let c = try decode(VivoWorkflowElectronicSystem.self, cfg)
                guard c.conformerIndex >= 0, c.alphaElectrons >= 0, c.betaElectrons >= 0 else { throw VivoChemistryError.invalid("electronic-system conformer or electron sector") }
            }, calculate: { cfg, inputs, budget in
                let c = try decode(VivoWorkflowElectronicSystem.self, cfg), source = try input(VivoMolecularStructureDocument.self, "structure", inputs)
                try source.validate()
                let structure = source.structure
                guard structure.periodicCell == nil else { throw VivoChemistryError.unsupported("periodic structure needs an explicitly prepared finite QM/MM region") }
                guard structure.atoms.count <= budget.maximumBasisFunctions, c.conformerIndex < structure.conformers.count else {
                    throw VivoChemistryError.invalid("structure lacks the selected explicit conformer or exceeds the nucleus budget")
                }
                let positions = structure.conformers[c.conformerIndex].positionsNM
                let system = VivoElectronicSystem(nuclei: structure.atoms.map { atom in
                    let p = positions[Int(atom.index)]
                    return .init(atomicNumber: Int(atom.element.atomicNumber),
                        positionBohr: SIMD3<Double>(p.x, p.y, p.z) / VivoAtomicUnits.bohrInNM, structureAtomIndex: atom.index)
                }, alphaElectrons: c.alphaElectrons, betaElectrons: c.betaElectrons)
                try system.validate()
                return ["system": try VivoCanonicalJSON.encode(system)]
            }))
        return definitions
    }
}
