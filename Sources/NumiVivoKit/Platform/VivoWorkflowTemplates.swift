import Foundation

public enum VivoWorkflowTemplates {
    /// Small reusable examples, not a specific publication's molecular inputs.
    public static func make(_ name: String) throws -> VivoWorkflowRecipe {
        switch name {
        case "molecular-analysis": return try molecularAnalysis()
        case "md-segments": return try molecularDynamics()
        case "md-electronic-analysis": return try mdElectronicAnalysis()
        default: throw VivoChemistryError.invalid("unknown workflow template: \(name)")
        }
    }
    public static func molecularAnalysis() throws -> VivoWorkflowRecipe {
        let hydrogen = VivoElement.from(symbol: "H")!
        let structure = VivoMolecularStructure(identifier: "hydrogen-example", atoms: [
            .init(index: 0, name: "H1", element: hydrogen), .init(index: 1, name: "H2", element: hydrogen)
        ], bonds: [.init(atomA: 0, atomB: 1)], conformers: [.init(positionsNM: [
            .init(0, 0, -0.7 * VivoAtomicUnits.bohrInNM), .init(0, 0, 0.7 * VivoAtomicUnits.bohrInNM)
        ])])
        let source = try VivoMolecularStructureDocument(structure: structure)
        let hf = try VivoPlatformOperations.json(VivoSCFConfiguration())
        return try .init(identifier: "molecular-analysis", artifacts: [
            .init(identifier: "molecule", source: .json(kind: "vivo.molecular-structure-document", payload: VivoPlatformOperations.json(source))),
            .init(identifier: "basis", source: .json(kind: "vivo.gaussian-basis", payload: VivoPlatformOperations.json(VivoGaussianBasis.hydrogenSTO3G(nucleusIndices: [0, 1]))))
        ], nodes: [
            .init(identifier: "system", operation: "vivo.platform.structure-electronic-system", inputs: ["structure": .artifact(identifier: "molecule")],
                configuration: VivoPlatformOperations.json(VivoWorkflowElectronicSystem(alphaElectrons: 1, betaElectrons: 1))),
            .init(identifier: "integrals", operation: "vivo.native.ao-integrals", inputs: ["system": .output(node: "system", port: "system"), "basis": .artifact(identifier: "basis")]),
            .init(identifier: "reference", operation: "vivo.native.hf-from-integrals", inputs: ["integrals": .output(node: "integrals", port: "integrals")], configuration: hf),
            .init(identifier: "hamiltonian", operation: "vivo.native.rhf-hamiltonian", inputs: [
                "integrals": .output(node: "integrals", port: "integrals"), "scf": .output(node: "reference", port: "scf")], configuration: hf),
            .init(identifier: "fci", operation: "vivo.native.many-body", inputs: ["hamiltonian": .output(node: "hamiltonian", port: "hamiltonian")],
                configuration: VivoPlatformOperations.json(VivoManyBodySolverRequest.configurationInteraction(method: .fci))),
            .init(identifier: "mp2", operation: "vivo.native.many-body", inputs: ["hamiltonian": .output(node: "hamiltonian", port: "hamiltonian")],
                configuration: VivoPlatformOperations.json(VivoManyBodySolverRequest.mp2(minimumGapHartree: 1e-8)))
        ], outputs: [.init(name: "fci", node: "fci", port: "manyBody"), .init(name: "mp2", node: "mp2", port: "manyBody"),
                     .init(name: "hamiltonian", node: "hamiltonian", port: "hamiltonian")])
    }
    public static func molecularDynamics() throws -> VivoWorkflowRecipe {
        let hydrogen = VivoElement.from(symbol: "H")!
        let structure = VivoMolecularStructure(identifier: "harmonic-md-example", atoms: [
            .init(index: 0, name: "H1", element: hydrogen), .init(index: 1, name: "H2", element: hydrogen)
        ], bonds: [.init(atomA: 0, atomB: 1)], conformers: [.init(positionsNM: [.zero, .init(0.11, 0, 0)])])
        let system = VivoClassicalSystem(identifier: "harmonic-md-example", structureFingerprint: try VivoStructureCodec.fingerprint(structure),
            particles: (0..<2).map { .init(index: UInt32($0), atomIndex: UInt32($0), typeIdentifier: "H", massDa: 1, chargeE: 0, sigmaNM: 0, epsilonKJPerMol: 0) },
            bonds: [.init(a: 0, b: 1, lengthNM: 0.1, forceConstant: 100)])
        let state = try VivoClassicalInitialState(systemFingerprint: system.fingerprint(), positionsNM: structure.conformers[0].positionsNM)
        let dynamics = VivoMDConfiguration(timeStepPS: 0.0002, cutoffNM: 0.8, electrostatics: .cutoff,
            ensemble: .nve, thermostat: .none, targetTemperatureK: nil, frictionPerPS: nil, neighborListEnabled: false)
        let resources = VivoChemistryResourceContract(numericalBackend: "metal-fp32")
        return try .init(identifier: "md-segments", artifacts: [
            .init(identifier: "structure", source: .json(kind: "vivo.molecular-structure-document", payload: VivoPlatformOperations.json(VivoMolecularStructureDocument(structure: structure)))),
            .init(identifier: "system", source: .json(kind: "vivo.classical-system", payload: VivoPlatformOperations.json(system))),
            .init(identifier: "initial", source: .json(kind: "vivo.classical-initial-state", payload: VivoPlatformOperations.json(state)))
        ], nodes: [
            .init(identifier: "segment-1", operation: "vivo.platform.md-start", inputs: ["system": .artifact(identifier: "system"), "state": .artifact(identifier: "initial")],
                configuration: VivoPlatformOperations.json(VivoWorkflowMDConfiguration(dynamics: dynamics, steps: 8)), resources: resources),
            .init(identifier: "segment-2", operation: "vivo.platform.md-continue", inputs: ["system": .artifact(identifier: "system"), "checkpoint": .output(node: "segment-1", port: "checkpoint")],
                configuration: VivoPlatformOperations.json(VivoWorkflowMDConfiguration(dynamics: dynamics, steps: 12)), resources: resources)
        ], outputs: [.init(name: "checkpoint", node: "segment-2", port: "checkpoint"), .init(name: "observables", node: "segment-2", port: "result")])
    }
    /// A finite-system trajectory snapshot feeds the standard molecular branch.
    /// Periodic snapshots deliberately require a separately prepared QM/MM region.
    public static func mdElectronicAnalysis() throws -> VivoWorkflowRecipe {
        var recipe = try molecularDynamics()
        var chemistry = try molecularAnalysis()
        recipe.identifier = "md-electronic-analysis"
        recipe.artifacts.append(contentsOf: chemistry.artifacts.filter { $0.identifier == "basis" })
        recipe.nodes.append(.init(identifier: "snapshot", operation: "vivo.platform.md-snapshot", inputs: [
            "structure": .artifact(identifier: "structure"), "system": .artifact(identifier: "system"),
            "checkpoint": .output(node: "segment-2", port: "checkpoint")]))
        chemistry.nodes[0].inputs["structure"] = .output(node: "snapshot", port: "structure")
        recipe.nodes += chemistry.nodes
        recipe.outputs += chemistry.outputs
        recipe.outputs.append(.init(name: "snapshot-mapping", node: "snapshot", port: "mapping"))
        return recipe
    }
}
