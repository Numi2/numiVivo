import Foundation

public struct VivoRefinedHamiltonianSelection: Codable, Sendable, Equatable {
    public let pointIdentifier: String
    public init(pointIdentifier: String) { self.pointIdentifier = pointIdentifier }
}
public extension VivoCorrelationRefinementOperations {
    /// Materialize one declared anchor, not a Hamiltonian for arbitrary new
    /// geometries. Its full refinement evidence and orbital frame are hashed
    /// into provenance. Existing native CAS/CI/ECC solvers consume the result.
    static func materializeHamiltonian(_ result: VivoPropertyDirectedSpaceResult,
                                       selection: VivoRefinedHamiltonianSelection,
                                       budget: VivoChemistryBudget) throws -> VivoEmbeddedHamiltonian {
        guard result.request.budget == budget, result.sensitivityEstablishedWithinDeclaredPool,
              let index = result.request.points.firstIndex(where: { $0.identifier == selection.pointIdentifier }) else {
            throw VivoChemistryError.invalid("refined Hamiltonian requires established sensitivity, a declared anchor and the same budget")
        }
        try VivoPropertyDirectedSpace.validate(result,request: result.request)
        guard result.transportRotations.count == result.request.points.count else {
            throw VivoChemistryError.invalid("refined Hamiltonian orbital-frame count")
        }
        let parent = result.request.points[index].hamiltonian
        let rotation = result.transportRotations[index]
        let reduced = try parent.rotated(by: rotation,budget: budget).frozenCore(active: result.finalSpace.active,
            doublyOccupiedCore: result.finalSpace.doublyOccupiedCore,budget: budget)
        var provenance = reduced.provenance
        provenance["refinementEvidenceSHA256"] = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(result)).hex
        provenance["refinementRequestSHA256"] = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(result.request)).hex
        provenance["refinementParentHamiltonianSHA256"] = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(parent)).hex
        provenance["refinementOrbitalFrameSHA256"] = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(rotation)).hex
        provenance["refinementPartitionSHA256"] = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(result.finalSpace)).hex
        provenance["refinementPointIdentifier"] = selection.pointIdentifier
        provenance["refinementScope"] = "fixed declared geometry; sensitivity within declared candidate pool only; not a free-energy/rate or arbitrary-geometry production qualification"
        let h = VivoEmbeddedHamiltonian(orbitalIdentifiers: reduced.orbitalIdentifiers,
            alphaElectrons: reduced.alphaElectrons,betaElectrons: reduced.betaElectrons,
            oneElectron: reduced.oneElectron,twoElectron: reduced.twoElectron,
            constantEnergyHartree: reduced.constantEnergyHartree,energyReference: reduced.energyReference,
            provenance: provenance)
        try h.validate(budget: budget)
        return h
    }
    static func refinedHamiltonian(implementationFingerprint id: VivoFingerprint) -> VivoChemistryOperation {
        .init(identifier: "vivo.native.refined-anchor-hamiltonian",version: "1",implementationFingerprint: id,
            outputs: [.init(name: "hamiltonian",kind: "vivo.embedded-hamiltonian")],execute: { cfg,inputs,budget in
                guard Set(inputs.keys) == ["refinement"], let data = inputs["refinement"] else {
                    throw VivoChemistryError.invalid("refined Hamiltonian input slots")
                }
                let result = try VivoCanonicalJSON.decode(VivoPropertyDirectedSpaceResult.self,from: data)
                let selection = try VivoCanonicalJSON.decode(VivoRefinedHamiltonianSelection.self,from: VivoCanonicalJSON.encode(cfg))
                return ["hamiltonian": try VivoCanonicalJSON.encode(materializeHamiltonian(result,selection: selection,budget: budget))]
            },validateOutputs: { cfg,inputs,outputs,budget in
                guard Set(inputs.keys) == ["refinement"], Set(outputs.keys) == ["hamiltonian"],
                      let input = inputs["refinement"], let output = outputs["hamiltonian"] else {
                    throw VivoChemistryError.invalid("refined Hamiltonian output slots")
                }
                let result = try VivoCanonicalJSON.decode(VivoPropertyDirectedSpaceResult.self,from: input)
                let selection = try VivoCanonicalJSON.decode(VivoRefinedHamiltonianSelection.self,from: VivoCanonicalJSON.encode(cfg))
                let h = try VivoCanonicalJSON.decode(VivoEmbeddedHamiltonian.self,from: output)
                guard h == (try materializeHamiltonian(result,selection: selection,budget: budget)) else {
                    throw VivoChemistryError.invalid("refined Hamiltonian differs from its physical source and refinement evidence")
                }
            })
    }
}
