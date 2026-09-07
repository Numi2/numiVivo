import Foundation

/// Molecular preparation reuses the existing AO/HF task identities. Its output
/// is a proposal request, not a completed refinement or chemical qualification.
public enum VivoMolecularSpacePreparationWorkflow {
    public static let inputKind = "vivo.molecular-space-preparation-request"
    public static let outputKind = "vivo.molecular-space-preparation-result"
    private static func decode<T:Decodable>(_ type: T.Type, _ name: String, _ inputs: [String:Data]) throws -> T {
        guard let data = inputs[name] else { throw VivoChemistryError.invalid("missing molecular preparation slot \(name)") }
        return try VivoCanonicalJSON.decode(type,from: data)
    }
    private static func calculate(_ config: VivoJSONValue, _ inputs: [String:Data], _ budget: VivoChemistryBudget,
                                  primitives: Bool) throws -> [String:Data] {
        guard config == .object([:]) else { throw VivoChemistryError.invalid("molecular preparation configuration") }
        let request = try decode(VivoMolecularSpacePreparationRequest.self,"request",inputs)
        try request.validate()
        let slots = ["request"] + (primitives ? request.snapshots.indices.flatMap { ["ao-\($0)","scf-\($0)"] } : [])
        guard Set(inputs.keys) == Set(slots), request.budget == budget else {
            throw VivoChemistryError.invalid("molecular preparation source slots or budget differ from request")
        }
        let result: VivoMolecularSpacePreparationResult
        if primitives {
            result = try VivoMolecularSpacePreparation.prepare(request,
                integrals: request.snapshots.indices.map { try decode(VivoAOIntegrals.self,"ao-\($0)",inputs) },
                references: request.snapshots.indices.map { try decode(VivoHartreeFockResult.self,"scf-\($0)",inputs) })
        } else { result = try VivoMolecularSpacePreparation.prepare(request) }
        return ["preparation": try VivoCanonicalJSON.encode(result),"request": try VivoCanonicalJSON.encode(result.refinementRequest)]
    }
    public static func operation(implementationFingerprint id: VivoFingerprint, suppliedPrimitives: Bool = true) -> VivoChemistryOperation {
        .init(identifier: suppliedPrimitives ? "vivo.native.molecular-space-from-primitives" : "vivo.native.molecular-space-preparation",
            version: "1",implementationFingerprint: id,
            outputs: [.init(name: "preparation",kind: outputKind),.init(name: "request",kind: "vivo.property-directed-space-request")],
            execute: { cfg,inputs,budget in try calculate(cfg,inputs,budget,primitives: suppliedPrimitives) },
            validateOutputs: { cfg,inputs,outputs,budget in
                guard outputs == (try calculate(cfg,inputs,budget,primitives: suppliedPrimitives)) else {
                    throw VivoChemistryError.invalid("molecular preparation artifact differs from native source reconstruction")
                }
            })
    }
    public static func plan(_ request: VivoMolecularSpacePreparationRequest, requestArtifact: VivoFingerprint,
                            systemArtifacts: [VivoFingerprint], basisArtifact: VivoFingerprint,
                            implementationFingerprint id: VivoFingerprint) throws -> VivoElectronicWorkflowPlan {
        try request.validate()
        guard systemArtifacts.count == request.snapshots.count else { throw VivoChemistryError.invalid("molecular source artifact count") }
        let resources = VivoChemistryResourceContract(budget: request.budget,
            maximumInputBytes: request.budget.maximumBytes,maximumOutputBytes: request.budget.maximumBytes)
        let cfg = try VivoCanonicalJSON.decode(VivoJSONValue.self,from: VivoCanonicalJSON.encode(request.reference))
        var nodes: [VivoChemistryDAGNode] = [], inputs: [VivoChemistryDAGInput] = [
            .artifact(name: "request",fingerprint: requestArtifact,kind: inputKind)]
        for i in request.snapshots.indices {
            let ao = "space-\(i)-integrals", hf = "space-\(i)-reference"
            nodes.append(.init(identifier: ao,operation: VivoElectronicWorkflowOperations.integrals(implementationFingerprint: id),
                inputs: [.artifact(name: "system",fingerprint: systemArtifacts[i],kind: "vivo.electronic-system"),
                         .artifact(name: "basis",fingerprint: basisArtifact,kind: "vivo.gaussian-basis")],
                configuration: .object([:]),resources: resources))
            nodes.append(.init(identifier: hf,operation: VivoElectronicWorkflowOperations.hartreeFock(implementationFingerprint: id),
                inputs: [.output(name: "integrals",node: ao,output: "integrals",kind: "vivo.ao-integrals")],
                configuration: cfg,resources: resources))
            inputs.append(.output(name: "ao-\(i)",node: ao,output: "integrals",kind: "vivo.ao-integrals"))
            inputs.append(.output(name: "scf-\(i)",node: hf,output: "scf",kind: "vivo.hartree-fock"))
        }
        nodes.append(.init(identifier: "prepare-space",operation: operation(implementationFingerprint: id),
            inputs: inputs,configuration: .object([:]),resources: resources))
        return .init(nodes: nodes,resultNode: "prepare-space",resultOutput: "request",resultKind: "vivo.property-directed-space-request")
    }
}
