import Foundation

public struct VivoMolecularECCPathOutput: Codable, Sendable, Equatable {
    public static let schema="numivivo.org/molecular-ecc-path-output/v1"
    public let schema: String
    public let requestFingerprint: VivoFingerprint
    public let snapshotFingerprints: [VivoFingerprint]
    public let coordinates: [Double]
    public let coordinateUnit: String
    public let path: VivoECCPathResult
    /// Explicit absence of stationary-point and thermodynamic qualification.
    public let energyMeaning: String
}
public enum VivoMolecularECCPathWorkflow {
    public static let outputKind="vivo.molecular-ecc-path-output"
    public static let energyMeaning="electronic-profile; no transition-state characterization, thermal correction, standard-state correction or kinetic export"
    private static func decode<T:Decodable>(_ type: T.Type, _ name: String, _ inputs: [String:Data]) throws -> T {
        guard let data=inputs[name] else { throw VivoChemistryError.invalid("missing path artifact slot \(name)") }
        return try VivoCanonicalJSON.decode(type,from:data)
    }
    private static func source(_ inputs: [String:Data], _ budget: VivoChemistryBudget) throws
        -> (VivoMolecularECCPathRequest,[VivoECCPathPoint]) {
        let request=try decode(VivoMolecularECCPathRequest.self,"request",inputs)
        try request.validate()
        let expected=["request"]+request.snapshots.indices.flatMap { ["ao-\($0)","scf-\($0)"] }
        guard Set(inputs.keys)==Set(expected),request.budget==budget else {
            throw VivoChemistryError.invalid("path input slots or execution budget differ from the request")
        }
        let integrals=try request.snapshots.indices.map { try decode(VivoAOIntegrals.self,"ao-\($0)",inputs) }
        let references=try request.snapshots.indices.map { try decode(VivoHartreeFockResult.self,"scf-\($0)",inputs) }
        return (request,try VivoMolecularECCPath.prepare(request,integrals:integrals,references:references))
    }
    private static func fingerprint<T:Encodable>(_ value: T) throws -> VivoFingerprint {
        try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(value))
    }
    public static func operation(implementationFingerprint id: VivoFingerprint) -> VivoChemistryOperation {
        .init(identifier:"vivo.native.molecular-ecc-path",version:"1",implementationFingerprint:id,
            outputs:[.init(name:"profile",kind:outputKind)],execute:{ config,inputs,budget in
                guard config == .object([:]) else { throw VivoChemistryError.invalid("unrecognized path-operation configuration") }
                let (request,points)=try source(inputs,budget)
                let path=try VivoECCReactionPath.solve(points:points,configuration:request.configuration,
                    initialSharedRotation:request.initialSharedRotation,budget:budget)
                guard path.converged else { throw VivoChemistryError.convergence("shared ECC path: \(path.termination)") }
                let output=VivoMolecularECCPathOutput(schema:VivoMolecularECCPathOutput.schema,
                    requestFingerprint:try fingerprint(request),snapshotFingerprints:try request.snapshots.map { try fingerprint($0) },
                    coordinates:request.snapshots.map(\.coordinate),coordinateUnit:request.coordinateUnit,path:path,energyMeaning:energyMeaning)
                return ["profile":try VivoCanonicalJSON.encode(output)]
            },validateOutputs:{ config,inputs,outputs,budget in
                guard config == .object([:]),Set(outputs.keys)==Set(["profile"]) else {
                    throw VivoChemistryError.invalid("path output slots/configuration")
                }
                let (request,points)=try source(inputs,budget)
                let result=try decode(VivoMolecularECCPathOutput.self,"profile",outputs)
                guard result.schema==VivoMolecularECCPathOutput.schema,
                      result.requestFingerprint == (try fingerprint(request)),
                      result.snapshotFingerprints == (try request.snapshots.map { try fingerprint($0) }),
                      result.coordinates==request.snapshots.map(\.coordinate),result.coordinateUnit==request.coordinateUnit,
                      result.energyMeaning==energyMeaning else { throw VivoChemistryError.invalid("path request/geometry/energy-meaning binding") }
                try VivoECCReactionPath.validate(result.path,points:points,configuration:request.configuration,budget:budget)
            })
    }
    /// AO and RHF primitives use the existing operation/version/configuration
    /// identities. Node names are not scientific identities, so isolated-point
    /// and path workflows can reuse exactly the same primitive results.
    public static func plan(_ request: VivoMolecularECCPathRequest, requestArtifact: VivoFingerprint,
                            systemArtifacts: [VivoFingerprint], basisArtifact: VivoFingerprint,
                            implementationFingerprint id: VivoFingerprint) throws -> VivoElectronicWorkflowPlan {
        try request.validate()
        guard systemArtifacts.count==request.snapshots.count else { throw VivoChemistryError.invalid("path system artifact count") }
        let resources=VivoChemistryResourceContract(budget:request.budget,
            maximumInputBytes:request.budget.maximumBytes,maximumOutputBytes:request.budget.maximumBytes)
        let scfConfig=try VivoCanonicalJSON.decode(VivoJSONValue.self,from:VivoCanonicalJSON.encode(request.reference))
        var nodes:[VivoChemistryDAGNode]=[],finalInputs:[VivoChemistryDAGInput]=[
            .artifact(name:"request",fingerprint:requestArtifact,kind:"vivo.molecular-ecc-path-request")]
        for i in request.snapshots.indices {
            let ao="path-\(i)-integrals",ref="path-\(i)-reference"
            nodes.append(.init(identifier:ao,operation:VivoElectronicWorkflowOperations.integrals(implementationFingerprint:id),
                inputs:[.artifact(name:"system",fingerprint:systemArtifacts[i],kind:"vivo.electronic-system"),
                        .artifact(name:"basis",fingerprint:basisArtifact,kind:"vivo.gaussian-basis")],
                configuration:.object([:]),resources:resources))
            nodes.append(.init(identifier:ref,operation:VivoElectronicWorkflowOperations.hartreeFock(implementationFingerprint:id),
                inputs:[.output(name:"integrals",node:ao,output:"integrals",kind:"vivo.ao-integrals")],
                configuration:scfConfig,resources:resources))
            finalInputs.append(.output(name:"ao-\(i)",node:ao,output:"integrals",kind:"vivo.ao-integrals"))
            finalInputs.append(.output(name:"scf-\(i)",node:ref,output:"scf",kind:"vivo.hartree-fock"))
        }
        nodes.append(.init(identifier:"path-profile",operation:operation(implementationFingerprint:id),
            inputs:finalInputs,configuration:.object([:]),resources:resources))
        return .init(nodes:nodes,resultNode:"path-profile",resultOutput:"profile",resultKind:outputKind)
    }
}
