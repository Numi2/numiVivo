import Foundation

public struct VivoSolvatedECCPathRequest: Codable, Sendable, Equatable {
    public static let schema="numivivo.org/reference-polarized-ecc-path/v1"
    public var schema:String
    public var path:VivoMolecularECCPathRequest
    public var solvent:VivoSmoothCPCMConfiguration
    public init(path:VivoMolecularECCPathRequest,solvent:VivoSmoothCPCMConfiguration = .init()) {
        schema=Self.schema;self.path=path;self.solvent=solvent
    }
    public func validate() throws {
        guard schema==Self.schema else {throw VivoChemistryError.invalid("reference-polarized path schema")}
        try path.validate()
    }
}
public struct VivoSolvatedECCPathResult: Codable, Sendable, Equatable {
    public static let schema="numivivo.org/reference-polarized-ecc-path-result/v1"
    public let schema:String
    public let request:VivoSolvatedECCPathRequest
    public let path:VivoECCPathResult
    public let referenceSolventFields:[VivoSmoothCPCMResult]
    public let frozenFieldConstantsHartree:[Double]
    public let convention:String
}
public enum VivoSolvatedECCPath {
    public static let convention="smooth C-PCM equilibrated to RHF at each geometry; correlated ECC solved in the frozen reference field; nuclear/surface-self contribution counted once; not correlated self-consistent PCM or an activation Gibbs energy"
    private static func prepare(_ request:VivoSolvatedECCPathRequest) throws
        -> (points:[VivoECCPathPoint],fields:[VivoSmoothCPCMResult],constants:[Double]) {
        try request.validate()
        let path=request.path
        var prepared:[VivoReferencePolarizedHamiltonian]=[],points:[VivoECCPathPoint]=[]
        for (i,snapshot) in path.snapshots.enumerated() {
            let ao=try VivoGaussianIntegralEngine.compute(system:snapshot.system,basis:path.basis,budget:path.budget)
            let current=try VivoReferencePolarization.prepare(system:snapshot.system,ao:ao,scf:path.reference,solvent:request.solvent,budget:path.budget)
            var overlap:VivoQMMatrix?
            if i>0 {
                let cross=try VivoGaussianIntegralEngine.crossOverlap(leftSystem:path.snapshots[i-1].system,leftBasis:path.basis,
                    rightSystem:snapshot.system,rightBasis:path.basis,budget:path.budget)
                overlap=try prepared[i-1].coefficients.transposed.multiplied(by:cross).multiplied(by:current.coefficients)
            }
            points.append(.init(identifier:snapshot.identifier,hamiltonian:current.hamiltonian,overlapWithPrevious:overlap))
            prepared.append(current)
        }
        try path.configuration.validate(points:points,budget:path.budget)
        return (points,prepared.map{$0.solvent!},prepared.map(\.frozenFieldConstantHartree))
    }
    public static func solve(_ request:VivoSolvatedECCPathRequest) throws -> VivoSolvatedECCPathResult {
        let prepared=try prepare(request),r=request.path
        let result=try VivoECCReactionPath.solve(points:prepared.points,configuration:r.configuration,initialSharedRotation:r.initialSharedRotation,budget:r.budget)
        guard result.converged else {throw VivoChemistryError.convergence("reference-polarized ECC path: \(result.termination)")}
        return .init(schema:VivoSolvatedECCPathResult.schema,request:request,path:result,referenceSolventFields:prepared.fields,
            frozenFieldConstantsHartree:prepared.constants,convention:convention)
    }
    public static func validate(_ result:VivoSolvatedECCPathResult,request:VivoSolvatedECCPathRequest) throws {
        guard result.schema==VivoSolvatedECCPathResult.schema,result.request==request,result.convention==convention else {
            throw VivoChemistryError.invalid("reference-polarized path output method binding")
        }
        let prepared=try prepare(request),r=request.path
        guard result.referenceSolventFields.count==prepared.fields.count,result.frozenFieldConstantsHartree.count==prepared.constants.count else {
            throw VivoChemistryError.invalid("reference-polarized path field count")
        }
        for i in prepared.fields.indices {
            guard result.frozenFieldConstantsHartree[i].isFinite,abs(prepared.constants[i]-result.frozenFieldConstantsHartree[i])<1e-9,
                  result.referenceSolventFields[i].configuration==prepared.fields[i].configuration,
                  try result.referenceSolventFields[i].reactionPotentialMatrix.adding(prepared.fields[i].reactionPotentialMatrix,scale:-1).frobeniusNorm<1e-8,
                  abs(result.referenceSolventFields[i].polarizationEnergyHartree-prepared.fields[i].polarizationEnergyHartree)<1e-9 else {
                throw VivoChemistryError.invalid("reference-polarized path frozen-field energy reconstruction")
            }
        }
        try VivoECCReactionPath.validate(result.path,points:prepared.points,configuration:r.configuration,budget:r.budget)
    }
    /// Freeze a verified shared orbital frame as a local nuclear surface. This
    /// does not differentiate a changing global QIO optimum independently at
    /// each geometry. Nuclear finite differences transport this one anchor.
    public static func nuclearModel(_ result:VivoSolvatedECCPathResult,point:Int) throws -> VivoNuclearElectronicModel {
        try validate(result,request:result.request)
        let path=result.request.path
        guard path.snapshots.indices.contains(point) else {throw VivoChemistryError.invalid("nuclear anchor point")}
        var embedding=path.configuration.embedding;embedding.localityGroups=[]
        let rotation=try result.path.transportRotations[point].multiplied(by:result.path.sharedRotation)
        let frame=VivoAnchoredECCFrame(anchorSystem:path.snapshots[point].system,anchorRotation:rotation,
            transportGroups:path.configuration.transportGroups,expectedBathOrbitals:path.configuration.expectedBathOrbitals,
            embedding:embedding,minimumTransportSingularValue:path.configuration.minimumTransportSingularValue)
        return .init(system:path.snapshots[point].system,basis:path.basis,solver:.anchoredECC,scf:path.reference,
            solvent:result.request.solvent,eccFrame:frame,budget:path.budget)
    }
}
