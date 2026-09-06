import Foundation
import CryptoKit

public enum VivoReproductionTarget:String,Codable,Sendable { case acrylamideMethanethiolate, btkSnapshot }
public enum VivoReproductionAssetRole:String,Codable,Sendable,CaseIterable {
    case precomplex,transitionState,product,quantumProtocol,referenceEnergies,preparedSnapshot,topology,environment
}
public enum VivoReproductionOrigin:String,Codable,Sendable { case authorSupplied, publicAuthorArchive, independentReconstruction }
public struct VivoReproductionAsset:Codable,Sendable,Equatable {
    public let role:VivoReproductionAssetRole
    public let origin:VivoReproductionOrigin
    public let sourceIdentifier:String
    /// Exact source bytes are verified before semantic JSON decoding. This is
    /// content integrity, NOT authentication of the asserted author/origin.
    public let sha256:String
    public let payload:Data
    public init(role:VivoReproductionAssetRole,origin:VivoReproductionOrigin,sourceIdentifier:String,sha256:String,payload:Data) {
        self.role=role;self.origin=origin;self.sourceIdentifier=sourceIdentifier;self.sha256=sha256;self.payload=payload
    }
}
public struct VivoReproductionGeometry:Codable,Sendable,Equatable {
    public let atomIdentifiers:[String]
    public let system:VivoElectronicSystem
}
public struct VivoReproductionProtocol:Codable,Sendable,Equatable {
    public let publicationIdentifier:String
    public let resultIdentifier:String
    public let basis:VivoGaussianBasis
    public let methodIdentifier:String
    public let methodConfiguration:VivoJSONValue
    public let solventDefinition:VivoJSONValue
    public let energyReference:String
    public let energyUnit:String
}
public struct VivoReproductionReferenceEnergies:Codable,Sendable,Equatable {
    public let publicationIdentifier:String
    public let resultIdentifier:String
    public let methodIdentifier:String
    public let energyReference:String
    public let energyUnit:String
    public let precomplexHartree:Double
    public let transitionStateHartree:Double
    public let productHartree:Double
}
public struct VivoReproductionTopology:Codable,Sendable,Equatable {
    public let atomIdentifiers:[String]
    public let bonds:[[Int]]
    public let forceFieldConfiguration:VivoJSONValue
}
public struct VivoReproductionEnvironment:Codable,Sendable,Equatable {
    public let quantumAtomIdentifiers:[String]
    public let classicalAtomIdentifiers:[String]
    public let pointCharges:[VivoQMPointCharge]
    public let boundaryTreatment:String
}
public struct VivoReproductionPackage:Codable,Sendable,Equatable {
    public static let schema="numivivo.org/reproduction-input-package/v1"
    public let schema:String
    public let target:VivoReproductionTarget
    public let publicationIdentifier:String
    public let resultIdentifier:String
    public let exactAuthorModel:Bool
    public let assets:[VivoReproductionAsset]
    public let budget:VivoChemistryBudget
    public init(target:VivoReproductionTarget,publicationIdentifier:String="arxiv:2604.10487v1",
                resultIdentifier:String,exactAuthorModel:Bool=true,assets:[VivoReproductionAsset]=[],budget:VivoChemistryBudget = .init()) {
        schema=Self.schema;self.target=target;self.publicationIdentifier=publicationIdentifier
        self.resultIdentifier=resultIdentifier;self.exactAuthorModel=exactAuthorModel;self.assets=assets;self.budget=budget
    }
}
public struct VivoReproductionReadiness:Codable,Sendable,Equatable {
    public let target:VivoReproductionTarget
    public let missingRoles:[VivoReproductionAssetRole]
    public let blockers:[String]
    public let verifiedContentDigests:[String:String]
    public let suppliedInputPackageConsistent:Bool
    public let reproductionExecuted:Bool
    public let sourceAuthenticityVerified:Bool
    public let meaning:String
}
public enum VivoReproductionPreflight {
    public static let meaning="content-integrity and declared molecular-input consistency only; author origin is an attestation, not authenticated; no scientific calculation, paper agreement or kinetic-rate qualification"
    public static func inspect(_ package:VivoReproductionPackage) throws -> VivoReproductionReadiness {
        try package.budget.validate()
        guard package.schema==VivoReproductionPackage.schema,!package.publicationIdentifier.isEmpty,
              !package.resultIdentifier.isEmpty,package.publicationIdentifier.utf8.count<=4096,
              package.resultIdentifier.utf8.count<=4096,Set(package.assets.map(\.role)).count==package.assets.count else {
            throw VivoChemistryError.invalid("reproduction package schema, publication/result or duplicate roles")
        }
        let required:[VivoReproductionAssetRole]=package.target == .btkSnapshot ? VivoReproductionAssetRole.allCases :
            [.precomplex,.transitionState,.product,.quantumProtocol,.referenceEnergies]
        var bytes=0,items:[VivoReproductionAssetRole:VivoReproductionAsset]=[:],digests:[String:String]=[:],blockers:[String]=[]
        for asset in package.assets {
            guard !asset.sourceIdentifier.isEmpty,asset.sourceIdentifier.utf8.count<=16384,!asset.payload.isEmpty,
                  asset.payload.count<=package.budget.maximumBytes-bytes,
                  asset.sha256.count==64,asset.sha256.allSatisfy({"0123456789abcdef".contains($0)}) else {
                throw VivoChemistryError.invalid("reproduction asset origin, data size or digest format")
            }
            bytes+=asset.payload.count
            let actual=SHA256.hash(data:asset.payload).map { String(format:"%02x",$0) }.joined()
            guard actual==asset.sha256 else { throw VivoChemistryError.invalid("reproduction asset digest mismatch: \(asset.role.rawValue)") }
            if package.exactAuthorModel && asset.origin == .independentReconstruction {
                blockers.append("independent input cannot establish the author model: \(asset.role.rawValue)")
            }
            items[asset.role]=asset;digests[asset.role.rawValue]=actual
        }
        func decode<T:Decodable>(_ type:T.Type,_ role:VivoReproductionAssetRole) throws -> T? {
            guard let asset=items[role] else { return nil };return try JSONDecoder().decode(type,from:asset.payload)
        }
        var geometries:[VivoReproductionAssetRole:VivoReproductionGeometry]=[:]
        for role in [VivoReproductionAssetRole.precomplex,.transitionState,.product,.preparedSnapshot] {
            guard let geometry=try decode(VivoReproductionGeometry.self,role) else { continue }
            try geometry.system.validate()
            guard geometry.atomIdentifiers.count==geometry.system.nuclei.count,
                  Set(geometry.atomIdentifiers).count==geometry.atomIdentifiers.count,
                  geometry.atomIdentifiers.allSatisfy({!$0.isEmpty && $0.utf8.count<=1024}) else {
                throw VivoChemistryError.invalid("reproduction mapped geometry identifiers")
            }
            geometries[role]=geometry
        }
        if let first=geometries[.precomplex] {
            for role in [VivoReproductionAssetRole.transitionState,.product] {
                guard let other=geometries[role] else { continue }
                guard first.atomIdentifiers==other.atomIdentifiers,
                      first.system.nuclei.map(\.atomicNumber)==other.system.nuclei.map(\.atomicNumber),
                      first.system.nuclei.map(\.structureAtomIndex)==other.system.nuclei.map(\.structureAtomIndex),
                      first.system.alphaElectrons==other.system.alphaElectrons,first.system.betaElectrons==other.system.betaElectrons else {
                    throw VivoChemistryError.invalid("reproduction geometries change atom order, mapping or electron sector")
                }
            }
        }
        let protocolSpec=try decode(VivoReproductionProtocol.self,.quantumProtocol)
        let energies=try decode(VivoReproductionReferenceEnergies.self,.referenceEnergies)
        if let q=protocolSpec {
            guard q.publicationIdentifier==package.publicationIdentifier,q.resultIdentifier==package.resultIdentifier,
                  q.energyUnit=="Hartree",!q.energyReference.isEmpty,!q.methodIdentifier.isEmpty,
                  q.methodConfiguration != .null,q.solventDefinition != .null else {
                throw VivoChemistryError.invalid("reproduction quantum protocol identity, units or omitted method/solvent definition")
            }
            if let g=geometries[.precomplex] { try q.basis.validate(nucleusCount:g.system.nuclei.count) }
        }
        if let e=energies {
            guard e.publicationIdentifier==package.publicationIdentifier,e.resultIdentifier==package.resultIdentifier,
                  e.energyUnit=="Hartree",[e.precomplexHartree,e.transitionStateHartree,e.productHartree].allSatisfy(\.isFinite) else {
                throw VivoChemistryError.invalid("reference energy result identity or units")
            }
            if let q=protocolSpec {
                guard q.methodIdentifier==e.methodIdentifier,q.energyReference==e.energyReference else {
                    throw VivoChemistryError.invalid("paper reference energies and protocol use different method or energy reference")
                }
            }
        }
        if let topology=try decode(VivoReproductionTopology.self,.topology) {
            guard Set(topology.atomIdentifiers).count==topology.atomIdentifiers.count,!topology.atomIdentifiers.isEmpty,
                  topology.forceFieldConfiguration != .null,
                  topology.bonds.allSatisfy({$0.count==2 && $0[0] != $0[1] && $0.allSatisfy({$0>=0 && $0<topology.atomIdentifiers.count})}) else {
                throw VivoChemistryError.invalid("reproduction topology indices or force-field definition")
            }
            if let snapshot=geometries[.preparedSnapshot],snapshot.atomIdentifiers != topology.atomIdentifiers {
                throw VivoChemistryError.invalid("prepared snapshot and topology atom mapping differ")
            }
        }
        if let env=try decode(VivoReproductionEnvironment.self,.environment) {
            let all=env.quantumAtomIdentifiers+env.classicalAtomIdentifiers
            guard !env.boundaryTreatment.isEmpty,Set(all).count==all.count,env.pointCharges.count==env.classicalAtomIdentifiers.count,
                  env.pointCharges.allSatisfy({$0.chargeE.isFinite && vivoQMFinite($0.positionBohr)}) else {
                throw VivoChemistryError.invalid("reproduction QM/MM partition or point-charge definition")
            }
            if let qm=geometries[.precomplex],qm.atomIdentifiers != env.quantumAtomIdentifiers {
                throw VivoChemistryError.invalid("QM geometry and environment partition mapping differ")
            }
            if let snapshot=geometries[.preparedSnapshot],Set(all) != Set(snapshot.atomIdentifiers) {
                throw VivoChemistryError.invalid("environment partition does not cover the prepared snapshot")
            }
        }
        let missing=required.filter { items[$0]==nil }
        return .init(target:package.target,missingRoles:missing,blockers:blockers,verifiedContentDigests:digests,
            suppliedInputPackageConsistent:missing.isEmpty && blockers.isEmpty,reproductionExecuted:false,
            sourceAuthenticityVerified:false,meaning:meaning)
    }
}
