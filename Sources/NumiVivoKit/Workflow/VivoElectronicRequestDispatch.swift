import Foundation

/// Explicit schema dispatch. A malformed advanced request is never retried as
/// an older method, and the public CLI uses the same planners as library callers.
public enum VivoElectronicRequestDocument: Sendable {
    case original(VivoElectronicWorkflowRequest)
    case advanced(VivoAdvancedElectronicWorkflowRequest)
    public static func decode(_ data:Data) throws -> Self {
        struct Header:Decodable {let schema:String}
        let header=try VivoCanonicalJSON.decode(Header.self,from:data)
        let document:Self
        switch header.schema {
        case VivoElectronicWorkflowRequest.schema:document = .original(try VivoCanonicalJSON.decode(VivoElectronicWorkflowRequest.self,from:data))
        case VivoAdvancedElectronicWorkflowRequest.schema:document = .advanced(try VivoCanonicalJSON.decode(VivoAdvancedElectronicWorkflowRequest.self,from:data))
        default:throw VivoChemistryError.unsupported("electronic request schema \(header.schema)")
        }
        switch document {case .original(let r):try r.validate();case .advanced(let r):try r.validate()}
        return document
    }
    public var system:VivoElectronicSystem {switch self{case .original(let r):return r.system;case .advanced(let r):return r.system}}
    public var basis:VivoGaussianBasis {switch self{case .original(let r):return r.basis;case .advanced(let r):return r.basis}}
    public var budget:VivoChemistryBudget {switch self{case .original(let r):return r.budget;case .advanced(let r):return r.budget}}
    public var artifactKind:String {switch self{case .original:return "vivo.electronic-workflow-request";case .advanced:return "vivo.advanced-electronic-workflow-request"}}
    public func canonicalData() throws -> Data {switch self{case .original(let r):return try VivoCanonicalJSON.encode(r);case .advanced(let r):return try VivoCanonicalJSON.encode(r)}}
    public func plan(systemArtifact:VivoFingerprint,basisArtifact:VivoFingerprint,implementationFingerprint:VivoFingerprint) throws -> VivoElectronicWorkflowPlan {
        switch self {
        case .original(let r):return try VivoElectronicWorkflowPlanner.plan(r,systemArtifact:systemArtifact,basisArtifact:basisArtifact,implementationFingerprint:implementationFingerprint)
        case .advanced(let r):return try VivoAdvancedElectronicWorkflowPlanner.plan(r,systemArtifact:systemArtifact,basisArtifact:basisArtifact,implementationFingerprint:implementationFingerprint)
        }
    }
}
public struct VivoECCDMETSolverRequest: Codable, Sendable, Equatable {
    public static let schema="numivivo.org/ecc-dmet-solver-request/v1"
    public let schema:String
    public let configuration:VivoECCDMETConfiguration
    public init(configuration:VivoECCDMETConfiguration){schema=Self.schema;self.configuration=configuration}
}
public struct VivoNativeSolverDispatch: Sendable {
    public let operation:VivoChemistryOperation
    public let configuration:VivoJSONValue
    public static func decode(_ data:Data,implementationFingerprint id:VivoFingerprint) throws -> Self {
        let json=try VivoCanonicalJSON.decode(VivoJSONValue.self,from:data)
        guard case .object(let object)=json else{throw VivoChemistryError.invalid("solver document must be an object")}
        func encoded<T:Encodable>(_ value:T) throws -> VivoJSONValue {try VivoCanonicalJSON.decode(VivoJSONValue.self,from:VivoCanonicalJSON.encode(value))}
        if object["schema"] != nil {
            let request=try VivoCanonicalJSON.decode(VivoECCDMETSolverRequest.self,from:data)
            guard request.schema==VivoECCDMETSolverRequest.schema,Set(object.keys)==["schema","configuration"] else {
                throw VivoChemistryError.unsupported("unknown ECC solver schema or fields")
            }
            return .init(operation:VivoAdvancedChemistryOperations.eccDMET(implementationFingerprint:id),configuration:try encoded(request.configuration))
        }
        guard object.count==1,let key=object.keys.first else{throw VivoChemistryError.invalid("exactly one explicit solver method is required")}
        if ["mp2","configurationInteraction","ccsd","casci","casscf"].contains(key) {
            return .init(operation:VivoElectronicWorkflowOperations.manyBody(implementationFingerprint:id),
                configuration:try encoded(VivoCanonicalJSON.decode(VivoManyBodySolverRequest.self,from:data)))
        }
        if ["directCI","tensorCCSD","multistateCASSCF"].contains(key) {
            return .init(operation:VivoAdvancedChemistryOperations.manyBody(implementationFingerprint:id),
                configuration:try encoded(VivoCanonicalJSON.decode(VivoAdvancedManyBodyRequest.self,from:data)))
        }
        throw VivoChemistryError.unsupported("unrecognized solver method \(key)")
    }
}
