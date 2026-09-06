import Foundation

public enum VivoQMMMFreeEnergyEnvironment: String, Codable, Sendable {
    case explicitSolution
    case proteinEnvironment
}

public struct VivoQMMMFreeEnergyProvenance: Codable, Sendable, Equatable {
    public static let schema="numivivo.org/qmmm-free-energy-provenance/v2"
    public var schema:String
    public var structureFingerprint:VivoFingerprint
    public var systemFingerprint:VivoFingerprint
    public var baseProviderFingerprint:VivoFingerprint
    public var dynamicsFingerprint:VivoFingerprint
    public var chemicalState:String
    public var environment:VivoQMMMFreeEnergyEnvironment
    public var environmentIdentifier:String
    public var methodDescription:String
    public init(structureFingerprint:VivoFingerprint,systemFingerprint:VivoFingerprint,baseProviderFingerprint:VivoFingerprint,
                dynamicsFingerprint:VivoFingerprint,chemicalState:String,environment:VivoQMMMFreeEnergyEnvironment,
                environmentIdentifier:String,methodDescription:String) {
        schema=Self.schema;self.structureFingerprint=structureFingerprint;self.systemFingerprint=systemFingerprint
        self.baseProviderFingerprint=baseProviderFingerprint;self.dynamicsFingerprint=dynamicsFingerprint
        self.chemicalState=chemicalState;self.environment=environment;self.environmentIdentifier=environmentIdentifier
        self.methodDescription=methodDescription
    }
    public func validate() throws {
        guard schema==Self.schema,
              !chemicalState.trimmingCharacters(in:.whitespacesAndNewlines).isEmpty,chemicalState.utf8.count<=1024,
              !environmentIdentifier.trimmingCharacters(in:.whitespacesAndNewlines).isEmpty,environmentIdentifier.utf8.count<=4096,
              !methodDescription.trimmingCharacters(in:.whitespacesAndNewlines).isEmpty,methodDescription.utf8.count<=16384 else {
            throw VivoChemistryError.invalid("QM/MM free-energy provenance identity")
        }
    }
}

public struct VivoQMMMQualifiedActivationFreeEnergy: Codable, Sendable, Equatable {
    public static let schema="numivivo.org/qmmm-qualified-activation-free-energy/v2"
    public let schema:String
    public let analysis:VivoQMMMActivationFreeEnergyResult
    public let provenance:VivoQMMMFreeEnergyProvenance
    public let evidenceFingerprint:VivoFingerprint
    public init(analysis:VivoQMMMActivationFreeEnergyResult,provenance:VivoQMMMFreeEnergyProvenance) throws {
        try provenance.validate();try VivoQMMMFreeEnergy.validate(analysis)
        guard analysis.converged else { throw VivoChemistryError.convergence("only a converged PMF can be qualified") }
        self.schema=Self.schema;self.analysis=analysis;self.provenance=provenance
        struct Evidence:Codable {let analysis:VivoQMMMActivationFreeEnergyResult;let provenance:VivoQMMMFreeEnergyProvenance}
        self.evidenceFingerprint=try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(Evidence(analysis:analysis,provenance:provenance)))
    }
    public func validate() throws {
        guard schema==Self.schema else { throw VivoChemistryError.invalid("qualified activation-free-energy schema") }
        try provenance.validate();try VivoQMMMFreeEnergy.validate(analysis)
        guard analysis.converged else { throw VivoChemistryError.convergence("qualified PMF is not converged") }
        struct Evidence:Codable {let analysis:VivoQMMMActivationFreeEnergyResult;let provenance:VivoQMMMFreeEnergyProvenance}
        let expected=try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(Evidence(analysis:analysis,provenance:provenance)))
        guard evidenceFingerprint==expected else { throw VivoChemistryError.invalid("qualified PMF evidence fingerprint mismatch") }
    }
}

public enum VivoQMMMFreeEnergyQualification {
    /// Builds provenance from the same system/provider/configuration used by the
    /// production runner. A caller cannot relabel an isolated or different
    /// Hamiltonian as a protein environment after sampling.
    public static func qualify(_ result:VivoQMMMActivationFreeEnergyResult,
                               sampling:VivoQMMMFreeEnergyRunRequest,
                               system:VivoClassicalSystem,
                               baseProvider:VivoMDCandidateForceProvider,
                               chemicalState:String,
                               environment:VivoQMMMFreeEnergyEnvironment,
                               environmentIdentifier:String,
                               methodDescription:String) throws -> VivoQMMMQualifiedActivationFreeEnergy {
        let systemID=try system.fingerprint()
        guard result.coordinate==sampling.coordinate,result.temperatureK==sampling.dynamics.targetTemperatureK,
              result.traces.map(\.window)==sampling.windows,
              result.traces.map(\.randomSeed)==sampling.randomSeeds,
              baseProvider.retainedSystemFingerprint==systemID else {
            throw VivoChemistryError.invalid("PMF result differs from sampling/provider identity")
        }
        let dynamicsID=try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(sampling.dynamics))
        let provenance=VivoQMMMFreeEnergyProvenance(structureFingerprint:system.structureFingerprint,systemFingerprint:systemID,
            baseProviderFingerprint:baseProvider.fingerprint,dynamicsFingerprint:dynamicsID,chemicalState:chemicalState,
            environment:environment,environmentIdentifier:environmentIdentifier,methodDescription:methodDescription)
        return try .init(analysis:result,provenance:provenance)
    }
}
