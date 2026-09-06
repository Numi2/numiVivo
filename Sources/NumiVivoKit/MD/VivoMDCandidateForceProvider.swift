import Foundation

public struct VivoMDCandidateGeometry: Codable, Sendable, Equatable {
    public let particlePositionsNM: [VivoVector3D]
    public let periodicCell: VivoPeriodicCell?
    public init(particlePositionsNM: [VivoVector3D],periodicCell: VivoPeriodicCell?) throws {
        guard !particlePositionsNM.isEmpty,particlePositionsNM.allSatisfy(\.isFinite),periodicCell?.isValid != false else {
            throw VivoChemistryError.invalid("candidate force geometry")
        }
        self.particlePositionsNM = particlePositionsNM;self.periodicCell = periodicCell
    }
    public func fingerprint() throws -> VivoFingerprint { try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(self)) }
}

public struct VivoMDCandidateForceEvaluation: Codable, Sendable, Equatable {
    public let providerFingerprint: VivoFingerprint
    public let geometryFingerprint: VivoFingerprint
    /// Kept in FP64 host accounting, never accumulated into a FP32 per-particle
    /// energy slot: large electronic offsets must not destroy small energy changes.
    public let additionalEnergyKJPerMol: Double
    public let physicalParticleForcesKJPerMolNM: [VivoVector3D]
    public let derivativeMethod: String
    public let convergenceResidual: Double
    public let requiredResidual: Double
    /// Additional contribution only, not the complete classical+electronic stress.
    public let additionalAffineStrainDerivativeKJPerMol: VivoQMMatrix?
    public init(providerFingerprint: VivoFingerprint,geometry: VivoMDCandidateGeometry,
                additionalEnergyKJPerMol: Double,physicalParticleForcesKJPerMolNM: [VivoVector3D],
                derivativeMethod: String,convergenceResidual: Double,requiredResidual: Double,
                additionalAffineStrainDerivativeKJPerMol: VivoQMMatrix? = nil) throws {
        self.providerFingerprint = providerFingerprint;geometryFingerprint = try geometry.fingerprint()
        self.additionalEnergyKJPerMol = additionalEnergyKJPerMol;self.physicalParticleForcesKJPerMolNM = physicalParticleForcesKJPerMolNM
        self.derivativeMethod = derivativeMethod;self.convergenceResidual = convergenceResidual;self.requiredResidual = requiredResidual
        self.additionalAffineStrainDerivativeKJPerMol = additionalAffineStrainDerivativeKJPerMol
    }
    public func validate(geometry: VivoMDCandidateGeometry,provider: VivoMDCandidateForceProvider,
                         system: VivoClassicalSystem) throws {
        guard providerFingerprint == provider.fingerprint,geometryFingerprint == (try geometry.fingerprint()),
              additionalEnergyKJPerMol.isFinite,!derivativeMethod.isEmpty,
              physicalParticleForcesKJPerMolNM.count == system.particles.count,
              physicalParticleForcesKJPerMolNM.allSatisfy(\.isFinite),
              convergenceResidual.isFinite,convergenceResidual >= 0,requiredResidual.isFinite,requiredResidual > 0,
              convergenceResidual <= min(requiredResidual,provider.maximumAcceptedResidual) else {
            throw VivoChemistryError.convergence("candidate force identity, completeness or convergence")
        }
        for particle in system.particles where particle.role == .virtualSite {
            guard physicalParticleForcesKJPerMolNM[Int(particle.index)] == .zero else {
                throw VivoChemistryError.invalid("candidate force provider must return already-mapped physical forces with zero virtual slots")
            }
        }
        if let strain = additionalAffineStrainDerivativeKJPerMol {
            guard strain.rows == 3,strain.columns == 3,strain.values.allSatisfy(\.isFinite) else {
                throw VivoChemistryError.invalid("candidate affine-strain derivative")
            }
        }
    }
}

/// A stateless Born–Oppenheimer force contribution to an explicitly retained
/// classical system. Providers must be deterministic functions of geometry and
/// their fingerprinted configuration. No clock-dependent force, mutable accepted
/// SCF history or implicit fallback is permitted by this contract.
public struct VivoMDCandidateForceProvider: Sendable {
    public let fingerprint: VivoFingerprint
    public let retainedSystemFingerprint: VivoFingerprint
    public let boundary: VivoQMMMEmbeddingBoundary
    public let supportsCellMoves: Bool
    public let polarizationModelFingerprint: VivoFingerprint?
    public let maximumAcceptedResidual: Double
    /// Original molecular connectivity for molecular-center NPT proposals. Removing
    /// QM bonded energies must not split the barostat's molecular components.
    public let molecularConnectivitySystem: VivoClassicalSystem
    public let evaluate: @Sendable (VivoMDCandidateGeometry) async throws -> VivoMDCandidateForceEvaluation
    public init(fingerprint: VivoFingerprint,retainedSystemFingerprint: VivoFingerprint,boundary: VivoQMMMEmbeddingBoundary,
                supportsCellMoves: Bool,maximumAcceptedResidual: Double,molecularConnectivitySystem: VivoClassicalSystem,
                polarizationModelFingerprint: VivoFingerprint? = nil,
                evaluate: @escaping @Sendable (VivoMDCandidateGeometry) async throws -> VivoMDCandidateForceEvaluation) throws {
        try VivoClassicalSystemValidator.validate(molecularConnectivitySystem)
        guard maximumAcceptedResidual.isFinite,maximumAcceptedResidual > 0 else { throw VivoChemistryError.invalid("candidate provider convergence tolerance") }
        self.fingerprint = fingerprint;self.retainedSystemFingerprint = retainedSystemFingerprint;self.boundary = boundary
        self.supportsCellMoves = supportsCellMoves;self.maximumAcceptedResidual = maximumAcceptedResidual
        self.polarizationModelFingerprint=polarizationModelFingerprint
        self.molecularConnectivitySystem = molecularConnectivitySystem;self.evaluate = evaluate
    }
    public func validate(system: VivoClassicalSystem,configuration: VivoMDConfiguration,cell: VivoPeriodicCell?) throws {
        guard retainedSystemFingerprint == (try system.fingerprint()),
              molecularConnectivitySystem.particles.count == system.particles.count,
              molecularConnectivitySystem.constraints == system.constraints,
              molecularConnectivitySystem.linearVirtualSites == system.linearVirtualSites,
              molecularConnectivitySystem.virtualSiteDefinitions == system.virtualSiteDefinitions,
              boundary == .finiteCluster ? cell == nil : cell != nil else {
            throw VivoChemistryError.invalid("candidate provider retained system, molecular manifold or boundary mismatch")
        }
        let polarizationID = try system.polarization.map { try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode($0)) }
        guard polarizationModelFingerprint == polarizationID else {
            throw VivoChemistryError.invalid("force provider does not own the canonical induced-dipole model")
        }
        for (source,retained) in zip(molecularConnectivitySystem.particles,system.particles) {
            guard source.index == retained.index,source.atomIndex == retained.atomIndex,source.role == retained.role,
                  source.massDa == retained.massDa else { throw VivoChemistryError.invalid("candidate provider changes particle or mass ownership") }
        }
        let hasLJ = system.particles.contains(where: { $0.epsilonKJPerMol > 0 })
            || (system.nonbondedTypePairs ?? []).contains(where: { $0.c6KJNM6PerMol > 0 || $0.c12KJNM12PerMol > 0 })
            || system.nonbondedExceptions.contains(where: { ($0.epsilonOverrideKJPerMol ?? 0) > 0 })
        if hasLJ,configuration.lennardJonesSwitchOnNM == nil {
            throw VivoChemistryError.unsupported("BO dynamics with LJ interactions requires an explicit smooth LJ switching radius")
        }
        if boundary == .periodicElectrostatic {
            guard configuration.electrostatics == .pme,configuration.relativeDielectric == 1 else {
                throw VivoChemistryError.unsupported("periodic electronic embedding requires retained PME in vacuum atomic-unit dielectric convention")
            }
        }
        if configuration.ensemble == .npt,boundary == .periodicElectrostatic,configuration.pmeGridDimensions == nil {
            throw VivoChemistryError.unsupported("periodic BO NPT requires fixed pmeGridDimensions; mesh replanning is a different numerical Hamiltonian")
        }
        if configuration.ensemble == .npt && !supportsCellMoves {
            throw VivoChemistryError.unsupported("candidate force provider has no full trial-cell energy support")
        }
    }
    static func executionFingerprint(configuration: VivoMDConfiguration,provider: VivoMDCandidateForceProvider?) throws -> VivoFingerprint {
        guard let provider else { return try configuration.fingerprint() }
        struct Identity: Codable { let configuration: VivoFingerprint; let hamiltonian: VivoFingerprint; let numericalProfile: String }
        return try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(Identity(configuration: configuration.fingerprint(),
            hamiltonian: provider.fingerprint,numericalProfile: "numivivo.org/md-bo-force-provider/v1")))
    }
}

/// Complete potential, not an electronic correction; virtual-site forces already
/// map to physical particles. Used by workflows and adaptive partition evaluation.
public struct VivoMDHamiltonianEvaluation: Codable, Sendable, Equatable {
    public let systemFingerprint: VivoFingerprint
    public let configurationFingerprint: VivoFingerprint
    public let evaluatedGeometry: VivoMDCandidateGeometry
    public let energyKJPerMol: Double
    public let physicalParticleForcesKJPerMolNM: [VivoVector3D]
    public var normalizedForceResidual: Double? = nil
}
