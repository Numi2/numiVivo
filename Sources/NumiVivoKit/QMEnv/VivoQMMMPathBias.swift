import Foundation

public struct VivoQMMMPathBiasConfiguration: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/qmmm-path-bias/v1"
    public var schema: String
    public var progressCenter: Double
    public var progressForceConstantKJPerMol: Double
    public var distanceCenter: Double
    public var distanceForceConstantKJPerMol: Double

    public init(progressCenter: Double,
                progressForceConstantKJPerMol: Double,
                distanceCenter: Double = 0,
                distanceForceConstantKJPerMol: Double = 0) {
        schema = Self.schema
        self.progressCenter = progressCenter
        self.progressForceConstantKJPerMol = progressForceConstantKJPerMol
        self.distanceCenter = distanceCenter
        self.distanceForceConstantKJPerMol = distanceForceConstantKJPerMol
    }

    public func validate() throws {
        guard schema == Self.schema,
              progressCenter.isFinite, progressCenter >= -1, progressCenter <= 2,
              progressForceConstantKJPerMol.isFinite, progressForceConstantKJPerMol >= 0,
              distanceCenter.isFinite,
              distanceForceConstantKJPerMol.isFinite, distanceForceConstantKJPerMol >= 0,
              progressForceConstantKJPerMol + distanceForceConstantKJPerMol > 0 else {
            throw VivoChemistryError.invalid("PathCV bias configuration")
        }
    }

    public func energyKJPerMol(progress: Double, distance: Double) -> Double {
        let ds = progress - progressCenter, dz = distance - distanceCenter
        return 0.5 * progressForceConstantKJPerMol * ds * ds
            + 0.5 * distanceForceConstantKJPerMol * dz * dz
    }
}

public enum VivoQMMMPathBias {
    /// Composes a differentiable path restraint around an existing complete BO
    /// provider. Both PathCV outputs are dimensionless; force constants therefore
    /// have units of kJ/mol. Cell derivatives are deliberately unsupported here.
    public static func provider(base: VivoMDCandidateForceProvider,
                                system: VivoClassicalSystem,
                                path: VivoQMMMPathCoordinate,
                                bias: VivoQMMMPathBiasConfiguration) throws -> VivoMDCandidateForceProvider {
        try bias.validate()
        let resolved = try VivoQMMMResolvedPathCoordinate(source: path, system: system)
        let systemID = try system.fingerprint()
        guard base.retainedSystemFingerprint == systemID else {
            throw VivoChemistryError.invalid("PathCV bias/base retained-system mismatch")
        }
        struct Identity: Codable {
            let schema: String
            let base: VivoFingerprint
            let system: VivoFingerprint
            let path: VivoQMMMPathCoordinate
            let bias: VivoQMMMPathBiasConfiguration
        }
        let id = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(Identity(
            schema: "numivivo.org/qmmm-path-bias-provider/v1",
            base: base.fingerprint, system: systemID, path: path, bias: bias)))
        return try VivoMDCandidateForceProvider(fingerprint: id,
            retainedSystemFingerprint: systemID,
            boundary: base.boundary,
            supportsCellMoves: false,
            maximumAcceptedResidual: base.maximumAcceptedResidual,
            molecularConnectivitySystem: base.molecularConnectivitySystem,
            polarizationModelFingerprint: base.polarizationModelFingerprint) { geometry in
            let raw = try await base.evaluate(geometry)
            try raw.validate(geometry: geometry, provider: base, system: system)
            let cv = try resolved.evaluate(geometry)
            let ds = cv.progress - bias.progressCenter
            let dz = cv.distanceFromPath - bias.distanceCenter
            let slopeS = bias.progressForceConstantKJPerMol * ds
            let slopeZ = bias.distanceForceConstantKJPerMol * dz
            var force = raw.physicalParticleForcesKJPerMolNM
            for (particle, gradient) in cv.progressGradients {
                force[Int(particle)] = force[Int(particle)] - gradient * slopeS
            }
            for (particle, gradient) in cv.distanceGradients {
                force[Int(particle)] = force[Int(particle)] - gradient * slopeZ
            }
            return try VivoMDCandidateForceEvaluation(providerFingerprint: id,
                geometry: geometry,
                additionalEnergyKJPerMol: raw.additionalEnergyKJPerMol + bias.energyKJPerMol(progress: cv.progress, distance: cv.distanceFromPath),
                physicalParticleForcesKJPerMolNM: force,
                derivativeMethod: raw.derivativeMethod + "; analytic multidimensional PathCV restraint",
                convergenceResidual: raw.convergenceResidual,
                requiredResidual: raw.requiredResidual,
                additionalAffineStrainDerivativeKJPerMol: nil)
        }
    }
}
