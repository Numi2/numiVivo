import Foundation

public struct VivoConstantPHMetalStateSpecification: Sendable {
    public let identifier: String
    public let system: VivoClassicalSystem
    public let configuration: VivoMDConfiguration
    public let forceProvider: VivoMDCandidateForceProvider?

    public init(identifier: String, system: VivoClassicalSystem,
                configuration: VivoMDConfiguration,
                forceProvider: VivoMDCandidateForceProvider? = nil) {
        self.identifier = identifier; self.system = system
        self.configuration = configuration; self.forceProvider = forceProvider
    }
}

private struct VivoConstantPHManifoldParticle: Codable, Sendable {
    let index: UInt32
    let atomIndex: UInt32?
    let role: VivoParticleRole
    let massDa: Double
}
private struct VivoConstantPHPhysicalManifold: Codable, Sendable {
    let schema: String
    let particles: [VivoConstantPHManifoldParticle]
    let constraints: [VivoDistanceConstraint]
    let linearVirtualSites: [VivoLinearVirtualSite]
    let dependentSites: [VivoDependentSite]
}

/// Production adapter for the generic constant-pH samplers. Every state uses the
/// existing transaction-owning Metal runtime, so potential energies and propagation
/// include the complete retained classical terms plus the supplied QM/MM provider.
public enum VivoConstantPHMetalStateFactory {
    public static let interpretation = "Native Apple-Metal constant-pH state adapter. Each state restores the existing MD runtime from the exact common physical snapshot and its own system/provider execution identity. Complete accepted-state potential energy is read from the authoritative MD evaluator; no second classical or QM/MM force path is introduced."

    public static func make(specifications: [VivoConstantPHMetalStateSpecification],
                            samplingTemperatureK: Double) throws -> [VivoConstantPHExecutableState] {
        guard !specifications.isEmpty, specifications.count <= 4096,
              Set(specifications.map(\.identifier)).count == specifications.count,
              samplingTemperatureK.isFinite, samplingTemperatureK > 0 else {
            throw VivoChemistryError.invalid("constant-pH Metal state specifications")
        }
        guard let firstConfiguration = specifications.first?.configuration else {
            throw VivoChemistryError.invalid("empty constant-pH Metal specification")
        }
        try firstConfiguration.validate()
        guard firstConfiguration.ensemble != .npt else {
            throw VivoChemistryError.unsupported("constant-pH v1 does not combine discrete state moves with barostat cell moves")
        }
        if firstConfiguration.ensemble == .nvt {
            guard firstConfiguration.targetTemperatureK == samplingTemperatureK else {
                throw VivoChemistryError.invalid("constant-pH sampling and NVT temperatures differ")
            }
        }
        let firstManifold = try manifoldFingerprint(specifications[0].system)
        var result: [VivoConstantPHExecutableState] = []
        result.reserveCapacity(specifications.count)
        for specification in specifications {
            guard !specification.identifier.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty,
                  specification.configuration == firstConfiguration else {
                throw VivoChemistryError.invalid("constant-pH Metal states require one MD numerical configuration")
            }
            try VivoClassicalSystemValidator.validate(specification.system)
            guard try manifoldFingerprint(specification.system) == firstManifold else {
                throw VivoChemistryError.invalid("constant-pH Metal states differ in particle masses, identities, constraints or dependent-site manifold")
            }
            let system = specification.system
            let configuration = specification.configuration
            let provider = specification.forceProvider
            let systemID = try system.fingerprint()
            let executionID = try VivoMDCandidateForceProvider.executionFingerprint(configuration: configuration, provider: provider)
            let hamiltonianID = try hamiltonianFingerprint(system: system, configuration: configuration, provider: provider)
            let identifier = specification.identifier
            let executable = VivoConstantPHExecutableState(identifier: identifier,
                physicalManifoldFingerprint: firstManifold,
                hamiltonianFingerprint: hamiltonianID,
                potentialEnergyKJPerMol: { physical in
                    try physical.validate()
                    guard physical.positionsNM.count == system.particles.count else {
                        throw VivoChemistryError.invalid("constant-pH physical snapshot does not match Metal particle count")
                    }
                    let checkpoint = VivoMDCheckpoint(systemFingerprint: systemID,
                        configurationFingerprint: executionID, acceptedStep: physical.stepIndex,
                        timePS: physical.timePS, positionsNM: physical.positionsNM,
                        velocitiesNMPerPS: physical.velocitiesNMPerPS, periodicCell: physical.periodicCell)
                    let runtime = try await VivoMDMetalRuntime.restore(system: system, configuration: configuration,
                        checkpoint: checkpoint, forceProvider: provider)
                    let observable = try await runtime.observables()
                    guard observable.potentialEnergyKJPerMol.isFinite else {
                        throw VivoChemistryError.convergence("nonfinite constant-pH Metal potential energy")
                    }
                    return observable.potentialEnergyKJPerMol
                },
                propagate: { physical, steps in
                    try physical.validate()
                    guard physical.positionsNM.count == system.particles.count, steps > 0 else {
                        throw VivoChemistryError.invalid("constant-pH Metal propagation snapshot or step count")
                    }
                    let checkpoint = VivoMDCheckpoint(systemFingerprint: systemID,
                        configurationFingerprint: executionID, acceptedStep: physical.stepIndex,
                        timePS: physical.timePS, positionsNM: physical.positionsNM,
                        velocitiesNMPerPS: physical.velocitiesNMPerPS, periodicCell: physical.periodicCell)
                    let runtime = try await VivoMDMetalRuntime.restore(system: system, configuration: configuration,
                        checkpoint: checkpoint, forceProvider: provider)
                    var remaining = steps
                    while remaining > 0 {
                        try Task.checkCancellation()
                        let certificate = try await runtime.step()
                        guard certificate.committed else {
                            throw VivoMDRuntimeError.candidateRejected(certificate.statusFlags)
                        }
                        remaining -= 1
                    }
                    let snapshot = try await runtime.snapshot()
                    return try VivoConstantPHPhysicalState(stepIndex: snapshot.stepIndex,
                        timePS: snapshot.timePS, positionsNM: snapshot.positionsNM,
                        velocitiesNMPerPS: snapshot.velocitiesNMPerPS, periodicCell: snapshot.periodicCell)
                })
            result.append(executable)
        }
        return result
    }

    private static func manifoldFingerprint(_ system: VivoClassicalSystem) throws -> VivoFingerprint {
        try VivoClassicalSystemValidator.validate(system)
        let signature = VivoConstantPHPhysicalManifold(
            schema: "numivivo.org/constant-ph-physical-manifold/v1",
            particles: system.particles.map { .init(index: $0.index, atomIndex: $0.atomIndex, role: $0.role, massDa: $0.massDa) },
            constraints: system.constraints,
            linearVirtualSites: system.linearVirtualSites ?? [],
            dependentSites: system.virtualSiteDefinitions ?? [])
        return try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(signature))
    }

    private static func hamiltonianFingerprint(system: VivoClassicalSystem,
                                               configuration: VivoMDConfiguration,
                                               provider: VivoMDCandidateForceProvider?) throws -> VivoFingerprint {
        struct Identity: Codable {
            let schema: String
            let system: VivoFingerprint
            let execution: VivoFingerprint
            let provider: VivoFingerprint?
        }
        let systemFingerprint: VivoFingerprint = try system.fingerprint()
        let executionFingerprint: VivoFingerprint = try VivoMDCandidateForceProvider.executionFingerprint(
            configuration: configuration, provider: provider)
        let providerFingerprint: VivoFingerprint? = provider?.fingerprint
        let identity: Identity = .init(
            schema: "numivivo.org/constant-ph-metal-hamiltonian/v1",
            system: systemFingerprint,
            execution: executionFingerprint,
            provider: providerFingerprint)
        let encoded: Data = try VivoCanonicalJSON.encode(identity)
        return try VivoCanonicalJSON.fingerprint(encoded)
    }
}