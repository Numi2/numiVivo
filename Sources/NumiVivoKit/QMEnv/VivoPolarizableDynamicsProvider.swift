import Foundation
@preconcurrency import Metal

public extension VivoMDCandidateForceProvider {
    /// Complete stationary polarization correction to the existing permanent-
    /// charge MD kernels. The nonlinear site Jacobian is applied exactly once.
    static func inducedDipoles(system: VivoClassicalSystem,ewald: VivoPeriodicElectrostaticConfiguration,
                               reciprocalOperator: VivoReciprocalElectrostaticOperator? = nil,
                               budget: VivoChemistryBudget = .init()) throws -> Self {
        try VivoClassicalSystemValidator.validate(system);try budget.validate();try ewald.validate()
        guard let model=system.polarization else { throw VivoChemistryError.invalid("classical system has no induced-dipole model") }
        guard abs(system.particles.reduce(0) { $0+$1.chargeE })<1e-6 else {
            throw VivoChemistryError.unsupported("polarizable dynamics v1 requires an explicitly neutral cell")
        }
        let graph=try system.resolvedVirtualSiteGraph(),systemID=try system.fingerprint()
        let modelID=try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(model))
        struct Identity: Encodable {
            let model: String;let system: VivoFingerprint;let ewald: VivoPeriodicElectrostaticConfiguration
            let mesh: VivoMultipoleMeshConfiguration?;let budget: VivoChemistryBudget
        }
        let id=try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(Identity(model:"mutual-induced-md/v1",
            system:systemID,ewald:ewald,mesh:reciprocalOperator?.configuration,budget:budget)))
        return try .init(fingerprint:id,retainedSystemFingerprint:systemID,boundary:.periodicElectrostatic,
            supportsCellMoves:true,maximumAcceptedResidual:1,molecularConnectivitySystem:system,
            polarizationModelFingerprint:modelID,evaluate:{ geometry in
                try Task.checkCancellation()
                guard let cell=geometry.periodicCell,geometry.particlePositionsNM.count==system.particles.count else {
                    throw VivoChemistryError.invalid("polarizable MD candidate geometry or cell")
                }
                let state=try graph.construct(positionsNM:geometry.particlePositionsNM,periodicCell:cell)
                let charges=system.particles.map { p -> VivoQMPointCharge in
                    let r=state.positionsNM[Int(p.index)]/VivoAtomicUnits.bohrInNM
                    return .init(chargeE:p.chargeE,positionBohr:.init(r.x,r.y,r.z),classicalParticleIndex:p.index)
                }
                let evaluated=try VivoInducedDipoles.evaluate(pointCharges:charges,cell:cell,configuration:model,
                    ewald:ewald,reciprocalOperator:reciprocalOperator,budget:budget)
                let physical=try graph.redistribute(rawForces:evaluated.centerForcesHartreePerBohr,state:state)
                let siteStrain=try graph.affineStrainCorrection(rawForces:evaluated.centerForcesHartreePerBohr,
                    state:state,periodicCell:cell,positionUnitScale:1/VivoAtomicUnits.bohrInNM)
                let strain=try evaluated.affineStrainDerivativeHartree.adding(siteStrain)
                let forceUnit=VivoAtomicUnits.hartreeInKJPerMol/VivoAtomicUnits.bohrInNM
                return try .init(providerFingerprint:id,geometry:geometry,
                    additionalEnergyKJPerMol:evaluated.energyHartree*VivoAtomicUnits.hartreeInKJPerMol,
                    physicalParticleForcesKJPerMolNM:physical.map { $0*forceUnit },
                    derivativeMethod:"variational-mutual-induced-dipoles; analytic-local-frame-and-Thole-pair; normalized-response-residual",
                    convergenceResidual:evaluated.response.maximumResponseResidual/model.residualToleranceHartreePerEBohr,
                    requiredResidual:1,additionalAffineStrainDerivativeKJPerMol:strain.scaled(VivoAtomicUnits.hartreeInKJPerMol))
            })
    }
}
