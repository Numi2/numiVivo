import Foundation

public struct VivoQMMMDynamicsElectronicConfiguration: Codable, Sendable, Equatable {
    public var basis: VivoGaussianBasis
    public var periodic: VivoPeriodicQMMMConfiguration
    public var cluster: VivoQMMMClusterConfiguration
    public var budget: VivoChemistryBudget
    public init(basis: VivoGaussianBasis,periodic: VivoPeriodicQMMMConfiguration = .init(),
                cluster: VivoQMMMClusterConfiguration = .init(),budget: VivoChemistryBudget = .init()) {
        self.basis = basis;self.periodic = periodic;self.cluster = cluster;self.budget = budget
    }
}

public extension VivoMDCandidateForceProvider {
    /// Reuses the existing periodic/finite compiler, Gaussian integral engine,
    /// shared SCF iterations, link Jacobian and dependent-site graph. The returned
    /// contribution contains NO LJ or retained classical terms: those are evaluated
    /// by VivoMDMetalRuntime using plan.retainedSystem.
    static func hartreeFock(document: VivoMolecularStructureDocument,sourceSystem: VivoClassicalSystem,
                           plan: VivoQMMMHamiltonianPlan,configuration cfg: VivoQMMMDynamicsElectronicConfiguration) throws -> Self {
        try plan.validate(document: document,source: sourceSystem,budget: cfg.budget)
        try cfg.periodic.validate();try cfg.cluster.validate();try cfg.budget.validate()
        guard !plan.configuration.region.qmAtomIndices.isEmpty else { throw VivoChemistryError.invalid("all-MM plans use the ordinary classical runtime") }
        struct Identity: Encodable {
            let implementation: String;let plan: VivoFingerprint;let electronic: VivoQMMMDynamicsElectronicConfiguration
        }
        let id = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(Identity(
            implementation: "numivivo.org/analytic-periodic-hf-md/v1",plan: plan.fingerprint(),electronic: cfg)))
        let retainedID = try plan.retainedSystem.fingerprint()
        let topology = try VivoQMMMTopology(document: document,system: sourceSystem)
        return try .init(fingerprint: id,retainedSystemFingerprint: retainedID,boundary: plan.configuration.boundary,
            supportsCellMoves: true,maximumAcceptedResidual: cfg.periodic.scf.commutatorTolerance,
            molecularConnectivitySystem: sourceSystem,evaluate: { geometry in
                try Task.checkCancellation()
                let cluster = try VivoQMMMClusterBuilder.reconstruct(document: document,system: sourceSystem,
                    particlePositionsNM: geometry.particlePositionsNM,periodicCell: geometry.periodicCell,
                    sourceFrameFingerprint: geometry.fingerprint(),qmAtomIndices: plan.configuration.region.qmAtomIndices,
                    configuration: cfg.cluster)
                var request = plan.configuration.region;request.coordinatesAreUnwrappedFiniteCluster = true
                let region = try VivoQMMMCompiler.prepare(document: document,system: sourceSystem,
                    particlePositionsNM: cluster.particlePositionsNM,request: request)
                var electronic = region.electronicSystem
                electronic.pointCharges = sourceSystem.particles.compactMap { particle in
                    let q = plan.embeddingChargesE[Int(particle.index)]
                    guard q != 0 else { return nil }
                    let p = cluster.particlePositionsNM[Int(particle.index)]/VivoAtomicUnits.bohrInNM
                    return .init(chargeE: q,positionBohr: .init(p.x,p.y,p.z),classicalParticleIndex: particle.index)
                }
                let energy: Double,nuclearForces: [VivoVector3D],chargeForces: [VivoVector3D],residual: Double
                let derivativeMethod: String
                var affine: VivoQMMatrix?
                switch plan.configuration.boundary {
                case .periodicElectrostatic:
                    guard let cell = geometry.periodicCell else { throw VivoChemistryError.invalid("periodic BO candidate has no cell") }
                    let result = try VivoPeriodicHartreeFock.evaluate(system: electronic,basis: cfg.basis,cell: cell,
                        configuration: cfg.periodic,budget: cfg.budget)
                    energy = result.reference.energyHartree;nuclearForces = result.nucleusForcesHartreePerBohr
                    chargeForces = result.pointChargeForcesHartreePerBohr;residual = result.reference.finalCommutatorNorm
                    affine = result.affineStrainDerivativeHartree
                    derivativeMethod = "analytic-gaussian-HF-Pulay+variational-quadrupolar-Ewald+exact-near-C2; reference-backend"
                case .finiteCluster:
                    guard geometry.periodicCell == nil else { throw VivoChemistryError.invalid("finite BO candidate unexpectedly has a cell") }
                    let ao = try VivoGaussianIntegralEngine.compute(system: electronic,basis: cfg.basis,budget: cfg.budget)
                    let reference = try VivoHartreeFock.solve(system: electronic,integrals: ao,configuration: cfg.periodic.scf,budget: cfg.budget)
                    let gradient = try VivoHartreeFockNuclearDerivatives.evaluate(integrals: ao,reference: reference,budget: cfg.budget)
                    energy = reference.energyHartree;nuclearForces = gradient.nuclearGradientsHartreePerBohr.map { $0 * -1 }
                    chargeForces = gradient.pointChargeGradientsHartreePerBohr.map { $0 * -1 };residual = reference.finalCommutatorNorm
                    derivativeMethod = "analytic-gaussian-HF-Pulay; finite-embedding"
                }
                let raw = try VivoQMMMForceMapper.projectCenters(electronicSystem: electronic,links: region.links,
                    atomToParticle: plan.atomToParticle,particlePositionsNM: cluster.particlePositionsNM,
                    nucleusForces: nuclearForces,pointChargeForces: chargeForces)
                let physical = try topology.siteGraph.redistribute(rawForces: raw,state: topology.siteGraph.construct(positionsNM: cluster.particlePositionsNM))
                if var strain = affine {
                    // Replace independent-center affine motion with the physical
                    // link/site Jacobian. The explicit lattice contribution is kept.
                    func virial(_ force: VivoVector3D,_ position: VivoVector3D,_ sign: Double) {
                        let f = [force.x,force.y,force.z],r = [position.x,position.y,position.z]
                        for a in 0..<3 { for b in 0..<3 { strain[a,b] += sign*f[a]*r[b] } }
                    }
                    for i in electronic.nuclei.indices {
                        let p = electronic.nuclei[i].positionBohr;virial(nuclearForces[i],.init(p.x,p.y,p.z),1)
                    }
                    for i in electronic.pointCharges.indices {
                        let p = electronic.pointCharges[i].positionBohr;virial(chargeForces[i],.init(p.x,p.y,p.z),1)
                    }
                    for i in physical.indices { virial(physical[i],cluster.particlePositionsNM[i]/VivoAtomicUnits.bohrInNM,-1) }
                    affine = try strain.scaled(VivoAtomicUnits.hartreeInKJPerMol)
                }
                let factor = VivoAtomicUnits.hartreeInKJPerMol/VivoAtomicUnits.bohrInNM
                try Task.checkCancellation()
                return try .init(providerFingerprint: id,geometry: geometry,additionalEnergyKJPerMol: energy*VivoAtomicUnits.hartreeInKJPerMol,
                    physicalParticleForcesKJPerMolNM: physical.map { $0*factor },derivativeMethod: derivativeMethod,
                    convergenceResidual: residual,requiredResidual: cfg.periodic.scf.commutatorTolerance,
                    additionalAffineStrainDerivativeKJPerMol: affine)
            })
    }
}
