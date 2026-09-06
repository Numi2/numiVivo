import Foundation

/// A finite-cluster sample, not a periodic QM/MM Hamiltonian and not an MD force
/// replacement. Original MD state and all source atom/particle identities survive.
public struct VivoMDQMMMPreparedFrame: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/md-qmmm-prepared-frame/v1"
    public let schema: String
    public let sourceStep: UInt64
    public let sourceTimePS: Double
    public let sourceConfigurationFingerprint: VivoFingerprint
    /// Present when prepared from an accepted checkpoint rather than a snapshot.
    public let sourceCheckpointFingerprint: VivoFingerprint?
    public let cluster: VivoQMMMFiniteCluster
    public let solventPromotion: VivoQMMMSolventPromotionResult?
    public let request: VivoQMMMRegionRequest
    public let region: VivoQMMMPreparedRegion

    /// Validate decoded artifacts before consuming forces. Preparation is replayed
    /// on the stored finite geometry, while the MD source identity remains separate.
    public func validate(document: VivoMolecularStructureDocument, system: VivoClassicalSystem) throws {
        let topology = try VivoQMMMTopology(document: document, system: system)
        guard schema == Self.schema, cluster.schema == VivoQMMMFiniteCluster.schema,
              sourceTimePS.isFinite, sourceTimePS >= 0,
              cluster.structureFingerprint == document.structureFingerprint,
              cluster.classicalSystemFingerprint == (try system.fingerprint()),
              cluster.atomToParticle == topology.atomToParticle,
              cluster.moleculeAtomIndices == topology.molecules,
              cluster.linearVirtualSites == topology.sites,
              cluster.particlePositionsNM.count == system.particles.count,
              cluster.particleImages.count == system.particles.count,
              request.coordinatesAreUnwrappedFiniteCluster else {
            throw VivoChemistryError.invalid("MD/QM/MM prepared-frame identity or mapping")
        }
        try cluster.configuration.validate()
        guard cluster.sourcePeriodicCell?.isValid != false,
              request.qmAtomIndices.contains(cluster.anchorAtomIndex), cluster.imagePairEvaluations >= 0 else {
            throw VivoChemistryError.invalid("MD/QM/MM source cell or anchor")
        }
        let siteByParticle = Dictionary(uniqueKeysWithValues: topology.sites.map { ($0.siteParticle, $0) })
        for (i, image) in cluster.particleImages.enumerated() {
            guard image.particleIndex == UInt32(i), image.structureAtomIndex == topology.particleToAtom[i] else {
                throw VivoChemistryError.invalid("MD/QM/MM particle identity changed")
            }
            if let atom = image.structureAtomIndex {
                guard image.moleculeIndex == topology.atomMolecule[Int(atom)], let coefficients = image.latticeImage,
                      coefficients.count == 3, coefficients.allSatisfy({ $0 > -(1 << 50) && $0 < (1 << 50) }),
                      cluster.sourcePeriodicCell != nil || coefficients == [0, 0, 0] else {
                    throw VivoChemistryError.invalid("MD/QM/MM physical image mapping")
                }
            } else {
                guard image.latticeImage == nil, let site = siteByParticle[UInt32(i)],
                      let parent = topology.particleToAtom[Int(site.parentParticles[0])],
                      image.moleculeIndex == topology.atomMolecule[Int(parent)] else {
                    throw VivoChemistryError.invalid("MD/QM/MM virtual-site image mapping")
                }
            }
        }
        guard try VivoQMMMCompiler.prepare(document: document, system: system,
                particlePositionsNM: cluster.particlePositionsNM, request: request) == region else {
            throw VivoChemistryError.invalid("MD/QM/MM region does not reconstruct from its finite frame")
        }
        if let promotion = solventPromotion {
            guard promotion.request == request, promotion.qmAtomIndices == request.qmAtomIndices,
                  !promotion.requiresUnwrappedFiniteClusterForQMMM else {
                throw VivoChemistryError.invalid("MD/QM/MM solvent selection and prepared region disagree")
            }
        }
    }
}

public extension VivoMDQMMMPreparedFrame {
    /// Optional source audit for decoded or transported artifacts. Requires the
    /// original immutable MD sample, not just the prepared geometry.
    func validateSource(document: VivoMolecularStructureDocument, system: VivoClassicalSystem,
                        snapshot: VivoMDStateSnapshot) throws {
        try validate(document: document, system: system)
        try snapshot.validate(particleCount: system.particles.count)
        guard cluster.sourceFrameFingerprint == (try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(snapshot))),
              snapshot.systemFingerprint == cluster.classicalSystemFingerprint,
              snapshot.configurationFingerprint == sourceConfigurationFingerprint,
              snapshot.stepIndex == sourceStep, snapshot.timePS == sourceTimePS,
              snapshot.periodicCell == cluster.sourcePeriodicCell else {
            throw VivoChemistryError.invalid("QM/MM source snapshot identity changed")
        }
        for image in cluster.particleImages {
            let slot = Int(image.particleIndex)
            if let coefficients = image.latticeImage {
                var position = snapshot.positionsNM[slot]
                if let cell = snapshot.periodicCell {
                    position = position + cell.a * Double(coefficients[0])
                        + cell.b * Double(coefficients[1]) + cell.c * Double(coefficients[2])
                }
                guard (position - cluster.particlePositionsNM[slot]).norm <= cluster.configuration.cycleToleranceNM else {
                    throw VivoChemistryError.invalid("QM/MM lattice mapping does not reconstruct its source particle")
                }
            } else {
                let residual = try VivoMDPreparationGeometry.minimumImage(
                    snapshot.positionsNM[slot] - cluster.particlePositionsNM[slot], cell: snapshot.periodicCell)
                guard residual.norm <= cluster.configuration.virtualSiteToleranceNM else {
                    throw VivoChemistryError.invalid("QM/MM virtual site differs from its sampled source")
                }
            }
        }
    }
}

public extension VivoQMMMCompiler {
    /// Accepted-checkpoint adapter preserves the checkpoint payload digest as well
    /// as the sampled state digest. Numerical-contract validation is not bypassed.
    static func prepareMD(document: VivoMolecularStructureDocument, system: VivoClassicalSystem,
                          checkpoint: VivoMDCheckpoint, request: VivoQMMMRegionRequest,
                          solventPolicy: VivoQMMMSolventPromotionPolicy? = nil,
                          configuration: VivoQMMMClusterConfiguration = .init(),
                          budget: VivoChemistryBudget = .init()) throws -> VivoMDQMMMPreparedFrame {
        try checkpoint.validate(particleCount: system.particles.count)
        let snapshot = VivoMDStateSnapshot(systemFingerprint: checkpoint.systemFingerprint,
            configurationFingerprint: checkpoint.configurationFingerprint, stepIndex: checkpoint.acceptedStep,
            timePS: checkpoint.timePS, positionsNM: checkpoint.positionsNM,
            velocitiesNMPerPS: checkpoint.velocitiesNMPerPS, periodicCell: checkpoint.periodicCell)
        let prepared = try prepareMD(document: document, system: system, snapshot: snapshot, request: request,
            solventPolicy: solventPolicy, configuration: configuration, budget: budget)
        return .init(schema: prepared.schema, sourceStep: prepared.sourceStep, sourceTimePS: prepared.sourceTimePS,
            sourceConfigurationFingerprint: prepared.sourceConfigurationFingerprint,
            sourceCheckpointFingerprint: try checkpoint.fingerprint(), cluster: prepared.cluster,
            solventPromotion: prepared.solventPromotion, request: prepared.request, region: prepared.region)
    }

    /// The sampled cell is authoritative, including nil for a nonperiodic sample.
    /// The input snapshot is immutable; this is an explicit host-side sampling step.
    static func prepareMD(document: VivoMolecularStructureDocument, system: VivoClassicalSystem,
                          snapshot: VivoMDStateSnapshot, request base: VivoQMMMRegionRequest,
                          solventPolicy: VivoQMMMSolventPromotionPolicy? = nil,
                          configuration: VivoQMMMClusterConfiguration = .init(),
                          budget: VivoChemistryBudget = .init()) throws -> VivoMDQMMMPreparedFrame {
        try snapshot.validate(particleCount: system.particles.count); try budget.validate()
        guard snapshot.systemFingerprint == (try system.fingerprint()) else {
            throw VivoChemistryError.invalid("MD snapshot does not identify the supplied classical system")
        }
        _ = try budget.elements([system.particles.count, 3], simultaneousArrays: 12)
        var limits = configuration
        limits.maximumImagePairs = min(limits.maximumImagePairs, budget.maximumOperatorApplications)
        let frameID = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(snapshot))
        let cluster = try VivoQMMMClusterBuilder.reconstruct(document: document, system: system,
            particlePositionsNM: snapshot.positionsNM, periodicCell: snapshot.periodicCell,
            sourceFrameFingerprint: frameID, qmAtomIndices: base.qmAtomIndices,
            solventPolicy: solventPolicy, configuration: limits)
        var request = base
        // Set only AFTER actual reconstruction, never to bypass wrapped input.
        request.coordinatesAreUnwrappedFiniteCluster = true
        var promotion: VivoQMMMSolventPromotionResult?
        if let policy = solventPolicy {
            let remaining = limits.maximumImagePairs - cluster.imagePairEvaluations
            guard remaining > 0 else { throw VivoChemistryError.resourceLimit("QM/MM solvent selection has no remaining geometry budget") }
            // Use Euclidean distances in the exact coordinates submitted to QM/MM.
            // A second atomwise minimum-image selection would be inconsistent.
            let selected = try VivoQMMMSolventPromotion.select(document: document, system: system,
                positions: cluster.particlePositionsNM, base: request, policy: policy,
                distanceCell: nil, requiresUnwrapping: false, maximumDistancePairs: remaining)
            request = selected.request
            promotion = .init(schema: selected.schema, policy: policy,
                originalQMAtomIndices: selected.originalQMAtomIndices, promotedResidues: selected.promotedResidues,
                qmAtomIndices: selected.qmAtomIndices, addedAlphaElectrons: selected.addedAlphaElectrons,
                addedBetaElectrons: selected.addedBetaElectrons, periodicCellUsed: snapshot.periodicCell != nil,
                requiresUnwrappedFiniteClusterForQMMM: false, request: request,
                interpretation: "whole connected solvent molecules selected in the reconstructed finite cluster; current MD cell used for whole-molecule imaging; physical atoms and virtual sites retain source mappings; closed-shell promotion assumption retained")
        }
        let region = try prepare(document: document, system: system,
            particlePositionsNM: cluster.particlePositionsNM, request: request)
        return .init(schema: VivoMDQMMMPreparedFrame.schema, sourceStep: snapshot.stepIndex,
            sourceTimePS: snapshot.timePS, sourceConfigurationFingerprint: snapshot.configurationFingerprint,
            sourceCheckpointFingerprint: nil, cluster: cluster, solventPromotion: promotion, request: request, region: region)
    }
}
