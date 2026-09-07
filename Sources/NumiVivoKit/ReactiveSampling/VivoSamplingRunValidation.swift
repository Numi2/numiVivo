import Foundation

public extension VivoRingPolymerRun {
    /// Validate stored state/energy accounting, not a fresh electronic replay.
    func validate(recordEvery: Int = 1) throws {
        try start.validate(); try end.validate()
        let cfg = start.configuration, n = start.definition.atomIndices.count, p = cfg.beadCount
        guard start.definition == end.definition, cfg == end.configuration, end.sweep > start.sweep,
              recordEvery > 0, end.sweep-start.sweep <= 1_000_000,
              interpretation == VivoRingPolymerSampling.interpretation else { throw VivoChemistryError.invalid("ring run source or scope") }
        let sweeps = Int(end.sweep-start.sweep)
        guard forceEvaluations == p*(1+sweeps*cfg.integrationSteps), forceEvaluations <= cfg.maximumForceEvaluations else {
            throw VivoChemistryError.invalid("ring evaluation accounting")
        }
        let indices = (1...sweeps).map { start.sweep+UInt64($0) }.filter { $0%UInt64(recordEvery) == 0 }
        guard observations.map(\.sweep) == indices,
              Double(sweeps)*Double(p)*Double(n)*3 <= Double(cfg.maximumPrimitiveWork) else {
            throw VivoChemistryError.invalid("ring stored sweep sequence or replay budget")
        }
        // Each HMC sweep draws 3*N*P normal variates and one acceptance uniform,
        // independently of acceptance. Thin output still binds the final RNG state.
        var replay = start.randomState, observationIndex = 0
        for offset in 1...sweeps {
            for _ in 0..<(3*n*p) { _ = replay.normal() }
            let logU = log(max(replay.unitInterval(),Double.leastNonzeroMagnitude))
            let sweep = start.sweep+UInt64(offset)
            if observationIndex < observations.count, observations[observationIndex].sweep == sweep {
                let record = observations[observationIndex]
                guard record.accepted == (logU < record.logAcceptanceProbability) else {
                    throw VivoChemistryError.invalid("ring acceptance RNG replay differs from stored decision")
                }
                observationIndex += 1
            }
        }
        guard replay == end.randomState else { throw VivoChemistryError.invalid("ring checkpoint RNG history mismatch") }
        let rt = VivoAtomicUnits.gasConstantJPerMolK*cfg.temperatureK/1000
        var previous = start.beadPositionsNM
        for record in observations {
            let q = record.beadPositionsNM
            guard q.count == p, q.allSatisfy({ $0.count == n && $0.allSatisfy(\.isFinite) }),
                  [record.logAcceptanceProbability,record.dimensionlessIntegrationWork,record.meanPotentialEnergyKJPerMol,
                   record.springEnergyKJPerMol,record.primitiveTotalEnergyKJPerMol].allSatisfy(\.isFinite),
                  record.logAcceptanceProbability == min(0,-record.dimensionlessIntegrationWork) else {
                throw VivoChemistryError.invalid("ring observation geometry or finite work")
            }
            if recordEvery == 1, !record.accepted, q != previous { throw VivoChemistryError.invalid("rejected ring state changed") }
            let spring = try VivoRingPolymerSampling.springEnergy(positions: q,massesDa: start.definition.massesDa,temperatureK: cfg.temperatureK)
            let centroid = (0..<n).map { i in q.reduce(VivoVector3D.zero) { $0+$1[i] }/Double(p) }
            let estimator = 1.5*Double(n*p)*rt-spring/Double(p)+record.meanPotentialEnergyKJPerMol
            guard record.springEnergyKJPerMol == spring, record.centroidPositionsNM == centroid,
                  record.primitiveTotalEnergyKJPerMol == estimator else { throw VivoChemistryError.invalid("ring energy/centroid reconstruction") }
            previous = q
        }
        if observations.last?.sweep == end.sweep, observations.last?.beadPositionsNM != end.beadPositionsNM {
            throw VivoChemistryError.invalid("ring final checkpoint differs from last retained sample")
        }
    }
}

public extension VivoReactiveSamplingRun {
    /// Replay RNG/Metropolis bookkeeping using stored authoritative endpoint
    /// labels. This does not reevaluate the force authority or a GPU trajectory.
    func validate(model: VivoReactiveSurrogateModel) throws {
        try start.validate(model: model); try end.validate(model: model)
        let cfg = start.configuration, definition = model.payload.authorityDefinition
        try initialAuthority.validate(definition: definition,positionsNM: start.positionsNM)
        guard end.sweep > start.sweep, end.sweep-start.sweep <= 1_000_000,
              end.modelFingerprint == start.modelFingerprint, end.authorityFingerprint == start.authorityFingerprint,
              end.configuration == cfg, end.backendProfile == start.backendProfile, end.chainIdentifier == start.chainIdentifier,
              observations.count == Int(end.sweep-start.sweep), authorityEvaluations == authorityLabels.count+1,
              authorityEvaluations <= cfg.maximumAuthorityEvaluations, baselineEvaluations <= cfg.maximumBaselineEvaluations,
              interpretation == VivoReactiveSurrogateSampling.interpretation else { throw VivoChemistryError.invalid("reactive run binding or evaluation count") }
        var rng = start.randomState, current = initialAuthority, q = start.positionsNM, labelIndex = 0
        let rt = VivoAtomicUnits.gasConstantJPerMolK*cfg.temperatureK/1000
        for (offset,record) in observations.enumerated() {
            let refresh = rng.unitInterval() < cfg.authorityRefreshProbability
            var kinetic = 0.0, refreshPositions: [VivoVector3D]?
            if refresh {
                refreshPositions = try definition.canonicalPositions(q.map { $0+VivoVector3D(rng.normal(),rng.normal(),rng.normal())*cfg.authorityRefreshLengthNM })
            } else {
                for mass in definition.massesDa {
                    let s = sqrt(rt/mass), v = VivoVector3D(rng.normal()*s,rng.normal()*s,rng.normal()*s)
                    kinetic += 0.5*mass*v.squaredNorm
                }
            }
            guard record.proposalKind == (refresh ? .authorityRandomWalk : .surrogateHMC), !refresh || record.rejectionReasons.isEmpty else {
                throw VivoChemistryError.invalid("reactive state-independent proposal mixture replay")
            }
            guard record.sweep == start.sweep+UInt64(offset+1), record.initialKineticEnergyKJPerMol == kinetic else {
                throw VivoChemistryError.invalid("reactive momentum RNG or sweep replay")
            }
            if record.rejectionReasons.isEmpty {
                guard labelIndex < authorityLabels.count, let finalKinetic = record.proposedKineticEnergyKJPerMol,
                      finalKinetic.isFinite, finalKinetic >= 0, let storedLogA = record.logAcceptanceProbability,
                      let storedUniform = record.logUniform else { throw VivoChemistryError.invalid("missing reactive endpoint/work evidence") }
                let label = authorityLabels[labelIndex]; labelIndex += 1
                let proposed = label.authority.requestedPositionsNM
                if refresh, (proposed != refreshPositions || finalKinetic != 0) { throw VivoChemistryError.invalid("authority refresh displacement replay") }
                try label.authority.validate(definition: definition,positionsNM: proposed)
                try label.baseline.validate(definition: model.payload.baselineDefinition,positionsNM: proposed)
                guard label.identifier == "\(start.chainIdentifier):sweep-\(record.sweep)", label.sourceGroup == start.chainIdentifier else {
                    throw VivoChemistryError.invalid("reactive acquisition source group")
                }
                let logA = min(0,-((label.authority.energyKJPerMol+finalKinetic)-(current.energyKJPerMol+kinetic))/rt)
                let logU = log(max(rng.unitInterval(),Double.leastNonzeroMagnitude))
                guard storedLogA == logA, storedUniform == logU, record.accepted == (logU < logA),
                      (refresh ? record.authorityMinusPredictedEnergyKJPerMol?.isFinite != false : record.authorityMinusPredictedEnergyKJPerMol?.isFinite == true) else {
                    throw VivoChemistryError.invalid("reactive acceptance accounting differs from authoritative endpoint energies")
                }
                if record.accepted { q = proposed; current = label.authority }
            } else {
                guard !record.accepted, record.logUniform == nil, record.logAcceptanceProbability == nil,
                      record.proposedKineticEnergyKJPerMol == nil, record.authorityMinusPredictedEnergyKJPerMol == nil else {
                    throw VivoChemistryError.invalid("domain rejection fabricated an authority acceptance")
                }
            }
            guard record.positionsNM == q, record.authorityEnergyKJPerMol == current.energyKJPerMol else {
                throw VivoChemistryError.invalid("reactive retained sample mismatch")
            }
        }
        guard end.randomState == rng, end.positionsNM == q, labelIndex == authorityLabels.count,
              acquisitions.allSatisfy({ $0.modelFingerprint == model.fingerprint && !$0.reason.isEmpty
                  && $0.positionsNM.count == definition.atomIndices.count && $0.positionsNM.allSatisfy(\.isFinite) }) else {
            throw VivoChemistryError.invalid("reactive checkpoint/RNG/acquisition replay mismatch")
        }
    }
}
