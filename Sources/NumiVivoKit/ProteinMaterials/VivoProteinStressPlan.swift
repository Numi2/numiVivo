import Foundation

public struct VivoProteinStressStage: Codable, Sendable, Equatable {
    public let name: String
    public let steps: UInt64
    public let temperatureK: Double
    /// A fixed harmonic reference for this stage; nil for an unbiased thermal run.
    public let pullReferenceNM: Double?
    public init(name: String, steps: UInt64, temperatureK: Double, pullReferenceNM: Double? = nil) {
        self.name = name; self.steps = steps; self.temperatureK = temperatureK; self.pullReferenceNM = pullReferenceNM
    }
}

/// An executable, finite protocol over a prepared, exact MD checkpoint. Preparation
/// and equilibration are separate existing MD workflows; no duration is invented.
public struct VivoProteinStressRequest: Codable, Sendable, Equatable {
    public static let schemaID = "numivivo.org/protein-stress-request/v1"
    public var schema: String = Self.schemaID
    public let system: VivoClassicalSystem
    public let sourceConfiguration: VivoMDConfiguration
    public let sourceCheckpoint: VivoMDCheckpoint
    public let sourceDescription: String
    public let replicaID: String
    public let randomSeed: UInt64
    public let pull: VivoProteinPullDefinition?
    public let stages: [VivoProteinStressStage]
    public let sampleEvery: UInt64
    public let selection: [UInt32]
    public let hydrogenBonds: [VivoProteinHydrogenBond]
    public let hydrogenBondCriteria: VivoProteinHydrogenBondCriteria
    public let nativeContacts: [VivoProteinNativeContact]
    public let maximumContactDistanceRatio: Double
    public init(system: VivoClassicalSystem, sourceConfiguration: VivoMDConfiguration,
                sourceCheckpoint: VivoMDCheckpoint, sourceDescription: String, replicaID: String,
                randomSeed: UInt64, pull: VivoProteinPullDefinition? = nil,
                stages: [VivoProteinStressStage], sampleEvery: UInt64, selection: [UInt32],
                hydrogenBonds: [VivoProteinHydrogenBond] = [],
                hydrogenBondCriteria: VivoProteinHydrogenBondCriteria = .init(),
                nativeContacts: [VivoProteinNativeContact] = [], maximumContactDistanceRatio: Double = 1.2) {
        self.system = system; self.sourceConfiguration = sourceConfiguration; self.sourceCheckpoint = sourceCheckpoint
        self.sourceDescription = sourceDescription; self.replicaID = replicaID; self.randomSeed = randomSeed
        self.pull = pull; self.stages = stages; self.sampleEvery = sampleEvery; self.selection = selection
        self.hydrogenBonds = hydrogenBonds; self.hydrogenBondCriteria = hydrogenBondCriteria
        self.nativeContacts = nativeContacts; self.maximumContactDistanceRatio = maximumContactDistanceRatio
    }
    public func fingerprint() throws -> VivoFingerprint { try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(self)) }
    public func configuration(for stage: Int) throws -> VivoMDConfiguration {
        guard stages.indices.contains(stage) else { throw VivoProteinStressError.invalid("stage index") }
        var result = sourceConfiguration
        result.targetTemperatureK = stages[stage].temperatureK; result.randomSeed = randomSeed
        try result.validate(); return result
    }
}

public struct VivoProteinStressFrame: Codable, Sendable, Equatable {
    public let coordinateNM: Double?
    public let tensileForcePN: Double?
    public let restraintEnergyKJPerMol: Double
    public let hydrogenBondsPresent: [Bool]
    public let structure: VivoProteinStructuralRetention
}

/// Native physical-particle and periodic-image owner shared by force evaluation,
/// analysis, stage transfer, and receipt verification.
public struct VivoProteinStressCompilation: Sendable {
    public let request: VivoProteinStressRequest
    public let layout: VivoProteinWholeMoleculeLayout
    public let massesDa: [Double]
    public let maximumJournalEntries: Int

    public init(_ request: VivoProteinStressRequest) throws {
        try VivoClassicalSystemValidator.validate(request.system)
        try request.sourceConfiguration.validate()
        try request.sourceCheckpoint.validate(particleCount: request.system.particles.count)
        guard request.schema == VivoProteinStressRequest.schemaID,
              !request.sourceDescription.isEmpty, request.sourceDescription.utf8.count <= 4096,
              !request.replicaID.isEmpty, request.replicaID.utf8.count <= 128,
              request.system.particles.count <= 200_000,
              request.sourceCheckpoint.systemFingerprint == (try request.system.fingerprint()),
              request.sourceCheckpoint.configurationFingerprint == (try request.sourceConfiguration.fingerprint()),
              request.sourceConfiguration.ensemble == .nvt,
              request.sourceConfiguration.thermostat == .langevinMiddle,
              request.sourceConfiguration.barostat == .none,
              request.sourceConfiguration.resolvedPositionPrecision == .fp32,
              (request.sourceCheckpoint.positionPrecision ?? .fp32) == .fp32,
              request.system.polarization == nil,
              !request.system.particles.contains(where: { $0.role == .drude }),
              !request.stages.isEmpty, request.stages.count <= 4096, request.sampleEvery > 0,
              Set(request.stages.map(\.name)).count == request.stages.count,
              !request.selection.isEmpty, Set(request.selection).count == request.selection.count,
              request.hydrogenBonds.count <= 10_000, request.nativeContacts.count <= 100_000 else {
            throw VivoProteinStressError.invalid("unsupported source identity, NVT profile, selection, or protocol size")
        }
        var steps: UInt64 = 0, entries: UInt64 = 1
        for stage in request.stages {
            guard !stage.name.isEmpty, stage.name.utf8.count <= 128, stage.steps > 0, stage.steps <= 1_000_000,
                  stage.temperatureK.isFinite, stage.temperatureK > 0, Float(stage.temperatureK).isFinite,
                  (request.pull == nil) == (stage.pullReferenceNM == nil), stage.pullReferenceNM?.isFinite != false else {
                throw VivoProteinStressError.invalid("invalid stress stage")
            }
            if request.pull?.projectionAxis == nil, let value = stage.pullReferenceNM, value < 0 {
                throw VivoProteinStressError.invalid("negative radial restraint reference")
            }
            steps += stage.steps; entries += 1 + (stage.steps - 1) / request.sampleEvery + 1
        }
        guard steps <= 1_000_000, request.sourceCheckpoint.acceptedStep <= UInt64.max - steps,
              entries <= 20_000,
              entries * UInt64(request.system.particles.count) <= 20_000_000,
              entries * UInt64(max(1, request.hydrogenBonds.count)) <= 5_000_000 else {
            throw VivoProteinStressError.invalid("declared steps, checkpoint storage, or analysis budget exceeded")
        }
        let finalTime = request.sourceCheckpoint.timePS + Double(steps) * request.sourceConfiguration.timeStepPS
        guard finalTime.isFinite, request.sourceCheckpoint.timePS + request.sourceConfiguration.timeStepPS > request.sourceCheckpoint.timePS else {
            throw VivoProteinStressError.invalid("stress clock would overflow or fail to advance")
        }
        let masses = request.system.particles.map(\.massDa)
        try request.pull?.validate(massesDa: masses)
        let selected = request.selection + (request.pull?.reference.particles ?? []) + (request.pull?.moving.particles ?? [])
            + request.hydrogenBonds.flatMap { [$0.donor, $0.hydrogen, $0.acceptor] }
        guard selected.allSatisfy({ Int($0) < masses.count && request.system.particles[Int($0)].role == .atom && masses[Int($0)] > 0 }) else {
            throw VivoProteinStressError.invalid("stress analysis requires physical atoms, not virtual sites")
        }
        let connectivity = try VivoMDBarostatPlan.make(system: request.system)
        let component = connectivity.componentIndexByParticle[Int(request.selection[0])]
        guard selected.allSatisfy({ connectivity.componentIndexByParticle[Int($0)] == component }) else {
            throw VivoProteinStressError.invalid("all pull and analysis atoms must belong to one connected molecule")
        }
        let range = Int(connectivity.componentOffsets[Int(component)])..<Int(connectivity.componentOffsets[Int(component) + 1])
        let members = Array(connectivity.componentParticles[range])
        let edgeRange = Int(connectivity.edgeOffsets[Int(component)])..<Int(connectivity.edgeOffsets[Int(component) + 1])
        let layout = VivoProteinWholeMoleculeLayout(orderedMembers: members,
            parents: members.map { connectivity.parentByParticle[Int($0)] },
            closureEdges: connectivity.edges[edgeRange].map { [$0.x, $0.y] })
        try layout.validate(particleCount: masses.count)
        self.request = request; self.massesDa = masses; self.layout = layout; self.maximumJournalEntries = Int(entries) + 1
        _ = try frame(checkpoint: request.sourceCheckpoint, stage: nil)
        for index in request.stages.indices {
            let config = try request.configuration(for: index)
            let provider = try forceProvider(stage: index)
            let initial = VivoClassicalInitialState(systemFingerprint: request.sourceCheckpoint.systemFingerprint,
                positionsNM: request.sourceCheckpoint.positionsNM, periodicCell: request.sourceCheckpoint.periodicCell,
                sourceTimePS: request.sourceCheckpoint.timePS)
            let capabilities = try VivoMDCapabilityAnalyzer.analyze(system: request.system, initialState: initial,
                configuration: config, forceProvider: provider)
            guard capabilities.executable else { throw VivoProteinStressError.invalid(capabilities.blockers.joined(separator: "; ")) }
        }
    }

    public func wholePositions(_ geometry: VivoMDCandidateGeometry) throws -> [VivoVector3D] {
        guard geometry.particlePositionsNM.count == massesDa.count,
              geometry.periodicCell == request.sourceCheckpoint.periodicCell else {
            throw VivoProteinStressError.invalid("stress geometry changed particle count or fixed cell")
        }
        return try layout.reconstruct(positionsNM: geometry.particlePositionsNM) { delta in
            guard let cell = geometry.periodicCell else { return delta }
            return try cell.minimumImage(delta)
        }
    }
    public func forceProvider(stage: Int) throws -> VivoMDCandidateForceProvider? {
        guard request.stages.indices.contains(stage) else { throw VivoProteinStressError.invalid("force stage index") }
        guard let pull = request.pull, let reference = request.stages[stage].pullReferenceNM else { return nil }
        struct Identity: Codable {
            let method: String
            let system: VivoFingerprint
            let pull: VivoProteinPullDefinition
            let referenceNM: Double
            let layout: VivoProteinWholeMoleculeLayout
            let fixedPeriodicCell: VivoPeriodicCell?
        }
        let systemID = try request.system.fingerprint()
        let identity = Identity(method: "numivivo.org/whole-molecule-harmonic-restraint/v1", system: systemID,
                                pull: pull, referenceNM: reference, layout: layout, fixedPeriodicCell: request.sourceCheckpoint.periodicCell)
        let fingerprint = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(identity))
        return try .init(fingerprint: fingerprint, retainedSystemFingerprint: systemID,
            boundary: request.sourceCheckpoint.periodicCell == nil ? .finiteCluster : .periodicElectrostatic,
            supportsCellMoves: false, maximumAcceptedResidual: 1e-12, molecularConnectivitySystem: request.system) { geometry in
                let whole = try self.wholePositions(geometry)
                let evaluation = try VivoProteinPullMath.evaluate(pull, wholePositionsNM: whole,
                    massesDa: self.massesDa, referenceNM: reference)
                return try .init(providerFingerprint: fingerprint, geometry: geometry,
                    additionalEnergyKJPerMol: evaluation.energyKJPerMol,
                    physicalParticleForcesKJPerMolNM: evaluation.forcesKJPerMolNM,
                    derivativeMethod: "analytic mass-weighted COM harmonic gradient; fixed periodic images",
                    convergenceResidual: 0, requiredResidual: 1e-12)
            }
    }
    public func configurationFingerprint(stage: Int) throws -> VivoFingerprint {
        try VivoMDCandidateForceProvider.executionFingerprint(configuration: request.configuration(for: stage), provider: forceProvider(stage: stage))
    }
    public func frame(checkpoint: VivoMDCheckpoint, stage: Int?) throws -> VivoProteinStressFrame {
        if let stage, !request.stages.indices.contains(stage) { throw VivoProteinStressError.invalid("frame stage index") }
        let geometry = try VivoMDCandidateGeometry(particlePositionsNM: checkpoint.positionsNM, periodicCell: checkpoint.periodicCell)
        let whole = try wholePositions(geometry)
        let pullEvaluation: VivoProteinPullEvaluation?
        if let stage, let pull = request.pull, let reference = request.stages[stage].pullReferenceNM {
            pullEvaluation = try VivoProteinPullMath.evaluate(pull, wholePositionsNM: whole, massesDa: massesDa, referenceNM: reference)
        } else { pullEvaluation = nil }
        let hydrogen = try VivoProteinHydrogenBondAnalysis.evaluate(request.hydrogenBonds, positionsNM: whole,
                                                                    criteria: request.hydrogenBondCriteria)
        let structure = try VivoProteinStructuralAnalysis.evaluate(positionsNM: whole, massesDa: massesDa,
            selection: request.selection, contacts: request.nativeContacts, maximumDistanceRatio: request.maximumContactDistanceRatio)
        return .init(coordinateNM: pullEvaluation?.coordinateNM, tensileForcePN: pullEvaluation?.tensileForcePN,
                     restraintEnergyKJPerMol: pullEvaluation?.energyKJPerMol ?? 0,
                     hydrogenBondsPresent: hydrogen.map(\.present), structure: structure)
    }

    /// A declared Hamiltonian/thermostat change at one exact accepted state.
    /// No state import, velocity regeneration, clock reset, or RNG-step reset.
    public func transfer(_ checkpoint: VivoMDCheckpoint, from sourceStage: Int?, to destination: Int) throws
    -> (checkpoint: VivoMDCheckpoint, workKJPerMol: Double) {
        try checkpoint.validate(particleCount: massesDa.count)
        let sourceID = try sourceStage.map { try configurationFingerprint(stage: $0) } ?? request.sourceConfiguration.fingerprint()
        guard destination == (sourceStage.map { $0 + 1 } ?? 0), request.stages.indices.contains(destination),
              checkpoint.systemFingerprint == (try request.system.fingerprint()), checkpoint.configurationFingerprint == sourceID,
              checkpoint.periodicCell == request.sourceCheckpoint.periodicCell else {
            throw VivoProteinStressError.invalid("stage transfer source identity or order")
        }
        let before = try frame(checkpoint: checkpoint, stage: sourceStage)
        let after = try frame(checkpoint: checkpoint, stage: destination)
        let work = after.restraintEnergyKJPerMol - before.restraintEnergyKJPerMol
        guard work.isFinite else { throw VivoProteinStressError.invalid("stage work overflow") }
        var next = checkpoint
        next.configurationFingerprint = try configurationFingerprint(stage: destination)
        try next.validate(particleCount: massesDa.count)
        return (next, work)
    }
}
