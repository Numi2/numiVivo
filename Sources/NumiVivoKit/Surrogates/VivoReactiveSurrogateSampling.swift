import Foundation

/// Backend approximation affects proposal efficiency; full authority energy is
/// still used in the Metropolis test. The backend is frozen for the entire run.
public struct VivoReactiveEnergyForceBackend: Sendable {
    public let modelFingerprint: VivoFingerprint
    public let numericalProfile: String
    public let predict: @Sendable ([[VivoVector3D]]) async throws -> [VivoReactiveSurrogatePrediction]
    public init(modelFingerprint: VivoFingerprint, numericalProfile: String,
                predict: @escaping @Sendable ([[VivoVector3D]]) async throws -> [VivoReactiveSurrogatePrediction]) {
        self.modelFingerprint = modelFingerprint; self.numericalProfile = numericalProfile; self.predict = predict
    }
}
public struct VivoReactiveSamplingConfiguration: Codable, Sendable, Equatable {
    public var temperatureK: Double
    public var timeStepPS: Double
    public var integrationSteps: Int
    /// State-independent exact-energy random-walk moves retain access outside the learned domain.
    public var authorityRefreshProbability: Double
    public var authorityRefreshLengthNM: Double
    public var maximumAuthorityEvaluations: Int
    public var maximumBaselineEvaluations: Int
    public var maximumStoredCoordinateElements: Int
    public init(temperatureK: Double, timeStepPS: Double = 0.0005, integrationSteps: Int = 16,
                authorityRefreshProbability: Double = 0.1, authorityRefreshLengthNM: Double = 0.02,
                maximumAuthorityEvaluations: Int = 100_000, maximumBaselineEvaluations: Int = 1_000_000,
                maximumStoredCoordinateElements: Int = 10_000_000) {
        self.temperatureK = temperatureK; self.timeStepPS = timeStepPS; self.integrationSteps = integrationSteps
        self.authorityRefreshProbability = authorityRefreshProbability; self.authorityRefreshLengthNM = authorityRefreshLengthNM
        self.maximumAuthorityEvaluations = maximumAuthorityEvaluations; self.maximumBaselineEvaluations = maximumBaselineEvaluations
        self.maximumStoredCoordinateElements = maximumStoredCoordinateElements
    }
    public func validate() throws {
        guard temperatureK.isFinite, temperatureK > 0, temperatureK <= 100_000,
              timeStepPS.isFinite, timeStepPS > 0, timeStepPS <= 0.1,
              (1...100_000).contains(integrationSteps), authorityRefreshProbability.isFinite, authorityRefreshProbability > 0, authorityRefreshProbability <= 1,
              authorityRefreshLengthNM.isFinite, authorityRefreshLengthNM > 0, authorityRefreshLengthNM <= 10, maximumAuthorityEvaluations > 0,
              maximumBaselineEvaluations > 0, maximumStoredCoordinateElements > 0 else {
            throw VivoChemistryError.invalid("reactive equilibrium sampling settings")
        }
    }
}
public struct VivoReactiveSamplingCheckpoint: Codable, Sendable, Equatable {
    public let schema: String
    public let modelFingerprint: VivoFingerprint
    public let authorityFingerprint: VivoFingerprint
    public let configuration: VivoReactiveSamplingConfiguration
    public let backendProfile: String
    public let chainIdentifier: String
    public let sweep: UInt64
    public let randomState: VivoSplitMix64
    public let positionsNM: [VivoVector3D]
    public let fingerprint: VivoFingerprint
    public init(model: VivoReactiveSurrogateModel, configuration: VivoReactiveSamplingConfiguration,
                chainIdentifier: String, seed: UInt64, positionsNM: [VivoVector3D], backendProfile: String = "native-fp64-rbf/v1") throws {
        self = try Self.make(model: model,configuration: configuration,chainIdentifier: chainIdentifier,
            sweep: 0,randomState: .init(state: seed),positionsNM: model.payload.authorityDefinition.canonicalPositions(positionsNM),backendProfile: backendProfile)
    }
    static func make(model: VivoReactiveSurrogateModel, configuration: VivoReactiveSamplingConfiguration,
                     chainIdentifier: String, sweep: UInt64, randomState: VivoSplitMix64,
                     positionsNM: [VivoVector3D], backendProfile: String) throws -> Self {
        try model.validate(); try configuration.validate()
        guard !chainIdentifier.isEmpty, chainIdentifier.utf8.count <= 512, !backendProfile.isEmpty,
              positionsNM.count == model.payload.authorityDefinition.atomIndices.count,
              positionsNM.allSatisfy(\.isFinite) else { throw VivoChemistryError.invalid("reactive checkpoint geometry or identity") }
        let authorityID = try model.payload.authorityDefinition.fingerprint()
        struct Identity: Encodable {
            let schema: String; let modelFingerprint: VivoFingerprint; let authorityFingerprint: VivoFingerprint
            let configuration: VivoReactiveSamplingConfiguration; let backendProfile: String; let chainIdentifier: String
            let sweep: UInt64; let randomState: VivoSplitMix64; let positionsNM: [VivoVector3D]
        }
        let id = Identity(schema: "numivivo.org/reactive-equilibrium-checkpoint/v1",modelFingerprint: model.fingerprint,
            authorityFingerprint: authorityID,configuration: configuration,backendProfile: backendProfile,
            chainIdentifier: chainIdentifier,sweep: sweep,randomState: randomState,positionsNM: positionsNM)
        let fp = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(id))
        return .init(schema: id.schema,modelFingerprint: model.fingerprint,authorityFingerprint: authorityID,
            configuration: configuration,backendProfile: backendProfile,chainIdentifier: chainIdentifier,
            sweep: sweep,randomState: randomState,positionsNM: positionsNM,fingerprint: fp)
    }
    private init(schema: String, modelFingerprint: VivoFingerprint, authorityFingerprint: VivoFingerprint,
                 configuration: VivoReactiveSamplingConfiguration, backendProfile: String, chainIdentifier: String,
                 sweep: UInt64, randomState: VivoSplitMix64, positionsNM: [VivoVector3D], fingerprint: VivoFingerprint) {
        self.schema = schema; self.modelFingerprint = modelFingerprint; self.authorityFingerprint = authorityFingerprint
        self.configuration = configuration; self.backendProfile = backendProfile; self.chainIdentifier = chainIdentifier
        self.sweep = sweep; self.randomState = randomState; self.positionsNM = positionsNM; self.fingerprint = fingerprint
    }
    public func validate(model: VivoReactiveSurrogateModel) throws {
        guard self == (try Self.make(model: model,configuration: configuration,chainIdentifier: chainIdentifier,
            sweep: sweep,randomState: randomState,positionsNM: positionsNM,backendProfile: backendProfile)) else {
            throw VivoChemistryError.invalid("reactive sampler checkpoint payload/model mismatch")
        }
    }
}
public enum VivoReactiveProposalKind: String, Codable, Sendable { case surrogateHMC, authorityRandomWalk }
public struct VivoReactiveSamplingObservation: Codable, Sendable, Equatable {
    public let sweep: UInt64
    public let proposalKind: VivoReactiveProposalKind
    public let positionsNM: [VivoVector3D]
    public let authorityEnergyKJPerMol: Double
    public let initialKineticEnergyKJPerMol: Double
    public let proposedKineticEnergyKJPerMol: Double?
    public let logUniform: Double?
    public let accepted: Bool
    public let logAcceptanceProbability: Double?
    public let authorityMinusPredictedEnergyKJPerMol: Double?
    public let rejectionReasons: [String]
}
public struct VivoReactiveAcquisition: Codable, Sendable, Equatable {
    public let positionsNM: [VivoVector3D]
    public let reason: [String]
    public let modelFingerprint: VivoFingerprint
}
public struct VivoReactiveSamplingRun: Codable, Sendable, Equatable {
    public let start: VivoReactiveSamplingCheckpoint
    public let end: VivoReactiveSamplingCheckpoint
    public let initialAuthority: VivoNuclearPotentialEvaluation
    public let observations: [VivoReactiveSamplingObservation]
    public let authorityLabels: [VivoReactiveTrainingLabel]
    public let acquisitions: [VivoReactiveAcquisition]
    public let authorityEvaluations: Int
    public let baselineEvaluations: Int
    public let interpretation: String
}
public enum VivoReactiveSurrogateSampling {
    public static let interpretation = "State-independent positive-probability exact-energy random-walk refreshes retain access outside surrogate support; surrogate domain rejections are self-transitions. Frozen energy-gradient delta-surrogate HMC proposals with full authoritative endpoint Metropolis correction, invariant physical masses and fixed cell. Labels from rejected endpoints remain usable; rejected chain states remain in observations. Model fitting occurs only between separately equilibrated epochs, not during production. These are correlated equilibrium Markov samples, not physical MD time, first-passage evidence or recrossing trajectories. Descriptor-domain rejection may reduce accessibility; qualification and coverage are not implied by exact-energy correction. Floating-point reversibility and numerical-potential errors require convergence checks."
    public static func run(model: VivoReactiveSurrogateModel, authority: VivoNuclearPotential, baseline: VivoNuclearPotential,
                           checkpoint: VivoReactiveSamplingCheckpoint, sweeps: Int,
                           backend: VivoReactiveEnergyForceBackend? = nil) async throws -> VivoReactiveSamplingRun {
        try model.validate(); try checkpoint.validate(model: model)
        let cfg = checkpoint.configuration, n = checkpoint.positionsNM.count
        guard model.payload.authorityDefinition == authority.definition, model.payload.baselineDefinition == baseline.definition,
              (1...100_000).contains(sweeps), checkpoint.sweep <= UInt64.max-UInt64(sweeps),
              model.payload.qualification.passed,
              backend?.modelFingerprint == nil || backend!.modelFingerprint == model.fingerprint,
              checkpoint.backendProfile == (backend?.numericalProfile ?? "native-fp64-rbf/v1") else {
            throw VivoChemistryError.invalid("reactive production requires a qualified frozen model and exact authority/baseline binding")
        }
        guard sweeps+1 <= cfg.maximumAuthorityEvaluations,
              Double(sweeps)*Double(cfg.integrationSteps+1)+1 <= Double(cfg.maximumBaselineEvaluations),
              Double(sweeps+1)*Double(n)*18 <= Double(cfg.maximumStoredCoordinateElements) else {
            throw VivoChemistryError.resourceLimit("reactive sampling aggregate call or coordinate storage budget")
        }
        var authorities = 0, baselines = 0
        func evaluateAuthority(_ q: [VivoVector3D]) async throws -> VivoNuclearPotentialEvaluation {
            guard authorities < cfg.maximumAuthorityEvaluations else { throw VivoChemistryError.resourceLimit("reactive authority budget") }
            authorities += 1; return try await authority.checked(q)
        }
        func evaluateBaseline(_ q: [VivoVector3D]) async throws -> VivoNuclearPotentialEvaluation {
            guard baselines < cfg.maximumBaselineEvaluations else { throw VivoChemistryError.resourceLimit("reactive baseline budget") }
            baselines += 1; return try await baseline.checked(q)
        }
        func prediction(_ q: [VivoVector3D]) async throws -> VivoReactiveSurrogatePrediction {
            let result: VivoReactiveSurrogatePrediction
            if let backend {
                let batch = try await backend.predict([q])
                guard batch.count == 1 else { throw VivoChemistryError.invalid("reactive backend batch shape") }
                result = batch[0]
            } else { result = try VivoReactiveDeltaSurrogate.predict(model,positionsNM: q) }
            guard result.modelFingerprint == model.fingerprint, result.positionsNM == q,
                  result.deltaEnergyKJPerMol.isFinite, result.deltaForcesKJPerMolNM.count == n,
                  result.deltaForcesKJPerMolNM.allSatisfy(\.isFinite), result.eligible == result.reasons.isEmpty else {
                throw VivoChemistryError.invalid("reactive backend identity, geometry or finite derivative contract")
            }
            return result
        }
        let masses = authority.definition.massesDa, rt = VivoAtomicUnits.gasConstantJPerMolK*cfg.temperatureK/1000
        func kinetic(_ velocity: [VivoVector3D]) -> Double {
            zip(velocity,masses).reduce(0) { $0+0.5*$1.1*$1.0.squaredNorm }
        }
        var q = checkpoint.positionsNM, rng = checkpoint.randomState
        var currentAuthority = try await evaluateAuthority(q)
        let initialAuthority = currentAuthority
        var records: [VivoReactiveSamplingObservation] = [], labels: [VivoReactiveTrainingLabel] = [], acquisitions: [VivoReactiveAcquisition] = []
        for step in 1...sweeps {
            try Task.checkCancellation()
            let sweep = checkpoint.sweep+UInt64(step)
            if rng.unitInterval() < cfg.authorityRefreshProbability {
                let displaced = q.map { position in position+VivoVector3D(rng.normal(),rng.normal(),rng.normal())*cfg.authorityRefreshLengthNM }
                let candidate = try authority.definition.canonicalPositions(displaced)
                let target = try await evaluateAuthority(candidate), base = try await evaluateBaseline(candidate)
                let work = (target.energyKJPerMol-currentAuthority.energyKJPerMol)/rt
                guard work.isFinite else { throw VivoChemistryError.convergence("authority refresh energy overflow") }
                let logA = min(0,-work), logU = log(max(rng.unitInterval(),Double.leastNonzeroMagnitude))
                let accepted = logU < logA
                let approximate = try await prediction(candidate)
                let error = target.energyKJPerMol-base.energyKJPerMol-approximate.deltaEnergyKJPerMol
                labels.append(.init(identifier: "\(checkpoint.chainIdentifier):sweep-\(sweep)",sourceGroup: checkpoint.chainIdentifier,
                    authority: target,baseline: base))
                acquisitions.append(.init(positionsNM: candidate,reason: ["state-independent authority refresh; eligible for a later frozen training epoch"],modelFingerprint: model.fingerprint))
                if accepted { q = candidate; currentAuthority = target }
                records.append(.init(sweep: sweep,proposalKind: .authorityRandomWalk,positionsNM: q,authorityEnergyKJPerMol: currentAuthority.energyKJPerMol,
                    initialKineticEnergyKJPerMol: 0,proposedKineticEnergyKJPerMol: 0,logUniform: logU,accepted: accepted,
                    logAcceptanceProbability: logA,authorityMinusPredictedEnergyKJPerMol: error,rejectionReasons: []))
                continue
            }
            var candidate = q
            var velocity = masses.map { m -> VivoVector3D in
                let s = sqrt(rt/m); return .init(rng.normal()*s,rng.normal()*s,rng.normal()*s)
            }
            let initialKinetic = kinetic(velocity)
            let initialH = currentAuthority.energyKJPerMol+initialKinetic
            var approximate = try await prediction(candidate)
            var candidateBase: VivoNuclearPotentialEvaluation?
            var reasons: [String] = []
            if !approximate.eligible { reasons = approximate.reasons }
            if reasons.isEmpty {
                candidateBase = try await evaluateBaseline(candidate)
                // Refreshing every attempt makes proposals independent of any
                // prior rejected momentum. The final reversal defines an involution.
                for _ in 0..<cfg.integrationSteps {
                    let force = zip(candidateBase!.forcesKJPerMolNM,approximate.deltaForcesKJPerMolNM).map(+)
                    for i in 0..<n {
                        velocity[i] = velocity[i]+force[i]*(0.5*cfg.timeStepPS/masses[i])
                        candidate[i] = candidate[i]+velocity[i]*cfg.timeStepPS
                    }
                    candidate = try authority.definition.canonicalPositions(candidate)
                    approximate = try await prediction(candidate)
                    if !approximate.eligible { reasons = approximate.reasons; break }
                    candidateBase = try await evaluateBaseline(candidate)
                    let nextForce = zip(candidateBase!.forcesKJPerMolNM,approximate.deltaForcesKJPerMolNM).map(+)
                    for i in 0..<n { velocity[i] = velocity[i]+nextForce[i]*(0.5*cfg.timeStepPS/masses[i]) }
                }
            }
            var accepted = false, logA: Double?, error: Double?, logUniform: Double?, finalKinetic: Double?
            if reasons.isEmpty {
                let target = try await evaluateAuthority(candidate)
                velocity = velocity.map { $0 * -1 }
                finalKinetic = kinetic(velocity)
                let finalH = target.energyKJPerMol+finalKinetic!
                let work = (finalH-initialH)/rt
                guard work.isFinite else { throw VivoChemistryError.convergence("reactive proposal Hamiltonian overflow") }
                logA = min(0,-work)
                logUniform = log(max(rng.unitInterval(),Double.leastNonzeroMagnitude))
                accepted = logUniform! < logA!
                error = target.energyKJPerMol-candidateBase!.energyKJPerMol-approximate.deltaEnergyKJPerMol
                labels.append(.init(identifier: "\(checkpoint.chainIdentifier):sweep-\(sweep)",sourceGroup: checkpoint.chainIdentifier,
                    authority: target,baseline: candidateBase!))
                if abs(error!) > model.payload.configuration.maximumHeldOutEnergyErrorKJPerMol {
                    acquisitions.append(.init(positionsNM: candidate,reason: ["observed authority energy error exceeds training qualification tolerance"],modelFingerprint: model.fingerprint))
                }
                if accepted { q = candidate; currentAuthority = target }
            } else {
                // A domain exit is a self-transition, never authority-free MD.
                acquisitions.append(.init(positionsNM: candidate,reason: reasons,modelFingerprint: model.fingerprint))
            }
            records.append(.init(sweep: sweep,proposalKind: .surrogateHMC,positionsNM: q,authorityEnergyKJPerMol: currentAuthority.energyKJPerMol,
                initialKineticEnergyKJPerMol: initialKinetic,proposedKineticEnergyKJPerMol: finalKinetic,logUniform: logUniform,
                accepted: accepted,logAcceptanceProbability: logA,authorityMinusPredictedEnergyKJPerMol: error,rejectionReasons: reasons))
        }
        let end = try VivoReactiveSamplingCheckpoint.make(model: model,configuration: cfg,chainIdentifier: checkpoint.chainIdentifier,
            sweep: checkpoint.sweep+UInt64(sweeps),randomState: rng,positionsNM: q,backendProfile: checkpoint.backendProfile)
        return .init(start: checkpoint,end: end,initialAuthority: initialAuthority,observations: records,authorityLabels: labels,acquisitions: acquisitions,
            authorityEvaluations: authorities,baselineEvaluations: baselines,interpretation: interpretation)
    }
}
