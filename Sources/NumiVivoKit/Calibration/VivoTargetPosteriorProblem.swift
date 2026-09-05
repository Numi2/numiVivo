import Foundation

public enum VivoKineticFitField: String, Codable, Sendable {
    case association, dissociation, inactivation, targetTurnover
    case competitorAssociation, competitorDissociation, competitorConcentration, exposureScale
    public var unit: String {
        switch self {
        case .association, .competitorAssociation: return "M^-1 s^-1"
        case .dissociation, .inactivation, .targetTurnover, .competitorDissociation: return "1/s"
        case .competitorConcentration: return "M"
        case .exposureScale: return "1"
        }
    }
}

public struct VivoTargetPosteriorBinding: Codable, Equatable, Sendable {
    public let parameter: VivoPosteriorParameter
    public let field: VivoKineticFitField
    /// Explicit sharing. Disjoint sets can have separate parameters; overlapping
    /// writers for the same field are rejected. No context transfer is implicit.
    public let caseIdentifiers: [String]
    public let priorEvidence: VivoKineticEvidence
    public init(parameter: VivoPosteriorParameter, field: VivoKineticFitField,
                caseIdentifiers: [String], priorEvidence: VivoKineticEvidence) {
        self.parameter = parameter; self.field = field; self.caseIdentifiers = caseIdentifiers
        self.priorEvidence = priorEvidence
    }
}

public struct VivoTargetPosteriorProblem: Codable, Equatable, Sendable {
    public let schemaVersion: UInt32
    public let study: VivoTargetEngagementStudy
    public let bindings: [VivoTargetPosteriorBinding]
    public let sampler: VivoSMCConfiguration
    public let numerics: VivoTargetEngagementNumerics
    public let parallelEvaluations: Int
    public init(study: VivoTargetEngagementStudy, bindings: [VivoTargetPosteriorBinding],
                sampler: VivoSMCConfiguration = .init(), numerics: VivoTargetEngagementNumerics = .init(),
                parallelEvaluations: Int = 4) {
        schemaVersion = 1; self.study = study; self.bindings = bindings; self.sampler = sampler
        self.numerics = numerics; self.parallelEvaluations = parallelEvaluations
    }
    public func validate() throws {
        try study.validate(); try numerics.validate(); try sampler.validate(dimension: bindings.count)
        guard schemaVersion == 1, (1...32).contains(parallelEvaluations),
              Set(bindings.map { $0.parameter.identifier }).count == bindings.count,
              study.cases.contains(where: { $0.partition == .calibration }) else {
            throw VivoPosteriorError.invalid("inference schema, binding identity, parallelism or missing calibration data")
        }
        let cases = Dictionary(uniqueKeysWithValues: study.cases.map { ($0.identifier, $0) })
        var writers = Set<String>()
        for binding in bindings {
            try binding.parameter.validate(); try binding.priorEvidence.validate(origin: .assumed)
            guard binding.parameter.unit == binding.field.unit, binding.parameter.lower >= 0,
                  !binding.caseIdentifiers.isEmpty, binding.caseIdentifiers.count <= study.cases.count,
                  Set(binding.caseIdentifiers).count == binding.caseIdentifiers.count else {
                throw VivoPosteriorError.invalid("binding units, prior range or case set")
            }
            var contexts: [VivoKineticContext] = [], hasTrainingCase = false
            for id in binding.caseIdentifiers {
                guard let item = cases[id], writers.insert(id + "\u{0}" + binding.field.rawValue).inserted else {
                    throw VivoPosteriorError.invalid("missing case or overlapping parameter writers")
                }
                contexts.append(item.experiment.kinetics.context)
                hasTrainingCase = hasTrainingCase || item.partition == .calibration
                switch binding.field {
                case .competitorAssociation, .competitorDissociation, .competitorConcentration:
                    guard item.experiment.kinetics.competitor != nil else {
                        throw VivoPosteriorError.invalid("competitor parameter applied to a model without competitor")
                    }
                case .exposureScale:
                    for knot in item.experiment.exposure.knots {
                        let value = knot.unboundDrugM * binding.parameter.upper
                        guard value.isFinite, value <= item.experiment.kinetics.maximumUnboundDrugM else {
                            throw VivoPosteriorError.invalid("exposure prior exceeds declared kinetic applicability domain")
                        }
                    }
                default: break
                }
            }
            guard hasTrainingCase, let context = contexts.first, contexts.allSatisfy({ $0 == context }) else {
                throw VivoPosteriorError.invalid("parameter has no calibration use or crosses incompatible kinetic contexts")
            }
        }
        var observations = 0, sampleCount = 0, knots = 0
        for item in study.cases where item.partition == .calibration {
            observations += item.observations.count; sampleCount += item.experiment.sampleTimesSeconds.count
            knots += item.experiment.exposure.knots.count
            for observation in item.observations {
                guard let sd = observation.standardDeviation, sd.isFinite, sd > 0 else {
                    throw VivoPosteriorError.invalid("calibration requires explicit measurement SD; unknown is not zero")
                }
            }
        }
        // Prevent a valid small parameter file from expanding into unbounded
        // per-likelihood time series. Sampler evaluation counts are a separate cap.
        guard observations <= 4096, sampleCount <= 8192, knots <= 8192 else {
            throw VivoPosteriorError.budget("per-likelihood calibration observation or boundary capacity")
        }
    }
}

/// Validated immutable study plus training-only likelihood identity. Independent
/// Gaussian measurement errors with supplied SD are explicit modelling assumptions.
/// No assay noise is inferred or silently assigned; fitting never reads held-out values.
public struct VivoPreparedTargetPosterior: Sendable {
    public let problem: VivoTargetPosteriorProblem
    public let plan: VivoPosteriorPlan
    public let calibrationCases: [VivoTargetEngagementStudyCase]

    private struct TrainingBinding: Encodable {
        let parameter: VivoPosteriorParameter
        let field: VivoKineticFitField
        let calibrationCaseIdentifiers: [String]
        let priorEvidence: VivoKineticEvidence
    }
    private struct TrainingIdentity: Encodable {
        let method: String
        let cases: [VivoTargetEngagementStudyCase]
        let bindings: [TrainingBinding]
        let numerics: VivoTargetEngagementNumerics
    }

    public init(_ problem: VivoTargetPosteriorProblem) throws {
        try problem.validate()
        self.problem = problem
        calibrationCases = problem.study.cases.filter { $0.partition == .calibration }.sorted { $0.identifier < $1.identifier }
        let ids = Set(calibrationCases.map(\.identifier))
        let bindings = problem.bindings.map {
            TrainingBinding(parameter: $0.parameter, field: $0.field,
                calibrationCaseIdentifiers: $0.caseIdentifiers.filter(ids.contains).sorted(), priorEvidence: $0.priorEvidence)
        }
        let identity = TrainingIdentity(method: "numivivo.occupancy-independent-gaussian-fp64.v1",
            cases: calibrationCases, bindings: bindings, numerics: problem.numerics)
        plan = try .init(likelihoodFingerprint: VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(identity)),
                         parameters: problem.bindings.map(\.parameter), configuration: problem.sampler)
    }

    public func physicalValues(_ particle: VivoPosteriorParticle) throws -> [Double] {
        guard particle.coordinates.count == plan.parameters.count else { throw VivoPosteriorError.invalid("particle dimensions") }
        return try zip(plan.parameters, particle.coordinates).map { try $0.value(at: $1) }
    }

    /// Applies shared parameters without altering the original experiment or
    /// claiming a proposal is a measured rate. The run record retains original
    /// evidence, priors and fitted joint particle coordinates.
    public func experiment(for item: VivoTargetEngagementStudyCase, values: [Double]) throws -> VivoTargetEngagementExperiment {
        guard values.count == problem.bindings.count,
              problem.study.cases.contains(where: { $0 == item }) else {
            throw VivoPosteriorError.invalid("candidate parameter layout or foreign study case")
        }
        let base = item.experiment.kinetics
        var association = base.association, dissociation = base.dissociation, inactivation = base.inactivation
        var turnover = base.targetTurnover, competitor = base.competitor, exposureScale = 1.0
        let evidence = VivoKineticEvidence(source: "NumiVivo posterior proposal",
            locator: "parameter proposal, not a measured or independently validated estimate",
            sourceFingerprint: plan.likelihoodFingerprint.hex)
        func proposed(_ value: Double, unit: VivoKineticUnit) -> VivoKineticParameter {
            .init(value: value, unit: unit, origin: .assumed, uncertainty: .unknown, evidence: evidence)
        }
        for (index, binding) in problem.bindings.enumerated() {
            let value = values[index]
            guard value.isFinite, value >= binding.parameter.lower, value <= binding.parameter.upper else {
                throw VivoPosteriorError.invalid("candidate outside prior support")
            }
            if !binding.caseIdentifiers.contains(item.identifier) { continue }
            switch binding.field {
            case .association: association = proposed(value, unit: .perMolarSecond)
            case .dissociation: dissociation = proposed(value, unit: .perSecond)
            case .inactivation: inactivation = proposed(value, unit: .perSecond)
            case .targetTurnover: turnover = proposed(value, unit: .perSecond)
            case .exposureScale: exposureScale = value
            case .competitorAssociation, .competitorDissociation, .competitorConcentration:
                guard let old = competitor else { throw VivoPosteriorError.invalid("missing competitor") }
                competitor = .init(identifier: old.identifier,
                    association: binding.field == .competitorAssociation ? proposed(value, unit: .perMolarSecond) : old.association,
                    dissociation: binding.field == .competitorDissociation ? proposed(value, unit: .perSecond) : old.dissociation,
                    unboundConcentration: binding.field == .competitorConcentration ? proposed(value, unit: .molar) : old.unboundConcentration)
            }
        }
        let model = VivoCovalentKineticPack(identifier: base.identifier, context: base.context,
            association: association, dissociation: dissociation, inactivation: inactivation,
            targetTurnover: turnover, baselineTarget: base.baselineTarget, competitor: competitor,
            maximumUnboundDrugM: base.maximumUnboundDrugM)
        let knots = try item.experiment.exposure.knots.map { knot -> VivoExposureKnot in
            let concentration = knot.unboundDrugM * exposureScale
            guard concentration.isFinite,
                  concentration != 0 || knot.unboundDrugM == 0 || exposureScale == 0 else {
                throw VivoPosteriorError.numerical("scaled exposure overflow or underflow")
            }
            return .init(timeSeconds: knot.timeSeconds, unboundDrugM: concentration)
        }
        let trace = VivoUnboundExposureTrace(context: base.context, knots: knots, origin: .assumed, evidence: evidence)
        let experiment = VivoTargetEngagementExperiment(kinetics: model, exposure: trace,
            initial: item.experiment.initial, sampleTimesSeconds: item.experiment.sampleTimesSeconds)
        try experiment.validate()
        return experiment
    }

    public func logLikelihood(values: [Double]) throws -> Double {
        var sum = 0.0, compensation = 0.0
        for item in calibrationCases {
            try Task.checkCancellation()
            let modified = try experiment(for: item, values: values)
            let result = try VivoTargetEngagementReference.run(modified, numerics: problem.numerics)
            let samples = Dictionary(uniqueKeysWithValues: result.samples.map { ($0.timeSeconds, $0) })
            for observation in item.observations {
                guard let sample = samples[observation.timeSeconds], let sd = observation.standardDeviation else {
                    throw VivoPosteriorError.invalid("missing calibration sample or measurement SD")
                }
                let residual = (observation.observable.value(sample) - observation.value) / sd
                let contribution = -0.5 * residual * residual - log(sd) - 0.5 * log(2 * Double.pi)
                guard contribution.isFinite, abs(contribution) <= 1e250 else {
                    throw VivoPosteriorError.numerical("Gaussian likelihood overflow; no finite penalty substitution")
                }
                let adjusted = contribution - compensation, next = sum + adjusted
                compensation = (next - sum) - adjusted; sum = next
            }
        }
        guard sum.isFinite, abs(sum) <= 1e250 else { throw VivoPosteriorError.numerical("likelihood accumulation overflow") }
        return sum
    }

    public func evaluate(_ candidates: [VivoPosteriorCandidate]) async throws -> [VivoPosteriorEvaluation] {
        guard candidates.count <= plan.configuration.particleCount else { throw VivoPosteriorError.budget("likelihood batch size") }
        return try await withThrowingTaskGroup(of: VivoPosteriorEvaluation.self) { group in
            func add(_ candidate: VivoPosteriorCandidate) {
                group.addTask {
                    do {
                        try Task.checkCancellation()
                        let reconstructed = try VivoPosteriorCandidate(ordinal: candidate.ordinal,
                            coordinates: candidate.coordinates, parameters: plan.parameters)
                        guard reconstructed == candidate else { throw VivoPosteriorError.invalid("candidate value binding") }
                        return try .init(candidate: candidate, logLikelihood: logLikelihood(values: candidate.values))
                    } catch is CancellationError { throw CancellationError() }
                    catch { return .init(candidate: candidate, failure: String(error.localizedDescription.prefix(4096))) }
                }
            }
            var next = 0, results: [VivoPosteriorEvaluation] = []
            for _ in 0..<min(problem.parallelEvaluations, candidates.count) { add(candidates[next]); next += 1 }
            while let result = try await group.next() {
                results.append(result)
                if next < candidates.count { add(candidates[next]); next += 1 }
            }
            return results.sorted { $0.candidate.ordinal < $1.candidate.ordinal }
        }
    }
}

public struct VivoTargetPosteriorRecord: Codable, Equatable, Sendable {
    public let schemaVersion: UInt32
    public let problem: VivoTargetPosteriorProblem
    public let posterior: VivoPosteriorRun
    public func validate(requireComplete: Bool = false) throws {
        let prepared = try VivoPreparedTargetPosterior(problem)
        guard schemaVersion == 1, posterior.plan == prepared.plan else {
            throw VivoPosteriorError.invalid("inference result does not bind its declared calibration problem")
        }
        try posterior.validate(requireComplete: requireComplete)
    }
}

public enum VivoTargetPosteriorFitter {
    public static func run(_ problem: VivoTargetPosteriorProblem, checkpoint: VivoPosteriorCheckpoint? = nil,
                           progress: VivoPosteriorCheckpointSink? = nil) async throws -> VivoTargetPosteriorRecord {
        let prepared = try VivoPreparedTargetPosterior(problem)
        let sampler = try VivoTemperedPosteriorSampler(plan: prepared.plan, checkpoint: checkpoint)
        let result = try await sampler.run(evaluate: { try await prepared.evaluate($0) }, progress: progress)
        return .init(schemaVersion: 1, problem: problem, posterior: result)
    }
}
