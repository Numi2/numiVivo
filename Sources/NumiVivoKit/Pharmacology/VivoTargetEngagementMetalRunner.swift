import Foundation
import Metal

public struct VivoTargetEngagementMetalPolicy: Codable, Equatable, Sendable {
    public let maximumTimeStepSeconds: Double
    public let maximumSteps: Int
    public let fractionBalanceTolerance: Double
    public init(maximumTimeStepSeconds: Double = 0.1, maximumSteps: Int = 1_000_000,
                fractionBalanceTolerance: Double = 1e-4) {
        self.maximumTimeStepSeconds = maximumTimeStepSeconds; self.maximumSteps = maximumSteps
        self.fractionBalanceTolerance = fractionBalanceTolerance
    }
}

public struct VivoTargetEngagementMetalBatchResult: Codable, Sendable {
    public let executionFingerprint: VivoFingerprint
    public let programFingerprint: VivoFingerprint
    public let deviceName: String
    public let experiments: [VivoTargetEngagementExperiment]
    public let policy: VivoTargetEngagementMetalPolicy
    public let timeStepSeconds: Double
    public let committedSteps: Int
    public let results: [VivoTargetEngagementResult]
}

/// Executes common-topology exposure cohorts in the EXISTING transactional F1
/// Metal runtime. No parallel physiology engine or alternate reaction kernels.
/// The RK2 step bound limits positivity/stability risk, not biological or global
/// truncation error; qualify accuracy against the FP64 reference and step refinement.
public enum VivoTargetEngagementMetalRunner {
    private struct Identity: Encodable {
        let version: UInt32
        let experiments: [VivoTargetEngagementExperiment]
        let policy: VivoTargetEngagementMetalPolicy
        let configuration: VivoRuntimeConfiguration
        let programFingerprint: VivoFingerprint
    }

    public static func run(_ experiments: [VivoTargetEngagementExperiment],
                           policy: VivoTargetEngagementMetalPolicy = .init(),
                           device: MTLDevice? = nil) async throws -> VivoTargetEngagementMetalBatchResult {
        guard (1...4096).contains(experiments.count), let first = experiments.first,
              policy.maximumTimeStepSeconds.isFinite, policy.maximumTimeStepSeconds > 0,
              (1...1_000_000).contains(policy.maximumSteps),
              policy.fractionBalanceTolerance.isFinite, policy.fractionBalanceTolerance >= 1e-7,
              policy.fractionBalanceTolerance <= 1e-4 else { throw VivoKineticsError.invalid("Metal batch/policy") }
        let sampleCount = first.sampleTimesSeconds.count.multipliedReportingOverflow(by: experiments.count)
        guard !sampleCount.overflow, sampleCount.partialValue <= 262_144 else {
            throw VivoKineticsError.capacity("Metal cohort output sample capacity")
        }
        var maximumScheduledDrugM = 0.0
        var totalKnots = 0
        var allBoundaries = Set(first.sampleTimesSeconds)
        for experiment in experiments {
            try experiment.validate()
            guard experiment.kinetics == first.kinetics, experiment.initial == first.initial,
                  experiment.sampleTimesSeconds == first.sampleTimesSeconds else {
                throw VivoKineticsError.unsupported("one Metal cohort requires identical kinetics, initial state and observation times")
            }
            totalKnots += experiment.exposure.knots.count
            guard totalKnots <= 262_144 else { throw VivoKineticsError.capacity("Metal exposure table") }
            for knot in experiment.exposure.knots {
                guard Double(Float(knot.timeSeconds)) == knot.timeSeconds,
                      knot.unboundDrugM <= Double(Float.greatestFiniteMagnitude),
                      knot.unboundDrugM == 0 || Float(knot.unboundDrugM) != 0 else {
                    throw VivoKineticsError.unsupported("Metal cohort requires FP32-exact time boundaries and representable exposure")
                }
                maximumScheduledDrugM = max(maximumScheduledDrugM, knot.unboundDrugM)
                allBoundaries.insert(knot.timeSeconds)
            }
        }
        guard first.sampleTimesSeconds.allSatisfy({ Double(Float($0)) == $0 }),
              let endTime = first.sampleTimesSeconds.last else { throw VivoKineticsError.invalid("Metal sample times") }
        let boundaries = allBoundaries.filter { $0 > 0 && $0 <= endTime }.sorted()
        let m = first.kinetics, d = m.targetTurnover.value
        let competitorRate = (m.competitor?.association.value ?? 0) * (m.competitor?.unboundConcentration.value ?? 0)
        // The full immutable schedule is validated above. Use its actual maximum,
        // not the looser applicability-domain ceiling, to size a common step.
        let maxRate = max(m.association.value * maximumScheduledDrugM + competitorRate + d,
                          m.dissociation.value + m.inactivation.value + d, (m.competitor?.dissociation.value ?? 0) + d)
        let limit = min(policy.maximumTimeStepSeconds, maxRate > 0 ? 0.1 / maxRate : policy.maximumTimeStepSeconds)
        guard limit.isFinite, limit > 0 else { throw VivoKineticsError.numerical("Metal step bound underflow") }
        let step = pow(2, floor(log2(limit)))
        guard step >= Double(Float.leastNormalMagnitude), step <= Double(Float.greatestFiniteMagnitude),
              endTime / step + Double(boundaries.count) <= Double(policy.maximumSteps) else {
            throw VivoKineticsError.capacity("explicit F1 step budget; use the FP64 reference for stiff/small systems")
        }
        let compiled = try VivoTargetEngagementCompiler.compile(first)
        guard let fidelity = VivoFidelity(rawValue: 1) else { throw VivoKineticsError.unsupported("F1 unavailable") }
        let configuration = VivoRuntimeConfiguration(fidelity: fidelity, environmentCount: UInt32(experiments.count),
            timeStep: Float(step), minimumTimeStep: Float.leastNormalMagnitude, maximumTimeStep: Float(step),
            maximumSubsteps: 1, eventCapacity: UInt32(max(1024, experiments.count * 4)))
        let runtime = try await VivoTransactionalMolecularRuntime.make(pack: compiled.pack, configuration: configuration, device: device)
        let identity = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(Identity(version: 1,
            experiments: experiments, policy: policy, configuration: configuration,
            programFingerprint: compiled.pack.header.contentFingerprint)))
        let requests = experiments.indices.flatMap { lane in
            compiled.targetSpeciesIndices.map { VivoPublicationRequest(speciesIndex: $0, laneIndex: UInt32(lane)) }
        }
        let observationTimes = Set(first.sampleTimesSeconds)
        var samples = [[VivoTargetEngagementSample]](repeating: [], count: experiments.count)
        var cursors = [Int](repeating: 0, count: experiments.count)
        var maximumErrors = [Double](repeating: abs(first.initial.total - 1), count: experiments.count)
        var time = 0.0, steps = 0
        if observationTimes.contains(0) {
            for lane in experiments.indices {
                samples[lane].append(.init(timeSeconds: 0, unboundDrugM: experiments[lane].exposure.knots[0].unboundDrugM,
                    fractions: first.initial, totalTargetM: first.initial.total * m.baselineTarget.value))
            }
        }
        for boundary in boundaries {
            while time < boundary {
                try Task.checkCancellation()
                guard steps < policy.maximumSteps else { throw VivoKineticsError.capacity("Metal step count") }
                let dt = min(step, boundary - time)
                guard Double(Float(dt)) == dt, dt >= Double(Float.leastNormalMagnitude), time + dt > time else {
                    throw VivoKineticsError.unsupported("boundary cannot be reached with the current FP32 step ABI")
                }
                let updates = experiments.indices.map { lane -> VivoCouplingUpdate in
                    let knots = experiments[lane].exposure.knots
                    while cursors[lane] + 1 < knots.count, knots[cursors[lane] + 1].timeSeconds <= time { cursors[lane] += 1 }
                    return .init(speciesIndex: compiled.drugSpeciesIndex, laneIndex: UInt32(lane),
                                 mode: .replace, value: Float(knots[cursors[lane]].unboundDrugM))
                }
                let publish = time + dt == boundary && observationTimes.contains(boundary)
                let candidate = try await runtime.prepareStep(.init(timeStep: Float(dt), coupling: updates,
                    publications: publish ? requests : [], permitAdaptiveReduction: false))
                guard candidate.canCommit, Double(candidate.candidateTimeStep) == dt else {
                    if candidate.canCommit { try await runtime.discardPreparedStep(transactionID: candidate.transactionID) }
                    throw VivoKineticsError.numerical("Metal candidate rejected; no failed run is discarded from an inferred ensemble")
                }
                if publish {
                    guard candidate.publications.count == experiments.count * 4 else {
                        try await runtime.discardPreparedStep(transactionID: candidate.transactionID)
                        throw VivoKineticsError.numerical("Metal publication shape")
                    }
                    for lane in experiments.indices {
                        let base = lane * 4
                        let f = VivoTargetFractions(free: Double(candidate.publications[base]),
                            reversible: Double(candidate.publications[base + 1]), covalent: Double(candidate.publications[base + 2]),
                            competitor: Double(candidate.publications[base + 3]))
                        do { try f.validate(hasCompetitor: m.competitor != nil, tolerance: policy.fractionBalanceTolerance) }
                        catch {
                            try await runtime.discardPreparedStep(transactionID: candidate.transactionID)
                            throw error
                        }
                        maximumErrors[lane] = max(maximumErrors[lane], abs(f.total - 1))
                    }
                }
                _ = try await runtime.commitPreparedStep(transactionID: candidate.transactionID)
                let after = await runtime.timeSeconds()
                guard after == time + dt else { throw VivoKineticsError.numerical("Metal accepted clock mismatch") }
                time = after; steps += 1
                if publish {
                    for lane in experiments.indices {
                        let knots = experiments[lane].exposure.knots
                        while cursors[lane] + 1 < knots.count, knots[cursors[lane] + 1].timeSeconds <= time { cursors[lane] += 1 }
                        let base = lane * 4
                        let f = VivoTargetFractions(free: Double(candidate.publications[base]),
                            reversible: Double(candidate.publications[base + 1]), covalent: Double(candidate.publications[base + 2]),
                            competitor: Double(candidate.publications[base + 3]))
                        samples[lane].append(.init(timeSeconds: time, unboundDrugM: knots[cursors[lane]].unboundDrugM,
                                                  fractions: f, totalTargetM: f.total * m.baselineTarget.value))
                    }
                }
            }
        }
        let results = experiments.indices.map { lane in
            VivoTargetEngagementResult(schemaVersion: 1, backend: "numivivo.programpack.f1-metal.v1", samples: samples[lane],
                propagationCount: steps, maximumScalingDepth: 0, maximumFractionMassError: maximumErrors[lane],
                containsAssumedParameters: m.containsAssumptions, hasUnknownParameterUncertainty: m.hasUnknownParameterUncertainty,
                limitations: ["FP32 F1 conditional simulation; global accuracy requires step refinement and reference comparison.",
                    "Maximum reported balance error is over published samples; internal bounds are enforced by ProgramPack monitors.",
                    "Externally maintained free drug; finite ligand depletion and downstream biological efficacy are not represented.",
                    "This result is not a predictive uncertainty interval or clinical recommendation."])
        }
        return .init(executionFingerprint: identity, programFingerprint: compiled.pack.header.contentFingerprint,
                     deviceName: runtime.capabilities.deviceName, experiments: experiments, policy: policy,
                     timeStepSeconds: step, committedSteps: steps, results: results)
    }
}
