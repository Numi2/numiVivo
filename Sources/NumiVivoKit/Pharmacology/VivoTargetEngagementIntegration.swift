import Foundation

public struct PreparedVivoTargetEngagement: Sendable {
    public let experiment: VivoTargetEngagementExperiment
    public let experimentFingerprint: VivoFingerprint
    public let kineticsFingerprint: VivoFingerprint
    public let source: Data
    public let pack: VivoProgramPack
    public let drugSpeciesIndex: UInt32
    public let targetSpeciesIndices: [UInt32]
}

public enum VivoTargetEngagementCompiler {
    public static func compile(_ experiment: VivoTargetEngagementExperiment) throws -> PreparedVivoTargetEngagement {
        try experiment.validate()
        let experimentID = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(experiment))
        let kineticsID = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(experiment.kinetics))
        let source = try VivoTargetEngagementProgramSource.make(experiment,
            kineticsFingerprint: kineticsID.hex, experimentFingerprint: experimentID.hex)
        guard let fidelity = VivoFidelity(rawValue: 1) else { throw VivoKineticsError.unsupported("F1 is unavailable") }
        let compilation = try VivoNativeCompilerBridge.compileProgram(source, fidelity: fidelity,
            strictUnits: true, strictSafety: true, requireTermination: true)
        guard compilation.succeeded else { throw VivoKineticsError.invalid(compilation.diagnosticsString()) }
        let pack = try VivoProgramPack(data: compilation.primary)
        let metadata = try pack.speciesMetadata()
        func index(_ id: String) throws -> UInt32 {
            guard let index = metadata.firstIndex(where: { $0.identifier == id }), index <= Int(UInt32.max) else {
                throw VivoKineticsError.invalid("compiled species is missing: \(id)")
            }
            return UInt32(index)
        }
        let drug = try index(VivoTargetEngagementProgramSource.drugInput)
        guard metadata[Int(drug)].isExternallyOwned else { throw VivoKineticsError.invalid("drug reservoir ownership was lost") }
        return .init(experiment: experiment, experimentFingerprint: experimentID, kineticsFingerprint: kineticsID,
                     source: source, pack: pack, drugSpeciesIndex: drug,
                     targetSpeciesIndices: try VivoTargetEngagementProgramSource.targetSpecies.map(index))
    }

    /// Uses the existing production physiology/reaction coordinator. The caller
    /// must identify an analyte that ALREADY represents unbound concentration.
    /// For total concentration use an explicitly qualified partition model first.
    public static func physiologyBridge(_ compiled: PreparedVivoTargetEngagement,
                                         physiology: PreparedVivoPhysiologyModel,
                                         sourceContext: VivoKineticContext,
                                         unboundAnalyte: String, compartment: String,
                                         sourceConcentrationIsUnbound: Bool) throws -> PreparedVivoMolecularPhysiologyBridge {
        guard sourceContext == compiled.experiment.kinetics.context, sourceConcentrationIsUnbound else {
            throw VivoKineticsError.unsupported("physiology context or unbound-concentration meaning is unresolved")
        }
        let coupling = VivoMolecularPhysiologyCouplingPack(
            programFingerprint: compiled.pack.header.contentFingerprint,
            physiologyFingerprint: physiology.fingerprint,
            physiologyToMolecular: [.init(identifier: "unbound-target-exposure", analyte: unboundAnalyte,
                compartment: compartment, molecularSpecies: VivoTargetEngagementProgramSource.drugInput,
                mapping: .oneToOne, molecularMode: .replace, transfer: .init(gain: 1), required: true)],
            labels: ["kinetics-sha256": compiled.kineticsFingerprint.hex,
                     "exposure-semantics": "externally-maintained-unbound-concentration",
                     "drug-depletion": "not-modelled"])
        return try VivoMolecularPhysiologyCouplingCompiler().compile(coupling,
            programPack: compiled.pack, physiology: physiology, molecularLaneCount: physiology.environmentCount)
    }
}

public enum VivoPhysiologyConcentrationMeaning: Codable, Equatable, Sendable {
    case alreadyUnbound
    /// A constant-fraction approximation. Saturable/time-varying protein binding
    /// must be computed upstream rather than silently approximated by this case.
    case totalWithConstantUnboundFraction(Double)
    public func factor() throws -> Double {
        switch self {
        case .alreadyUnbound: return 1
        case .totalWithConstantUnboundFraction(let value):
            guard value.isFinite, value >= 0, value <= 1 else { throw VivoKineticsError.invalid("unbound fraction outside [0,1]") }
            return value
        }
    }
}

public struct VivoPhysiologyExposureProjection: Sendable {
    public let trace: VivoUnboundExposureTrace
    /// Persist these bytes in the shared artifact store: trace evidence identifies
    /// this exact model/snapshot/projection document, not a fabricated source hash.
    public let evidenceData: Data
    public let evidenceFingerprint: VivoFingerprint
}

public enum VivoPhysiologyExposureAdapter {
    private struct Source: Encodable {
        let schemaVersion: UInt32
        let model: PreparedVivoPhysiologyModel
        let snapshots: [VivoPhysiologySnapshot]
        let context: VivoKineticContext
        let analyte: String
        let compartment: String
        let environment: UInt32
        let meaning: VivoPhysiologyConcentrationMeaning
        let interpretationEvidence: VivoKineticEvidence
        let interpretationOrigin: VivoKineticOrigin
    }
    public static func project(model: PreparedVivoPhysiologyModel, snapshots: [VivoPhysiologySnapshot],
                               kinetics: VivoCovalentKineticPack, sourceContext: VivoKineticContext,
                               analyte: String, compartment: String, environment: UInt32,
                               meaning: VivoPhysiologyConcentrationMeaning,
                               interpretationEvidence: VivoKineticEvidence,
                               interpretationOrigin: VivoKineticOrigin) throws -> VivoPhysiologyExposureProjection {
        try VivoPhysiologyModelValidator.validate(model); try kinetics.validate()
        try interpretationEvidence.validate(origin: interpretationOrigin)
        guard sourceContext == kinetics.context, environment < model.environmentCount,
              (2...65_536).contains(snapshots.count), snapshots.first?.absoluteTimeSeconds == 0,
              let analyteMetadata = model.analytes.first(where: { $0.identifier == analyte }) else {
            throw VivoKineticsError.invalid("projection context, analyte, environment or snapshot support")
        }
        let pair = try model.pairIndex(analyte: analyte, compartment: compartment)
        let elementCount = UInt64(model.pairCount) * UInt64(model.environmentCount)
        guard elementCount <= UInt64(Int.max) else { throw VivoKineticsError.capacity("snapshot shape") }
        let factor = try meaning.factor()
        var previous = -Double.infinity
        var knots: [VivoExposureKnot] = []
        for snapshot in snapshots {
            guard snapshot.modelFingerprint == model.fingerprint, snapshot.pairCount == model.pairCount,
                  snapshot.environmentCount == model.environmentCount, snapshot.values.count == Int(elementCount),
                  snapshot.absoluteTimeSeconds.isFinite, snapshot.absoluteTimeSeconds > previous else {
                throw VivoKineticsError.invalid("physiology snapshot identity, shape or clock mismatch")
            }
            let index = Int(pair) * Int(model.environmentCount) + Int(environment)
            let raw = Double(snapshot.values[index])
            guard raw.isFinite, raw >= 0 else { throw VivoKineticsError.invalid("negative/nonfinite exposure publication") }
            let molar = try VivoUnitSystem.standard.convert(raw, from: analyteMetadata.concentrationUnit, to: "M") * factor
            knots.append(.init(timeSeconds: snapshot.absoluteTimeSeconds, unboundDrugM: molar))
            previous = snapshot.absoluteTimeSeconds
        }
        let data = try VivoCanonicalJSON.encode(Source(schemaVersion: 1, model: model, snapshots: snapshots,
            context: sourceContext, analyte: analyte, compartment: compartment, environment: environment,
            meaning: meaning, interpretationEvidence: interpretationEvidence, interpretationOrigin: interpretationOrigin))
        let identity = try VivoCanonicalJSON.fingerprint(data)
        let trace = VivoUnboundExposureTrace(context: sourceContext, knots: knots, origin: .calculated,
            evidence: .init(source: "NumiVivo physiology exposure projection", locator: "snapshots + concentration meaning",
                            sourceFingerprint: identity.hex))
        try trace.validate(for: kinetics)
        return .init(trace: trace, evidenceData: data, evidenceFingerprint: identity)
    }
}

/// Self-contained input/output provenance. The fingerprints detect changed
/// content; they do not authenticate a researcher or validate biological claims.
public struct VivoTargetEngagementRunRecord: Codable, Equatable, Sendable {
    public let schemaVersion: UInt32
    public let experiment: VivoTargetEngagementExperiment
    public let numerics: VivoTargetEngagementNumerics
    public let experimentFingerprint: VivoFingerprint
    public let executionFingerprint: VivoFingerprint
    public let resultFingerprint: VivoFingerprint
    public let result: VivoTargetEngagementResult

    private struct Execution: Encodable {
        let numericalVersion: UInt32
        let experiment: VivoTargetEngagementExperiment
        let numerics: VivoTargetEngagementNumerics
    }
    public static func reference(_ experiment: VivoTargetEngagementExperiment,
                                 numerics: VivoTargetEngagementNumerics = .init()) throws -> Self {
        let result = try VivoTargetEngagementReference.run(experiment, numerics: numerics)
        return try .init(schemaVersion: 1, experiment: experiment, numerics: numerics,
            experimentFingerprint: VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(experiment)),
            executionFingerprint: VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(
                Execution(numericalVersion: 1, experiment: experiment, numerics: numerics))),
            resultFingerprint: VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(result)), result: result)
    }
    public func validate() throws {
        try experiment.validate(); try numerics.validate()
        guard schemaVersion == 1, result.schemaVersion == 1,
              result.backend == "numivivo.target-engagement.fp64-uniformization.v1",
              result.containsAssumedParameters == experiment.kinetics.containsAssumptions,
              result.hasUnknownParameterUncertainty == experiment.kinetics.hasUnknownParameterUncertainty,
              result.propagationCount >= 0, result.maximumScalingDepth >= 0,
              result.maximumFractionMassError >= 0,
              experimentFingerprint == (try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(experiment))),
              executionFingerprint == (try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(
                Execution(numericalVersion: 1, experiment: experiment, numerics: numerics)))),
              resultFingerprint == (try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(result))),
              result.samples.count == experiment.sampleTimesSeconds.count,
              result.propagationCount <= numerics.maximumPropagations,
              result.maximumScalingDepth <= numerics.maximumMatrixSquarings,
              result.maximumFractionMassError.isFinite,
              result.maximumFractionMassError <= numerics.fractionConservationTolerance else {
            throw VivoKineticsError.invalid("run record identity, shape or numerical policy mismatch")
        }
        var exposureIndex = 0
        for (sample, time) in zip(result.samples, experiment.sampleTimesSeconds) {
            while exposureIndex + 1 < experiment.exposure.knots.count,
                  experiment.exposure.knots[exposureIndex + 1].timeSeconds <= time { exposureIndex += 1 }
            guard sample.timeSeconds == time, sample.totalTargetM.isFinite, sample.totalTargetM > 0,
                  sample.unboundDrugM.isFinite, sample.unboundDrugM >= 0,
                  sample.unboundDrugM == experiment.exposure.knots[exposureIndex].unboundDrugM,
                  sample.totalTargetM == sample.fractions.total * experiment.kinetics.baselineTarget.value else {
                throw VivoKineticsError.invalid("run record sample support")
            }
            try sample.fractions.validate(hasCompetitor: experiment.kinetics.competitor != nil,
                                          tolerance: numerics.fractionConservationTolerance)
        }
    }
}
