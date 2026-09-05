import Foundation

public struct VivoMetalTargetScreenConfiguration: Codable, Equatable, Sendable {
    public let maximumTimeStepSeconds: Double
    public let maximumStepsPerCase: Int
    public let maximumTotalPlannedSteps: Int
    public let operationsPerCommand: Int
    public let maximumRetainedBufferBytes: Int
    public let maximumScreeningEvaluations: Int
    public init(maximumTimeStepSeconds: Double = 0.25, maximumStepsPerCase: Int = 262_144,
                maximumTotalPlannedSteps: Int = 1_048_576, operationsPerCommand: Int = 128,
                maximumRetainedBufferBytes: Int = 268_435_456, maximumScreeningEvaluations: Int = 2_000_000) {
        self.maximumTimeStepSeconds = maximumTimeStepSeconds; self.maximumStepsPerCase = maximumStepsPerCase
        self.maximumTotalPlannedSteps = maximumTotalPlannedSteps; self.operationsPerCommand = operationsPerCommand
        self.maximumRetainedBufferBytes = maximumRetainedBufferBytes
        self.maximumScreeningEvaluations = maximumScreeningEvaluations
    }
    public func validate() throws {
        guard maximumTimeStepSeconds.isFinite, maximumTimeStepSeconds > 0,
              (1...1_048_576).contains(maximumStepsPerCase),
              (1...4_194_304).contains(maximumTotalPlannedSteps),
              (1...256).contains(operationsPerCommand),
              (1_048_576...1_073_741_824).contains(maximumRetainedBufferBytes),
              (1...20_000_000).contains(maximumScreeningEvaluations) else {
            throw VivoPosteriorError.invalid("Metal screening policy or capacity")
        }
    }
}

public struct VivoMetalTargetCaseDescription: Codable, Equatable, Sendable {
    public let caseIdentifier: String
    public let programFingerprint: VivoFingerprint
    public let scheduleFingerprint: VivoFingerprint
    public let timeStepSeconds: Double
    public let steps: Int
    public let operations: Int
    public let capacity: Int
}

internal enum VivoMetalTargetOperation: Codable, Sendable {
    case exposure(Float)
    case step(dt: Float, after: Float)
    case observe(Int)
}

/// Finite FP32 inputs must survive normal-range execution. Subnormal values can
/// be flushed by GPU arithmetic and must not silently become zero parameters.
internal func vivoScreenFloat(_ value: Double, label: String) throws -> Float {
    guard value.isFinite, abs(value) <= Double(Float.greatestFiniteMagnitude),
          value == 0 || abs(value) >= Double(Float.leastNormalMagnitude) else {
        throw VivoPosteriorError.invalid("\(label) is outside the supported FP32 normal range")
    }
    return Float(value)
}

internal struct VivoMetalTargetCasePlan: Sendable {
    let item: VivoTargetEngagementStudyCase
    let pack: VivoProgramPack
    let configuration: VivoRuntimeConfiguration
    let description: VivoMetalTargetCaseDescription
    let operations: [VivoMetalTargetOperation]
    let drugIndex: UInt32
    let competitorIndex: UInt32?
    let targetIndices: SIMD4<UInt32>
    let parameterNames: [String]

    static func make(_ item: VivoTargetEngagementStudyCase, prepared: VivoPreparedTargetPosterior,
                     policy: VivoMetalTargetScreenConfiguration) throws -> Self {
        try policy.validate()
        let problem = prepared.problem, model = item.experiment.kinetics
        guard item.partition == .calibration, item.observations.allSatisfy({
            problem.assaySupport(caseIdentifier: item.identifier, observationIdentifier: $0.identifier) == .exact
        }) else { throw VivoPosteriorError.invalid("Metal screen v1 requires exact-valued calibration observations") }
        // Correlation, when present, remains in the exact CPU likelihood. This
        // screen deliberately uses independent marginal Gaussian observations.
        func upper(_ field: VivoKineticFitField, _ fallback: Double) -> Double {
            problem.bindings.first { $0.field == field && $0.caseIdentifiers.contains(item.identifier) }?.parameter.upper ?? fallback
        }
        func lower(_ field: VivoKineticFitField, _ fallback: Double) -> Double {
            problem.bindings.first { $0.field == field && $0.caseIdentifiers.contains(item.identifier) }?.parameter.lower ?? fallback
        }
        let kon = upper(.association, model.association.value)
        let koff = upper(.dissociation, model.dissociation.value)
        let kinact = upper(.inactivation, model.inactivation.value)
        let turnover = upper(.targetTurnover, model.targetTurnover.value)
        let xon = upper(.competitorAssociation, model.competitor?.association.value ?? 0)
        let xoff = upper(.competitorDissociation, model.competitor?.dissociation.value ?? 0)
        let x = upper(.competitorConcentration, model.competitor?.unboundConcentration.value ?? 0)
        let scale = upper(.exposureScale, 1)
        let exposure = (item.experiment.exposure.knots.map(\.unboundDrugM).max() ?? 0) * scale
        let rate = max(kon * exposure + xon * x + turnover, koff + kinact + turnover, xoff + turnover)
        guard rate.isFinite, rate >= 0 else { throw VivoPosteriorError.invalid("Metal prior-rate bound overflow") }
        let bound = min(policy.maximumTimeStepSeconds, rate > 0 ? 0.1 / rate : policy.maximumTimeStepSeconds)
        let dt = pow(2, floor(log2(bound)))
        _ = try vivoScreenFloat(dt, label: "fixed integration step")
        guard dt > 0, let end = item.observations.map(\.timeSeconds).max(),
              end / dt <= Double(policy.maximumStepsPerCase) else {
            throw VivoPosteriorError.budget("prior-bounded F1 screen is too stiff for its step budget")
        }
        let kineticsID = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(model))
        let experimentID = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(item.experiment))
        let original = try VivoTargetEngagementProgramSource.make(item.experiment,
            kineticsFingerprint: kineticsID.hex, experimentFingerprint: experimentID.hex)
        guard var source = try JSONSerialization.jsonObject(with: original) as? [String: Any],
              var spec = source["spec"] as? [String: Any], var parameters = spec["parameters"] as? [[String: Any]],
              var inputs = spec["inputs"] as? [[String: Any]], var metadata = source["metadata"] as? [String: Any] else {
            throw VivoPosteriorError.invalid("target compiler source structure")
        }
        let fields: [String: VivoKineticFitField] = ["kon": .association, "koff": .dissociation,
            "kinact": .inactivation, "kturnover": .targetTurnover,
            "competitor-kon": .competitorAssociation, "competitor-koff": .competitorDissociation]
        for index in parameters.indices {
            guard let name = parameters[index]["id"] as? String, let field = fields[name],
                  let value = parameters[index]["value"] as? Double else {
                throw VivoPosteriorError.invalid("unrecognized target parameter")
            }
            // Computational bounds are not statistical uncertainty intervals.
            // Cover both the initialization value and the declared prior support.
            let lo = min(value, lower(field, value)), hi = max(value, upper(field, value))
            _ = try vivoScreenFloat(lo, label: name + " lower bound")
            _ = try vivoScreenFloat(hi, label: name + " upper bound")
            // Metadata stores FP64 bounds but runtime values are FP32. Include
            // their outward rounding so an endpoint is not spuriously rejected.
            parameters[index]["bounds"] = ["min": min(lo, Double(Float(lo))), "max": max(hi, Double(Float(hi)))]
        }
        for index in inputs.indices {
            if inputs[index]["id"] as? String == VivoTargetEngagementProgramSource.competitorInput {
                let maximum = max(x, model.competitor?.unboundConcentration.value ?? 0)
                inputs[index]["bounds"] = ["min": 0, "max": max(maximum, Double(try vivoScreenFloat(maximum, label: "competitor")))]
            } else if inputs[index]["id"] as? String == VivoTargetEngagementProgramSource.drugInput {
                let maximum = model.maximumUnboundDrugM
                inputs[index]["bounds"] = ["min": 0, "max": max(maximum, Double(try vivoScreenFloat(maximum, label: "exposure")))]
            }
        }
        var labels = metadata["labels"] as? [String: String] ?? [:]
        labels["execution-profile"] = "independent-f1-screen-not-authoritative-likelihood"
        labels["parameter-bounds"] = "numerical-prior-support-not-confidence-intervals"
        metadata["labels"] = labels; source["metadata"] = metadata
        spec["parameters"] = parameters; spec["inputs"] = inputs; source["spec"] = spec
        let bytes = try JSONSerialization.data(withJSONObject: source, options: [.sortedKeys, .withoutEscapingSlashes])
        let compiled = try VivoNativeCompilerBridge.compileProgram(bytes, fidelity: .deterministic,
            strictUnits: true, strictSafety: true, requireTermination: true)
        guard compiled.succeeded else { throw VivoPosteriorError.invalid(compiled.diagnosticsString()) }
        let pack = try VivoProgramPack(data: compiled.primary)
        let species = try pack.speciesMetadata()
        func index(_ id: String) throws -> UInt32 {
            guard let i = species.firstIndex(where: { $0.identifier == id }) else {
                throw VivoPosteriorError.invalid("missing compiled target state \(id)")
            }
            return UInt32(i)
        }
        let indices = try VivoTargetEngagementProgramSource.targetSpecies.map(index)
        let drugIndex = try index(VivoTargetEngagementProgramSource.drugInput)
        let competitorIndex = try model.competitor.map { _ in try index(VivoTargetEngagementProgramSource.competitorInput) }
        guard species[Int(drugIndex)].isExternallyOwned,
              competitorIndex.map({ species[Int($0)].isExternallyOwned }) ?? true,
              pack.runtimeContract.temporalStateCount == 0, pack.runtimeContract.ruleCount == 0,
              !species.contains(where: \.isCountValued) else {
            throw VivoPosteriorError.invalid("screen compiler produced unsupported ownership or temporal state")
        }
        let capacity = problem.sampler.particleCount
        let config = VivoRuntimeConfiguration(fidelity: .deterministic, environmentCount: UInt32(capacity),
            timeStep: Float(dt), minimumTimeStep: .leastNormalMagnitude, maximumTimeStep: Float(dt),
            maximumSubsteps: 1, eventCapacity: UInt32(max(1024, capacity * 4)))
        try VivoProgramExecutionContract.validate(pack: pack, configuration: config)
        var operations: [VivoMetalTargetOperation] = [], time = 0.0, steps = 0
        let knots = item.experiment.exposure.knots
        var boundaries = Set(item.observations.map(\.timeSeconds))
        for knot in knots where knot.timeSeconds <= end { boundaries.insert(knot.timeSeconds) }
        var knotIndex = 0
        operations.append(.exposure(try vivoScreenFloat(knots[0].unboundDrugM, label: "exposure knot")))
        for boundary in boundaries.sorted() {
            guard Double(Float(boundary)) == boundary else {
                throw VivoPosteriorError.invalid("screen observation/exposure clock is not exactly FP32-representable")
            }
            while time < boundary {
                let step = min(dt, boundary - time), after = time + step
                guard steps < policy.maximumStepsPerCase, Double(Float(step)) == step,
                      Float(after) > Float(time), Float(after).isFinite else {
                    throw VivoPosteriorError.budget("screen clock resolution or fixed-step count")
                }
                operations.append(.step(dt: Float(step), after: Float(after)))
                time = after; steps += 1
            }
            while knotIndex + 1 < knots.count, knots[knotIndex + 1].timeSeconds <= boundary {
                knotIndex += 1
                operations.append(.exposure(try vivoScreenFloat(knots[knotIndex].unboundDrugM, label: "exposure knot")))
            }
            for index in item.observations.indices where item.observations[index].timeSeconds == boundary {
                let observation = item.observations[index]
                _ = try vivoScreenFloat(observation.value, label: "observation")
                _ = try vivoScreenFloat(observation.standardDeviation ?? 0, label: "reported SD")
                operations.append(.observe(index))
            }
        }
        let description = try VivoMetalTargetCaseDescription(caseIdentifier: item.identifier,
            programFingerprint: pack.header.contentFingerprint,
            scheduleFingerprint: VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(operations)),
            timeStepSeconds: dt, steps: steps, operations: operations.count, capacity: capacity)
        return .init(item: item, pack: pack, configuration: config, description: description, operations: operations,
            drugIndex: drugIndex, competitorIndex: competitorIndex,
            targetIndices: SIMD4(indices[0], indices[1], indices[2], indices[3]),
            parameterNames: try pack.parameterMetadata().map(\.identifier))
    }
}
