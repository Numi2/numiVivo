import Foundation

public enum VivoPlatformAdaptiveOperations {
    public static let metricKind = "vivo.refinement-metric-evidence"
    /// Only these adapters may supply campaign acceptance metrics. Their native
    /// validators run on fresh execution and on every cached workflow reuse.
    static let metricOperations: Set<String> = ["vivo.platform.mbar-target-refinement",
        "vivo.platform.metric-electronic-refinement","vivo.platform.metric-replicated-rate",
        "vivo.platform.metric-chemical-sensitivity","vivo.platform.kinetic-observable-covariance"]
    private static func evidence<T: Encodable>(_ value: T) throws -> String {
        try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(value)).hex
    }
    static func rateSources(_ request: VivoQMMMReplicatedFreeEnergyRateRequest) -> [String] {
        var values = request.replicas.map { $0.freeEnergy.provenance.samplingExecution.requestFingerprint.hex }
        values += request.replicas.compactMap { $0.transmissionEvidence.sourceFingerprint }
        return Set(values).sorted()
    }
    public static func definitions(implementationFingerprint id: VivoFingerprint) -> [VivoWorkflowDefinition] {
        let resultPort = VivoChemistryTaskOutput(name: "metric",kind: metricKind)
        return [
            VivoPlatformOperations.pure(identifier: "vivo.platform.mbar-target-refinement",id: id,
                inputs: ["request":"vivo.mbar-target-refinement-request"],
                outputs: [.init(name: "result",kind: "vivo.mbar-target-refinement-result"),resultPort],
                summary: "Shared native MBAR with paired target/profile covariance, per-bin support checks and explicit exploratory targets.",
                configure: VivoPlatformOperations.empty,calculate: { _,data,budget in
                    let request = try VivoPlatformOperations.input(VivoMBARTargetRefinementRequest.self,"request",data)
                    guard request.configuration.maximumPrimitiveElements <= budget.maximumOperatorApplications else {
                        throw VivoChemistryError.resourceLimit("target MBAR request exceeds workflow numerical budget")
                    }
                    let result = try VivoMBARTargetRefinement.calculate(request)
                    let sd = result.corrections.compactMap(\.bootstrapStandardDeviation).max()
                    let metric = try VivoRefinementMetricEvidence(metricUnit: "dimensionless-bin-free-energy-standard-deviation",
                        observableFingerprint: VivoAdaptiveMetricContext.targetProfile(request).hex,value: sd,nativeChecksPassed: result.qualification == .sampledTargetsPassDiagnostics && sd != nil,
                        evidenceIdentifier: evidence(["request":try VivoPlatformOperations.json(request),"result":try VivoPlatformOperations.json(result)]),
                        executionSourceIdentifiers: Set(request.samples.map { "declared-block:"+$0.independentBlockIdentifier }).sorted(),
                        interpretation: "Conditional equilibrium histogram uncertainty; independence is declared by the input block partition, not measured here. Not log-rate uncertainty or dynamical evidence.")
                    return ["result":try VivoCanonicalJSON.encode(result),"metric":try VivoCanonicalJSON.encode(metric)]
                }),
            VivoPlatformOperations.pure(identifier: "vivo.platform.chemical-exchange-sensitivity",id: id,
                inputs: ["request":"vivo.chemical-exchange-sensitivity-request"],
                outputs: [.init(name: "result",kind: "vivo.chemical-exchange-sensitivity-result")],
                summary: "Evidence-bound chemical/exchange log-rate derivatives through the shared transient generator.",
                configure: VivoPlatformOperations.empty,calculate: { _,data,budget in
                    let request = try VivoPlatformOperations.input(VivoQMMMChemicalExchangeSensitivityRequest.self,"request",data)
                    guard request.configuration.maximumPrimitiveWork <= budget.maximumOperatorApplications else {
                        throw VivoChemistryError.resourceLimit("kinetic sensitivity exceeds workflow budget")
                    }
                    return ["result":try VivoCanonicalJSON.encode(VivoQMMMChemicalExchangeSensitivity.calculate(request))]
                }),
            VivoPlatformOperations.pure(identifier: "vivo.platform.metric-electronic-refinement",id: id,
                inputs: ["result":"vivo.property-directed-space-result"],outputs: [resultPort],
                summary: "Reconstruct held-out electronic-space sensitivity before producing a campaign metric.",
                configure: VivoPlatformOperations.empty,calculate: { _,data,_ in
                    let result = try VivoPlatformOperations.input(VivoPropertyDirectedSpaceResult.self,"result",data)
                    try VivoPropertyDirectedSpace.validate(result,request: result.request)
                    let metric = try VivoRefinementMetricEvidence(metricUnit: "normalized-electronic-profile-sensitivity",
                        observableFingerprint: VivoAdaptiveMetricContext.electronic(request: result.request).hex,value: result.confirmation?.shift.normalizedPropertyImpact,
                        nativeChecksPassed: result.sensitivityEstablishedWithinDeclaredPool,evidenceIdentifier: evidence(result),
                        executionSourceIdentifiers: result.request.points.filter { result.request.confirmationPointIdentifiers.contains($0.identifier) }.map {
                            try "held-out-hamiltonian:"+evidence($0.hamiltonian)
                        }.sorted(),interpretation: result.meaning)
                    return ["metric":try VivoCanonicalJSON.encode(metric)]
                }),
            VivoPlatformOperations.pure(identifier: "vivo.platform.metric-replicated-rate",id: id,
                inputs: ["request":"vivo.qmmm-replicated-free-energy-rate-request","result":"vivo.qmmm-replicated-free-energy-rate-result"],
                outputs: [resultPort],summary: "Validated replicated conditional log-rate uncertainty; no invented missing uncertainty.",
                configure: VivoPlatformOperations.empty,calculate: { _,data,_ in
                    let request = try VivoPlatformOperations.input(VivoQMMMReplicatedFreeEnergyRateRequest.self,"request",data)
                    let result = try VivoPlatformOperations.input(VivoQMMMReplicatedFreeEnergyRateResult.self,"result",data)
                    try VivoQMMMReplicatedFreeEnergyRate.validate(result,request: request)
                    let metric = try VivoRefinementMetricEvidence(metricUnit: "conditional-log-rate-standard-deviation",
                        observableFingerprint: VivoAdaptiveMetricContext.rate(request).hex,
                        value: result.combinedConditionalLogRateStandardDeviation,
                        nativeChecksPassed: result.converged && result.combinedConditionalLogRateStandardDeviation != nil,
                        evidenceIdentifier: result.evidenceFingerprint.hex,executionSourceIdentifiers: rateSources(request),
                        interpretation: result.interpretation)
                    return ["metric":try VivoCanonicalJSON.encode(metric)]
                }),
            VivoPlatformOperations.pure(identifier: "vivo.platform.metric-chemical-sensitivity",id: id,
                inputs: ["request":"vivo.qmmm-chemical-qualification-request","result":"vivo.qmmm-chemical-qualification-result"],
                outputs: [resultPort],summary: "Existing chemical sensitivity gate plus conservative correlated statistical guard.",
                configure: { cfg in
                    let factor = try VivoPlatformOperations.decode(VivoChemicalSensitivityMetricConfiguration.self,cfg)
                    guard factor.standardDeviationMultiplier.isFinite, factor.standardDeviationMultiplier >= 0,
                          factor.standardDeviationMultiplier <= 10 else { throw VivoChemistryError.invalid("statistical sensitivity guard") }
                },calculate: { cfg,data,_ in
                    let request = try VivoPlatformOperations.input(VivoQMMMChemicalQualificationRequest.self,"request",data)
                    let result = try VivoPlatformOperations.input(VivoQMMMChemicalQualificationResult.self,"result",data)
                    try VivoQMMMChemicalQualification.validate(result,request: request)
                    let factor = try VivoPlatformOperations.decode(VivoChemicalSensitivityMetricConfiguration.self,cfg).standardDeviationMultiplier
                    var guarded: [Double] = [], known = true
                    for variant in request.variants {
                        guard let a = request.baselineResult.combinedConditionalLogRateStandardDeviation,
                              let b = variant.result.combinedConditionalLogRateStandardDeviation else { known = false; continue }
                        guarded.append(abs(variant.result.meanLogRatePerSecond-request.baselineResult.meanLogRatePerSecond)+factor*(a+b))
                    }
                    let value = known ? guarded.max() : nil
                    let sources = Set(rateSources(request.baselineRequest)+request.variants.flatMap { rateSources($0.request) }).sorted()
                    let metric = try VivoRefinementMetricEvidence(metricUnit: "statistically-guarded-absolute-log-rate-shift",
                        observableFingerprint: VivoAdaptiveMetricContext.rate(request.baselineRequest).hex,
                        value: value,nativeChecksPassed: result.converged && value != nil,evidenceIdentifier: result.evidenceFingerprint.hex,
                        executionSourceIdentifiers: sources,
                        interpretation: "Absolute conditional log-rate shift plus a declared multiple of sigma_baseline + sigma_variant. This conservative standard-deviation inequality does not assume independent errors and is not a confidence interval or model-error bound.")
                    return ["metric":try VivoCanonicalJSON.encode(metric)]
                }),
            kineticCovariance(id: id)
        ]
    }
}
public struct VivoChemicalSensitivityMetricConfiguration: Codable, Sendable, Equatable {
    public let standardDeviationMultiplier: Double
    public init(standardDeviationMultiplier: Double = 2) { self.standardDeviationMultiplier = standardDeviationMultiplier }
}
public enum VivoKineticSensitivityObservable: Codable, Sendable, Equatable {
    case survivalProbability, reactedProbability, hazardPerSecond, stateProbability(index: Int)
}
public struct VivoKineticObservableCovarianceRequest: Codable, Sendable, Equatable {
    public let parameterIdentifiers: [String]
    public let meanRatesPerSecond: [Double]
    public let logRateCovariance: VivoQMMatrix
    public let observationIndex: Int
    public let observable: VivoKineticSensitivityObservable
    public let covarianceEvidence: VivoKineticEvidence
    public let covarianceOrigin: VivoKineticOrigin
    public init(parameterIdentifiers: [String],meanRatesPerSecond: [Double],logRateCovariance: VivoQMMatrix,
                observationIndex: Int,observable: VivoKineticSensitivityObservable,covarianceEvidence: VivoKineticEvidence,
                covarianceOrigin: VivoKineticOrigin = .fitted) {
        self.parameterIdentifiers = parameterIdentifiers; self.meanRatesPerSecond = meanRatesPerSecond
        self.logRateCovariance = logRateCovariance; self.observationIndex = observationIndex
        self.observable = observable; self.covarianceEvidence = covarianceEvidence; self.covarianceOrigin = covarianceOrigin
    }
}
private extension VivoPlatformAdaptiveOperations {
    static func kineticCovariance(id: VivoFingerprint) -> VivoWorkflowDefinition {
        VivoPlatformOperations.pure(identifier: "vivo.platform.kinetic-observable-covariance",id: id,
            inputs: ["request":"vivo.chemical-exchange-sensitivity-request","result":"vivo.chemical-exchange-sensitivity-result",
                     "covariance":"vivo.kinetic-observable-covariance-request"],
            outputs: [.init(name: "result",kind: "vivo.local-observable-covariance"),.init(name: "metric",kind: metricKind)],
            summary: "Local observable uncertainty with an explicitly sourced joint log-rate covariance, not inferred independence.",
            configure: VivoPlatformOperations.empty,calculate: { _,data,_ in
                let request = try VivoPlatformOperations.input(VivoQMMMChemicalExchangeSensitivityRequest.self,"request",data)
                let result = try VivoPlatformOperations.input(VivoQMMMChemicalExchangeSensitivityResult.self,"result",data)
                let covariance = try VivoPlatformOperations.input(VivoKineticObservableCovarianceRequest.self,"covariance",data)
                try VivoQMMMChemicalExchangeSensitivity.validate(result,request: request)
                try covariance.covarianceEvidence.validate(origin: covariance.covarianceOrigin)
                guard covariance.parameterIdentifiers == result.parameters.map(\.identifier),
                      covariance.meanRatesPerSecond == result.parameters.map(\.ratePerSecond),
                      result.kinetics.observations.indices.contains(covariance.observationIndex) else {
                    throw VivoChemistryError.invalid("covariance parameter/expansion-point binding")
                }
                let point = result.kinetics.observations[covariance.observationIndex]
                let gradient: [Double], unit: String
                switch covariance.observable {
                case .survivalProbability: gradient = point.derivatives.map(\.survivalProbability); unit = "linearized-survival-probability-standard-deviation"
                case .reactedProbability: gradient = point.derivatives.map(\.reactedProbability); unit = "linearized-reacted-probability-standard-deviation"
                case .hazardPerSecond:
                    guard point.derivatives.allSatisfy({ $0.hazardPerSecond != nil }) else { throw VivoChemistryError.convergence("hazard derivative unavailable after survival underflow") }
                    gradient = point.derivatives.map { $0.hazardPerSecond! }; unit = "linearized-hazard-standard-deviation-per-second"
                case .stateProbability(let index):
                    guard point.probabilityByState.indices.contains(index) else { throw VivoChemistryError.invalid("kinetic covariance state index") }
                    gradient = point.derivatives.map { $0.probabilityByState[index] }; unit = "linearized-state-probability-standard-deviation"
                }
                let projection = try VivoObservableCovariance.project(derivatives: gradient,covariance: covariance.logRateCovariance)
                let metric = try VivoRefinementMetricEvidence(metricUnit: unit,observableFingerprint: VivoAdaptiveMetricContext.kinetic(request: request,covariance: covariance).hex,value: projection.standardDeviation,nativeChecksPassed: covariance.covarianceOrigin != .assumed,
                    evidenceIdentifier: evidence(["sensitivity":try VivoPlatformOperations.json(result),"covariance":try VivoPlatformOperations.json(covariance)]),
                    executionSourceIdentifiers: [result.requestFingerprint.hex],
                    interpretation: "Local Jacobian propagation of the supplied evidence-bound joint covariance at the bound rate values; not global nonlinear predictive uncertainty or experimental agreement.")
                return ["result":try VivoCanonicalJSON.encode(projection),"metric":try VivoCanonicalJSON.encode(metric)]
            })
    }
}

public enum VivoAdaptiveMetricContext {
    public static func electronic(request: VivoPropertyDirectedSpaceRequest) throws -> VivoFingerprint {
        struct Context: Encodable { let points: [VivoECCPathPoint]; let target: VivoElectronicProfileTarget }
        return try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(Context(points: request.points,target: request.target)))
    }
    public static func targetProfile(_ request: VivoMBARTargetRefinementRequest) throws -> VivoFingerprint {
        struct Context: Encodable {
            let sampledStates: [String]; let targets: [String]; let edges: [Double]
            let unit: String; let referenceTarget: Int; let referenceBin: Int
        }
        return try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(Context(sampledStates: request.sampledStateIdentifiers,
            targets: request.targetIdentifiers,edges: request.binEdges,unit: request.coordinateUnit,
            referenceTarget: request.referenceTargetIndex,referenceBin: request.referenceBinIndex)))
    }
    public static func rate(_ request: VivoQMMMReplicatedFreeEnergyRateRequest) throws -> VivoFingerprint {
        guard let first = request.replicas.first else { throw VivoChemistryError.invalid("empty rate context") }
        struct Context: Encodable { let context: VivoKineticContext; let reaction: VivoQMMMReactionConnectivityBinding }
        return try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(Context(context: first.context,
            reaction: first.freeEnergy.provenance.reactionConnectivity)))
    }
    public static func kinetic(request: VivoQMMMChemicalExchangeSensitivityRequest,
                               covariance: VivoKineticObservableCovarianceRequest) throws -> VivoFingerprint {
        struct Context: Encodable { let request: VivoQMMMChemicalExchangeSensitivityRequest; let observable: VivoKineticSensitivityObservable; let observationIndex: Int }
        return try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(Context(request: request,observable: covariance.observable,
            observationIndex: covariance.observationIndex)))
    }
}
