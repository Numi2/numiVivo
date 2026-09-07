import Foundation

/// Additional adaptive evidence adapters kept separate from the core campaign
/// policy so dimension-specific QM/MM validation remains owned by the existing
/// chemical qualification implementation.
public enum VivoPlatformAdaptiveVariantOperations {
    public static let metricOperationIdentifier = "vivo.platform.metric-qmmm-variant"
    public static func definitions(implementationFingerprint id: VivoFingerprint) -> [VivoWorkflowDefinition] {
        [VivoPlatformOperations.pure(identifier: metricOperationIdentifier,id: id,
            inputs: ["request":"vivo.qmmm-variant-sensitivity-request"],
            outputs: [.init(name:"result",kind:"vivo.qmmm-variant-sensitivity-result"),
                      .init(name:"metric",kind:VivoPlatformAdaptiveOperations.metricKind)],
            summary: "Reconstruct one independently replicated electronic/QM-region or other QM/MM sensitivity variant and expose its statistically guarded log-rate shift.",
            configure: VivoPlatformOperations.empty,calculate: { _,data,_ in
                let request = try VivoPlatformOperations.input(VivoQMMMVariantSensitivityRequest.self,"request",data)
                let result = try VivoQMMMVariantSensitivity.calculate(request)
                try VivoQMMMVariantSensitivity.validate(result,request:request)
                let source = Set(VivoPlatformAdaptiveOperations.rateSources(request.baselineRequest)
                    + VivoPlatformAdaptiveOperations.rateSources(request.variant.request)).sorted()
                let evidence = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(result)).hex
                let metric = try VivoRefinementMetricEvidence(
                    metricUnit:"statistically-guarded-absolute-log-rate-shift",
                    observableFingerprint:VivoAdaptiveMetricContext.rate(request.baselineRequest).hex,
                    value:result.statisticallyGuardedAbsoluteLogRateShift,
                    nativeChecksPassed:result.statisticallyGuardedAbsoluteLogRateShift != nil,
                    evidenceIdentifier:evidence,executionSourceIdentifiers:source,
                    interpretation:result.interpretation)
                return ["result":try VivoCanonicalJSON.encode(result),"metric":try VivoCanonicalJSON.encode(metric)]
            })]
    }
}
