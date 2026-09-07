import Foundation

/// Deterministic artifact adapters for bounded reaction-mechanism hypotheses.
/// These operations enumerate mapped graph candidates only; they do not certify
/// electronic feasibility, a saddle, connectivity, a barrier, or a kinetic rate.
public enum VivoPlatformMechanismOperations {
    public static func definitions(implementationFingerprint id: VivoFingerprint) -> [VivoWorkflowDefinition] {
        [VivoPlatformOperations.pure(identifier: "vivo.platform.mechanism-candidates", id: id,
            inputs: ["request": "vivo.mechanism-candidate-request"],
            outputs: [.init(name: "candidates", kind: "vivo.mechanism-candidate-result")],
            summary: "Bounded mapped bond-change hypotheses for downstream electronic/path/saddle qualification; not mechanism validation.",
            configure: VivoPlatformOperations.empty,
            calculate: { _, inputs, _ in
                let request = try VivoPlatformOperations.input(VivoMechanismCandidateRequest.self, "request", inputs)
                let result = try VivoMechanismCandidateEnumeration.enumerate(request)
                return ["candidates": try VivoCanonicalJSON.encode(result)]
            })]
    }
}
