import Foundation

/// Uses the same artifact scheduler as electronic and nuclear calculations.
/// The request embeds exact mapped endpoint records; its fingerprint therefore
/// binds their geometries, methods, masses and qualification data as well.
public enum VivoReactionConnectivityWorkflow {
    public static let requestKind = "vivo.mapped-reaction-connectivity-request"
    public static let resultKind = "vivo.mapped-reaction-connectivity-result"

    private static func request(_ inputs: [String: Data], budget: VivoChemistryBudget) throws -> VivoReactionConnectivityRequest {
        guard Set(inputs.keys) == Set(["request"]), let data = inputs["request"] else {
            throw VivoChemistryError.invalid("reaction connectivity operation input slots")
        }
        let value = try VivoCanonicalJSON.decode(VivoReactionConnectivityRequest.self, from: data)
        try VivoReactionConnectivity.validateRequest(value)
        guard value.saddle.request.model.budget == budget,
              value.endpoints.allSatisfy({ endpoint in
                  endpoint.components.allSatisfy { $0.point.request.model.budget == budget }
              }) else {
            throw VivoChemistryError.invalid("reaction endpoint or saddle numerical resource contract differs from the task")
        }
        return value
    }

    public static func operation(implementationFingerprint id: VivoFingerprint) -> VivoChemistryOperation {
        .init(identifier: "vivo.native.mapped-reaction-connectivity", version: "1", implementationFingerprint: id,
              outputs: [.init(name: "connection", kind: resultKind)], execute: { configuration, inputs, budget in
            guard configuration == .object([:]) else {
                throw VivoChemistryError.invalid("unrecognized reaction connectivity operation configuration")
            }
            let source = try request(inputs, budget: budget)
            let result = try VivoReactionConnectivity.run(source)
            guard result.converged else {
                throw VivoChemistryError.convergence("reaction connectivity failed displacement or step refinement; no qualified artifact published")
            }
            return ["connection": try VivoCanonicalJSON.encode(result)]
        }, validateOutputs: { configuration, inputs, outputs, budget in
            guard configuration == .object([:]), Set(outputs.keys) == Set(["connection"]),
                  let data = outputs["connection"] else {
                throw VivoChemistryError.invalid("reaction connectivity output slot or configuration contract")
            }
            let source = try request(inputs, budget: budget)
            let result = try VivoCanonicalJSON.decode(VivoReactionConnectivityResult.self, from: data)
            try VivoReactionConnectivity.validate(result, request: source)
        })
    }

    public static func plan(_ source: VivoReactionConnectivityRequest, requestArtifact: VivoFingerprint,
                            implementationFingerprint id: VivoFingerprint) throws -> VivoElectronicWorkflowPlan {
        try VivoReactionConnectivity.validateRequest(source)
        let budget = source.saddle.request.model.budget
        guard source.endpoints.allSatisfy({ $0.components.allSatisfy { $0.point.request.model.budget == budget } }) else {
            throw VivoChemistryError.invalid("mapped endpoint resource contract differs from the saddle")
        }
        let node = VivoChemistryDAGNode(identifier: "reaction-connectivity", operation: operation(implementationFingerprint: id),
            inputs: [.artifact(name: "request", fingerprint: requestArtifact, kind: requestKind)],
            configuration: .object([:]), resources: .init(budget: budget,
                maximumInputBytes: budget.maximumBytes, maximumOutputBytes: budget.maximumBytes))
        return .init(nodes: [node], resultNode: "reaction-connectivity", resultOutput: "connection", resultKind: resultKind)
    }
}
