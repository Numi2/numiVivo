import Foundation

public struct VivoConditionalEncounterWorkflowResult: Codable, Sendable, Equatable {
    public static let schemaID = "numivivo.org/conditional-encounter-workflow-result/v1"
    public let schema: String
    /// Digest of the exact decoded reaction input payload bytes. The chemistry
    /// task/receipt retains the input artifact (including its workflow envelope).
    public let sourceReactionPayloadFingerprint: VivoFingerprint
    public let encounter: VivoConditionalEncounterResult
}

/// A pure consumer of the existing reaction operation's full TST output. Input
/// artifact identity remains owned by the existing chemistry workflow receipts.
public enum VivoConditionalEncounterWorkflow {
    public static let operationIdentifier = "vivo.native.conditional-encounter"
    public static let requestKind = "vivo.conditional-encounter-request"
    public static let resultKind = "vivo.conditional-encounter-result"

    public static func operation(implementationFingerprint id: VivoFingerprint) -> VivoChemistryOperation {
        .init(identifier: operationIdentifier, version: "1", implementationFingerprint: id,
            outputs: [.init(name: "result", kind: resultKind)], execute: { cfg, inputs, budget in
                ["result": try VivoCanonicalJSON.encode(calculate(cfg, inputs: inputs, budget: budget))]
            }, validateOutputs: { cfg, inputs, outputs, budget in
                try Task.checkCancellation()
                guard Set(outputs.keys) == ["result"], let data = outputs["result"] else {
                    throw VivoChemistryError.invalid("conditional encounter output slots")
                }
                let result = try VivoCanonicalJSON.decode(VivoConditionalEncounterWorkflowResult.self, from: data)
                guard result == (try calculate(cfg, inputs: inputs, budget: budget)) else {
                    throw VivoChemistryError.invalid("conditional encounter output/source reconstruction")
                }
            })
    }

    private static func calculate(_ cfg: VivoJSONValue, inputs: [String: Data],
                                  budget: VivoChemistryBudget) throws -> VivoConditionalEncounterWorkflowResult {
        try Task.checkCancellation(); try budget.validate()
        guard cfg == .object([:]), Set(inputs.keys) == ["request", "reaction"],
              let requestData = inputs["request"], let reactionData = inputs["reaction"],
              requestData.count <= budget.maximumBytes, reactionData.count <= budget.maximumBytes - requestData.count else {
            throw VivoChemistryError.invalid("conditional encounter input slots or aggregate byte budget")
        }
        let request = try VivoCanonicalJSON.decode(VivoConditionalEncounterRequest.self, from: requestData)
        let wrapper = try VivoCanonicalJSON.decode(VivoReactionCalculationResult.self, from: reactionData)
        let source: VivoTransitionStateTheoryResult
        switch wrapper {
        case .connectedReaction(let result), .transitionStateTheory(let result): source = result
        default: throw VivoChemistryError.invalid("conditional encounter needs a complete molecular TST reaction result")
        }
        guard source.request.connectivity.request.saddle.request.model.budget == budget,
              source.request.connectivity.request.endpoints.allSatisfy({ $0.components.allSatisfy { $0.point.request.model.budget == budget } }),
              request.numerics.maximumPrimitiveWork <= budget.maximumOperatorApplications else {
            throw VivoChemistryError.invalid("conditional encounter source/numerical resource contracts differ from the task")
        }
        _ = try budget.elements([request.observationTimesSeconds.count, 1 + request.reservoirs.count, 12], simultaneousArrays: 2)
        let result = try VivoConditionalEncounter.calculate(request, source: source)
        return .init(schema: VivoConditionalEncounterWorkflowResult.schemaID,
                     sourceReactionPayloadFingerprint: try VivoCanonicalJSON.fingerprint(reactionData), encounter: result)
    }
}
