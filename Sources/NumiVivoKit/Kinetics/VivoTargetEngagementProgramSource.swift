import Foundation

/// Lowers a covalent target model to the existing VivoProgram language. Target
/// states use normalized fractions; unbound input uses M. This avoids tiny target
/// concentrations being compared to an absolute FP32 concentration tolerance.
public enum VivoTargetEngagementProgramSource {
    public static let targetSpecies = ["target-free", "target-reversible", "target-covalent", "target-competitor"]
    public static let drugInput = "drug-unbound"
    public static let competitorInput = "competitor-unbound"

    public static func make(_ experiment: VivoTargetEngagementExperiment,
                            kineticsFingerprint: String, experimentFingerprint: String) throws -> Data {
        try experiment.validate()
        guard VivoKineticEvidence.isSHA256(kineticsFingerprint), VivoKineticEvidence.isSHA256(experimentFingerprint) else {
            throw VivoKineticsError.invalid("lowering requires immutable kinetics and experiment identities")
        }
        let m = experiment.kinetics
        let duration = experiment.exposure.knots.last!.timeSeconds
        let stopTime = duration + max(1, duration * 1e-5)
        let represented = m.parameters.map(\.value) + experiment.initial.values
            + [m.maximumUnboundDrugM, stopTime] + experiment.exposure.knots.map(\.unboundDrugM)
        for value in represented {
            guard abs(value) <= Double(Float.greatestFiniteMagnitude), value == 0 || Float(value) != 0 else {
                throw VivoKineticsError.unsupported("a kinetic value cannot be lowered to FP32 without overflow/underflow")
            }
        }
        func evidence(_ p: VivoKineticParameter) -> [String: Any] {
            let classification: String
            switch p.origin {
            case .assumed: classification = "assumed"
            case .measured: classification = "observed"
            case .fitted: classification = "calibrated"
            case .calculated: classification = "derived"
            }
            return ["class": classification, "source": p.evidence.source,
                    "dataset": p.evidence.sourceFingerprint ?? "", "context": m.context.hostContext,
                    "note": p.evidence.locator]
        }
        func parameter(_ id: String, _ p: VivoKineticParameter) -> [String: Any] {
            var result: [String: Any] = ["id": id, "value": p.value, "unit": p.unit.rawValue, "evidence": evidence(p)]
            if case .bounded(let lower, let upper) = p.uncertainty { result["bounds"] = ["min": lower, "max": upper] }
            return result
        }
        func reaction(_ id: String, _ from: [String], _ to: [String], _ rate: String, law: String = "mass-action") -> [String: Any] {
            ["id": id, "compartment": "cell", "reactants": from, "products": to,
             "rate": ["law": law, "parameters": [rate]]]
        }
        var inputs: [[String: Any]] = [["id": drugInput, "source": "experiment", "unit": "M",
            "default": experiment.exposure.knots[0].unboundDrugM,
            "bounds": ["min": 0, "max": m.maximumUnboundDrugM]]]
        let species: [[String: Any]] = zip(targetSpecies, experiment.initial.values).map { id, value in
            ["id": id, "kind": "occupancy", "compartment": "cell", "unit": "fraction", "initial": value,
             "bounds": ["min": 0, "max": 1.0001], "externallyOwned": false, "conserved": false]
        }
        var parameters = [parameter("kon", m.association), parameter("koff", m.dissociation),
                          parameter("kinact", m.inactivation), parameter("kturnover", m.targetTurnover)]
        var reactions = [reaction("association", [drugInput, targetSpecies[0]], [targetSpecies[1]], "kon"),
                         reaction("dissociation", [targetSpecies[1]], [targetSpecies[0]], "koff"),
                         reaction("inactivation", [targetSpecies[1]], [targetSpecies[2]], "kinact"),
                         reaction("target-synthesis", [], [targetSpecies[0]], "kturnover", law: "zero-order")]
        for speciesID in targetSpecies { reactions.append(reaction("loss-" + speciesID, [speciesID], [], "kturnover")) }
        if let competitor = m.competitor {
            inputs.append(["id": competitorInput, "source": "experiment", "unit": "M",
                           "default": competitor.unboundConcentration.value,
                           "bounds": ["min": competitor.unboundConcentration.value, "max": competitor.unboundConcentration.value]])
            parameters += [parameter("competitor-kon", competitor.association), parameter("competitor-koff", competitor.dissociation)]
            reactions += [reaction("competitor-association", [competitorInput, targetSpecies[0]], [targetSpecies[3]], "competitor-kon"),
                          reaction("competitor-dissociation", [targetSpecies[3]], [targetSpecies[0]], "competitor-koff")]
        }
        var total: [String: Any] = ["species": targetSpecies[0]]
        for id in targetSpecies.dropFirst() { total = ["add": [total, ["species": id]]] }
        let source: [String: Any] = [
            "apiVersion": "numivivo.org/v1alpha1", "kind": "VivoProgram",
            "metadata": ["name": "target-engagement", "version": "1.0.0",
                "description": "Conditional covalent target occupancy; not a clinical efficacy model.",
                "labels": ["kinetics-sha256": kineticsFingerprint, "experiment-sha256": experimentFingerprint,
                           "exposure": m.exposureTreatment.rawValue, "turnover": m.turnoverModel.rawValue,
                           "host-context": m.context.hostContext, "compound": m.context.compound,
                           "target": m.context.target, "target-variant": m.context.targetVariant,
                           "chemical-state": m.context.chemicalState]],
            "spec": ["minimumFidelity": "F1", "target": ["cellType": "declared-kinetic-context"],
                "inputs": inputs, "species": species, "parameters": parameters, "reactions": reactions,
                "constraints": [["id": "normalized-target-balance", "severity": "error", "response": "reject-step",
                    "expression": ["any": [
                        ["gt": [total, ["literal": ["value": 1.0001, "unit": "fraction"]]]],
                        ["lt": [total, ["literal": ["value": 0.9999, "unit": "fraction"]]]]
                    ]], "message": "Target fraction balance exceeded the declared numerical tolerance."]],
                "termination": [["id": "bounded-computational-horizon", "when": ["gt": [
                    ["time": true], ["literal": ["value": stopTime, "unit": "s"]]]],
                    "action": "permanent-shutdown", "reason": "Declared computational horizon ended."]]
            ]
        ]
        return try JSONSerialization.data(withJSONObject: source, options: [.sortedKeys, .withoutEscapingSlashes])
    }
}
