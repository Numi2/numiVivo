import Foundation

public struct VivoSBML3ExportOptions: Sendable, Codable, Equatable {
    public var modelIdentifier: String
    public var modelName: String?
    public var includeProvenance: Bool
    public var includeNumiVivoAnnotations: Bool

    public init(
        modelIdentifier: String = "numivivo_model",
        modelName: String? = nil,
        includeProvenance: Bool = true,
        includeNumiVivoAnnotations: Bool = true
    ) {
        self.modelIdentifier = modelIdentifier
        self.modelName = modelName
        self.includeProvenance = includeProvenance
        self.includeNumiVivoAnnotations = includeNumiVivoAnnotations
    }
}

public struct VivoSBML3Export: Sendable {
    public let xml: Data
    public let fingerprint: VivoFingerprint
    public let identifiers: VivoIdentifierMap
    public let diagnostics: [VivoStandardsDiagnostic]

    public init(
        xml: Data,
        fingerprint: VivoFingerprint,
        identifiers: VivoIdentifierMap,
        diagnostics: [VivoStandardsDiagnostic]
    ) {
        self.xml = xml
        self.fingerprint = fingerprint
        self.identifiers = identifiers
        self.diagnostics = diagnostics
    }
}

public struct VivoSBML3Exporter: Sendable {
    public init() {}

    public func export(
        _ pack: VivoProgramPack,
        options: VivoSBML3ExportOptions = .init()
    ) throws -> VivoSBML3Export {
        let species = try pack.speciesMetadata()
        let parameters = try pack.parameterMetadata()
        let reactions = try pack.reactionMetadata()
        let rules = try pack.ruleMetadata()
        let monitors = try pack.monitorMetadata()

        let compartmentNames = Array(Set(
            species.map(\.compartment) + reactions.map(\.compartment)
        )).sorted()
        let compartmentIDs = VivoXMLCodec.stableIdentifiers(compartmentNames, prefix: "compartment")
        let speciesIDs = VivoXMLCodec.stableIdentifiers(species.map(\.identifier), prefix: "species")
        let parameterIDs = VivoXMLCodec.stableIdentifiers(parameters.map(\.identifier), prefix: "parameter")
        let reactionIDs = VivoXMLCodec.stableIdentifiers(reactions.map(\.identifier), prefix: "reaction")
        let modelID = VivoXMLCodec.identifier(options.modelIdentifier, fallback: "numivivo_model")
        let identifiers = VivoIdentifierMap(
            model: modelID,
            compartments: compartmentIDs,
            species: speciesIDs,
            parameters: parameterIDs,
            reactions: reactionIDs
        )

        var diagnostics: [VivoStandardsDiagnostic] = []
        var lines: [String] = []
        lines.reserveCapacity(128 + species.count * 8 + reactions.count * 20)
        lines.append(#"<?xml version="1.0" encoding="UTF-8"?>"#)
        lines.append(#"<sbml xmlns="http://www.sbml.org/sbml/level3/version2/core" level="3" version="2">"#)
        var modelAttributes = #"id="\#(VivoXMLCodec.escapedAttribute(modelID))""#
        if let name = options.modelName, !name.isEmpty {
            modelAttributes += #" name="\#(VivoXMLCodec.escapedAttribute(name))""#
        }
        lines.append("  <model \(modelAttributes)>")
        lines.append("    <notes>")
        lines.append(#"      <body xmlns="http://www.w3.org/1999/xhtml">"#)
        lines.append("        <p>Generated from a NumiVivo ProgramPack. This document represents an abstract computational reaction model and does not assert sequence-level realizability.</p>")
        lines.append("      </body>")
        lines.append("    </notes>")

        if options.includeNumiVivoAnnotations {
            lines.append("    <annotation>")
            lines.append(#"      <numivivo:program xmlns:numivivo="https://numivivo.org/ns#" sourceFingerprint="\#(pack.header.sourceFingerprint.hex)" programPackFingerprint="\#(pack.header.contentFingerprint.hex)" fidelity="\#(pack.header.fidelity.label)" rules="\#(rules.count)" monitors="\#(monitors.count)"/>"#)
            lines.append("    </annotation>")
        }

        if !compartmentNames.isEmpty {
            lines.append("    <listOfCompartments>")
            for name in compartmentNames {
                guard let identifier = compartmentIDs[name] else { continue }
                lines.append(#"      <compartment id="\#(VivoXMLCodec.escapedAttribute(identifier))" name="\#(VivoXMLCodec.escapedAttribute(name))" spatialDimensions="3" size="1" constant="true"/>"#)
            }
            lines.append("    </listOfCompartments>")
        }

        if !species.isEmpty {
            lines.append("    <listOfSpecies>")
            for value in species {
                guard let identifier = speciesIDs[value.identifier],
                      let compartment = compartmentIDs[value.compartment] else {
                    throw VivoArtifactValidationError.unresolved("SBML identifier map is incomplete for species \(value.identifier)")
                }
                var attributes = #"id="\#(VivoXMLCodec.escapedAttribute(identifier))" name="\#(VivoXMLCodec.escapedAttribute(value.identifier))" compartment="\#(VivoXMLCodec.escapedAttribute(compartment))""#
                if value.isCountValued {
                    attributes += #" initialAmount="\#(VivoXMLCodec.finite(value.initialValue))" hasOnlySubstanceUnits="true""#
                } else {
                    attributes += #" initialConcentration="\#(VivoXMLCodec.finite(value.initialValue))" hasOnlySubstanceUnits="false""#
                }
                attributes += #" boundaryCondition="\#(value.isExternallyOwned ? "true" : "false")" constant="false""#
                if options.includeNumiVivoAnnotations {
                    lines.append("      <species \(attributes)>")
                    lines.append("        <annotation>")
                    lines.append(#"          <numivivo:state xmlns:numivivo="https://numivivo.org/ns#" unit="\#(VivoXMLCodec.escapedAttribute(value.unit))" minimum="\#(VivoXMLCodec.finite(value.minimum))" maximum="\#(VivoXMLCodec.finite(value.maximum))" flags="\#(value.flags)"/>"#)
                    lines.append("        </annotation>")
                    lines.append("      </species>")
                } else {
                    lines.append("      <species \(attributes)/>")
                }
            }
            lines.append("    </listOfSpecies>")
        }

        if !parameters.isEmpty {
            lines.append("    <listOfParameters>")
            for value in parameters {
                guard let identifier = parameterIDs[value.identifier] else {
                    throw VivoArtifactValidationError.unresolved("SBML identifier map is incomplete for parameter \(value.identifier)")
                }
                let attributes = #"id="\#(VivoXMLCodec.escapedAttribute(identifier))" name="\#(VivoXMLCodec.escapedAttribute(value.identifier))" value="\#(VivoXMLCodec.finite(value.value))" constant="true""#
                if options.includeNumiVivoAnnotations {
                    lines.append("      <parameter \(attributes)>")
                    lines.append("        <annotation>")
                    var annotation = #"          <numivivo:parameter xmlns:numivivo="https://numivivo.org/ns#" unit="\#(VivoXMLCodec.escapedAttribute(value.unit))" minimum="\#(VivoXMLCodec.finite(value.minimum))" maximum="\#(VivoXMLCodec.finite(value.maximum))" evidenceClass="\#(value.evidenceClass)""#
                    if options.includeProvenance, !value.evidenceSource.isEmpty {
                        annotation += #" evidenceSource="\#(VivoXMLCodec.escapedAttribute(value.evidenceSource))""#
                    }
                    annotation += "/>"
                    lines.append(annotation)
                    lines.append("        </annotation>")
                    lines.append("      </parameter>")
                } else {
                    lines.append("      <parameter \(attributes)/>")
                }
            }
            lines.append("    </listOfParameters>")
        }

        if !reactions.isEmpty {
            lines.append("    <listOfReactions>")
            for reaction in reactions {
                guard let identifier = reactionIDs[reaction.identifier],
                      let compartment = compartmentIDs[reaction.compartment] else {
                    throw VivoArtifactValidationError.unresolved("SBML identifier map is incomplete for reaction \(reaction.identifier)")
                }
                let reversible = reaction.rateLaw == .reversibleMassAction
                lines.append(#"      <reaction id="\#(VivoXMLCodec.escapedAttribute(identifier))" name="\#(VivoXMLCodec.escapedAttribute(reaction.identifier))" compartment="\#(VivoXMLCodec.escapedAttribute(compartment))" reversible="\#(reversible ? "true" : "false")" fast="false">"#)
                appendSpeciesReferences(
                    reaction.reactants,
                    element: "listOfReactants",
                    speciesIDs: speciesIDs,
                    indentation: "        ",
                    into: &lines
                )
                appendSpeciesReferences(
                    reaction.products,
                    element: "listOfProducts",
                    speciesIDs: speciesIDs,
                    indentation: "        ",
                    into: &lines
                )

                if let math = kineticMath(
                    reaction: reaction,
                    speciesIDs: speciesIDs,
                    parameterIDs: parameterIDs
                ) {
                    lines.append("        <kineticLaw>")
                    lines.append(#"          <math xmlns="http://www.w3.org/1998/Math/MathML">"#)
                    lines.append("            \(math)")
                    lines.append("          </math>")
                    lines.append("        </kineticLaw>")
                } else {
                    diagnostics.append(.init(
                        severity: .warning,
                        code: "NVSBML001",
                        subject: reaction.identifier,
                        message: "The \(reaction.rateLaw) rate law is retained as a NumiVivo annotation because SBML Core export cannot reproduce its runtime semantics without additional context."
                    ))
                }

                if options.includeNumiVivoAnnotations {
                    lines.append("        <annotation>")
                    lines.append(#"          <numivivo:reaction xmlns:numivivo="https://numivivo.org/ns#" rateLaw="\#(reaction.rateLaw.rawValue)" flags="\#(reaction.flags)" delaySeconds="\#(VivoXMLCodec.finite(reaction.delaySeconds))" characteristicRate="\#(VivoXMLCodec.finite(reaction.characteristicRate))" cohortIndex="\#(reaction.cohortIndex)" expressionOffset="\#(reaction.expressionOffset)" expressionCount="\#(reaction.expressionCount)"\#(reaction.gateExpressionOffset.map { " gateExpressionOffset=\"\($0)\"" } ?? "")/>"#)
                    lines.append("        </annotation>")
                }
                lines.append("      </reaction>")

                if reaction.hasGate {
                    diagnostics.append(.init(
                        severity: .warning,
                        code: "NVSBML002",
                        subject: reaction.identifier,
                        message: "Reaction gate bytecode is not represented as an SBML Core construct; the ProgramPack offset is preserved in annotation."
                    ))
                }
                if reaction.isDelayed {
                    diagnostics.append(.init(
                        severity: .warning,
                        code: "NVSBML003",
                        subject: reaction.identifier,
                        message: "Transactional delayed-product scheduling is not equivalent to a plain SBML delay and is preserved only as annotation."
                    ))
                }
            }
            lines.append("    </listOfReactions>")
        }

        if !rules.isEmpty || !monitors.isEmpty {
            diagnostics.append(.init(
                severity: .warning,
                code: "NVSBML004",
                subject: modelID,
                message: "NumiVivo behavioral rules, temporal state, safety monitors, and transactional responses are summarized in model annotations rather than converted to SBML rules."
            ))
        }
        diagnostics.append(.init(
            severity: .note,
            code: "NVSBML005",
            subject: modelID,
            message: "Units are retained as NumiVivo annotations. They are not emitted as inferred SBML unit definitions because context-dependent count and concentration conversions require explicit compartment volumes."
        ))

        lines.append("  </model>")
        lines.append("</sbml>")
        let text = lines.joined(separator: "\n") + "\n"
        guard let data = text.data(using: .utf8) else {
            throw VivoArtifactValidationError.invalid("SBML serialization did not produce UTF-8")
        }
        return VivoSBML3Export(
            xml: data,
            fingerprint: try VivoCanonicalJSON.fingerprint(data),
            identifiers: identifiers,
            diagnostics: diagnostics
        )
    }

    private func appendSpeciesReferences(
        _ terms: [VivoProgramPack.StoichiometryMetadata],
        element: String,
        speciesIDs: [String: String],
        indentation: String,
        into lines: inout [String]
    ) {
        guard !terms.isEmpty else { return }
        lines.append("\(indentation)<\(element)>")
        for term in terms {
            guard let identifier = speciesIDs[term.speciesIdentifier] else { continue }
            let coefficient = max(Int(term.coefficient), 0)
            lines.append(#"\#(indentation)  <speciesReference species="\#(VivoXMLCodec.escapedAttribute(identifier))" stoichiometry="\#(coefficient)" constant="true"/>"#)
        }
        lines.append("\(indentation)</\(element)>")
    }

    private func kineticMath(
        reaction: VivoProgramPack.ReactionMetadata,
        speciesIDs: [String: String],
        parameterIDs: [String: String]
    ) -> String? {
        let parameters = reaction.parameterIdentifiers.compactMap { parameterIDs[$0] }
        let reactants = reaction.reactants.compactMap { term -> (String, Int)? in
            guard let identifier = speciesIDs[term.speciesIdentifier] else { return nil }
            return (identifier, max(Int(term.coefficient), 1))
        }
        let products = reaction.products.compactMap { term -> (String, Int)? in
            guard let identifier = speciesIDs[term.speciesIdentifier] else { return nil }
            return (identifier, max(Int(term.coefficient), 1))
        }

        switch reaction.rateLaw {
        case .zeroOrder:
            guard let rate = parameters.first else { return nil }
            return ci(rate)
        case .massAction, .degradation:
            guard let rate = parameters.first else { return nil }
            return times([ci(rate)] + reactants.map { poweredSpecies($0.0, exponent: $0.1) })
        case .reversibleMassAction:
            guard parameters.count >= 2 else { return nil }
            let forward = times([ci(parameters[0])] + reactants.map { poweredSpecies($0.0, exponent: $0.1) })
            let reverse = times([ci(parameters[1])] + products.map { poweredSpecies($0.0, exponent: $0.1) })
            return apply("minus", [forward, reverse])
        case .michaelisMenten:
            guard parameters.count >= 2, let substrate = reactants.first else { return nil }
            let s = ci(substrate.0)
            return apply("divide", [times([ci(parameters[0]), s]), apply("plus", [ci(parameters[1]), s])])
        case .hillActivation:
            guard parameters.count >= 3, let substrate = reactants.first else { return nil }
            let signal = ci(substrate.0)
            let signalPower = apply("power", [signal, ci(parameters[2])])
            let thresholdPower = apply("power", [ci(parameters[1]), ci(parameters[2])])
            return apply("divide", [times([ci(parameters[0]), signalPower]), apply("plus", [thresholdPower, signalPower])])
        case .hillRepression:
            guard parameters.count >= 3, let substrate = reactants.first else { return nil }
            let signalPower = apply("power", [ci(substrate.0), ci(parameters[2])])
            let thresholdPower = apply("power", [ci(parameters[1]), ci(parameters[2])])
            return apply("divide", [times([ci(parameters[0]), thresholdPower]), apply("plus", [thresholdPower, signalPower])])
        case .passiveTransport, .saturableTransport, .customBytecode:
            return nil
        }
    }

    private func poweredSpecies(_ identifier: String, exponent: Int) -> String {
        exponent == 1 ? ci(identifier) : apply("power", [ci(identifier), cn(Double(exponent))])
    }

    private func times(_ operands: [String]) -> String {
        switch operands.count {
        case 0: cn(1)
        case 1: operands[0]
        default: apply("times", operands)
        }
    }

    private func apply(_ operatorName: String, _ operands: [String]) -> String {
        "<apply><\(operatorName)/>\(operands.joined())</apply>"
    }

    private func ci(_ identifier: String) -> String {
        "<ci>\(VivoXMLCodec.escapedText(identifier))</ci>"
    }

    private func cn(_ value: Double) -> String {
        "<cn>\(VivoXMLCodec.finite(value))</cn>"
    }
}
