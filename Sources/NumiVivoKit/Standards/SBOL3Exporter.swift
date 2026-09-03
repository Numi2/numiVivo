import Foundation

public struct VivoSBOL3ExportOptions: Sendable, Codable, Equatable {
    public var baseURI: String
    public var componentIdentifier: String
    public var componentName: String?
    public var includeBehavioralLayer: Bool
    public var includeParameterProvenance: Bool

    public init(
        baseURI: String = "https://numivivo.org/designs",
        componentIdentifier: String = "numivivo_design",
        componentName: String? = nil,
        includeBehavioralLayer: Bool = true,
        includeParameterProvenance: Bool = true
    ) {
        self.baseURI = baseURI
        self.componentIdentifier = componentIdentifier
        self.componentName = componentName
        self.includeBehavioralLayer = includeBehavioralLayer
        self.includeParameterProvenance = includeParameterProvenance
    }
}

public struct VivoSBOL3Export: Sendable {
    public let turtle: Data
    public let fingerprint: VivoFingerprint
    public let diagnostics: [VivoStandardsDiagnostic]

    public init(turtle: Data, fingerprint: VivoFingerprint, diagnostics: [VivoStandardsDiagnostic]) {
        self.turtle = turtle
        self.fingerprint = fingerprint
        self.diagnostics = diagnostics
    }
}

public struct VivoSBOL3Exporter: Sendable {
    public init() {}

    public func export(
        _ pack: VivoProgramPack,
        options: VivoSBOL3ExportOptions = .init()
    ) throws -> VivoSBOL3Export {
        guard let baseURL = URL(string: options.baseURI),
              baseURL.scheme != nil,
              !options.componentIdentifier.isEmpty else {
            throw VivoArtifactValidationError.invalid("SBOL baseURI must be absolute and componentIdentifier must be non-empty")
        }

        let species = try pack.speciesMetadata()
        let parameters = try pack.parameterMetadata()
        let reactions = try pack.reactionMetadata()
        let rules = try pack.ruleMetadata()
        let actions = try pack.actionMetadata()
        let monitors = try pack.monitorMetadata()
        var diagnostics: [VivoStandardsDiagnostic] = []

        let base = options.baseURI.trimmingCharacters(in: CharacterSet(charactersIn: "/#"))
        let componentID = sbolToken(options.componentIdentifier)
        let componentIRI = "\(base)/\(componentID)"
        var writer = VivoTurtleWriter()
        writer.preamble()

        writer.subject(componentIRI, predicates: [
            ("a", ":Component"),
            (":displayId", writer.literal(componentID)),
            (":type", "<https://numivivo.org/type/InVivoMolecularProgram>"),
            (":name", writer.literal(options.componentName ?? options.componentIdentifier)),
            ("prov:wasDerivedFrom", "<urn:sha256:\(pack.header.sourceFingerprint.hex)>"),
            ("numivivo:programPackFingerprint", writer.literal(pack.header.contentFingerprint.hex)),
            ("numivivo:fidelity", writer.literal(pack.header.fidelity.label))
        ])

        var speciesIRIs: [String: String] = [:]
        for (index, value) in species.enumerated() {
            let iri = "\(componentIRI)/feature/\(sbolToken(value.identifier))_\(index)"
            speciesIRIs[value.identifier] = iri
            writer.add(componentIRI, predicate: ":hasFeature", object: "<\(iri)>")
            writer.subject(iri, predicates: [
                ("a", ":LocalSubComponent"),
                (":displayId", writer.literal(sbolToken(value.identifier))),
                (":name", writer.literal(value.identifier)),
                (":type", "<https://numivivo.org/type/MolecularState>"),
                ("numivivo:compartment", writer.literal(value.compartment)),
                ("numivivo:unit", writer.literal(value.unit)),
                ("numivivo:initialValue", writer.typedNumber(Double(value.initialValue))),
                ("numivivo:minimum", writer.typedNumber(Double(value.minimum))),
                ("numivivo:maximum", writer.typedNumber(Double(value.maximum))),
                ("numivivo:flags", writer.typedInteger(UInt64(value.flags)))
            ])
        }

        for reaction in reactions {
            let interactionIRI = "\(componentIRI)/interaction/\(sbolToken(reaction.identifier))_\(reaction.index)"
            writer.add(componentIRI, predicate: ":hasInteraction", object: "<\(interactionIRI)>")
            var predicates: [(String, String)] = [
                ("a", ":Interaction"),
                (":displayId", writer.literal(sbolToken(reaction.identifier))),
                (":name", writer.literal(reaction.identifier)),
                (":type", "<https://numivivo.org/interaction/\(reaction.rateLaw)>"),
                ("numivivo:rateLaw", writer.literal(String(describing: reaction.rateLaw))),
                ("numivivo:delaySeconds", writer.typedNumber(Double(reaction.delaySeconds))),
                ("numivivo:characteristicRate", writer.typedNumber(Double(reaction.characteristicRate))),
                ("numivivo:cohortIndex", writer.typedInteger(UInt64(reaction.cohortIndex))),
                ("numivivo:flags", writer.typedInteger(UInt64(reaction.flags)))
            ]
            if reaction.expressionCount > 0 {
                predicates.append(("numivivo:expressionOffset", writer.typedInteger(UInt64(reaction.expressionOffset))))
                predicates.append(("numivivo:expressionCount", writer.typedInteger(UInt64(reaction.expressionCount))))
            }
            if let gate = reaction.gateExpressionOffset {
                predicates.append(("numivivo:gateExpressionOffset", writer.typedInteger(UInt64(gate))))
            }
            writer.subject(interactionIRI, predicates: predicates)

            for (role, terms) in [("reactant", reaction.reactants), ("product", reaction.products)] {
                for (position, term) in terms.enumerated() {
                    guard let participant = speciesIRIs[term.speciesIdentifier] else {
                        throw VivoArtifactValidationError.unresolved("SBOL participant \(term.speciesIdentifier) is unresolved")
                    }
                    let participationIRI = "\(interactionIRI)/participation/\(role)_\(position)"
                    writer.add(interactionIRI, predicate: ":hasParticipation", object: "<\(participationIRI)>")
                    writer.subject(participationIRI, predicates: [
                        ("a", ":Participation"),
                        (":participant", "<\(participant)>"),
                        (":role", "<https://numivivo.org/role/\(role)>"),
                        ("numivivo:stoichiometry", writer.typedInteger(UInt64(max(term.coefficient, 0))))
                    ])
                }
            }
            for (position, identifier) in reaction.parameterIdentifiers.enumerated() {
                let parameterIRI = "\(componentIRI)/parameter/\(sbolToken(identifier))"
                writer.add(interactionIRI, predicate: "numivivo:parameter", object: "<\(parameterIRI)>")
                writer.add(interactionIRI, predicate: "numivivo:parameterPosition", object: writer.literal("\(position):\(identifier)"))
            }
        }

        if options.includeParameterProvenance {
            for parameter in parameters {
                let iri = "\(componentIRI)/parameter/\(sbolToken(parameter.identifier))"
                writer.subject(iri, predicates: [
                    ("a", "prov:Entity"),
                    ("dcterms:title", writer.literal(parameter.identifier)),
                    ("numivivo:unit", writer.literal(parameter.unit)),
                    ("numivivo:value", writer.typedNumber(parameter.value)),
                    ("numivivo:minimum", writer.typedNumber(parameter.minimum)),
                    ("numivivo:maximum", writer.typedNumber(parameter.maximum)),
                    ("numivivo:evidenceClass", writer.typedInteger(UInt64(parameter.evidenceClass))),
                    ("prov:wasDerivedFrom", parameter.evidenceSource.isEmpty
                        ? "<urn:numivivo:unspecified-evidence>"
                        : "<\(safeIRI(parameter.evidenceSource, fallback: "urn:numivivo:evidence:\(sbolToken(parameter.identifier))"))>")
                ])
            }
        }

        if options.includeBehavioralLayer {
            for rule in rules {
                let iri = "\(componentIRI)/rule/\(sbolToken(rule.identifier))_\(rule.index)"
                writer.add(componentIRI, predicate: "numivivo:hasRule", object: "<\(iri)>")
                writer.subject(iri, predicates: [
                    ("a", "numivivo:Rule"),
                    ("dcterms:title", writer.literal(rule.identifier)),
                    ("numivivo:conditionOffset", writer.typedInteger(UInt64(rule.conditionOffset))),
                    ("numivivo:conditionCount", writer.typedInteger(UInt64(rule.conditionCount))),
                    ("numivivo:actionOffset", writer.typedInteger(UInt64(rule.actionOffset))),
                    ("numivivo:actionCount", writer.typedInteger(UInt64(rule.actionCount))),
                    ("numivivo:priority", writer.typedSignedInteger(Int64(rule.priority))),
                    ("numivivo:refractorySeconds", writer.typedNumber(Double(rule.refractorySeconds)))
                ])
                let begin = Int(rule.actionOffset)
                let endResult = begin.addingReportingOverflow(Int(rule.actionCount))
                guard !endResult.overflow, begin >= 0, endResult.partialValue <= actions.count else {
                    throw VivoArtifactValidationError.invalid("rule \(rule.identifier) action range is out of bounds")
                }
                for action in actions[begin..<endResult.partialValue] {
                    let actionIRI = "\(iri)/action/\(action.index)"
                    writer.add(iri, predicate: "numivivo:hasAction", object: "<\(actionIRI)>")
                    writer.subject(actionIRI, predicates: [
                        ("a", "numivivo:Action"),
                        ("numivivo:kind", writer.literal(String(describing: action.kind))),
                        ("numivivo:targetIndex", writer.typedInteger(UInt64(action.targetIndex))),
                        ("numivivo:constantValue", writer.typedNumber(Double(action.constantValue))),
                        ("numivivo:maximumRate", writer.typedNumber(Double(action.maximumRate))),
                        ("numivivo:unit", writer.literal(action.unit))
                    ])
                }
            }

            for monitor in monitors {
                let iri = "\(componentIRI)/monitor/\(sbolToken(monitor.identifier))_\(monitor.index)"
                writer.add(componentIRI, predicate: "numivivo:hasMonitor", object: "<\(iri)>")
                writer.subject(iri, predicates: [
                    ("a", "numivivo:Monitor"),
                    ("dcterms:title", writer.literal(monitor.identifier)),
                    ("dcterms:description", writer.literal(monitor.message)),
                    ("numivivo:expressionOffset", writer.typedInteger(UInt64(monitor.expressionOffset))),
                    ("numivivo:expressionCount", writer.typedInteger(UInt64(monitor.expressionCount))),
                    ("numivivo:severity", writer.typedInteger(UInt64(monitor.severity))),
                    ("numivivo:response", writer.literal(String(describing: monitor.response))),
                    ("numivivo:isTermination", writer.typedBoolean(monitor.isTerminationRule))
                ])
            }
            diagnostics.append(.init(
                severity: .note,
                code: "NVSBOL001",
                subject: options.componentIdentifier,
                message: "Behavioral rules and safety monitors are represented as NumiVivo ontology extensions because SBOL does not prescribe their execution semantics."
            ))
        }

        diagnostics.append(.init(
            severity: .note,
            code: "NVSBOL002",
            subject: options.componentIdentifier,
            message: "This export is an abstract, sequence-free design graph. It does not assert DNA, RNA, or protein sequence realizability."
        ))

        let text = writer.render()
        guard let data = text.data(using: .utf8) else {
            throw VivoArtifactValidationError.invalid("SBOL Turtle serialization failed to produce UTF-8")
        }
        return VivoSBOL3Export(
            turtle: data,
            fingerprint: try VivoCanonicalJSON.fingerprint(data),
            diagnostics: diagnostics
        )
    }

    private func sbolToken(_ value: String) -> String {
        let scalars = value.unicodeScalars.map { scalar -> Character in
            if CharacterSet.alphanumerics.contains(scalar) || scalar.value == 95 || scalar.value == 45 {
                return Character(String(scalar))
            }
            return "_"
        }
        let token = String(scalars)
        return token.isEmpty ? "item" : token
    }

    private func safeIRI(_ candidate: String, fallback: String) -> String {
        guard let url = URL(string: candidate), url.scheme != nil else { return fallback }
        return candidate
            .replacingOccurrences(of: ">", with: "%3E")
            .replacingOccurrences(of: "<", with: "%3C")
            .replacingOccurrences(of: " ", with: "%20")
    }
}

private struct VivoTurtleWriter {
    private var lines: [String] = []

    mutating func preamble() {
        lines.append("@prefix : <http://sbols.org/v3#> .")
        lines.append("@prefix dcterms: <http://purl.org/dc/terms/> .")
        lines.append("@prefix prov: <http://www.w3.org/ns/prov#> .")
        lines.append("@prefix xsd: <http://www.w3.org/2001/XMLSchema#> .")
        lines.append("@prefix numivivo: <https://numivivo.org/ns#> .")
        lines.append("")
    }

    mutating func subject(_ iri: String, predicates: [(String, String)]) {
        guard !predicates.isEmpty else { return }
        lines.append("<\(iri)>")
        for (index, predicate) in predicates.enumerated() {
            let ending = index == predicates.count - 1 ? " ." : " ;"
            lines.append("    \(predicate.0) \(predicate.1)\(ending)")
        }
        lines.append("")
    }

    mutating func add(_ subject: String, predicate: String, object: String) {
        lines.append("<\(subject)> \(predicate) \(object) .")
    }

    func literal(_ value: String) -> String {
        let escaped = value
            .replacingOccurrences(of: "\\", with: "\\\\")
            .replacingOccurrences(of: "\"", with: "\\\"")
            .replacingOccurrences(of: "\n", with: "\\n")
            .replacingOccurrences(of: "\r", with: "\\r")
        return "\"\(escaped)\""
    }

    func typedNumber(_ value: Double) -> String {
        let finite = value.isFinite ? value : 0
        let formatted = String(format: "%.17g", locale: Locale(identifier: "en_US_POSIX"), finite)
        return "\"\(formatted)\"^^xsd:double"
    }

    func typedInteger(_ value: UInt64) -> String {
        "\"\(value)\"^^xsd:nonNegativeInteger"
    }

    func typedSignedInteger(_ value: Int64) -> String {
        "\"\(value)\"^^xsd:integer"
    }

    func typedBoolean(_ value: Bool) -> String {
        "\"\(value ? "true" : "false")\"^^xsd:boolean"
    }

    func render() -> String {
        lines.joined(separator: "\n") + "\n"
    }
}
