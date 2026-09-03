import CryptoKit
import Foundation

public enum VivoCompositionConflictPolicy: String, Codable, Sendable {
    case reject
    case requireSameValue
    case lastWriterWins
}

public struct VivoContextMetadata: Codable, Sendable, Hashable {
    public var id: String
    public var version: String
    public var description: String
    public var authors: [String]
    public var labels: [String: String]

    public init(
        id: String,
        version: String,
        description: String = "",
        authors: [String] = [],
        labels: [String: String] = [:]
    ) {
        self.id = id
        self.version = version
        self.description = description
        self.authors = authors
        self.labels = labels
    }
}

public struct VivoEvidenceReference: Codable, Sendable, Hashable {
    public var classification: String
    public var source: String
    public var detail: String
    public var sourceFingerprint: String?

    public init(
        classification: String,
        source: String,
        detail: String = "",
        sourceFingerprint: String? = nil
    ) {
        self.classification = classification
        self.source = source
        self.detail = detail
        self.sourceFingerprint = sourceFingerprint
    }
}

public struct VivoNumericOverride: Codable, Sendable, Hashable {
    public var id: String
    public var value: Double
    public var unit: String?
    public var minimum: Double?
    public var maximum: Double?
    public var evidence: VivoEvidenceReference?

    public init(
        id: String,
        value: Double,
        unit: String? = nil,
        minimum: Double? = nil,
        maximum: Double? = nil,
        evidence: VivoEvidenceReference? = nil
    ) {
        self.id = id
        self.value = value
        self.unit = unit
        self.minimum = minimum
        self.maximum = maximum
        self.evidence = evidence
    }
}

public struct VivoContextConstraint: Codable, Sendable, Hashable {
    public var id: String
    public var expression: VivoJSONValue
    public var severity: String
    public var response: String
    public var message: String

    public init(
        id: String,
        expression: VivoJSONValue,
        severity: String = "error",
        response: String = "reject-step",
        message: String
    ) {
        self.id = id
        self.expression = expression
        self.severity = severity
        self.response = response
        self.message = message
    }
}

public struct VivoContextPack: Codable, Sendable, Hashable {
    public var apiVersion: String
    public var kind: String
    public var metadata: VivoContextMetadata
    public var expectedProgramSourceFingerprint: String?
    public var targetMatch: [String: VivoJSONValue]
    public var parameterOverrides: [VivoNumericOverride]
    public var speciesInitialOverrides: [VivoNumericOverride]
    public var stateInitialOverrides: [VivoNumericOverride]
    public var inputDefaultOverrides: [VivoNumericOverride]
    public var addedConstraints: [VivoContextConstraint]
    public var targetAnnotations: [String: VivoJSONValue]

    public init(
        apiVersion: String = "numivivo.org/v1alpha1",
        kind: String = "VivoContextPack",
        metadata: VivoContextMetadata,
        expectedProgramSourceFingerprint: String? = nil,
        targetMatch: [String: VivoJSONValue] = [:],
        parameterOverrides: [VivoNumericOverride] = [],
        speciesInitialOverrides: [VivoNumericOverride] = [],
        stateInitialOverrides: [VivoNumericOverride] = [],
        inputDefaultOverrides: [VivoNumericOverride] = [],
        addedConstraints: [VivoContextConstraint] = [],
        targetAnnotations: [String: VivoJSONValue] = [:]
    ) {
        self.apiVersion = apiVersion
        self.kind = kind
        self.metadata = metadata
        self.expectedProgramSourceFingerprint = expectedProgramSourceFingerprint
        self.targetMatch = targetMatch
        self.parameterOverrides = parameterOverrides
        self.speciesInitialOverrides = speciesInitialOverrides
        self.stateInitialOverrides = stateInitialOverrides
        self.inputDefaultOverrides = inputDefaultOverrides
        self.addedConstraints = addedConstraints
        self.targetAnnotations = targetAnnotations
    }
}

public struct VivoAppliedContext: Codable, Sendable, Hashable {
    public var id: String
    public var version: String
    public var contentSHA256: String
}

public struct VivoProgramCompositionReceipt: Codable, Sendable, Hashable {
    public static let schema = "numivivo.org/composition-receipt/v1"

    public var schema: String
    public var baseSourceSHA256: String
    public var appliedContexts: [VivoAppliedContext]
    public var conflictPolicy: VivoCompositionConflictPolicy
    public var outputSourceSHA256: String
}

public struct VivoComposedProgram: Sendable, Hashable {
    public var sourceJSON: Data
    public var receipt: VivoProgramCompositionReceipt
}

public enum VivoCompositionError: Error, Sendable, CustomStringConvertible {
    case invalidProgram(String)
    case invalidContext(String)
    case fingerprintMismatch(context: String, expected: String, actual: String)
    case targetMismatch(context: String, key: String)
    case missingRecord(section: String, id: String)
    case unitMismatch(section: String, id: String, expected: String, actual: String)
    case conflict(path: String, firstContext: String, secondContext: String)
    case duplicateConstraint(String)

    public var description: String {
        switch self {
        case .invalidProgram(let message), .invalidContext(let message):
            return message
        case .fingerprintMismatch(let context, let expected, let actual):
            return "Context '\(context)' requires source \(expected), but the base source is \(actual)."
        case .targetMismatch(let context, let key):
            return "Context '\(context)' does not match target field '\(key)'."
        case .missingRecord(let section, let id):
            return "No record '\(id)' exists in program section '\(section)'."
        case .unitMismatch(let section, let id, let expected, let actual):
            return "Override \(section).\(id) uses unit '\(actual)', but the program declares '\(expected)'."
        case .conflict(let path, let first, let second):
            return "Contexts '\(first)' and '\(second)' both modify '\(path)' with conflicting values."
        case .duplicateConstraint(let id):
            return "Constraint '\(id)' already exists in the program or another context."
        }
    }
}

public enum VivoProgramComposer {
    public static func decodeContext(_ data: Data) throws -> VivoContextPack {
        let decoder = JSONDecoder()
        let context = try decoder.decode(VivoContextPack.self, from: data)
        guard context.apiVersion == "numivivo.org/v1alpha1" else {
            throw VivoCompositionError.invalidContext(
                "Unsupported context apiVersion '\(context.apiVersion)'."
            )
        }
        guard context.kind == "VivoContextPack" else {
            throw VivoCompositionError.invalidContext(
                "Expected kind VivoContextPack, found '\(context.kind)'."
            )
        }
        guard !context.metadata.id.isEmpty, !context.metadata.version.isEmpty else {
            throw VivoCompositionError.invalidContext(
                "Context metadata.id and metadata.version are required."
            )
        }
        return context
    }

    public static func compose(
        programJSON: Data,
        contextJSON: [Data],
        conflictPolicy: VivoCompositionConflictPolicy = .reject
    ) throws -> VivoComposedProgram {
        try compose(
            programJSON: programJSON,
            contexts: contextJSON.map(decodeContext),
            conflictPolicy: conflictPolicy
        )
    }

    public static func compose(
        programJSON: Data,
        contexts: [VivoContextPack],
        conflictPolicy: VivoCompositionConflictPolicy = .reject
    ) throws -> VivoComposedProgram {
        guard var root = try JSONSerialization.jsonObject(with: programJSON) as? [String: Any],
              var specification = root["spec"] as? [String: Any] else {
            throw VivoCompositionError.invalidProgram(
                "VivoProgram source must be a JSON object containing a spec object."
            )
        }
        guard root["kind"] as? String == "VivoProgram" else {
            throw VivoCompositionError.invalidProgram("Source kind must be VivoProgram.")
        }
        guard root["apiVersion"] as? String == "numivivo.org/v1alpha1" else {
            throw VivoCompositionError.invalidProgram(
                "Source apiVersion must be numivivo.org/v1alpha1."
            )
        }

        let baseHash = sha256(programJSON)
        var writes: [String: (context: String, value: Any)] = [:]
        var applied: [VivoAppliedContext] = []
        applied.reserveCapacity(contexts.count)

        for context in contexts {
            guard context.apiVersion == "numivivo.org/v1alpha1",
                  context.kind == "VivoContextPack",
                  !context.metadata.id.isEmpty,
                  !context.metadata.version.isEmpty else {
                throw VivoCompositionError.invalidContext(
                    "Context apiVersion, kind, id, or version is invalid."
                )
            }
            if let expected = context.expectedProgramSourceFingerprint,
               expected.lowercased() != baseHash {
                throw VivoCompositionError.fingerprintMismatch(
                    context: context.metadata.id,
                    expected: expected,
                    actual: baseHash
                )
            }

            var target = specification["target"] as? [String: Any] ?? [:]
            for (key, expected) in context.targetMatch {
                guard let actual = target[key], jsonEquivalent(actual, expected.foundationValue) else {
                    throw VivoCompositionError.targetMismatch(
                        context: context.metadata.id,
                        key: key
                    )
                }
            }
            for key in context.targetAnnotations.keys.sorted() {
                guard let value = context.targetAnnotations[key] else { continue }
                try registerWrite(
                    path: "spec.target.\(key)",
                    value: value.foundationValue,
                    context: context.metadata.id,
                    policy: conflictPolicy,
                    writes: &writes
                )
                target[key] = value.foundationValue
            }
            specification["target"] = target

            try apply(
                context.parameterOverrides,
                section: "parameters",
                valueField: "value",
                context: context.metadata.id,
                policy: conflictPolicy,
                specification: &specification,
                writes: &writes
            )
            try apply(
                context.speciesInitialOverrides,
                section: "species",
                valueField: "initial",
                context: context.metadata.id,
                policy: conflictPolicy,
                specification: &specification,
                writes: &writes
            )
            try apply(
                context.stateInitialOverrides,
                section: "state",
                valueField: "initial",
                context: context.metadata.id,
                policy: conflictPolicy,
                specification: &specification,
                writes: &writes
            )
            try apply(
                context.inputDefaultOverrides,
                section: "inputs",
                valueField: "default",
                context: context.metadata.id,
                policy: conflictPolicy,
                specification: &specification,
                writes: &writes
            )

            var constraints = specification["constraints"] as? [[String: Any]] ?? []
            var constraintIDs = Set(
                constraints.compactMap { $0["id"] as? String }
            )
            for constraint in context.addedConstraints {
                guard constraintIDs.insert(constraint.id).inserted else {
                    throw VivoCompositionError.duplicateConstraint(constraint.id)
                }
                constraints.append([
                    "id": constraint.id,
                    "expression": constraint.expression.foundationValue,
                    "severity": constraint.severity,
                    "response": constraint.response,
                    "message": constraint.message
                ])
            }
            specification["constraints"] = constraints

            let contextData = try JSONEncoder.nvivoCanonical.encode(context)
            applied.append(
                VivoAppliedContext(
                    id: context.metadata.id,
                    version: context.metadata.version,
                    contentSHA256: sha256(contextData)
                )
            )
        }

        root["spec"] = specification
        var metadata = root["metadata"] as? [String: Any] ?? [:]
        var labels = metadata["labels"] as? [String: String] ?? [:]
        labels["numivivo.contexts"] = applied.map { "\($0.id)@\($0.version)" }.joined(separator: ",")
        metadata["labels"] = labels
        root["metadata"] = metadata

        guard JSONSerialization.isValidJSONObject(root) else {
            throw VivoCompositionError.invalidProgram(
                "Context composition produced an invalid JSON object."
            )
        }
        let output = try JSONSerialization.data(
            withJSONObject: root,
            options: [.sortedKeys, .withoutEscapingSlashes]
        )
        return VivoComposedProgram(
            sourceJSON: output,
            receipt: VivoProgramCompositionReceipt(
                schema: VivoProgramCompositionReceipt.schema,
                baseSourceSHA256: baseHash,
                appliedContexts: applied,
                conflictPolicy: conflictPolicy,
                outputSourceSHA256: sha256(output)
            )
        )
    }

    private static func apply(
        _ overrides: [VivoNumericOverride],
        section: String,
        valueField: String,
        context: String,
        policy: VivoCompositionConflictPolicy,
        specification: inout [String: Any],
        writes: inout [String: (context: String, value: Any)]
    ) throws {
        guard !overrides.isEmpty else { return }
        guard var records = specification[section] as? [[String: Any]] else {
            throw VivoCompositionError.invalidProgram(
                "Program does not contain a '\(section)' array."
            )
        }
        var indices: [String: Int] = [:]
        for (index, record) in records.enumerated() {
            if let id = record["id"] as? String { indices[id] = index }
        }

        for override in overrides {
            guard override.value.isFinite,
                  override.minimum?.isFinite ?? true,
                  override.maximum?.isFinite ?? true else {
                throw VivoCompositionError.invalidContext(
                    "Override '\(override.id)' contains a non-finite value or bound."
                )
            }
            if let minimum = override.minimum, let maximum = override.maximum,
               minimum > override.value || override.value > maximum || minimum > maximum {
                throw VivoCompositionError.invalidContext(
                    "Override '\(override.id)' value and bounds are not ordered."
                )
            }
            guard let index = indices[override.id] else {
                throw VivoCompositionError.missingRecord(section: section, id: override.id)
            }
            var record = records[index]
            if let overrideUnit = override.unit,
               let declaredUnit = record["unit"] as? String,
               overrideUnit != declaredUnit {
                throw VivoCompositionError.unitMismatch(
                    section: section,
                    id: override.id,
                    expected: declaredUnit,
                    actual: overrideUnit
                )
            }

            try registerWrite(
                path: "spec.\(section).\(override.id).\(valueField)",
                value: override.value,
                context: context,
                policy: policy,
                writes: &writes
            )
            record[valueField] = override.value
            if override.minimum != nil || override.maximum != nil {
                var bounds = record["bounds"] as? [String: Any] ?? [:]
                if let minimum = override.minimum {
                    try registerWrite(
                        path: "spec.\(section).\(override.id).bounds.minimum",
                        value: minimum,
                        context: context,
                        policy: policy,
                        writes: &writes
                    )
                    bounds["minimum"] = minimum
                }
                if let maximum = override.maximum {
                    try registerWrite(
                        path: "spec.\(section).\(override.id).bounds.maximum",
                        value: maximum,
                        context: context,
                        policy: policy,
                        writes: &writes
                    )
                    bounds["maximum"] = maximum
                }
                record["bounds"] = bounds
            }
            if let evidence = override.evidence {
                var evidenceObject: [String: Any] = [
                    "class": evidence.classification,
                    "source": evidence.source,
                    "detail": evidence.detail
                ]
                if let sourceFingerprint = evidence.sourceFingerprint {
                    evidenceObject["sourceFingerprint"] = sourceFingerprint
                }
                record["evidence"] = evidenceObject
            }
            records[index] = record
        }
        specification[section] = records
    }

    private static func registerWrite(
        path: String,
        value: Any,
        context: String,
        policy: VivoCompositionConflictPolicy,
        writes: inout [String: (context: String, value: Any)]
    ) throws {
        guard let previous = writes[path] else {
            writes[path] = (context, value)
            return
        }
        switch policy {
        case .reject:
            throw VivoCompositionError.conflict(
                path: path,
                firstContext: previous.context,
                secondContext: context
            )
        case .requireSameValue:
            guard jsonEquivalent(previous.value, value) else {
                throw VivoCompositionError.conflict(
                    path: path,
                    firstContext: previous.context,
                    secondContext: context
                )
            }
        case .lastWriterWins:
            writes[path] = (context, value)
        }
    }

    private static func jsonEquivalent(_ left: Any, _ right: Any) -> Bool {
        let leftObject = ["value": left]
        let rightObject = ["value": right]
        guard JSONSerialization.isValidJSONObject(leftObject),
              JSONSerialization.isValidJSONObject(rightObject),
              let leftData = try? JSONSerialization.data(
                withJSONObject: leftObject,
                options: [.sortedKeys, .withoutEscapingSlashes]
              ),
              let rightData = try? JSONSerialization.data(
                withJSONObject: rightObject,
                options: [.sortedKeys, .withoutEscapingSlashes]
              ) else {
            return false
        }
        return leftData == rightData
    }

    private static func sha256(_ data: Data) -> String {
        let digest = SHA256.hash(data: data)
        let alphabet = Array("0123456789abcdef".utf8)
        var output = [UInt8]()
        output.reserveCapacity(64)
        for byte in digest {
            output.append(alphabet[Int(byte >> 4)])
            output.append(alphabet[Int(byte & 0x0f)])
        }
        return String(decoding: output, as: UTF8.self)
    }
}

private extension JSONEncoder {
    static var nvivoCanonical: JSONEncoder {
        let encoder = JSONEncoder()
        encoder.outputFormatting = [.sortedKeys, .withoutEscapingSlashes]
        return encoder
    }
}
