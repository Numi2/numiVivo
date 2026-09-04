import Foundation

public struct VivoStandardsDiagnostic: Codable, Sendable, Equatable, Hashable {
    public enum Severity: String, Codable, Sendable, CaseIterable {
        case note
        case warning
        case error
    }

    public let severity: Severity
    public let code: String
    public let subject: String
    public let message: String

    public init(
        severity: Severity,
        code: String,
        subject: String,
        message: String
    ) {
        self.severity = severity
        self.code = code
        self.subject = subject
        self.message = message
    }
}

public struct VivoStandardDocument: Codable, Sendable, Equatable {
    public enum Format: String, Codable, Sendable, CaseIterable {
        case sbol3Turtle
        case sbmlLevel3Version2
        case sedmlLevel1Version4
        case combineManifest
    }

    public let format: Format
    public let mediaType: String
    public let data: Data
    public let fingerprint: VivoFingerprint
    public let diagnostics: [VivoStandardsDiagnostic]

    public init(
        format: Format,
        mediaType: String,
        data: Data,
        fingerprint: VivoFingerprint,
        diagnostics: [VivoStandardsDiagnostic]
    ) {
        self.format = format
        self.mediaType = mediaType
        self.data = data
        self.fingerprint = fingerprint
        self.diagnostics = diagnostics
    }
}

public struct VivoIdentifierMap: Codable, Sendable, Equatable {
    public let model: String
    public let compartments: [String: String]
    public let species: [String: String]
    public let parameters: [String: String]
    public let reactions: [String: String]

    public init(
        model: String,
        compartments: [String: String],
        species: [String: String],
        parameters: [String: String],
        reactions: [String: String]
    ) {
        self.model = model
        self.compartments = compartments
        self.species = species
        self.parameters = parameters
        self.reactions = reactions
    }
}

enum VivoXMLCodec {
    static func escapedText(_ value: String) -> String {
        value
            .replacingOccurrences(of: "&", with: "&amp;")
            .replacingOccurrences(of: "<", with: "&lt;")
            .replacingOccurrences(of: ">", with: "&gt;")
    }

    static func escapedAttribute(_ value: String) -> String {
        escapedText(value)
            .replacingOccurrences(of: "\"", with: "&quot;")
            .replacingOccurrences(of: "'", with: "&apos;")
    }

    static func identifier(_ source: String, fallback: String = "item") -> String {
        let scalars = source.unicodeScalars
        var result = ""
        result.reserveCapacity(source.utf8.count)
        for scalar in scalars {
            let allowed = CharacterSet.alphanumerics.contains(scalar) || scalar.value == 95
            result.append(allowed ? Character(String(scalar)) : "_")
        }
        if result.isEmpty { result = fallback }
        if let first = result.unicodeScalars.first,
           !CharacterSet.letters.contains(first),
           first.value != 95 {
            result = "_" + result
        }
        return result
    }

    static func stableIdentifiers(_ source: [String], prefix: String) -> [String: String] {
        var occupied = Set<String>()
        var result: [String: String] = [:]
        result.reserveCapacity(source.count)
        for (index, value) in source.enumerated() {
            let base = identifier(value, fallback: "\(prefix)_\(index)")
            var candidate = base
            var suffix = 1
            while occupied.contains(candidate) {
                candidate = "\(base)_\(suffix)"
                suffix += 1
            }
            occupied.insert(candidate)
            result[value] = candidate
        }
        return result
    }

    static func finite(_ value: Double) -> String {
        let safe = value.isFinite ? value : 0
        return String(format: "%.17g", locale: Locale(identifier: "en_US_POSIX"), safe)
    }

    static func finite(_ value: Float) -> String {
        finite(Double(value))
    }
}
