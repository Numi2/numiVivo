import Foundation

/// Lossless reader for the sectioned AMBER7+ parm7/prmtop container. The parser
/// honors each section's fixed-width Fortran format instead of whitespace-splitting
/// character fields. Scientific interpretation is performed by `VivoAmberImporter`.
public struct VivoAmberPrmtop: Sendable {
    public struct Format: Sendable, Equatable {
        public enum Kind: Sendable, Equatable { case integer, real, character }
        public var repeatCount: Int
        public var kind: Kind
        public var width: Int
        public var precision: Int?
    }

    public let source: Data
    public let versionLine: String?
    private let sections: [String: (Format, [String])]

    public init(data: Data) throws {
        guard let text = String(data: data, encoding: .utf8) else {
            throw VivoArtifactValidationError.invalid("AMBER prmtop is not UTF-8/ASCII text")
        }
        source = data
        let lines = text.split(separator: "\n", omittingEmptySubsequences: false).map(String.init)
        versionLine = lines.first(where: { $0.hasPrefix("%VERSION") })
        var parsed: [String: (Format, [String])] = [:]
        var i = 0
        while i < lines.count {
            guard lines[i].hasPrefix("%FLAG ") else { i += 1; continue }
            let name = String(lines[i].dropFirst(6)).trimmingCharacters(in: .whitespacesAndNewlines)
            guard !name.isEmpty, parsed[name] == nil else {
                throw VivoArtifactValidationError.invalid("AMBER prmtop contains empty or duplicate %FLAG '\(name)'")
            }
            i += 1
            // The AMBER specification permits any number of comments between
            // %FLAG and %FORMAT. They are descriptive only and never data.
            while i < lines.count, lines[i].hasPrefix("%COMMENT") { i += 1 }
            guard i < lines.count, lines[i].hasPrefix("%FORMAT") else {
                throw VivoArtifactValidationError.invalid("AMBER prmtop section \(name) has no %FORMAT")
            }
            let format = try Self.parseFormat(lines[i])
            i += 1
            var fields: [String] = []
            while i < lines.count, !lines[i].hasPrefix("%FLAG "), !lines[i].hasPrefix("%VERSION") {
                let line = lines[i]
                if line.hasPrefix("%COMMENT") { i += 1; continue }
                if line.hasPrefix("%") {
                    throw VivoArtifactValidationError.incompatible("unexpected AMBER directive inside section \(name): \(line)")
                }
                if !line.isEmpty {
                    let bytes = Array(line.utf8)
                    var offset = 0
                    while offset < bytes.count {
                        let end = Swift.min(offset + format.width, bytes.count)
                        let raw = String(decoding: bytes[offset..<end], as: UTF8.self)
                        if format.kind == .character || !raw.trimmingCharacters(in: .whitespaces).isEmpty {
                            fields.append(raw)
                        }
                        offset += format.width
                    }
                }
                i += 1
            }
            parsed[name] = (format, fields)
        }
        guard parsed["POINTERS"] != nil, parsed["ATOM_NAME"] != nil else {
            throw VivoArtifactValidationError.invalid("AMBER prmtop lacks required POINTERS/ATOM_NAME sections")
        }
        sections = parsed
    }

    public func contains(_ flag: String) -> Bool { sections[flag] != nil }

    public func strings(_ flag: String, required: Bool = true) throws -> [String] {
        guard let (format, values) = sections[flag] else {
            if required { throw VivoArtifactValidationError.unresolved("AMBER prmtop lacks %FLAG \(flag)") }
            return []
        }
        guard format.kind == .character else {
            throw VivoArtifactValidationError.invalid("AMBER prmtop \(flag) is not a character section")
        }
        return values.map { $0.trimmingCharacters(in: .whitespacesAndNewlines) }
    }

    public func integers(_ flag: String, required: Bool = true) throws -> [Int] {
        guard let (format, values) = sections[flag] else {
            if required { throw VivoArtifactValidationError.unresolved("AMBER prmtop lacks %FLAG \(flag)") }
            return []
        }
        guard format.kind == .integer else {
            throw VivoArtifactValidationError.invalid("AMBER prmtop \(flag) is not an integer section")
        }
        return try values.map { raw in
            guard let value = Int(raw.trimmingCharacters(in: .whitespacesAndNewlines)) else {
                throw VivoArtifactValidationError.invalid("AMBER prmtop \(flag) contains invalid integer '\(raw)'")
            }
            return value
        }
    }

    public func reals(_ flag: String, required: Bool = true) throws -> [Double] {
        guard let (format, values) = sections[flag] else {
            if required { throw VivoArtifactValidationError.unresolved("AMBER prmtop lacks %FLAG \(flag)") }
            return []
        }
        guard format.kind == .real else {
            throw VivoArtifactValidationError.invalid("AMBER prmtop \(flag) is not a real section")
        }
        return try values.map { raw in
            let normalized = raw.trimmingCharacters(in: .whitespacesAndNewlines)
                .replacingOccurrences(of: "D", with: "E")
                .replacingOccurrences(of: "d", with: "e")
            guard let value = Double(normalized), value.isFinite else {
                throw VivoArtifactValidationError.invalid("AMBER prmtop \(flag) contains invalid real '\(raw)'")
            }
            return value
        }
    }

    public func fingerprint() throws -> VivoFingerprint {
        try VivoCanonicalJSON.fingerprint(source)
    }

    private static func parseFormat(_ raw: String) throws -> Format {
        guard let open = raw.firstIndex(of: "("), let close = raw.lastIndex(of: ")"), open < close else {
            throw VivoArtifactValidationError.invalid("invalid AMBER %FORMAT line '\(raw)'")
        }
        let body = raw[raw.index(after: open)..<close].replacingOccurrences(of: " ", with: "")
        var cursor = body.startIndex
        var repeatDigits = ""
        while cursor < body.endIndex, body[cursor].isNumber {
            repeatDigits.append(body[cursor]); cursor = body.index(after: cursor)
        }
        guard let repeatCount = Int(repeatDigits), repeatCount > 0, cursor < body.endIndex else {
            throw VivoArtifactValidationError.invalid("unsupported AMBER %FORMAT '\(raw)'")
        }
        let code = body[cursor].uppercased()
        cursor = body.index(after: cursor)
        var widthDigits = ""
        while cursor < body.endIndex, body[cursor].isNumber {
            widthDigits.append(body[cursor]); cursor = body.index(after: cursor)
        }
        guard let width = Int(widthDigits), width > 0 else {
            throw VivoArtifactValidationError.invalid("AMBER %FORMAT has invalid width '\(raw)'")
        }
        var precision: Int?
        if cursor < body.endIndex, body[cursor] == "." {
            cursor = body.index(after: cursor)
            var digits = ""
            while cursor < body.endIndex, body[cursor].isNumber {
                digits.append(body[cursor]); cursor = body.index(after: cursor)
            }
            precision = Int(digits)
        }
        guard cursor == body.endIndex else {
            throw VivoArtifactValidationError.incompatible("unsupported AMBER %FORMAT modifier '\(raw)'")
        }
        let kind: Format.Kind
        switch code {
        case "I": kind = .integer
        case "E", "F", "D": kind = .real
        case "A": kind = .character
        default: throw VivoArtifactValidationError.incompatible("unsupported AMBER %FORMAT type '\(code)'")
        }
        return .init(repeatCount: repeatCount, kind: kind, width: width, precision: precision)
    }
}
