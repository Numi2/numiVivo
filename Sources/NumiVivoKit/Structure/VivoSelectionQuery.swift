import Foundation

/// Small deterministic selection language used by CLI, QM-region definitions and
/// later MD analysis. Function syntax avoids the precedence ambiguities common in
/// ad-hoc textual atom selections.
///
/// Examples:
///   `element(C,N,O,S)`
///   `and(chain(A),resseq(450,500),not(hydrogen(true)))`
///   `within(0.35,resname(LIG),false)`
public enum VivoSelectionQuery {
    public static func parse(_ source: String) throws -> VivoAtomSelection {
        var parser = Parser(source)
        let selection = try parser.selection()
        guard parser.peek() == nil else {
            throw VivoArtifactValidationError.invalid("unexpected selection token '\(parser.peek()!)'")
        }
        return selection
    }

    private struct Parser {
        var tokens: [String]
        var cursor = 0

        init(_ source: String) {
            var values: [String] = []
            var current = ""
            for ch in source {
                if ch == "(" || ch == ")" || ch == "," {
                    if !current.isEmpty { values.append(current); current = "" }
                    values.append(String(ch))
                } else if ch.isWhitespace {
                    if !current.isEmpty { values.append(current); current = "" }
                } else {
                    current.append(ch)
                }
            }
            if !current.isEmpty { values.append(current) }
            tokens = values
        }

        func peek() -> String? { cursor < tokens.count ? tokens[cursor] : nil }

        mutating func take() throws -> String {
            guard cursor < tokens.count else { throw VivoArtifactValidationError.invalid("unexpected end of atom selection") }
            defer { cursor += 1 }
            return tokens[cursor]
        }

        mutating func expect(_ value: String) throws {
            let actual = try take()
            guard actual == value else { throw VivoArtifactValidationError.invalid("expected '\(value)' in atom selection, found '\(actual)'") }
        }

        mutating func selection() throws -> VivoAtomSelection {
            let function = try take().lowercased()
            if function == "all" { return .all }
            if function == "none" { return .none }
            if function == "protein" { return .residueNames(Self.proteinResidues) }
            if function == "water" { return .residueNames(["HOH", "WAT", "SOL", "TIP3", "TIP3P", "OPC"]) }
            try expect("(")
            switch function {
            case "index", "indices":
                let values = try scalarList().map { token -> UInt32 in
                    guard let value = UInt32(token) else { throw VivoArtifactValidationError.invalid("atom index '\(token)' is not UInt32") }
                    return value
                }
                return .atomIndices(values)
            case "element", "elements":
                let values = try scalarList().map { token -> UInt16 in
                    guard let element = VivoElement.from(symbol: token) else {
                        throw VivoArtifactValidationError.invalid("unknown element '\(token)' in selection")
                    }
                    return element.atomicNumber
                }
                return .elements(values)
            case "name", "atomname": return .atomNames(try scalarList())
            case "resname", "residue": return .residueNames(try scalarList())
            case "chain": return .chainIdentifiers(try scalarList())
            case "resseq":
                let values = try scalarList()
                guard values.count == 2, let lower = Int32(values[0]), let upper = Int32(values[1]) else {
                    throw VivoArtifactValidationError.invalid("resseq(lower,upper) requires two Int32 values")
                }
                return .residueSequenceRange(lower, upper)
            case "hetero": return .hetero(try booleanArgument())
            case "hydrogen": return .hydrogen(try booleanArgument())
            case "bonded":
                let inner = try selection(); try expect(")"); return .bonded(to: inner)
            case "not":
                let inner = try selection(); try expect(")"); return .not(inner)
            case "and": return .and(try selectionList())
            case "or": return .or(try selectionList())
            case "within":
                let distanceToken = try take()
                guard let distance = Double(distanceToken), distance.isFinite, distance >= 0 else {
                    throw VivoArtifactValidationError.invalid("within distance must be finite nonnegative nanometres")
                }
                try expect(",")
                let inner = try selection()
                var include = false
                if peek() == "," { _ = try take(); include = try parseBoolean(try take()) }
                try expect(")")
                return .within(distanceNM: distance, of: inner, includeSource: include)
            default:
                throw VivoArtifactValidationError.invalid("unknown atom-selection function '\(function)'")
            }
        }

        mutating func scalarList() throws -> [String] {
            var result: [String] = []
            if peek() == ")" { _ = try take(); return result }
            while true {
                let value = try take()
                guard value != "(" && value != ")" && value != "," else {
                    throw VivoArtifactValidationError.invalid("expected scalar atom-selection argument")
                }
                result.append(value)
                if peek() == "," { _ = try take(); continue }
                try expect(")")
                return result
            }
        }

        mutating func selectionList() throws -> [VivoAtomSelection] {
            var result: [VivoAtomSelection] = []
            if peek() == ")" { _ = try take(); return result }
            while true {
                result.append(try selection())
                if peek() == "," { _ = try take(); continue }
                try expect(")")
                return result
            }
        }

        mutating func booleanArgument() throws -> Bool {
            let value = try parseBoolean(try take())
            try expect(")")
            return value
        }

        func parseBoolean(_ raw: String) throws -> Bool {
            switch raw.lowercased() {
            case "true", "yes", "1": true
            case "false", "no", "0": false
            default: throw VivoArtifactValidationError.invalid("expected boolean atom-selection argument")
            }
        }

        static let proteinResidues = [
            "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
            "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
            "HID", "HIE", "HIP", "ASH", "GLH", "CYX", "CYM", "LYN"
        ]
    }
}
