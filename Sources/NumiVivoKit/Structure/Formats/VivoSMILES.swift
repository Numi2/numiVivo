import Foundation

/// Native graph-level SMILES ingestion for Wave A. It preserves atom identity,
/// isotope/formal charge, branching, disconnected components, aromatic bonds and
/// ring closures. It rejects semantic features that the current structure IR
/// cannot represent exactly instead of silently flattening them.
public enum VivoSMILES {
    private struct ParsedAtom {
        var element: VivoElement
        var isotope: UInt16?
        var charge: Int16
        var aromatic: Bool
    }
    private struct PendingBond {
        var order: VivoBondOrder?
    }
    private struct RingOpen {
        var atom: UInt32
        var bond: PendingBond
    }

    public static func read(_ smiles: String, identifier: String = "smiles") throws -> VivoMolecularStructure {
        let chars = Array(smiles.trimmingCharacters(in: .whitespacesAndNewlines))
        guard !chars.isEmpty else { throw VivoArtifactValidationError.invalid("SMILES is empty") }
        let builder = VivoStructureAssemblyBuilder()
        var bonds: [VivoMolecularBond] = []
        var aromatic: [Bool] = []
        var current: UInt32?
        var branches: [UInt32] = []
        var rings: [String: RingOpen] = [:]
        var pending = PendingBond()
        var i = 0

        func append(_ atom: ParsedAtom) throws {
            let index = try builder.appendAtom(name: "\(atom.element.symbol)\(builder.atoms.count + 1)",
                                               element: atom.element, positionNM: .zero,
                                               sourceSerial: Int32(exactly: builder.atoms.count + 1),
                                               formalCharge: atom.charge, isHetero: true)
            builder.atoms[Int(index)].isotopeMassNumber = atom.isotope
            aromatic.append(atom.aromatic)
            if let previous = current {
                let defaultOrder: VivoBondOrder = aromatic[Int(previous)] && atom.aromatic ? .aromatic : .single
                bonds.append(.init(atomA: previous, atomB: index,
                                   order: pending.order ?? defaultOrder))
            }
            current = index
            pending = PendingBond()
        }

        while i < chars.count {
            let ch = chars[i]
            if ch.isWhitespace { i += 1; continue }
            switch ch {
            case "(":
                guard let current else { throw VivoArtifactValidationError.invalid("SMILES branch starts without an atom") }
                branches.append(current); i += 1
            case ")":
                guard let resumed = branches.popLast() else { throw VivoArtifactValidationError.invalid("SMILES has unmatched ')' ") }
                current = resumed; pending = PendingBond(); i += 1
            case ".":
                current = nil; pending = PendingBond(); i += 1
            case "-", "=", "#", ":":
                guard pending.order == nil else {
                    throw VivoArtifactValidationError.invalid("SMILES contains consecutive bond tokens")
                }
                if ch == "-" { pending.order = .single }
                else if ch == "=" { pending.order = .double }
                else if ch == "#" { pending.order = .triple }
                else { pending.order = .aromatic }
                i += 1
            case "/", "\\":
                throw VivoArtifactValidationError.incompatible("SMILES directional double-bond stereochemistry is not yet represented; import is rejected rather than made non-stereospecific")
            case "@":
                throw VivoArtifactValidationError.incompatible("SMILES tetrahedral stereochemistry requires the forthcoming stereo IR; import is rejected rather than made achiral")
            case "[":
                var j = i + 1
                var content = ""
                while j < chars.count, chars[j] != "]" { content.append(chars[j]); j += 1 }
                guard j < chars.count else { throw VivoArtifactValidationError.invalid("unterminated bracket atom in SMILES") }
                let atom = try parseBracket(content)
                try append(atom)
                i = j + 1
            case "%":
                guard i + 2 < chars.count, chars[i + 1].isNumber, chars[i + 2].isNumber else {
                    throw VivoArtifactValidationError.invalid("SMILES '%' ring closure requires two digits")
                }
                let label = String(chars[(i + 1)...(i + 2)])
                try closeOrOpenRing(label, current: current, pending: &pending, rings: &rings,
                                    bonds: &bonds, aromatic: aromatic)
                i += 3
            default:
                if ch.isNumber {
                    try closeOrOpenRing(String(ch), current: current, pending: &pending, rings: &rings,
                                        bonds: &bonds, aromatic: aromatic)
                    i += 1
                } else {
                    let parsed = try parseOrganic(chars, index: i)
                    try append(parsed.atom)
                    i += parsed.consumed
                }
            }
        }
        guard branches.isEmpty else { throw VivoArtifactValidationError.invalid("SMILES has unmatched '('") }
        guard rings.isEmpty else { throw VivoArtifactValidationError.invalid("SMILES has unclosed ring labels: \(rings.keys.sorted().joined(separator: ","))") }
        guard !builder.atoms.isEmpty else { throw VivoArtifactValidationError.invalid("SMILES contains no atoms") }
        var structure = try builder.makeStructure(identifier: identifier, bonds: bonds,
                                                  conformers: [], metadata: ["sourceFormat": "SMILES", "smiles": smiles])
        // Zero coordinates used while assembling are not published as a conformer.
        structure.conformers = []
        _ = try VivoStructureValidator.validate(structure)
        return structure
    }

    private static func closeOrOpenRing(_ label: String, current: UInt32?, pending: inout PendingBond,
                                        rings: inout [String: RingOpen], bonds: inout [VivoMolecularBond],
                                        aromatic: [Bool]) throws {
        guard let current else { throw VivoArtifactValidationError.invalid("SMILES ring label appears without an atom") }
        if let open = rings.removeValue(forKey: label) {
            guard open.atom != current else { throw VivoArtifactValidationError.invalid("SMILES ring cannot close an atom to itself") }
            if let a = open.bond.order, let b = pending.order, a != b {
                throw VivoArtifactValidationError.invalid("SMILES ring closure specifies conflicting bond orders")
            }
            let defaultOrder: VivoBondOrder = aromatic[Int(open.atom)] && aromatic[Int(current)] ? .aromatic : .single
            bonds.append(.init(atomA: open.atom, atomB: current,
                               order: pending.order ?? open.bond.order ?? defaultOrder))
        } else {
            rings[label] = .init(atom: current, bond: pending)
        }
        pending = PendingBond()
    }

    private static func parseOrganic(_ chars: [Character], index: Int) throws -> (atom: ParsedAtom, consumed: Int) {
        let ch = chars[index]
        let aromaticSymbols: [Character: String] = ["b": "B", "c": "C", "n": "N", "o": "O", "p": "P", "s": "S"]
        if let symbol = aromaticSymbols[ch], let element = VivoElement.from(symbol: symbol) {
            return (.init(element: element, isotope: nil, charge: 0, aromatic: true), 1)
        }
        if index + 1 < chars.count {
            let pair = String([ch, chars[index + 1]])
            if ["Cl", "Br"].contains(pair), let element = VivoElement.from(symbol: pair) {
                return (.init(element: element, isotope: nil, charge: 0, aromatic: false), 2)
            }
        }
        let symbol = String(ch)
        guard ["B", "C", "N", "O", "P", "S", "F", "I"].contains(symbol),
              let element = VivoElement.from(symbol: symbol) else {
            throw VivoArtifactValidationError.incompatible("unsupported or malformed unbracketed SMILES atom near '\(symbol)'")
        }
        return (.init(element: element, isotope: nil, charge: 0, aromatic: false), 1)
    }

    private static func parseBracket(_ content: String) throws -> ParsedAtom {
        let chars = Array(content)
        guard !chars.isEmpty else { throw VivoArtifactValidationError.invalid("empty bracket atom in SMILES") }
        var i = 0
        var isotopeDigits = ""
        while i < chars.count, chars[i].isNumber { isotopeDigits.append(chars[i]); i += 1 }
        let isotope = isotopeDigits.isEmpty ? nil : UInt16(isotopeDigits)
        guard i < chars.count else { throw VivoArtifactValidationError.invalid("bracket atom contains isotope but no element") }

        var aromatic = false
        var symbol: String
        if ["b", "c", "n", "o", "p", "s"].contains(chars[i]) {
            aromatic = true; symbol = String(chars[i]).uppercased(); i += 1
        } else {
            symbol = String(chars[i]); i += 1
            if i < chars.count, chars[i].isLowercase {
                let candidate = symbol + String(chars[i])
                if VivoElement.from(symbol: candidate) != nil { symbol = candidate; i += 1 }
            }
        }
        guard let element = VivoElement.from(symbol: symbol) else {
            throw VivoArtifactValidationError.invalid("unknown bracket atom element \(symbol)")
        }
        var charge: Int16 = 0
        while i < chars.count {
            let ch = chars[i]
            if ch == "@" {
                throw VivoArtifactValidationError.incompatible("SMILES tetrahedral stereochemistry requires the forthcoming stereo IR")
            }
            if ch == "H" {
                throw VivoArtifactValidationError.incompatible("SMILES bracket hydrogen counts are chemically significant but not yet represented; provide explicit hydrogen atoms or a coordinate format")
            }
            if ch == "+" || ch == "-" {
                let sign: Int16 = ch == "+" ? 1 : -1
                i += 1
                var digits = ""
                while i < chars.count, chars[i].isNumber { digits.append(chars[i]); i += 1 }
                if digits.isEmpty {
                    var repeated: Int16 = 1
                    while i < chars.count, chars[i] == ch { repeated += 1; i += 1 }
                    charge += sign * repeated
                } else if let magnitude = Int16(digits) { charge += sign * magnitude }
                else { throw VivoArtifactValidationError.invalid("SMILES bracket charge is out of range") }
                continue
            }
            if ch == ":" {
                throw VivoArtifactValidationError.incompatible("SMILES atom-map labels are not yet represented in VivoMolecularStructure")
            }
            throw VivoArtifactValidationError.incompatible("unsupported bracket-atom annotation '\(ch)' in SMILES")
        }
        return .init(element: element, isotope: isotope, charge: charge, aromatic: aromatic)
    }
}
