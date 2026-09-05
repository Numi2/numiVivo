import Foundation

public enum VivoMMCIF {
    private struct Loop {
        var headers: [String]
        var rows: [[String]]
    }
    private struct AtomRow: Equatable {
        var model: Int
        var serial: Int32?
        var name: String
        var alt: String?
        var residue: String
        var chain: String
        var sequence: Int32?
        var insertion: String?
        var element: VivoElement
        var charge: Int16
        var hetero: Bool
        var occupancy: Double?
        var bFactor: Double?
        var position: VivoVector3D
    }

    public static func read(_ text: String, identifier: String = "mmcif") throws -> VivoMolecularStructure {
        let tokens = try tokenize(text)
        let parsed = try parse(tokens)
        guard let atomLoop = parsed.loops.first(where: { $0.headers.contains(where: { $0.hasPrefix("_atom_site.") }) }) else {
            throw VivoArtifactValidationError.invalid("mmCIF contains no _atom_site loop")
        }
        let atoms = try atomRows(atomLoop)
        guard !atoms.isEmpty else { throw VivoArtifactValidationError.invalid("mmCIF _atom_site loop is empty") }
        let grouped = Dictionary(grouping: atoms, by: \.model)
        let modelIDs = grouped.keys.sorted()
        guard let reference = grouped[modelIDs[0]] else { throw VivoArtifactValidationError.invalid("mmCIF has no model") }
        for model in modelIDs.dropFirst() {
            guard let candidate = grouped[model], candidate.count == reference.count else {
                throw VivoArtifactValidationError.invalid("mmCIF model \(model) atom count differs from first model")
            }
            for index in reference.indices where !sameIdentity(reference[index], candidate[index]) {
                throw VivoArtifactValidationError.invalid("mmCIF model \(model) atom identity/order differs at atom \(index)")
            }
        }

        let builder = VivoStructureAssemblyBuilder()
        for atom in reference {
            try builder.appendAtom(name: atom.name, element: atom.element, positionNM: atom.position,
                                   residueName: atom.residue, chainIdentifier: atom.chain,
                                   residueSequence: atom.sequence, insertionCode: atom.insertion,
                                   sourceSerial: atom.serial, alternateLocation: atom.alt,
                                   occupancy: atom.occupancy, bFactor: atom.bFactor,
                                   formalCharge: atom.charge, isHetero: atom.hetero)
        }
        let conformers = try modelIDs.map { model -> VivoMolecularConformer in
            guard let rows = grouped[model] else { throw VivoArtifactValidationError.invalid("mmCIF model disappeared") }
            return .init(identifier: "model-\(model)", positionsNM: rows.map(\.position))
        }
        let cell = try periodicCell(parsed.scalars)
        return try builder.makeStructure(identifier: identifier, conformers: conformers,
                                         periodicCell: cell, metadata: ["sourceFormat": "mmCIF"])
    }

    private static func atomRows(_ loop: Loop) throws -> [AtomRow] {
        let columns = Dictionary(uniqueKeysWithValues: loop.headers.enumerated().map { ($1, $0) })
        func column(_ names: [String]) -> Int? { names.compactMap { columns[$0] }.first }
        guard let x = column(["_atom_site.Cartn_x"]),
              let y = column(["_atom_site.Cartn_y"]),
              let z = column(["_atom_site.Cartn_z"]),
              let element = column(["_atom_site.type_symbol"]),
              let name = column(["_atom_site.auth_atom_id", "_atom_site.label_atom_id"]),
              let residue = column(["_atom_site.auth_comp_id", "_atom_site.label_comp_id"]),
              let chain = column(["_atom_site.auth_asym_id", "_atom_site.label_asym_id"]) else {
            throw VivoArtifactValidationError.invalid("mmCIF _atom_site lacks required coordinate/identity columns")
        }
        let serial = column(["_atom_site.id"])
        let sequence = column(["_atom_site.auth_seq_id", "_atom_site.label_seq_id"])
        let insertion = column(["_atom_site.pdbx_PDB_ins_code"])
        let alternate = column(["_atom_site.label_alt_id"])
        let occupancy = column(["_atom_site.occupancy"])
        let bFactor = column(["_atom_site.B_iso_or_equiv"])
        let charge = column(["_atom_site.pdbx_formal_charge"])
        let group = column(["_atom_site.group_PDB"])
        let model = column(["_atom_site.pdbx_PDB_model_num"])

        return try loop.rows.map { row in
            guard let xx = Double(row[x]), let yy = Double(row[y]), let zz = Double(row[z]),
                  let resolvedElement = VivoElement.from(symbol: row[element]) else {
                throw VivoArtifactValidationError.invalid("mmCIF _atom_site row has invalid coordinate or element")
            }
            return AtomRow(model: model.flatMap { Int(normalized(row[$0]) ?? "") } ?? 1,
                           serial: serial.flatMap { Int32(normalized(row[$0]) ?? "") },
                           name: normalized(row[name]) ?? row[name],
                           alt: alternate.flatMap { normalized(row[$0]) },
                           residue: normalized(row[residue]) ?? "UNK",
                           chain: normalized(row[chain]) ?? "",
                           sequence: sequence.flatMap { Int32(normalized(row[$0]) ?? "") },
                           insertion: insertion.flatMap { normalized(row[$0]) },
                           element: resolvedElement,
                           charge: charge.flatMap { Int16(normalized(row[$0]) ?? "") } ?? 0,
                           hetero: group.map { row[$0].uppercased() == "HETATM" } ?? false,
                           occupancy: occupancy.flatMap { Double(normalized(row[$0]) ?? "") },
                           bFactor: bFactor.flatMap { Double(normalized(row[$0]) ?? "") },
                           position: .init(xx * 0.1, yy * 0.1, zz * 0.1))
        }
    }

    private static func periodicCell(_ scalars: [String: String]) throws -> VivoPeriodicCell? {
        guard let aa = scalars["_cell.length_a"].flatMap(Double.init),
              let bb = scalars["_cell.length_b"].flatMap(Double.init),
              let cc = scalars["_cell.length_c"].flatMap(Double.init) else { return nil }
        let alpha = scalars["_cell.angle_alpha"].flatMap(Double.init) ?? 90
        let beta = scalars["_cell.angle_beta"].flatMap(Double.init) ?? 90
        let gamma = scalars["_cell.angle_gamma"].flatMap(Double.init) ?? 90
        return try VivoPDB.cellFromParameters(aA: aa, bA: bb, cA: cc, alpha: alpha, beta: beta, gamma: gamma)
    }

    private static func sameIdentity(_ a: AtomRow, _ b: AtomRow) -> Bool {
        a.serial == b.serial && a.name == b.name && a.alt == b.alt && a.residue == b.residue &&
        a.chain == b.chain && a.sequence == b.sequence && a.insertion == b.insertion &&
        a.element == b.element && a.charge == b.charge && a.hetero == b.hetero
    }

    private static func normalized(_ value: String) -> String? {
        value == "." || value == "?" || value.isEmpty ? nil : value
    }

    private static func parse(_ tokens: [String]) throws -> (scalars: [String: String], loops: [Loop]) {
        var scalars: [String: String] = [:]
        var loops: [Loop] = []
        var i = 0
        while i < tokens.count {
            let token = tokens[i]
            if token == "loop_" {
                i += 1
                var headers: [String] = []
                while i < tokens.count, tokens[i].hasPrefix("_") { headers.append(tokens[i]); i += 1 }
                guard !headers.isEmpty else { throw VivoArtifactValidationError.invalid("mmCIF loop has no headers") }
                var flat: [String] = []
                while i < tokens.count {
                    let next = tokens[i]
                    if flat.count % headers.count == 0,
                       next == "loop_" || next == "stop_" || next.hasPrefix("data_") || next.hasPrefix("save_") || next.hasPrefix("_") {
                        break
                    }
                    flat.append(next); i += 1
                }
                guard flat.count % headers.count == 0 else {
                    throw VivoArtifactValidationError.invalid("mmCIF loop value count is not divisible by header count")
                }
                var rows: [[String]] = []
                for start in stride(from: 0, to: flat.count, by: headers.count) {
                    rows.append(Array(flat[start..<start + headers.count]))
                }
                loops.append(.init(headers: headers, rows: rows))
                if i < tokens.count, tokens[i] == "stop_" { i += 1 }
            } else if token.hasPrefix("_") {
                guard i + 1 < tokens.count else { throw VivoArtifactValidationError.invalid("mmCIF scalar \(token) has no value") }
                scalars[token] = tokens[i + 1]
                i += 2
            } else {
                i += 1
            }
        }
        return (scalars, loops)
    }

    /// CIF tokenizer supporting comments, single/double quoted values and
    /// semicolon-delimited multiline text. It intentionally does not interpret
    /// chemistry while tokenizing.
    private static func tokenize(_ text: String) throws -> [String] {
        let lines = text.split(separator: "\n", omittingEmptySubsequences: false).map(String.init)
        var tokens: [String] = []
        var lineIndex = 0
        while lineIndex < lines.count {
            var line = lines[lineIndex]
            if line.hasPrefix(";") {
                var block = String(line.dropFirst())
                lineIndex += 1
                var closed = false
                while lineIndex < lines.count {
                    if lines[lineIndex].hasPrefix(";") { closed = true; break }
                    if !block.isEmpty { block += "\n" }
                    block += lines[lineIndex]
                    lineIndex += 1
                }
                guard closed else { throw VivoArtifactValidationError.invalid("unterminated mmCIF multiline value") }
                tokens.append(block)
                lineIndex += 1
                continue
            }
            var cursor = line.startIndex
            while cursor < line.endIndex {
                while cursor < line.endIndex, line[cursor].isWhitespace { cursor = line.index(after: cursor) }
                if cursor >= line.endIndex || line[cursor] == "#" { break }
                let first = line[cursor]
                if first == "'" || first == "\"" {
                    let quote = first
                    cursor = line.index(after: cursor)
                    var value = ""
                    var closed = false
                    while cursor < line.endIndex {
                        if line[cursor] == quote { closed = true; cursor = line.index(after: cursor); break }
                        value.append(line[cursor]); cursor = line.index(after: cursor)
                    }
                    guard closed else { throw VivoArtifactValidationError.invalid("unterminated quoted mmCIF value") }
                    tokens.append(value)
                } else {
                    let start = cursor
                    while cursor < line.endIndex, !line[cursor].isWhitespace, line[cursor] != "#" {
                        cursor = line.index(after: cursor)
                    }
                    tokens.append(String(line[start..<cursor]))
                    if cursor < line.endIndex, line[cursor] == "#" { break }
                }
            }
            lineIndex += 1
        }
        return tokens
    }
}
