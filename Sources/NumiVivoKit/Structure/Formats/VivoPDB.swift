import Foundation

public enum VivoPDB {
    public struct ReadOptions: Sendable {
        public var identifier: String
        public var keepAlternateLocations: Bool
        public init(identifier: String = "pdb", keepAlternateLocations: Bool = true) {
            self.identifier = identifier
            self.keepAlternateLocations = keepAlternateLocations
        }
    }

    private struct AtomRecord: Equatable {
        var serial: Int32?
        var name: String
        var altLoc: String?
        var residueName: String
        var chain: String
        var residueSequence: Int32?
        var insertion: String?
        var positionNM: VivoVector3D
        var occupancy: Double?
        var bFactor: Double?
        var element: VivoElement
        var formalCharge: Int16
        var hetero: Bool
    }

    public static func read(_ text: String, options: ReadOptions = .init()) throws -> VivoMolecularStructure {
        let lines = text.split(whereSeparator: \.isNewline).map(String.init)
        var models: [[AtomRecord]] = [[]]
        var inExplicitModel = false
        var periodicCell: VivoPeriodicCell?
        var conectPairs = Set<UInt64>()

        for line in lines {
            let record = line.vivoField(0..<6).trimmingCharacters(in: .whitespaces)
            switch record {
            case "MODEL":
                if !models[0].isEmpty || inExplicitModel { models.append([]) }
                inExplicitModel = true
            case "ENDMDL":
                inExplicitModel = false
            case "ATOM", "HETATM":
                let atom = try parseAtom(line, hetero: record == "HETATM")
                if !options.keepAlternateLocations, let alt = atom.altLoc, alt != "A" { continue }
                if models.isEmpty { models = [[]] }
                models[models.count - 1].append(atom)
            case "CRYST1":
                periodicCell = try parseCell(line)
            case "CONECT":
                let serials = parseConect(line)
                if let first = serials.first {
                    for other in serials.dropFirst() where other != first {
                        let low = UInt32(bitPattern: min(first, other)), high = UInt32(bitPattern: max(first, other))
                        conectPairs.insert(UInt64(low) << 32 | UInt64(high))
                    }
                }
            default:
                continue
            }
        }
        models.removeAll(where: { $0.isEmpty })
        guard let reference = models.first, !reference.isEmpty else {
            throw VivoArtifactValidationError.invalid("PDB contains no ATOM/HETATM records")
        }
        for (modelIndex, model) in models.enumerated() {
            guard model.count == reference.count else {
                throw VivoArtifactValidationError.invalid("PDB model \(modelIndex + 1) atom count differs from model 1")
            }
            for index in reference.indices {
                guard sameIdentity(reference[index], model[index]) else {
                    throw VivoArtifactValidationError.invalid("PDB model \(modelIndex + 1) atom ordering/identity differs at atom \(index)")
                }
            }
        }

        let builder = VivoStructureAssemblyBuilder()
        var serialToIndex: [Int32: UInt32] = [:]
        for atom in reference {
            let index = try builder.appendAtom(name: atom.name, element: atom.element, positionNM: atom.positionNM,
                                               residueName: atom.residueName, chainIdentifier: atom.chain,
                                               residueSequence: atom.residueSequence, insertionCode: atom.insertion,
                                               sourceSerial: atom.serial, alternateLocation: atom.altLoc,
                                               occupancy: atom.occupancy, bFactor: atom.bFactor,
                                               formalCharge: atom.formalCharge, isHetero: atom.hetero)
            if let serial = atom.serial {
                guard serialToIndex.updateValue(index, forKey: serial) == nil else {
                    throw VivoArtifactValidationError.invalid("PDB model 1 contains duplicate atom serial \(serial)")
                }
            }
        }
        var bonds: [VivoMolecularBond] = []
        bonds.reserveCapacity(conectPairs.count)
        for encoded in conectPairs.sorted() {
            let lowSerial = Int32(bitPattern: UInt32(encoded >> 32))
            let highSerial = Int32(bitPattern: UInt32(encoded & 0xffff_ffff))
            guard let a = serialToIndex[lowSerial], let b = serialToIndex[highSerial] else { continue }
            bonds.append(.init(atomA: a, atomB: b, order: .unknown))
        }
        let conformers = models.enumerated().map { offset, model in
            VivoMolecularConformer(identifier: "model-\(offset + 1)", positionsNM: model.map(\.positionNM))
        }
        return try builder.makeStructure(identifier: options.identifier, bonds: bonds,
                                         conformers: conformers, periodicCell: periodicCell,
                                         metadata: ["sourceFormat": "PDB"])
    }

    public static func write(_ structure: VivoMolecularStructure, conformerIndex: Int = 0) throws -> String {
        _ = try VivoStructureValidator.validate(structure)
        guard structure.conformers.indices.contains(conformerIndex) else {
            throw VivoArtifactValidationError.unresolved("PDB writer conformer index is absent")
        }
        var output = ""
        if let cell = structure.periodicCell {
            let values = try cellParameters(cell)
            output += String(format: "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1\n",
                             values.aA, values.bA, values.cA, values.alpha, values.beta, values.gamma)
        }
        let coordinates = structure.conformers[conformerIndex].positionsNM
        for atom in structure.atoms {
            let residue = atom.residueIndex.map { structure.residues[Int($0)] }
            let chain = residue?.chainIndex.map { structure.chains[Int($0)].identifier } ?? ""
            let p = coordinates[Int(atom.index)] * 10.0
            let serial = atom.sourceSerial ?? Int32(atom.index + 1)
            let record = atom.isHetero ? "HETATM" : "ATOM  "
            let name = paddedAtomName(atom.name, element: atom.element.symbol)
            let alt = atom.alternateLocation?.prefix(1) ?? " "
            let resName = String((residue?.name ?? "UNK").prefix(3)).padding(toLength: 3, withPad: " ", startingAt: 0)
            let chainID = String(chain.prefix(1))
            let seq = residue?.sequenceNumber ?? 0
            let ins = residue?.insertionCode?.prefix(1) ?? " "
            let occ = atom.occupancy ?? 1.0
            let b = atom.bFactor ?? 0.0
            let element = atom.element.symbol.padding(toLength: 2, withPad: " ", startingAt: 0)
            let charge = pdbCharge(atom.formalCharge)
            output += String(format: "%@%5d %@%@%@ %@%4d%@   %8.3f%8.3f%8.3f%6.2f%6.2f          %2@%2@\n",
                             record, serial, name, String(alt), resName, chainID, seq, String(ins),
                             p.x, p.y, p.z, occ, b, element, charge)
        }
        if !structure.bonds.isEmpty {
            let serials = structure.atoms.map { $0.sourceSerial ?? Int32($0.index + 1) }
            var neighbors: [UInt32: [UInt32]] = [:]
            for bond in structure.bonds {
                neighbors[bond.atomA, default: []].append(bond.atomB)
                neighbors[bond.atomB, default: []].append(bond.atomA)
            }
            for atom in neighbors.keys.sorted() {
                let list = neighbors[atom, default: []].sorted()
                var cursor = 0
                while cursor < list.count {
                    let chunk = list[cursor..<min(cursor + 4, list.count)]
                    output += String(format: "CONECT%5d", serials[Int(atom)])
                    for n in chunk { output += String(format: "%5d", serials[Int(n)]) }
                    output += "\n"
                    cursor += chunk.count
                }
            }
        }
        output += "END\n"
        return output
    }

    private static func parseAtom(_ line: String, hetero: Bool) throws -> AtomRecord {
        let rawName = line.vivoField(12..<16)
        let name = rawName.trimmingCharacters(in: .whitespaces)
        guard !name.isEmpty else { throw VivoArtifactValidationError.invalid("PDB atom has empty name") }
        let elementField = line.vivoField(76..<78).trimmingCharacters(in: .whitespaces)
        guard let element = VivoElement.from(symbol: elementField.isEmpty ? inferredElement(rawName) : elementField) else {
            throw VivoArtifactValidationError.invalid("PDB atom \(name) has unknown element")
        }
        guard let x = Double(line.vivoField(30..<38).trimmingCharacters(in: .whitespaces)),
              let y = Double(line.vivoField(38..<46).trimmingCharacters(in: .whitespaces)),
              let z = Double(line.vivoField(46..<54).trimmingCharacters(in: .whitespaces)) else {
            throw VivoArtifactValidationError.invalid("PDB atom \(name) has invalid coordinates")
        }
        return AtomRecord(serial: Int32(line.vivoField(6..<11).trimmingCharacters(in: .whitespaces)),
                          name: name,
                          altLoc: normalized(line.vivoField(16..<17)),
                          residueName: line.vivoField(17..<20).trimmingCharacters(in: .whitespaces),
                          chain: line.vivoField(21..<22).trimmingCharacters(in: .whitespaces),
                          residueSequence: Int32(line.vivoField(22..<26).trimmingCharacters(in: .whitespaces)),
                          insertion: normalized(line.vivoField(26..<27)),
                          positionNM: .init(x * 0.1, y * 0.1, z * 0.1),
                          occupancy: Double(line.vivoField(54..<60).trimmingCharacters(in: .whitespaces)),
                          bFactor: Double(line.vivoField(60..<66).trimmingCharacters(in: .whitespaces)),
                          element: element,
                          formalCharge: parseCharge(line.vivoField(78..<80)), hetero: hetero)
    }

    private static func parseConect(_ line: String) -> [Int32] {
        var result: [Int32] = []
        for start in stride(from: 6, to: min(line.utf8.count, 31), by: 5) {
            if let value = Int32(line.vivoField(start..<min(start + 5, line.utf8.count)).trimmingCharacters(in: .whitespaces)) {
                result.append(value)
            }
        }
        return result
    }

    private static func parseCell(_ line: String) throws -> VivoPeriodicCell {
        guard let a = Double(line.vivoField(6..<15).trimmingCharacters(in: .whitespaces)),
              let b = Double(line.vivoField(15..<24).trimmingCharacters(in: .whitespaces)),
              let c = Double(line.vivoField(24..<33).trimmingCharacters(in: .whitespaces)),
              let alpha = Double(line.vivoField(33..<40).trimmingCharacters(in: .whitespaces)),
              let beta = Double(line.vivoField(40..<47).trimmingCharacters(in: .whitespaces)),
              let gamma = Double(line.vivoField(47..<54).trimmingCharacters(in: .whitespaces)) else {
            throw VivoArtifactValidationError.invalid("PDB CRYST1 record is invalid")
        }
        return try cellFromParameters(aA: a, bA: b, cA: c, alpha: alpha, beta: beta, gamma: gamma)
    }

    static func cellFromParameters(aA: Double, bA: Double, cA: Double,
                                   alpha: Double, beta: Double, gamma: Double) throws -> VivoPeriodicCell {
        let d2r = Double.pi / 180
        let ca = cos(alpha * d2r), cb = cos(beta * d2r), cg = cos(gamma * d2r), sg = sin(gamma * d2r)
        guard aA > 0, bA > 0, cA > 0, abs(sg) > 1e-12 else {
            throw VivoArtifactValidationError.invalid("unit cell lengths/angles are degenerate")
        }
        let av = VivoVector3D(aA * 0.1, 0, 0)
        let bv = VivoVector3D(bA * cg * 0.1, bA * sg * 0.1, 0)
        let cx = cA * cb
        let cy = cA * (ca - cb * cg) / sg
        let cz2 = cA * cA - cx * cx - cy * cy
        guard cz2 > -1e-8 else { throw VivoArtifactValidationError.invalid("unit cell angles are inconsistent") }
        let cv = VivoVector3D(cx * 0.1, cy * 0.1, sqrt(max(0, cz2)) * 0.1)
        let cell = VivoPeriodicCell(a: av, b: bv, c: cv)
        guard cell.isValid else { throw VivoArtifactValidationError.invalid("unit cell is invalid") }
        return cell
    }

    private static func cellParameters(_ cell: VivoPeriodicCell) throws -> (aA: Double, bA: Double, cA: Double, alpha: Double, beta: Double, gamma: Double) {
        guard cell.isValid else { throw VivoArtifactValidationError.invalid("cannot write invalid unit cell") }
        let an = cell.a.norm, bn = cell.b.norm, cn = cell.c.norm
        func angle(_ u: VivoVector3D, _ v: VivoVector3D, _ un: Double, _ vn: Double) -> Double {
            acos(max(-1, min(1, u.dot(v) / (un * vn)))) * 180 / Double.pi
        }
        return (an * 10, bn * 10, cn * 10,
                angle(cell.b, cell.c, bn, cn), angle(cell.a, cell.c, an, cn), angle(cell.a, cell.b, an, bn))
    }

    private static func sameIdentity(_ lhs: AtomRecord, _ rhs: AtomRecord) -> Bool {
        lhs.serial == rhs.serial && lhs.name == rhs.name && lhs.altLoc == rhs.altLoc &&
        lhs.residueName == rhs.residueName && lhs.chain == rhs.chain &&
        lhs.residueSequence == rhs.residueSequence && lhs.insertion == rhs.insertion &&
        lhs.element == rhs.element && lhs.hetero == rhs.hetero
    }

    private static func inferredElement(_ rawName: String) -> String {
        let trimmed = rawName.trimmingCharacters(in: .whitespacesAndNewlines)
        guard !trimmed.isEmpty else { return "" }
        if rawName.first == " " { return String(trimmed.prefix(1)) }
        let letters = trimmed.drop(while: { $0.isNumber })
        if letters.count >= 2 {
            let pair = String(letters.prefix(2))
            if VivoElement.from(symbol: pair) != nil { return pair }
        }
        return String(letters.prefix(1))
    }

    private static func normalized(_ value: String) -> String? {
        let trimmed = value.trimmingCharacters(in: .whitespacesAndNewlines)
        return trimmed.isEmpty ? nil : trimmed
    }

    private static func parseCharge(_ raw: String) -> Int16 {
        let text = raw.trimmingCharacters(in: .whitespaces)
        guard text.count == 2, let magnitude = Int16(String(text.prefix(1))) else { return 0 }
        return text.last == "-" ? -magnitude : magnitude
    }

    private static func pdbCharge(_ charge: Int16) -> String {
        guard charge != 0, abs(charge) <= 9 else { return "  " }
        return "\(abs(charge))\(charge < 0 ? "-" : "+")"
    }

    private static func paddedAtomName(_ name: String, element: String) -> String {
        let clipped = String(name.prefix(4))
        if element.count == 1 && clipped.count < 4 { return " " + clipped.padding(toLength: 3, withPad: " ", startingAt: 0) }
        return clipped.padding(toLength: 4, withPad: " ", startingAt: 0)
    }
}
