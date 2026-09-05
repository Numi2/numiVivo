import Foundation

public enum VivoMOL2 {
    public static func read(_ text: String, identifier: String = "mol2") throws -> VivoMolecularStructure {
        enum Section { case none, molecule, atom, bond, substructure }
        var section: Section = .none
        var title = identifier
        var atomLines: [String] = []
        var bondLines: [String] = []
        var substructureNames: [Int: String] = [:]
        let lines = text.split(separator: "\n", omittingEmptySubsequences: false).map(String.init)
        var moleculeTitlePending = false
        for raw in lines {
            let line = raw.trimmingCharacters(in: .whitespacesAndNewlines)
            if line.hasPrefix("@<TRIPOS>") {
                switch line.uppercased() {
                case "@<TRIPOS>MOLECULE": section = .molecule; moleculeTitlePending = true
                case "@<TRIPOS>ATOM": section = .atom
                case "@<TRIPOS>BOND": section = .bond
                case "@<TRIPOS>SUBSTRUCTURE": section = .substructure
                default: section = .none
                }
                continue
            }
            guard !line.isEmpty && !line.hasPrefix("#") else { continue }
            if section == .molecule, moleculeTitlePending {
                title = line; moleculeTitlePending = false; continue
            }
            if section == .atom { atomLines.append(line) }
            if section == .bond { bondLines.append(line) }
            if section == .substructure {
                let parts = line.split(whereSeparator: \.isWhitespace)
                if parts.count >= 2, let id = Int(parts[0]) { substructureNames[id] = String(parts[1]) }
            }
        }
        guard !atomLines.isEmpty else { throw VivoArtifactValidationError.invalid("MOL2 contains no ATOM section") }

        struct ParsedAtom {
            var id: Int
            var name: String
            var element: VivoElement
            var position: VivoVector3D
            var substructureID: Int?
            var substructureName: String?
        }
        var parsedAtoms: [ParsedAtom] = []
        parsedAtoms.reserveCapacity(atomLines.count)
        var idToDense: [Int: UInt32] = [:]
        for (dense, line) in atomLines.enumerated() {
            let parts = line.split(whereSeparator: \.isWhitespace).map(String.init)
            guard parts.count >= 6, let id = Int(parts[0]),
                  let x = Double(parts[2]), let y = Double(parts[3]), let z = Double(parts[4]) else {
                throw VivoArtifactValidationError.invalid("MOL2 atom row \(dense + 1) is malformed")
            }
            let elementToken = parts[5].split(separator: ".").first.map(String.init) ?? parts[5]
            guard let element = VivoElement.from(symbol: elementToken) ?? inferredElement(parts[1]) else {
                throw VivoArtifactValidationError.invalid("MOL2 atom \(parts[1]) has unknown element/type \(parts[5])")
            }
            guard idToDense.updateValue(UInt32(dense), forKey: id) == nil else {
                throw VivoArtifactValidationError.invalid("MOL2 contains duplicate atom id \(id)")
            }
            let subID = parts.count > 6 ? Int(parts[6]) : nil
            let subName = parts.count > 7 ? parts[7] : subID.flatMap { substructureNames[$0] }
            parsedAtoms.append(.init(id: id, name: parts[1], element: element,
                                     position: .init(x * 0.1, y * 0.1, z * 0.1),
                                     substructureID: subID, substructureName: subName))
        }

        let builder = VivoStructureAssemblyBuilder()
        for atom in parsedAtoms {
            try builder.appendAtom(name: atom.name, element: atom.element, positionNM: atom.position,
                                   residueName: atom.substructureName,
                                   residueSequence: atom.substructureID.flatMap(Int32.init),
                                   sourceSerial: Int32(exactly: atom.id), isHetero: true)
        }
        var bonds: [VivoMolecularBond] = []
        bonds.reserveCapacity(bondLines.count)
        for (row, line) in bondLines.enumerated() {
            let parts = line.split(whereSeparator: \.isWhitespace).map(String.init)
            guard parts.count >= 4, let from = Int(parts[1]), let to = Int(parts[2]),
                  let a = idToDense[from], let b = idToDense[to] else {
                throw VivoArtifactValidationError.invalid("MOL2 bond row \(row + 1) is malformed")
            }
            let order: VivoBondOrder = switch parts[3].lowercased() {
            case "1", "am": .single
            case "2": .double
            case "3": .triple
            case "ar": .aromatic
            default: .unknown
            }
            bonds.append(.init(atomA: a, atomB: b, order: order))
        }
        return try builder.makeStructure(identifier: title, bonds: bonds,
                                         metadata: ["sourceFormat": "MOL2"])
    }

    public static func write(_ structure: VivoMolecularStructure, conformerIndex: Int = 0) throws -> String {
        _ = try VivoStructureValidator.validate(structure)
        guard structure.conformers.indices.contains(conformerIndex) else {
            throw VivoArtifactValidationError.unresolved("MOL2 writer conformer index is absent")
        }
        var output = "@<TRIPOS>MOLECULE\n\(structure.identifier)\n"
        output += "\(structure.atoms.count) \(structure.bonds.count) \(structure.residues.count) 0 0\nSMALL\nNO_CHARGES\n\n"
        output += "@<TRIPOS>ATOM\n"
        let positions = structure.conformers[conformerIndex].positionsNM
        for atom in structure.atoms {
            let p = positions[Int(atom.index)] * 10
            let residue = atom.residueIndex.map { structure.residues[Int($0)] }
            let sid = residue.map { Int($0.index) + 1 } ?? 1
            let sname = residue?.name ?? "MOL"
            output += String(format: "%7d %-8@ %10.5f %10.5f %10.5f %-6@ %4d %-8@ 0.0000\n",
                             atom.index + 1, atom.name, p.x, p.y, p.z, atom.element.symbol, sid, sname)
        }
        output += "@<TRIPOS>BOND\n"
        for (index, bond) in structure.bonds.enumerated() {
            let type: String = switch bond.order {
            case .single: "1"; case .double: "2"; case .triple: "3"; case .aromatic: "ar"; default: "un"
            }
            output += "\(index + 1) \(bond.atomA + 1) \(bond.atomB + 1) \(type)\n"
        }
        if !structure.residues.isEmpty {
            output += "@<TRIPOS>SUBSTRUCTURE\n"
            for residue in structure.residues {
                let root = residue.atomIndices.first.map { $0 + 1 } ?? 1
                output += "\(residue.index + 1) \(residue.name) \(root) RESIDUE\n"
            }
        }
        return output
    }

    private static func inferredElement(_ atomName: String) -> VivoElement? {
        let letters = atomName.drop(while: { $0.isNumber })
        if letters.count >= 2, let pair = VivoElement.from(symbol: String(letters.prefix(2))) { return pair }
        return VivoElement.from(symbol: String(letters.prefix(1)))
    }
}
