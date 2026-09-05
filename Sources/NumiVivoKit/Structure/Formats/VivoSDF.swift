import Foundation

public enum VivoSDF {
    public static func readAll(_ text: String, identifierPrefix: String = "sdf") throws -> [VivoMolecularStructure] {
        let records = text.components(separatedBy: "$$$$")
            .map { $0.trimmingCharacters(in: .whitespacesAndNewlines) }
            .filter { !$0.isEmpty }
        guard !records.isEmpty else { throw VivoArtifactValidationError.invalid("SDF contains no molecule records") }
        return try records.enumerated().map { index, record in
            try readRecord(record, fallbackIdentifier: "\(identifierPrefix)-\(index + 1)")
        }
    }

    public static func readFirst(_ text: String, identifier: String = "sdf") throws -> VivoMolecularStructure {
        guard let first = text.components(separatedBy: "$$$$").first else {
            throw VivoArtifactValidationError.invalid("SDF is empty")
        }
        return try readRecord(first, fallbackIdentifier: identifier)
    }

    private static func readRecord(_ text: String, fallbackIdentifier: String) throws -> VivoMolecularStructure {
        let lines = text.split(separator: "\n", omittingEmptySubsequences: false).map(String.init)
        guard lines.count >= 4 else { throw VivoArtifactValidationError.invalid("SDF record is shorter than V2000 header") }
        let title = lines[0].trimmingCharacters(in: .whitespacesAndNewlines)
        let counts = lines[3]
        if counts.contains("V3000") {
            throw VivoArtifactValidationError.incompatible("SDF V3000 is not yet supported; no lossy fallback is applied")
        }
        guard counts.contains("V2000"),
              let atomCount = Int(counts.vivoField(0..<3).trimmingCharacters(in: .whitespaces)),
              let bondCount = Int(counts.vivoField(3..<6).trimmingCharacters(in: .whitespaces)),
              atomCount >= 0, bondCount >= 0,
              lines.count >= 4 + atomCount + bondCount else {
            throw VivoArtifactValidationError.invalid("SDF V2000 counts line or record length is invalid")
        }
        let chiralFlag = Int(counts.vivoField(12..<15).trimmingCharacters(in: .whitespaces)) ?? 0
        guard chiralFlag == 0 else {
            throw VivoArtifactValidationError.incompatible("SDF V2000 chiral flag is set; authoritative tetrahedral stereo IR is not implemented yet")
        }

        let builder = VivoStructureAssemblyBuilder()
        var charges = [Int16](repeating: 0, count: atomCount)
        var isotopes = [UInt16?](repeating: nil, count: atomCount)
        for index in 0..<atomCount {
            let line = lines[4 + index]
            guard let x = Double(line.vivoField(0..<10).trimmingCharacters(in: .whitespaces)),
                  let y = Double(line.vivoField(10..<20).trimmingCharacters(in: .whitespaces)),
                  let z = Double(line.vivoField(20..<30).trimmingCharacters(in: .whitespaces)),
                  let element = VivoElement.from(symbol: line.vivoField(31..<34).trimmingCharacters(in: .whitespaces)) else {
                throw VivoArtifactValidationError.invalid("SDF atom \(index + 1) is malformed or uses a query/pseudo atom")
            }
            let massDifference = Int(line.vivoField(34..<36).trimmingCharacters(in: .whitespaces)) ?? 0
            guard massDifference == 0 else {
                throw VivoArtifactValidationError.incompatible("legacy SDF mass-difference isotope encoding is not imported; use an M  ISO record")
            }
            let chargeCode = Int(line.vivoField(36..<39).trimmingCharacters(in: .whitespaces)) ?? 0
            guard chargeCode != 4 else {
                throw VivoArtifactValidationError.incompatible("SDF radical atom encoding is not represented in VivoMolecularStructure")
            }
            let atomStereo = Int(line.vivoField(39..<42).trimmingCharacters(in: .whitespaces)) ?? 0
            guard atomStereo == 0 else {
                throw VivoArtifactValidationError.incompatible("SDF atom stereo parity is not represented; import is rejected rather than made achiral")
            }
            charges[index] = formalCharge(code: chargeCode)
            try builder.appendAtom(name: "\(element.symbol)\(index + 1)", element: element,
                                   positionNM: .init(x * 0.1, y * 0.1, z * 0.1),
                                   sourceSerial: Int32(index + 1), formalCharge: charges[index], isHetero: true)
        }
        var bonds: [VivoMolecularBond] = []
        for index in 0..<bondCount {
            let line = lines[4 + atomCount + index]
            guard let a1 = Int(line.vivoField(0..<3).trimmingCharacters(in: .whitespaces)),
                  let b1 = Int(line.vivoField(3..<6).trimmingCharacters(in: .whitespaces)),
                  a1 > 0, b1 > 0, a1 <= atomCount, b1 <= atomCount else {
                throw VivoArtifactValidationError.invalid("SDF bond \(index + 1) has invalid endpoints")
            }
            let type = Int(line.vivoField(6..<9).trimmingCharacters(in: .whitespaces)) ?? 0
            let order: VivoBondOrder
            switch type {
            case 1: order = .single
            case 2: order = .double
            case 3: order = .triple
            case 4: order = .aromatic
            default:
                throw VivoArtifactValidationError.incompatible("SDF query/unknown bond type \(type) is not a concrete molecular bond")
            }
            let stereo = Int(line.vivoField(9..<12).trimmingCharacters(in: .whitespaces)) ?? 0
            guard stereo == 0 else {
                throw VivoArtifactValidationError.incompatible("SDF bond stereo flag \(stereo) requires authoritative stereo semantics; import is rejected")
            }
            bonds.append(.init(atomA: UInt32(a1 - 1), atomB: UInt32(b1 - 1), order: order))
        }

        let propertyLines = lines.dropFirst(4 + atomCount + bondCount)
        for line in propertyLines {
            if line.hasPrefix("M  CHG") {
                let parts = line.split(whereSeparator: \.isWhitespace)
                guard parts.count >= 3, let pairCount = Int(parts[2]), parts.count >= 3 + pairCount * 2 else {
                    throw VivoArtifactValidationError.invalid("SDF M  CHG record is malformed")
                }
                for pair in 0..<pairCount {
                    guard let atom = Int(parts[3 + pair * 2]), let charge = Int16(parts[4 + pair * 2]),
                          atom > 0, atom <= atomCount else {
                        throw VivoArtifactValidationError.invalid("SDF M  CHG pair is invalid")
                    }
                    charges[atom - 1] = charge
                }
            } else if line.hasPrefix("M  ISO") {
                let parts = line.split(whereSeparator: \.isWhitespace)
                guard parts.count >= 3, let pairCount = Int(parts[2]), parts.count >= 3 + pairCount * 2 else {
                    throw VivoArtifactValidationError.invalid("SDF M  ISO record is malformed")
                }
                for pair in 0..<pairCount {
                    guard let atom = Int(parts[3 + pair * 2]), let isotope = UInt16(parts[4 + pair * 2]),
                          atom > 0, atom <= atomCount, isotope > 0 else {
                        throw VivoArtifactValidationError.invalid("SDF M  ISO pair is invalid")
                    }
                    isotopes[atom - 1] = isotope
                }
            } else if line.hasPrefix("M  RAD") || line.hasPrefix("M  ALS") {
                throw VivoArtifactValidationError.incompatible("SDF radical/query atom properties are not represented by the concrete structure IR")
            } else if line.hasPrefix("M  "), !line.hasPrefix("M  END") {
                // Unknown CTAB properties can affect chemistry. Data fields after
                // M END are application metadata and are intentionally separate.
                throw VivoArtifactValidationError.incompatible("unsupported SDF CTAB property '\(line.prefix(6))'")
            }
        }

        var structure = try builder.makeStructure(identifier: title.isEmpty ? fallbackIdentifier : title,
                                                  bonds: bonds, metadata: ["sourceFormat": "SDF-V2000"])
        for index in charges.indices {
            structure.atoms[index].formalCharge = charges[index]
            structure.atoms[index].isotopeMassNumber = isotopes[index]
        }
        _ = try VivoStructureValidator.validate(structure)
        return structure
    }

    public static func write(_ structure: VivoMolecularStructure, conformerIndex: Int = 0) throws -> String {
        _ = try VivoStructureValidator.validate(structure)
        guard structure.atoms.count <= 999, structure.bonds.count <= 999,
              structure.conformers.indices.contains(conformerIndex) else {
            throw VivoArtifactValidationError.incompatible("SDF V2000 writer requires <=999 atoms/bonds and an existing conformer")
        }
        guard structure.bonds.allSatisfy({ $0.stereo == .unspecified }) else {
            throw VivoArtifactValidationError.incompatible("SDF writer cannot yet emit authoritative stereochemistry")
        }
        var output = "\(structure.identifier)\nNumiVivo\n\n"
        output += String(format: "%3d%3d  0  0  0  0            999 V2000\n", structure.atoms.count, structure.bonds.count)
        let positions = structure.conformers[conformerIndex].positionsNM
        for atom in structure.atoms {
            let p = positions[Int(atom.index)] * 10
            let symbol = atom.element.symbol.padding(toLength: 3, withPad: " ", startingAt: 0)
            output += String(format: "%10.4f%10.4f%10.4f %@ 0  0  0  0  0  0  0  0  0  0  0  0\n",
                             p.x, p.y, p.z, symbol)
        }
        for bond in structure.bonds {
            let type: Int
            switch bond.order {
            case .single: type = 1
            case .double: type = 2
            case .triple: type = 3
            case .aromatic: type = 4
            default: throw VivoArtifactValidationError.incompatible("SDF writer cannot encode bond order \(bond.order.rawValue)")
            }
            output += String(format: "%3d%3d%3d  0  0  0  0\n", bond.atomA + 1, bond.atomB + 1, type)
        }
        let charged = structure.atoms.filter { $0.formalCharge != 0 }
        var cursor = 0
        while cursor < charged.count {
            let chunk = charged[cursor..<Swift.min(cursor + 8, charged.count)]
            output += String(format: "M  CHG%3d", chunk.count)
            for atom in chunk { output += String(format: "%4d%4d", atom.index + 1, atom.formalCharge) }
            output += "\n"
            cursor += chunk.count
        }
        let isotopic = structure.atoms.filter { $0.isotopeMassNumber != nil }
        cursor = 0
        while cursor < isotopic.count {
            let chunk = isotopic[cursor..<Swift.min(cursor + 8, isotopic.count)]
            output += String(format: "M  ISO%3d", chunk.count)
            for atom in chunk { output += String(format: "%4d%4d", atom.index + 1, atom.isotopeMassNumber!) }
            output += "\n"
            cursor += chunk.count
        }
        output += "M  END\n$$$$\n"
        return output
    }

    private static func formalCharge(code: Int) -> Int16 {
        switch code { case 1: 3; case 2: 2; case 3: 1; case 5: -1; case 6: -2; case 7: -3; default: 0 }
    }
}
