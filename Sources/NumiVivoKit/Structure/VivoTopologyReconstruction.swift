import Foundation

public struct VivoTopologyBuildReport: Codable, Sendable, Equatable {
    public var templateBondsAdded: UInt32
    public var peptideBondsAdded: UInt32
    public var disulfideBondsAdded: UInt32
    public var hydrogenBondsAdded: UInt32
    public var unresolvedResidues: [UInt32]
}

public enum VivoBiopolymerTopologyBuilder {
    /// Adds deterministic standard-amino-acid topology without changing atom
    /// identity or coordinates. Unknown residues are left untouched and reported.
    @discardableResult
    public static func apply(to structure: inout VivoMolecularStructure,
                             conformerIndex: Int = 0,
                             inferHydrogenAttachment: Bool = true,
                             inferDisulfides: Bool = true) throws -> VivoTopologyBuildReport {
        _ = try VivoStructureValidator.validate(structure)
        if inferHydrogenAttachment || inferDisulfides {
            guard structure.conformers.indices.contains(conformerIndex) else {
                throw VivoArtifactValidationError.unresolved("topology reconstruction requires conformer \(conformerIndex)")
            }
        }
        var existing = Set<UInt64>()
        for bond in structure.bonds { existing.insert(key(bond.atomA, bond.atomB)) }
        var templateAdded: UInt32 = 0, peptideAdded: UInt32 = 0, disulfideAdded: UInt32 = 0, hydrogenAdded: UInt32 = 0
        var unresolved: [UInt32] = []

        func add(_ a: UInt32, _ b: UInt32, order: VivoBondOrder = .single) {
            guard a != b, existing.insert(key(a, b)).inserted else { return }
            structure.bonds.append(.init(atomA: a, atomB: b, order: order))
        }

        var nameMaps: [[String: UInt32]] = []
        nameMaps.reserveCapacity(structure.residues.count)
        for residue in structure.residues {
            var map: [String: UInt32] = [:]
            for atomIndex in residue.atomIndices {
                let atom = structure.atoms[Int(atomIndex)]
                guard map.updateValue(atomIndex, forKey: atom.name.uppercased()) == nil else {
                    throw VivoArtifactValidationError.invalid("residue \(residue.index) has duplicate atom name '\(atom.name)'; resolve alternate locations before deterministic topology reconstruction")
                }
            }
            nameMaps.append(map)
            let canonical = canonicalResidue(residue.name)
            guard let template = templates[canonical] else { unresolved.append(residue.index); continue }
            for spec in template {
                guard let a = map[spec.a], let b = map[spec.b] else { continue }
                let before = existing.count
                add(a, b, order: spec.order)
                if existing.count > before { templateAdded += 1 }
            }
            // Standard peptide backbone shared by all amino acids.
            for spec in backbone {
                guard let a = map[spec.a], let b = map[spec.b] else { continue }
                let before = existing.count
                add(a, b, order: spec.order)
                if existing.count > before { templateAdded += 1 }
            }
        }

        // Peptide links follow chain residue ordering, never mere spatial proximity.
        for chain in structure.chains {
            for pair in zip(chain.residueIndices, chain.residueIndices.dropFirst()) {
                let left = structure.residues[Int(pair.0)], right = structure.residues[Int(pair.1)]
                guard templates[canonicalResidue(left.name)] != nil, templates[canonicalResidue(right.name)] != nil,
                      let c = nameMaps[Int(pair.0)]["C"], let n = nameMaps[Int(pair.1)]["N"] else { continue }
                let before = existing.count
                add(c, n)
                if existing.count > before { peptideAdded += 1 }
            }
        }

        let positions = structure.conformers.indices.contains(conformerIndex) ? structure.conformers[conformerIndex].positionsNM : []
        if inferDisulfides {
            var sulfurs: [UInt32] = []
            for residue in structure.residues where ["CYS", "CYX"].contains(canonicalResidue(residue.name)) {
                if let sg = nameMaps[Int(residue.index)]["SG"] { sulfurs.append(sg) }
            }
            for i in sulfurs.indices {
                for j in sulfurs.indices where j > i {
                    let a = sulfurs[i], b = sulfurs[j]
                    var delta = positions[Int(a)] - positions[Int(b)]
                    if let cell = structure.periodicCell { delta = try cell.minimumImage(delta) }
                    if delta.norm <= 0.235 {
                        let before = existing.count
                        add(a, b)
                        if existing.count > before { disulfideAdded += 1 }
                    }
                }
            }
        }

        if inferHydrogenAttachment {
            for residue in structure.residues {
                let atomSet = Set(residue.atomIndices)
                let heavy = residue.atomIndices.filter { structure.atoms[Int($0)].element.atomicNumber != 1 }
                for hydrogen in residue.atomIndices where structure.atoms[Int(hydrogen)].element.atomicNumber == 1 {
                    if structure.bonds.contains(where: { $0.atomA == hydrogen || $0.atomB == hydrogen }) { continue }
                    var best: (UInt32, Double)?
                    for candidate in heavy where atomSet.contains(candidate) {
                        var delta = positions[Int(hydrogen)] - positions[Int(candidate)]
                        if let cell = structure.periodicCell { delta = try cell.minimumImage(delta) }
                        let distance = delta.norm
                        let cutoff = hydrogenCutoffNM(for: structure.atoms[Int(candidate)].element.atomicNumber)
                        if distance <= cutoff, distance > 0.04, best == nil || distance < best!.1 { best = (candidate, distance) }
                    }
                    if let best {
                        let before = existing.count
                        add(hydrogen, best.0)
                        if existing.count > before { hydrogenAdded += 1 }
                    }
                }
            }
        }

        _ = try VivoStructureValidator.validate(structure)
        return .init(templateBondsAdded: templateAdded, peptideBondsAdded: peptideAdded,
                     disulfideBondsAdded: disulfideAdded, hydrogenBondsAdded: hydrogenAdded,
                     unresolvedResidues: unresolved)
    }

    private struct BondSpec {
        var a: String; var b: String; var order: VivoBondOrder = .single
    }
    private static let backbone: [BondSpec] = [
        .init(a: "N", b: "CA"), .init(a: "CA", b: "C"),
        .init(a: "C", b: "O", order: .double), .init(a: "C", b: "OXT")
    ]
    private static let templates: [String: [BondSpec]] = [
        "ALA": [.init(a:"CA",b:"CB")],
        "ARG": [.init(a:"CA",b:"CB"),.init(a:"CB",b:"CG"),.init(a:"CG",b:"CD"),.init(a:"CD",b:"NE"),.init(a:"NE",b:"CZ"),.init(a:"CZ",b:"NH1"),.init(a:"CZ",b:"NH2")],
        "ASN": [.init(a:"CA",b:"CB"),.init(a:"CB",b:"CG"),.init(a:"CG",b:"OD1",order:.double),.init(a:"CG",b:"ND2")],
        "ASP": [.init(a:"CA",b:"CB"),.init(a:"CB",b:"CG"),.init(a:"CG",b:"OD1",order:.double),.init(a:"CG",b:"OD2")],
        "CYS": [.init(a:"CA",b:"CB"),.init(a:"CB",b:"SG")],
        "GLN": [.init(a:"CA",b:"CB"),.init(a:"CB",b:"CG"),.init(a:"CG",b:"CD"),.init(a:"CD",b:"OE1",order:.double),.init(a:"CD",b:"NE2")],
        "GLU": [.init(a:"CA",b:"CB"),.init(a:"CB",b:"CG"),.init(a:"CG",b:"CD"),.init(a:"CD",b:"OE1",order:.double),.init(a:"CD",b:"OE2")],
        "GLY": [],
        "HIS": [.init(a:"CA",b:"CB"),.init(a:"CB",b:"CG"),.init(a:"CG",b:"ND1",order:.aromatic),.init(a:"ND1",b:"CE1",order:.aromatic),.init(a:"CE1",b:"NE2",order:.aromatic),.init(a:"NE2",b:"CD2",order:.aromatic),.init(a:"CD2",b:"CG",order:.aromatic)],
        "ILE": [.init(a:"CA",b:"CB"),.init(a:"CB",b:"CG1"),.init(a:"CG1",b:"CD1"),.init(a:"CB",b:"CG2")],
        "LEU": [.init(a:"CA",b:"CB"),.init(a:"CB",b:"CG"),.init(a:"CG",b:"CD1"),.init(a:"CG",b:"CD2")],
        "LYS": [.init(a:"CA",b:"CB"),.init(a:"CB",b:"CG"),.init(a:"CG",b:"CD"),.init(a:"CD",b:"CE"),.init(a:"CE",b:"NZ")],
        "MET": [.init(a:"CA",b:"CB"),.init(a:"CB",b:"CG"),.init(a:"CG",b:"SD"),.init(a:"SD",b:"CE")],
        "PHE": [.init(a:"CA",b:"CB"),.init(a:"CB",b:"CG"),.init(a:"CG",b:"CD1",order:.aromatic),.init(a:"CD1",b:"CE1",order:.aromatic),.init(a:"CE1",b:"CZ",order:.aromatic),.init(a:"CZ",b:"CE2",order:.aromatic),.init(a:"CE2",b:"CD2",order:.aromatic),.init(a:"CD2",b:"CG",order:.aromatic)],
        "PRO": [.init(a:"CA",b:"CB"),.init(a:"CB",b:"CG"),.init(a:"CG",b:"CD"),.init(a:"CD",b:"N")],
        "SER": [.init(a:"CA",b:"CB"),.init(a:"CB",b:"OG")],
        "THR": [.init(a:"CA",b:"CB"),.init(a:"CB",b:"OG1"),.init(a:"CB",b:"CG2")],
        "TRP": [.init(a:"CA",b:"CB"),.init(a:"CB",b:"CG"),.init(a:"CG",b:"CD1",order:.aromatic),.init(a:"CD1",b:"NE1",order:.aromatic),.init(a:"NE1",b:"CE2",order:.aromatic),.init(a:"CE2",b:"CD2",order:.aromatic),.init(a:"CD2",b:"CG",order:.aromatic),.init(a:"CE2",b:"CZ2",order:.aromatic),.init(a:"CZ2",b:"CH2",order:.aromatic),.init(a:"CH2",b:"CZ3",order:.aromatic),.init(a:"CZ3",b:"CE3",order:.aromatic),.init(a:"CE3",b:"CD2",order:.aromatic)],
        "TYR": [.init(a:"CA",b:"CB"),.init(a:"CB",b:"CG"),.init(a:"CG",b:"CD1",order:.aromatic),.init(a:"CD1",b:"CE1",order:.aromatic),.init(a:"CE1",b:"CZ",order:.aromatic),.init(a:"CZ",b:"CE2",order:.aromatic),.init(a:"CE2",b:"CD2",order:.aromatic),.init(a:"CD2",b:"CG",order:.aromatic),.init(a:"CZ",b:"OH")],
        "VAL": [.init(a:"CA",b:"CB"),.init(a:"CB",b:"CG1"),.init(a:"CB",b:"CG2")]
    ]

    private static func canonicalResidue(_ raw: String) -> String {
        switch raw.uppercased() {
        case "HID", "HIE", "HIP": "HIS"
        case "ASH": "ASP"
        case "GLH": "GLU"
        case "CYX", "CYM": "CYS"
        case "LYN": "LYS"
        default: raw.uppercased()
        }
    }
    private static func key(_ a: UInt32, _ b: UInt32) -> UInt64 {
        UInt64(min(a,b)) << 32 | UInt64(max(a,b))
    }
    private static func hydrogenCutoffNM(for atomicNumber: UInt16) -> Double {
        switch atomicNumber {
        case 7: 0.125
        case 8: 0.120
        case 9: 0.115
        case 15, 16: 0.150
        default: 0.130
        }
    }
}

public struct VivoCovalentPerceptionOptions: Codable, Sendable, Equatable {
    public var distanceScale: Double
    public var minimumDistanceNM: Double
    public var allowInterResidue: Bool
    public var includeMetals: Bool
    public init(distanceScale: Double = 1.18, minimumDistanceNM: Double = 0.04,
                allowInterResidue: Bool = false, includeMetals: Bool = false) {
        self.distanceScale = distanceScale; self.minimumDistanceNM = minimumDistanceNM
        self.allowInterResidue = allowInterResidue; self.includeMetals = includeMetals
    }
}

public enum VivoCovalentBondPerception {
    /// Conservative distance perception for ligands/unknown residues. Existing
    /// bonds are preserved. It is intentionally opt-in for inter-residue bonds.
    public static func apply(to structure: inout VivoMolecularStructure,
                             conformerIndex: Int = 0,
                             options: VivoCovalentPerceptionOptions = .init()) throws -> UInt32 {
        _ = try VivoStructureValidator.validate(structure)
        guard structure.conformers.indices.contains(conformerIndex),
              options.distanceScale.isFinite, options.distanceScale >= 1,
              options.minimumDistanceNM.isFinite, options.minimumDistanceNM > 0 else {
            throw VivoArtifactValidationError.invalid("invalid covalent-bond perception configuration or conformer")
        }
        let coordinates = structure.conformers[conformerIndex].positionsNM
        var existing = Set(structure.bonds.map { UInt64(min($0.atomA,$0.atomB)) << 32 | UInt64(max($0.atomA,$0.atomB)) })
        var degree = [Int](repeating: 0, count: structure.atoms.count)
        for bond in structure.bonds { degree[Int(bond.atomA)] += 1; degree[Int(bond.atomB)] += 1 }
        var candidates: [(Double, UInt32, UInt32)] = []
        for i in structure.atoms.indices {
            for j in (i + 1)..<structure.atoms.count {
                let a = structure.atoms[i], b = structure.atoms[j]
                if !options.allowInterResidue, a.residueIndex != b.residueIndex { continue }
                if !options.includeMetals, isMetal(a.element.atomicNumber) || isMetal(b.element.atomicNumber) { continue }
                let key = UInt64(UInt32(i)) << 32 | UInt64(UInt32(j))
                if existing.contains(key) { continue }
                guard let ra = covalentRadiusNM[a.element.atomicNumber], let rb = covalentRadiusNM[b.element.atomicNumber] else { continue }
                var delta = coordinates[i] - coordinates[j]
                if let cell = structure.periodicCell { delta = try cell.minimumImage(delta) }
                let distance = delta.norm
                if distance >= options.minimumDistanceNM, distance <= (ra + rb) * options.distanceScale {
                    candidates.append((distance, UInt32(i), UInt32(j)))
                }
            }
        }
        candidates.sort { $0.0 < $1.0 }
        var added: UInt32 = 0
        for (_, a, b) in candidates {
            guard degree[Int(a)] < maximumValence(structure.atoms[Int(a)].element.atomicNumber),
                  degree[Int(b)] < maximumValence(structure.atoms[Int(b)].element.atomicNumber) else { continue }
            let k = UInt64(min(a,b)) << 32 | UInt64(max(a,b))
            guard existing.insert(k).inserted else { continue }
            structure.bonds.append(.init(atomA: a, atomB: b, order: .unknown))
            degree[Int(a)] += 1; degree[Int(b)] += 1; added += 1
        }
        _ = try VivoStructureValidator.validate(structure)
        return added
    }

    private static func maximumValence(_ z: UInt16) -> Int {
        switch z { case 1,9,17,35,53: 1; case 8: 2; case 7,15: 4; case 6: 4; case 16: 6; default: 6 }
    }
    private static func isMetal(_ z: UInt16) -> Bool {
        switch z { case 3,4,11,12,13,19,20,21...31,37,38,39...50,55,56,57...84,87,88,89...112: true; default: false }
    }
    private static let covalentRadiusNM: [UInt16: Double] = [
        1:0.031, 5:0.084, 6:0.076, 7:0.071, 8:0.066, 9:0.057,
        14:0.111, 15:0.107, 16:0.105, 17:0.102, 35:0.120, 53:0.139
    ]
}
