import Foundation

public struct VivoClassicalInitialState: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/classical-initial-state/v1"
    public var schema: String
    public var systemFingerprint: VivoFingerprint
    public var positionsNM: [VivoVector3D]
    public var periodicCell: VivoPeriodicCell?
    public var sourceTimePS: Double?

    public init(systemFingerprint: VivoFingerprint, positionsNM: [VivoVector3D],
                periodicCell: VivoPeriodicCell? = nil, sourceTimePS: Double? = nil) {
        self.schema = Self.schema; self.systemFingerprint = systemFingerprint
        self.positionsNM = positionsNM; self.periodicCell = periodicCell; self.sourceTimePS = sourceTimePS
    }

    public func validate(particleCount: Int) throws {
        guard schema == Self.schema, positionsNM.count == particleCount,
              positionsNM.allSatisfy(\.isFinite), periodicCell?.isValid != false,
              sourceTimePS?.isFinite != false else {
            throw VivoArtifactValidationError.invalid("classical initial state shape, coordinates, cell, or time is invalid")
        }
    }
}

public struct VivoAmberImportResult: Codable, Sendable, Equatable {
    public var structure: VivoMolecularStructureDocument
    public var system: VivoClassicalSystem
    public var initialState: VivoClassicalInitialState?
    /// AMBER particle index -> physical structure atom index. nil is a virtual site.
    public var particleToStructureAtom: [UInt32?]
    public var amberAtomTypes: [String]
}

public enum VivoAmberImporter {
    private static let kcalToKJ = 4.184
    private static let amberChargeScale = 18.2223

    public static func importSystem(prmtopData: Data, restartData: Data? = nil,
                                    identifier: String = "amber") throws -> VivoAmberImportResult {
        let top = try VivoAmberPrmtop(data: prmtopData)
        let names = try top.strings("ATOM_NAME")
        let charges = try top.reals("CHARGE")
        let masses = try top.reals("MASS")
        let atomicNumbers = try top.integers("ATOMIC_NUMBER")
        let ljIndices = try top.integers("ATOM_TYPE_INDEX")
        let amberTypes = try top.strings("AMBER_ATOM_TYPE")
        let natom = names.count
        guard natom > 0, charges.count == natom, masses.count == natom,
              atomicNumbers.count == natom, ljIndices.count == natom, amberTypes.count == natom else {
            throw VivoArtifactValidationError.invalid("AMBER atom arrays disagree on NATOM")
        }
        let restart = try restartData.map(VivoAmberRestart.init(data:))
        if let restart, restart.atomCount != UInt32(natom) {
            throw VivoArtifactValidationError.incompatible("AMBER restart NATOM does not match prmtop")
        }

        let residueLabels = try top.strings("RESIDUE_LABEL")
        let residuePointers = try top.integers("RESIDUE_POINTER")
        guard !residueLabels.isEmpty, residueLabels.count == residuePointers.count,
              residuePointers.first == 1 else {
            throw VivoArtifactValidationError.invalid("AMBER residue arrays are invalid")
        }
        for i in residuePointers.indices {
            guard residuePointers[i] > 0, residuePointers[i] <= natom,
                  i == 0 || residuePointers[i] > residuePointers[i - 1] else {
                throw VivoArtifactValidationError.invalid("AMBER RESIDUE_POINTER must be strictly increasing and 1-based")
            }
        }
        let chainIDs = try top.strings("RESIDUE_CHAINID", required: false)
        if !chainIDs.isEmpty, chainIDs.count != residueLabels.count {
            throw VivoArtifactValidationError.invalid("AMBER RESIDUE_CHAINID count differs from residue count")
        }
        var residueForParticle = [Int](repeating: 0, count: natom)
        for residue in residuePointers.indices {
            let start = residuePointers[residue] - 1
            let end = residue + 1 < residuePointers.count ? residuePointers[residue + 1] - 1 : natom
            guard start < end else { throw VivoArtifactValidationError.invalid("AMBER residue has no particles") }
            for particle in start..<end { residueForParticle[particle] = residue }
        }

        let builder = VivoStructureAssemblyBuilder()
        var particleToAtom = [UInt32?](repeating: nil, count: natom)
        for particle in 0..<natom {
            let z = atomicNumbers[particle]
            if z <= 0 {
                guard masses[particle] >= 0, masses[particle] < 1e-6 else {
                    throw VivoArtifactValidationError.incompatible("AMBER particle \(particle) has no atomic number but nonzero mass; element semantics are unresolved")
                }
                continue
            }
            guard z <= 118, let element = VivoElement.from(symbol: elementSymbol(z)) else {
                throw VivoArtifactValidationError.invalid("AMBER particle \(particle) has invalid atomic number \(z)")
            }
            let residue = residueForParticle[particle]
            let position = restart?.positionsNM[particle] ?? .zero
            let atomIndex = try builder.appendAtom(
                name: names[particle], element: element, positionNM: position,
                residueName: residueLabels[residue],
                chainIdentifier: chainIDs.isEmpty ? nil : chainIDs[residue],
                residueSequence: Int32(exactly: residue + 1),
                sourceSerial: Int32(exactly: particle + 1), isHetero: false
            )
            particleToAtom[particle] = atomIndex
        }
        guard !builder.atoms.isEmpty else {
            throw VivoArtifactValidationError.invalid("AMBER topology contains no physical atoms")
        }

        let bondRecords = try bondedRecords(top, hydrogenFlag: "BONDS_INC_HYDROGEN",
                                            heavyFlag: "BONDS_WITHOUT_HYDROGEN", width: 3)
        var structureBonds: [VivoMolecularBond] = []
        for record in bondRecords {
            let a = try amberAtom(record[0], natom: natom)
            let b = try amberAtom(record[1], natom: natom)
            if let aa = particleToAtom[a], let bb = particleToAtom[b] {
                structureBonds.append(.init(atomA: aa, atomB: bb, order: .unknown))
            }
        }
        var structure = try builder.makeStructure(
            identifier: identifier, bonds: structureBonds,
            conformers: restart == nil ? [] : [.init(identifier: "amber-restart", positionsNM: builder.positions)],
            periodicCell: restart?.periodicCell,
            metadata: ["sourceFormat": "AMBER-prmtop", "particleCount": String(natom)]
        )
        // Source serials point to AMBER particle IDs, which intentionally include
        // gaps at virtual-site positions.
        _ = try VivoStructureValidator.validate(structure)
        let structureDocument = try VivoMolecularStructureDocument(
            structure: structure, sourceFingerprint: try top.fingerprint(), sourceFormat: nil
        )

        let usedLJIndices = Set(ljIndices)
        guard usedLJIndices.allSatisfy({ $0 > 0 }) else {
            throw VivoArtifactValidationError.invalid("AMBER ATOM_TYPE_INDEX values must be positive")
        }
        let ntypes = ljIndices.max() ?? 0
        let nbIndex = try top.integers("NONBONDED_PARM_INDEX")
        guard ntypes > 0, nbIndex.count >= ntypes * ntypes else {
            throw VivoArtifactValidationError.invalid("AMBER NONBONDED_PARM_INDEX is smaller than the active LJ type matrix")
        }
        let aCoefficients = try top.reals("LENNARD_JONES_ACOEF")
        let bCoefficients = try top.reals("LENNARD_JONES_BCOEF")
        var typePairs: [VivoNonbondedTypePair] = []
        let sortedTypes = usedLJIndices.sorted()
        for (offsetI, i) in sortedTypes.enumerated() {
            for j in sortedTypes[offsetI...] {
                let parameterIndex = nbIndex[(i - 1) * ntypes + (j - 1)]
                guard parameterIndex > 0, parameterIndex <= aCoefficients.count,
                      parameterIndex <= bCoefficients.count else {
                    throw VivoArtifactValidationError.incompatible("AMBER active LJ pair uses unsupported non-positive or absent parameter index")
                }
                let converted = convertLJ(a: aCoefficients[parameterIndex - 1], b: bCoefficients[parameterIndex - 1])
                typePairs.append(.init(typeA: typeID(i), typeB: typeID(j),
                                       c12KJNM12PerMol: converted.c12,
                                       c6KJNM6PerMol: converted.c6))
            }
        }
        var diagonal: [Int: (sigma: Double, epsilon: Double)] = [:]
        for i in sortedTypes {
            let parameterIndex = nbIndex[(i - 1) * ntypes + (i - 1)]
            guard parameterIndex > 0, parameterIndex <= aCoefficients.count,
                  parameterIndex <= bCoefficients.count else {
                throw VivoArtifactValidationError.invalid("AMBER diagonal LJ parameter is unavailable")
            }
            let converted = convertLJ(a: aCoefficients[parameterIndex - 1], b: bCoefficients[parameterIndex - 1])
            if converted.c12 > 0, converted.c6 > 0 {
                let sigma = pow(converted.c12 / converted.c6, 1.0 / 6.0)
                let epsilon = converted.c6 * converted.c6 / (4 * converted.c12)
                diagonal[i] = (sigma, epsilon)
            } else {
                diagonal[i] = (0, 0)
            }
        }

        var particles: [VivoClassicalParticle] = []
        particles.reserveCapacity(natom)
        for particle in 0..<natom {
            let type = ljIndices[particle]
            let diag = diagonal[type] ?? (0, 0)
            particles.append(.init(index: UInt32(particle), atomIndex: particleToAtom[particle],
                                   typeIdentifier: typeID(type),
                                   role: particleToAtom[particle] == nil ? .virtualSite : .atom,
                                   massDa: masses[particle], chargeE: charges[particle] / amberChargeScale,
                                   sigmaNM: diag.sigma, epsilonKJPerMol: diag.epsilon))
        }

        let bondK = try top.reals("BOND_FORCE_CONSTANT")
        let bondR0 = try top.reals("BOND_EQUIL_VALUE")
        var bonds: [VivoHarmonicBond] = []
        for record in bondRecords {
            let p = try parameterIndex(record[2], count: min(bondK.count, bondR0.count), label: "bond")
            bonds.append(.init(a: UInt32(try amberAtom(record[0], natom: natom)),
                               b: UInt32(try amberAtom(record[1], natom: natom)),
                               lengthNM: bondR0[p] * 0.1,
                               forceConstant: 2 * bondK[p] * kcalToKJ * 100))
        }

        let angleRecords = try bondedRecords(top, hydrogenFlag: "ANGLES_INC_HYDROGEN",
                                             heavyFlag: "ANGLES_WITHOUT_HYDROGEN", width: 4)
        let angleK = try top.reals("ANGLE_FORCE_CONSTANT")
        let angle0 = try top.reals("ANGLE_EQUIL_VALUE")
        var angles: [VivoHarmonicAngle] = []
        for record in angleRecords {
            let p = try parameterIndex(record[3], count: min(angleK.count, angle0.count), label: "angle")
            angles.append(.init(a: UInt32(try amberAtom(record[0], natom: natom)),
                                b: UInt32(try amberAtom(record[1], natom: natom)),
                                c: UInt32(try amberAtom(record[2], natom: natom)),
                                angleRadians: angle0[p],
                                forceConstant: 2 * angleK[p] * kcalToKJ))
        }

        let dihedralRecords = try bondedRecords(top, hydrogenFlag: "DIHEDRALS_INC_HYDROGEN",
                                                heavyFlag: "DIHEDRALS_WITHOUT_HYDROGEN", width: 5)
        let torsionK = try top.reals("DIHEDRAL_FORCE_CONSTANT")
        let periodicity = try top.reals("DIHEDRAL_PERIODICITY")
        let phase = try top.reals("DIHEDRAL_PHASE")
        let scee = try top.reals("SCEE_SCALE_FACTOR", required: false)
        let scnb = try top.reals("SCNB_SCALE_FACTOR", required: false)
        var torsions: [VivoPeriodicTorsion] = []
        var active14: [UInt64: (Double, Double)] = [:]
        for record in dihedralRecords {
            let p = try parameterIndex(record[4], count: min(torsionK.count, periodicity.count, phase.count), label: "dihedral")
            let a = try amberAtom(record[0], natom: natom)
            let b = try amberAtom(record[1], natom: natom)
            let c = try amberAtom(record[2], natom: natom)
            let d = try amberAtom(record[3], natom: natom)
            let improper = record[3] < 0
            let ignoreEnd = record[2] < 0 || improper
            let nRaw = abs(periodicity[p])
            let rounded = nRaw.rounded()
            guard nRaw > 0, abs(nRaw - rounded) < 1e-6, rounded <= Double(UInt16.max) else {
                throw VivoArtifactValidationError.incompatible("AMBER dihedral periodicity is not a positive integer")
            }
            torsions.append(.init(a: UInt32(a), b: UInt32(b), c: UInt32(c), d: UInt32(d),
                                  periodicity: UInt16(rounded), phaseRadians: phase[p],
                                  barrierKJPerMol: torsionK[p] * kcalToKJ, improper: improper))
            if !ignoreEnd {
                guard p < scee.count, p < scnb.count, scee[p] > 0, scnb[p] > 0 else {
                    throw VivoArtifactValidationError.unresolved("AMBER proper dihedral requires explicit SCEE/SCNB scale factors")
                }
                let key = pairKey(UInt32(a), UInt32(d))
                let scales = (1.0 / scee[p], 1.0 / scnb[p])
                if let existing = active14[key], abs(existing.0 - scales.0) > 1e-12 || abs(existing.1 - scales.1) > 1e-12 {
                    throw VivoArtifactValidationError.incompatible("AMBER gives conflicting 1-4 scales for the same particle pair")
                }
                active14[key] = scales
            }
        }

        let exclusions = try amberExclusions(top, natom: natom)
        var exceptionMap: [UInt64: VivoNonbondedException] = [:]
        for key in exclusions {
            let pair = decodePair(key)
            exceptionMap[key] = .init(a: pair.0, b: pair.1, coulombScale: 0, lennardJonesScale: 0)
        }
        for (key, scale) in active14 {
            let pair = decodePair(key)
            exceptionMap[key] = .init(a: pair.0, b: pair.1,
                                      coulombScale: scale.0, lennardJonesScale: scale.1)
        }

        let topFingerprint = try top.fingerprint()
        var system = VivoClassicalSystem(
            identifier: identifier + ":amber", structureFingerprint: structureDocument.structureFingerprint,
            parameterSourceFingerprints: [topFingerprint], mixingRule: .explicitPairTable,
            particles: particles, bonds: bonds, angles: angles, torsions: torsions,
            constraints: [], nonbondedTypePairs: typePairs,
            nonbondedExceptions: exceptionMap.values.sorted { lhs, rhs in pairKey(lhs.a, lhs.b) < pairKey(rhs.a, rhs.b) },
            metadata: ["sourceFormat": "AMBER-prmtop",
                       "amberVersion": top.versionLine ?? "unknown",
                       "chargeScale": String(amberChargeScale)]
        )
        try VivoClassicalSystemValidator.validate(system, atomCount: UInt32(structure.atoms.count))
        let systemFingerprint = try system.fingerprint()
        let initialState: VivoClassicalInitialState?
        if let restart {
            let value = VivoClassicalInitialState(systemFingerprint: systemFingerprint,
                                                  positionsNM: restart.positionsNM,
                                                  periodicCell: restart.periodicCell,
                                                  sourceTimePS: restart.timePS)
            try value.validate(particleCount: natom)
            initialState = value
            if restart.rawVelocityValues != nil {
                system.metadata["amberRestartVelocities"] = "present-but-not-converted-wave-B-boundary"
            }
        } else { initialState = nil }

        return .init(structure: structureDocument, system: system, initialState: initialState,
                     particleToStructureAtom: particleToAtom, amberAtomTypes: amberTypes)
    }

    private static func bondedRecords(_ top: VivoAmberPrmtop, hydrogenFlag: String,
                                      heavyFlag: String, width: Int) throws -> [[Int]] {
        let values = try top.integers(hydrogenFlag, required: false) + top.integers(heavyFlag, required: false)
        guard values.count % width == 0 else {
            throw VivoArtifactValidationError.invalid("AMBER bonded section width is inconsistent")
        }
        return stride(from: 0, to: values.count, by: width).map { Array(values[$0..<$0 + width]) }
    }

    private static func amberAtom(_ encoded: Int, natom: Int) throws -> Int {
        let magnitude = abs(encoded)
        guard magnitude % 3 == 0 else {
            throw VivoArtifactValidationError.invalid("AMBER bonded atom pointer \(encoded) is not divisible by 3")
        }
        let atom = magnitude / 3
        guard atom >= 0, atom < natom else {
            throw VivoArtifactValidationError.invalid("AMBER bonded atom pointer is out of range")
        }
        return atom
    }

    private static func parameterIndex(_ oneBased: Int, count: Int, label: String) throws -> Int {
        guard oneBased > 0, oneBased <= count else {
            throw VivoArtifactValidationError.invalid("AMBER \(label) parameter pointer \(oneBased) is out of range")
        }
        return oneBased - 1
    }

    private static func amberExclusions(_ top: VivoAmberPrmtop, natom: Int) throws -> Set<UInt64> {
        let counts = try top.integers("NUMBER_EXCLUDED_ATOMS")
        let list = try top.integers("EXCLUDED_ATOMS_LIST")
        guard counts.count == natom, counts.allSatisfy({ $0 >= 0 }) else {
            throw VivoArtifactValidationError.invalid("AMBER exclusion counts are invalid")
        }
        var cursor = 0
        var result = Set<UInt64>()
        for atom in 0..<natom {
            let amount = counts[atom]
            guard cursor + amount <= list.count else {
                throw VivoArtifactValidationError.invalid("AMBER excluded-atoms list is truncated")
            }
            for raw in list[cursor..<cursor + amount] where raw != 0 {
                guard raw > 0, raw <= natom else {
                    throw VivoArtifactValidationError.invalid("AMBER excluded atom index is out of range")
                }
                let other = raw - 1
                if other != atom { result.insert(pairKey(UInt32(atom), UInt32(other))) }
            }
            cursor += amount
        }
        guard cursor == list.count || list[cursor...].allSatisfy({ $0 == 0 }) else {
            throw VivoArtifactValidationError.invalid("AMBER excluded-atoms list has unexplained trailing entries")
        }
        return result
    }

    private static func convertLJ(a: Double, b: Double) -> (c12: Double, c6: Double) {
        (a * kcalToKJ * pow(0.1, 12), b * kcalToKJ * pow(0.1, 6))
    }

    private static func typeID(_ amberIndex: Int) -> String { "amber-lj-\(amberIndex)" }
    private static func pairKey(_ a: UInt32, _ b: UInt32) -> UInt64 { UInt64(min(a,b)) << 32 | UInt64(max(a,b)) }
    private static func decodePair(_ value: UInt64) -> (UInt32, UInt32) {
        (UInt32(value >> 32), UInt32(value & 0xffff_ffff))
    }

    private static func elementSymbol(_ atomicNumber: Int) -> String {
        let symbols = ["", "H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn","Nh","Fl","Mc","Lv","Ts","Og"]
        return atomicNumber > 0 && atomicNumber < symbols.count ? symbols[atomicNumber] : ""
    }
}
