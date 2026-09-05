import Foundation

// MARK: - Stable molecular identity and geometry

/// NumiVivo stores molecular coordinates in nanometres. File adapters are
/// responsible for explicit conversion at the boundary (PDB/mmCIF typically Å).
public struct VivoVector3D: Codable, Sendable, Equatable, Hashable {
    public var x: Double
    public var y: Double
    public var z: Double

    public init(_ x: Double, _ y: Double, _ z: Double) {
        self.x = x; self.y = y; self.z = z
    }

    public static let zero = VivoVector3D(0, 0, 0)
    public var isFinite: Bool { x.isFinite && y.isFinite && z.isFinite }

    public static func +(lhs: Self, rhs: Self) -> Self { .init(lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z) }
    public static func -(lhs: Self, rhs: Self) -> Self { .init(lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z) }
    public static func *(lhs: Self, rhs: Double) -> Self { .init(lhs.x * rhs, lhs.y * rhs, lhs.z * rhs) }
    public static func /(lhs: Self, rhs: Double) -> Self { .init(lhs.x / rhs, lhs.y / rhs, lhs.z / rhs) }
    public func dot(_ other: Self) -> Double { x * other.x + y * other.y + z * other.z }
    public func cross(_ other: Self) -> Self {
        .init(y * other.z - z * other.y,
              z * other.x - x * other.z,
              x * other.y - y * other.x)
    }
    public var squaredNorm: Double { dot(self) }
    public var norm: Double { sqrt(squaredNorm) }
}

public struct VivoPeriodicCell: Codable, Sendable, Equatable, Hashable {
    /// Triclinic lattice vectors in nanometres.
    public var a: VivoVector3D
    public var b: VivoVector3D
    public var c: VivoVector3D

    public init(a: VivoVector3D, b: VivoVector3D, c: VivoVector3D) {
        self.a = a; self.b = b; self.c = c
    }

    public var volumeNM3: Double { abs(a.dot(b.cross(c))) }
    public var isValid: Bool {
        a.isFinite && b.isFinite && c.isFinite && volumeNM3.isFinite && volumeNM3 > 1e-18
    }
}

public struct VivoElement: Codable, Sendable, Equatable, Hashable {
    public let atomicNumber: UInt16
    public let symbol: String

    public init(atomicNumber: UInt16, symbol: String) {
        self.atomicNumber = atomicNumber
        self.symbol = symbol
    }

    public static func from(symbol raw: String) -> VivoElement? {
        let key = raw.trimmingCharacters(in: .whitespacesAndNewlines).lowercased()
        guard let z = periodicNumbers[key] else { return nil }
        return .init(atomicNumber: z, symbol: canonicalSymbols[Int(z)])
    }

    private static let periodicNumbers: [String: UInt16] = {
        var result: [String: UInt16] = [:]
        for (index, symbol) in canonicalSymbols.enumerated() where index > 0 {
            result[symbol.lowercased()] = UInt16(index)
        }
        return result
    }()

    private static let canonicalSymbols: [String] = [
        "", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
        "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
        "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
        "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
        "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
        "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
        "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
        "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
        "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
        "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
        "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
        "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"
    ]
}

public enum VivoBondOrder: String, Codable, Sendable, CaseIterable {
    case unknown
    case single
    case double
    case triple
    case aromatic
    case dative
}

public enum VivoBondStereo: String, Codable, Sendable, CaseIterable {
    case unspecified
    case up
    case down
    case cis
    case trans
    case either
}

public struct VivoMolecularAtom: Codable, Sendable, Equatable {
    /// Dense immutable identity within one `VivoMolecularStructure`.
    public let index: UInt32
    public var name: String
    public var element: VivoElement
    public var isotopeMassNumber: UInt16?
    public var formalCharge: Int16
    public var residueIndex: UInt32?
    public var sourceSerial: Int32?
    public var alternateLocation: String?
    public var occupancy: Double?
    public var bFactor: Double?
    public var isHetero: Bool

    public init(index: UInt32, name: String, element: VivoElement,
                isotopeMassNumber: UInt16? = nil, formalCharge: Int16 = 0,
                residueIndex: UInt32? = nil, sourceSerial: Int32? = nil,
                alternateLocation: String? = nil, occupancy: Double? = nil,
                bFactor: Double? = nil, isHetero: Bool = false) {
        self.index = index
        self.name = name
        self.element = element
        self.isotopeMassNumber = isotopeMassNumber
        self.formalCharge = formalCharge
        self.residueIndex = residueIndex
        self.sourceSerial = sourceSerial
        self.alternateLocation = alternateLocation
        self.occupancy = occupancy
        self.bFactor = bFactor
        self.isHetero = isHetero
    }
}

public struct VivoMolecularBond: Codable, Sendable, Equatable, Hashable {
    public let atomA: UInt32
    public let atomB: UInt32
    public var order: VivoBondOrder
    public var stereo: VivoBondStereo

    public init(atomA: UInt32, atomB: UInt32, order: VivoBondOrder = .single,
                stereo: VivoBondStereo = .unspecified) {
        self.atomA = atomA; self.atomB = atomB; self.order = order; self.stereo = stereo
    }

    public var canonicalPair: (UInt32, UInt32) { atomA <= atomB ? (atomA, atomB) : (atomB, atomA) }
}

public struct VivoMolecularResidue: Codable, Sendable, Equatable {
    public let index: UInt32
    public var name: String
    public var chainIndex: UInt32?
    public var sequenceNumber: Int32?
    public var insertionCode: String?
    public var atomIndices: [UInt32]

    public init(index: UInt32, name: String, chainIndex: UInt32? = nil,
                sequenceNumber: Int32? = nil, insertionCode: String? = nil,
                atomIndices: [UInt32] = []) {
        self.index = index; self.name = name; self.chainIndex = chainIndex
        self.sequenceNumber = sequenceNumber; self.insertionCode = insertionCode
        self.atomIndices = atomIndices
    }
}

public struct VivoMolecularChain: Codable, Sendable, Equatable {
    public let index: UInt32
    public var identifier: String
    public var residueIndices: [UInt32]

    public init(index: UInt32, identifier: String, residueIndices: [UInt32] = []) {
        self.index = index; self.identifier = identifier; self.residueIndices = residueIndices
    }
}

public struct VivoMolecularConformer: Codable, Sendable, Equatable {
    public var identifier: String
    /// One position per atom, in nanometres and structure atom order.
    public var positionsNM: [VivoVector3D]
    public var potentialEnergyKJPerMol: Double?

    public init(identifier: String = "model-1", positionsNM: [VivoVector3D],
                potentialEnergyKJPerMol: Double? = nil) {
        self.identifier = identifier
        self.positionsNM = positionsNM
        self.potentialEnergyKJPerMol = potentialEnergyKJPerMol
    }
}

/// Authoritative topology + conformer container used by native MD, QM/MM and
/// structure adapters. Array position and public `index` are deliberately the
/// same: this gives stable O(1) indexing and an unambiguous GPU packing order.
public struct VivoMolecularStructure: Codable, Sendable, Equatable {
    public static let schemaVersion: UInt32 = 1

    public var schemaVersion: UInt32
    public var identifier: String
    public var atoms: [VivoMolecularAtom]
    public var bonds: [VivoMolecularBond]
    public var residues: [VivoMolecularResidue]
    public var chains: [VivoMolecularChain]
    public var conformers: [VivoMolecularConformer]
    public var periodicCell: VivoPeriodicCell?
    public var metadata: [String: String]

    public init(identifier: String,
                atoms: [VivoMolecularAtom], bonds: [VivoMolecularBond] = [],
                residues: [VivoMolecularResidue] = [], chains: [VivoMolecularChain] = [],
                conformers: [VivoMolecularConformer] = [], periodicCell: VivoPeriodicCell? = nil,
                metadata: [String: String] = [:]) {
        self.schemaVersion = Self.schemaVersion
        self.identifier = identifier
        self.atoms = atoms
        self.bonds = bonds
        self.residues = residues
        self.chains = chains
        self.conformers = conformers
        self.periodicCell = periodicCell
        self.metadata = metadata
    }
}
