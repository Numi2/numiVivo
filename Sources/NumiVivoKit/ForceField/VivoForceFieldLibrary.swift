import Foundation

public struct VivoForceFieldAtomType: Codable, Sendable, Equatable, Hashable {
    public var identifier: String
    public var elementAtomicNumber: UInt16?
    public var massDa: Double
    public var sigmaNM: Double
    public var epsilonKJPerMol: Double
    public init(identifier: String, elementAtomicNumber: UInt16? = nil,
                massDa: Double, sigmaNM: Double, epsilonKJPerMol: Double) {
        self.identifier = identifier; self.elementAtomicNumber = elementAtomicNumber
        self.massDa = massDa; self.sigmaNM = sigmaNM; self.epsilonKJPerMol = epsilonKJPerMol
    }
}

public struct VivoForceFieldTemplateAtom: Codable, Sendable, Equatable, Hashable {
    public var name: String
    public var typeIdentifier: String
    public var chargeE: Double
    public init(name: String, typeIdentifier: String, chargeE: Double) {
        self.name = name; self.typeIdentifier = typeIdentifier; self.chargeE = chargeE
    }
}

public struct VivoForceFieldTemplateBond: Codable, Sendable, Equatable, Hashable {
    public var atomA: String
    public var atomB: String
    public var order: VivoBondOrder
    public init(atomA: String, atomB: String, order: VivoBondOrder = .single) {
        self.atomA = atomA; self.atomB = atomB; self.order = order
    }
}

public struct VivoForceFieldResidueTemplate: Codable, Sendable, Equatable {
    public var residueNames: [String]
    public var atoms: [VivoForceFieldTemplateAtom]
    public var bonds: [VivoForceFieldTemplateBond]
    public init(residueNames: [String], atoms: [VivoForceFieldTemplateAtom], bonds: [VivoForceFieldTemplateBond] = []) {
        self.residueNames = residueNames; self.atoms = atoms; self.bonds = bonds
    }
}

public struct VivoForceFieldBondParameter: Codable, Sendable, Equatable, Hashable {
    public var typeA: String
    public var typeB: String
    public var lengthNM: Double
    public var forceConstant: Double
    public init(typeA: String, typeB: String, lengthNM: Double, forceConstant: Double) {
        self.typeA = typeA; self.typeB = typeB; self.lengthNM = lengthNM; self.forceConstant = forceConstant
    }
}

public struct VivoForceFieldAngleParameter: Codable, Sendable, Equatable, Hashable {
    public var typeA: String
    public var typeB: String
    public var typeC: String
    public var angleRadians: Double
    public var forceConstant: Double
    public init(typeA: String, typeB: String, typeC: String, angleRadians: Double, forceConstant: Double) {
        self.typeA = typeA; self.typeB = typeB; self.typeC = typeC
        self.angleRadians = angleRadians; self.forceConstant = forceConstant
    }
}

public struct VivoForceFieldTorsionParameter: Codable, Sendable, Equatable, Hashable {
    /// `*` may be used at terminal positions only. Exact matches always win.
    public var typeA: String
    public var typeB: String
    public var typeC: String
    public var typeD: String
    public var periodicity: UInt16
    public var phaseRadians: Double
    public var barrierKJPerMol: Double
    public var improper: Bool
    public init(typeA: String, typeB: String, typeC: String, typeD: String,
                periodicity: UInt16, phaseRadians: Double,
                barrierKJPerMol: Double, improper: Bool = false) {
        self.typeA = typeA; self.typeB = typeB; self.typeC = typeC; self.typeD = typeD
        self.periodicity = periodicity; self.phaseRadians = phaseRadians
        self.barrierKJPerMol = barrierKJPerMol; self.improper = improper
    }
}

public struct VivoForceFieldLibrary: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/force-field-library/v1"
    public var schema: String
    public var identifier: String
    public var version: String
    public var mixingRule: VivoMixingRule
    public var coulomb14Scale: Double
    public var lennardJones14Scale: Double
    public var atomTypes: [VivoForceFieldAtomType]
    public var residueTemplates: [VivoForceFieldResidueTemplate]
    public var bondParameters: [VivoForceFieldBondParameter]
    public var angleParameters: [VivoForceFieldAngleParameter]
    public var torsionParameters: [VivoForceFieldTorsionParameter]
    public var provenance: [String: String]

    public init(identifier: String, version: String, mixingRule: VivoMixingRule = .lorentzBerthelot,
                coulomb14Scale: Double = 1, lennardJones14Scale: Double = 1,
                atomTypes: [VivoForceFieldAtomType], residueTemplates: [VivoForceFieldResidueTemplate] = [],
                bondParameters: [VivoForceFieldBondParameter] = [], angleParameters: [VivoForceFieldAngleParameter] = [],
                torsionParameters: [VivoForceFieldTorsionParameter] = [], provenance: [String: String] = [:]) {
        self.schema = Self.schema; self.identifier = identifier; self.version = version
        self.mixingRule = mixingRule; self.coulomb14Scale = coulomb14Scale; self.lennardJones14Scale = lennardJones14Scale
        self.atomTypes = atomTypes; self.residueTemplates = residueTemplates
        self.bondParameters = bondParameters; self.angleParameters = angleParameters
        self.torsionParameters = torsionParameters; self.provenance = provenance
    }

    public func validate() throws {
        guard schema == Self.schema, !identifier.isEmpty, !version.isEmpty,
              coulomb14Scale.isFinite, coulomb14Scale >= 0,
              lennardJones14Scale.isFinite, lennardJones14Scale >= 0 else {
            throw VivoArtifactValidationError.invalid("force-field library header/scales are invalid")
        }
        var types = Set<String>()
        for type in atomTypes {
            guard !type.identifier.isEmpty, types.insert(type.identifier).inserted,
                  type.massDa.isFinite, type.massDa > 0,
                  type.sigmaNM.isFinite, type.sigmaNM >= 0,
                  type.epsilonKJPerMol.isFinite, type.epsilonKJPerMol >= 0 else {
                throw VivoArtifactValidationError.invalid("force-field atom type is duplicate or invalid")
            }
        }
        var residueAliases = Set<String>()
        for template in residueTemplates {
            guard !template.residueNames.isEmpty else { throw VivoArtifactValidationError.invalid("force-field residue template has no names") }
            for name in template.residueNames {
                guard residueAliases.insert(name.uppercased()).inserted else {
                    throw VivoArtifactValidationError.invalid("force-field residue alias '\(name)' is assigned to multiple templates")
                }
            }
            var names = Set<String>()
            for atom in template.atoms {
                guard names.insert(atom.name).inserted, types.contains(atom.typeIdentifier), atom.chargeE.isFinite else {
                    throw VivoArtifactValidationError.invalid("force-field template atom is duplicate, references unknown type, or has invalid charge")
                }
            }
            for bond in template.bonds where !names.contains(bond.atomA) || !names.contains(bond.atomB) || bond.atomA == bond.atomB {
                throw VivoArtifactValidationError.invalid("force-field template bond references absent/identical atoms")
            }
        }
        for bond in bondParameters {
            guard types.contains(bond.typeA), types.contains(bond.typeB),
                  bond.lengthNM.isFinite, bond.lengthNM > 0,
                  bond.forceConstant.isFinite, bond.forceConstant >= 0 else {
                throw VivoArtifactValidationError.invalid("force-field bond parameter is invalid")
            }
        }
        for angle in angleParameters {
            guard types.contains(angle.typeA), types.contains(angle.typeB), types.contains(angle.typeC),
                  angle.angleRadians.isFinite, angle.angleRadians > 0, angle.angleRadians < Double.pi,
                  angle.forceConstant.isFinite, angle.forceConstant >= 0 else {
                throw VivoArtifactValidationError.invalid("force-field angle parameter is invalid")
            }
        }
        for torsion in torsionParameters {
            let terminalsValid = (torsion.typeA == "*" || types.contains(torsion.typeA)) &&
                                 (torsion.typeD == "*" || types.contains(torsion.typeD))
            guard terminalsValid, types.contains(torsion.typeB), types.contains(torsion.typeC),
                  torsion.periodicity > 0, torsion.phaseRadians.isFinite, torsion.barrierKJPerMol.isFinite else {
                throw VivoArtifactValidationError.invalid("force-field torsion parameter is invalid")
            }
        }
    }

    public func fingerprint() throws -> VivoFingerprint {
        try validate()
        return try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(self))
    }
}

public struct VivoForceFieldAssignment: Codable, Sendable, Equatable {
    public var atomTypeByAtom: [String]
    public var chargeByAtomE: [Double]
    public var unresolvedAtoms: [UInt32]
    public init(atomTypeByAtom: [String], chargeByAtomE: [Double], unresolvedAtoms: [UInt32]) {
        self.atomTypeByAtom = atomTypeByAtom; self.chargeByAtomE = chargeByAtomE; self.unresolvedAtoms = unresolvedAtoms
    }
}

public enum VivoResidueTemplateAssigner {
    /// Deterministic assignment for residues represented by an imported force-field
    /// library. Small-molecule atom typing is deliberately a separate subsystem.
    public static func assign(structure: VivoMolecularStructure,
                              library: VivoForceFieldLibrary) throws -> VivoForceFieldAssignment {
        _ = try VivoStructureValidator.validate(structure)
        try library.validate()
        var templates: [String: VivoForceFieldResidueTemplate] = [:]
        for template in library.residueTemplates {
            for alias in template.residueNames { templates[alias.uppercased()] = template }
        }
        var types = [String](repeating: "", count: structure.atoms.count)
        var charges = [Double](repeating: 0, count: structure.atoms.count)
        var unresolved: [UInt32] = []
        for residue in structure.residues {
            guard let template = templates[residue.name.uppercased()] else {
                unresolved.append(contentsOf: residue.atomIndices); continue
            }
            let byName = Dictionary(uniqueKeysWithValues: template.atoms.map { ($0.name, $0) })
            for atomIndex in residue.atomIndices {
                let atom = structure.atoms[Int(atomIndex)]
                guard let parameter = byName[atom.name] else { unresolved.append(atomIndex); continue }
                types[Int(atomIndex)] = parameter.typeIdentifier
                charges[Int(atomIndex)] = parameter.chargeE
            }
        }
        for atom in structure.atoms where atom.residueIndex == nil { unresolved.append(atom.index) }
        return .init(atomTypeByAtom: types, chargeByAtomE: charges, unresolvedAtoms: Array(Set(unresolved)).sorted())
    }
}
