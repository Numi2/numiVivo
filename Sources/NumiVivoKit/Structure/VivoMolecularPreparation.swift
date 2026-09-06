import Foundation

/// References survive hydrogen removal and dense atom reindexing. Generated
/// hydrogens have explicit identities rather than pretending to be source atoms.
public enum VivoPreparationAtomReference: Codable, Sendable, Equatable, Hashable {
    case source(UInt32)
    case hydrogen(String)
}

public struct VivoHydrogenPlacement: Codable, Sendable, Equatable {
    public var identifier: String
    public var name: String
    public var parent: UInt32
    public var angleReference: UInt32
    public var planeReference: UInt32
    public var lengthNM: Double
    public var angleRadians: Double
    public var dihedralRadians: Double
    public init(identifier: String, name: String, parent: UInt32,
                angleReference: UInt32, planeReference: UInt32, lengthNM: Double,
                angleRadians: Double, dihedralRadians: Double) {
        self.identifier = identifier; self.name = name; self.parent = parent
        self.angleReference = angleReference; self.planeReference = planeReference
        self.lengthNM = lengthNM; self.angleRadians = angleRadians
        self.dihedralRadians = dihedralRadians
    }
}

/// Ordered-neighbor orientation is deliberately NOT labeled R/S: no incomplete
/// CIP implementation is used to invent an absolute stereochemical assignment.
public enum VivoPreparedStereochemistry: Codable, Sendable, Equatable {
    case tetrahedral(center: VivoPreparationAtomReference,
                     orderedNeighbors: [VivoPreparationAtomReference], positive: Bool)
    case doubleBond(a: VivoPreparationAtomReference, b: VivoPreparationAtomReference,
                    substituentA: VivoPreparationAtomReference,
                    substituentB: VivoPreparationAtomReference, sameSide: Bool)
}

public struct VivoPreparationChargeEdit: Codable, Sendable, Equatable {
    public var sourceAtom: UInt32
    public var formalCharge: Int16
    public init(sourceAtom: UInt32, formalCharge: Int16) {
        self.sourceAtom = sourceAtom; self.formalCharge = formalCharge
    }
}
public struct VivoPreparationResidueEdit: Codable, Sendable, Equatable {
    public var sourceResidue: UInt32
    public var name: String
    public init(sourceResidue: UInt32, name: String) {
        self.sourceResidue = sourceResidue; self.name = name
    }
}

public struct VivoMolecularPreparationRequest: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/molecular-preparation/v1"
    public var schema: String
    public var structure: VivoMolecularStructure
    public var microstateIdentifier: String
    public var protonationSourceIdentifier: String
    public var pH: Double
    public var expectedFormalCharge: Int
    public var removeHydrogens: [UInt32]
    public var addHydrogens: [VivoHydrogenPlacement]
    public var chargeEdits: [VivoPreparationChargeEdit]
    public var residueEdits: [VivoPreparationResidueEdit]
    public var stereochemistry: [VivoPreparedStereochemistry]
    public var forceField: VivoForceFieldLibrary?
    /// Additional impropers use prepared atom indices after hydrogen edits.
    public var compilationOptions: VivoForceFieldCompilationOptions?
    public init(structure: VivoMolecularStructure, microstateIdentifier: String,
                protonationSourceIdentifier: String, pH: Double, expectedFormalCharge: Int,
                removeHydrogens: [UInt32] = [], addHydrogens: [VivoHydrogenPlacement] = [],
                chargeEdits: [VivoPreparationChargeEdit] = [],
                residueEdits: [VivoPreparationResidueEdit] = [],
                stereochemistry: [VivoPreparedStereochemistry] = [],
                forceField: VivoForceFieldLibrary? = nil,
                compilationOptions: VivoForceFieldCompilationOptions? = nil) {
        schema = Self.schema; self.structure = structure
        self.microstateIdentifier = microstateIdentifier
        self.protonationSourceIdentifier = protonationSourceIdentifier; self.pH = pH
        self.expectedFormalCharge = expectedFormalCharge
        self.removeHydrogens = removeHydrogens; self.addHydrogens = addHydrogens
        self.chargeEdits = chargeEdits; self.residueEdits = residueEdits
        self.stereochemistry = stereochemistry; self.forceField = forceField
        self.compilationOptions = compilationOptions
    }
}

public struct VivoMolecularPreparationResult: Codable, Sendable, Equatable {
    public let requestFingerprint: VivoFingerprint
    public let sourceFingerprint: VivoFingerprint
    public let structure: VivoMolecularStructure
    public let sourceToPrepared: [UInt32?]
    public let preparedToSource: [UInt32?]
    public let hydrogenIndices: [String: UInt32]
    public let compiledForceField: VivoCompiledForceField?
    public let requiresCoordinateRelaxation: Bool
    public let interpretation: String
}

public enum VivoMolecularPreparation {
    public static let interpretation = "Explicit protonation microstate and mapped stereochemical constraints; pH records context, not an inferred pKa or population. Hydrogen coordinates are internal-coordinate initial guesses, not minimized conformers. Native parameters come exclusively from the supplied library."

    /// Pure transaction: validation failure never mutates the input structure.
    /// Heavy atoms retain their source identities and all conformers are edited.
    public static func prepare(_ request: VivoMolecularPreparationRequest) throws -> VivoMolecularPreparationResult {
        let source = request.structure
        _ = try VivoStructureValidator.validate(source)
        guard request.schema == VivoMolecularPreparationRequest.schema,
              request.pH.isFinite, !request.microstateIdentifier.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty,
              !request.protonationSourceIdentifier.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty,
              !source.conformers.isEmpty,
              request.addHydrogens.count <= 1_000_000,
              source.atoms.count + request.addHydrogens.count <= Int(UInt32.max) else {
            throw VivoArtifactValidationError.invalid("preparation schema, microstate, pH, conformers, alternate locations or capacity")
        }
        // The existing resolver preserves the selected alternate-location label.
        // Reject duplicate sites, not a resolved A/B label on a unique site.
        struct Site: Hashable { let residue: UInt32?; let name: String }
        var sites = Set<Site>()
        for atom in source.atoms {
            guard sites.insert(.init(residue: atom.residueIndex, name: atom.name)).inserted else {
                throw VivoArtifactValidationError.invalid("resolve duplicate atom sites before chemical preparation")
            }
        }
        let removed = Set(request.removeHydrogens)
        guard removed.count == request.removeHydrogens.count else {
            throw VivoArtifactValidationError.invalid("duplicate removed hydrogen")
        }
        let index = try VivoStructureIndex(validated: source)
        for atom in removed {
            guard Int(atom) < source.atoms.count, source.atoms[Int(atom)].element.atomicNumber == 1,
                  index.bondedAtoms[Int(atom)].count == 1,
                  source.atoms[Int(index.bondedAtoms[Int(atom)][0])].element.atomicNumber != 1,
                  !removed.contains(index.bondedAtoms[Int(atom)][0]) else {
                throw VivoArtifactValidationError.invalid("protonation may remove only singly attached source hydrogens")
            }
        }
        let rewrite = try VivoStructureSlicer.slice(source, atomIndices: source.atoms.map(\.index).filter { !removed.contains($0) })
        var output = rewrite.structure
        var reverse: [UInt32?] = rewrite.newToOld.map { Optional($0) }
        var hydrogenIndices: [String: UInt32] = [:]
        func retained(_ atom: UInt32) throws -> UInt32 {
            guard Int(atom) < rewrite.oldToNew.count, let target = rewrite.oldToNew[Int(atom)] else {
                throw VivoArtifactValidationError.unresolved("preparation references absent source atom \(atom)")
            }
            return target
        }
        var editedCharges = Set<UInt32>()
        for edit in request.chargeEdits {
            guard editedCharges.insert(edit.sourceAtom).inserted else {
                throw VivoArtifactValidationError.invalid("duplicate formal-charge edit")
            }
            let target = try retained(edit.sourceAtom)
            output.atoms[Int(target)].formalCharge = edit.formalCharge
        }
        var editedResidues = Set<UInt32>()
        for edit in request.residueEdits {
            guard Int(edit.sourceResidue) < source.residues.count,
                  editedResidues.insert(edit.sourceResidue).inserted, !edit.name.isEmpty,
                  let surviving = source.residues[Int(edit.sourceResidue)].atomIndices.first(where: { !removed.contains($0) }) else {
                throw VivoArtifactValidationError.invalid("duplicate, missing or empty residue-state edit")
            }
            let target = try retained(surviving)
            guard let residue = output.atoms[Int(target)].residueIndex else {
                throw VivoArtifactValidationError.invalid("residue-state mapping disappeared")
            }
            output.residues[Int(residue)].name = edit.name
        }
        for hydrogen in request.addHydrogens {
            guard !hydrogen.identifier.isEmpty, !hydrogen.name.isEmpty,
                  hydrogenIndices[hydrogen.identifier] == nil,
                  Set([hydrogen.parent, hydrogen.angleReference, hydrogen.planeReference]).count == 3,
                  hydrogen.lengthNM.isFinite, hydrogen.lengthNM >= 0.05, hydrogen.lengthNM <= 0.2,
                  hydrogen.angleRadians.isFinite, hydrogen.angleRadians > 0, hydrogen.angleRadians < .pi,
                  hydrogen.dihedralRadians.isFinite else {
                throw VivoArtifactValidationError.invalid("hydrogen identity or internal coordinates")
            }
            let parent = try retained(hydrogen.parent)
            let angle = try retained(hydrogen.angleReference), plane = try retained(hydrogen.planeReference)
            guard output.atoms[Int(parent)].element.atomicNumber != 1 else {
                throw VivoArtifactValidationError.invalid("hydrogen parent must be a heavy atom")
            }
            let residue = output.atoms[Int(parent)].residueIndex
            guard !output.atoms.contains(where: { $0.residueIndex == residue && $0.name == hydrogen.name }) else {
                throw VivoArtifactValidationError.invalid("hydrogen name collides with an existing atom")
            }
            let target = UInt32(output.atoms.count)
            for conformer in output.conformers.indices {
                let coordinates = output.conformers[conformer].positionsNM
                func displacement(_ atom: UInt32) throws -> VivoVector3D {
                    let d = coordinates[Int(atom)] - coordinates[Int(parent)]
                    return try output.periodicCell?.minimumImage(d) ?? d
                }
                let u0 = try displacement(angle), w = try displacement(plane)
                guard u0.norm > 1e-10 else { throw VivoArtifactValidationError.invalid("zero hydrogen angle reference") }
                let u = u0 / u0.norm, n0 = u.cross(w)
                guard n0.norm > 1e-10 else { throw VivoArtifactValidationError.invalid("collinear hydrogen plane references") }
                let n = n0 / n0.norm, v = n.cross(u)
                let transverse = v * cos(hydrogen.dihedralRadians) + n * sin(hydrogen.dihedralRadians)
                let direction = u * cos(hydrogen.angleRadians) + transverse * sin(hydrogen.angleRadians)
                let position = coordinates[Int(parent)] + direction * hydrogen.lengthNM
                guard position.isFinite else { throw VivoArtifactValidationError.invalid("hydrogen coordinate overflow") }
                output.conformers[conformer].positionsNM.append(position)
            }
            output.atoms.append(.init(index: target, name: hydrogen.name,
                                      element: .init(atomicNumber: 1, symbol: "H"),
                                      residueIndex: residue, isHetero: output.atoms[Int(parent)].isHetero))
            if let residue { output.residues[Int(residue)].atomIndices.append(target) }
            output.bonds.append(.init(atomA: parent, atomB: target))
            reverse.append(nil); hydrogenIndices[hydrogen.identifier] = target
        }
        let formalCharge = output.atoms.reduce(0) { $0 + Int($1.formalCharge) }
        guard formalCharge == request.expectedFormalCharge else {
            throw VivoArtifactValidationError.invalid("prepared formal charge differs from the declared microstate")
        }
        // Source potential energies are invalid after ANY topology/charge edit.
        if !removed.isEmpty || !request.addHydrogens.isEmpty || !request.chargeEdits.isEmpty || !request.residueEdits.isEmpty {
            for c in output.conformers.indices { output.conformers[c].potentialEnergyKJPerMol = nil }
        }
        output.metadata["preparationMicrostate"] = request.microstateIdentifier
        output.metadata["protonationSource"] = request.protonationSourceIdentifier
        output.metadata["preparationPH"] = String(request.pH)
        _ = try VivoStructureValidator.validate(output)
        func resolve(_ ref: VivoPreparationAtomReference) throws -> UInt32 {
            switch ref {
            case .source(let atom): return try retained(atom)
            case .hydrogen(let name):
                guard let atom = hydrogenIndices[name] else { throw VivoArtifactValidationError.unresolved("unknown generated hydrogen \(name)") }
                return atom
            }
        }
        let topology = try VivoStructureIndex(validated: output)
        for stereo in request.stereochemistry {
            for conformer in output.conformers {
                let p = conformer.positionsNM
                func delta(_ a: UInt32, _ b: UInt32) throws -> VivoVector3D {
                    let d = p[Int(a)] - p[Int(b)]
                    return try output.periodicCell?.minimumImage(d) ?? d
                }
                switch stereo {
                case .tetrahedral(let centerRef, let refs, let positive):
                    let center = try resolve(centerRef), neighbors = try refs.map(resolve)
                    guard neighbors.count == 4, Set(neighbors).count == 4,
                          Set(topology.bondedAtoms[Int(center)]) == Set(neighbors) else {
                        throw VivoArtifactValidationError.invalid("tetrahedral constraint requires all four bonded neighbors")
                    }
                    let vectors = try neighbors.map { try delta($0, center) }
                    let scale = vectors.map(\.norm).max() ?? 0
                    guard scale > 1e-10 else { throw VivoArtifactValidationError.invalid("degenerate tetrahedral geometry") }
                    let a = (vectors[0] - vectors[3]) / scale
                    let b = (vectors[1] - vectors[3]) / scale
                    let c = (vectors[2] - vectors[3]) / scale
                    let orientation = a.dot(b.cross(c))
                    guard orientation.isFinite, abs(orientation) > 1e-5, (orientation > 0) == positive else {
                        throw VivoArtifactValidationError.incompatible("tetrahedral stereochemistry is planar or inverted in \(conformer.identifier)")
                    }
                case .doubleBond(let ar, let br, let xr, let yr, let sameSide):
                    let a = try resolve(ar), b = try resolve(br), x = try resolve(xr), y = try resolve(yr)
                    guard Set([a,b,x,y]).count == 4,
                          output.bonds.contains(where: { (($0.atomA == a && $0.atomB == b) || ($0.atomA == b && $0.atomB == a)) && $0.order == .double }),
                          topology.bondedAtoms[Int(a)].contains(x), topology.bondedAtoms[Int(b)].contains(y) else {
                        throw VivoArtifactValidationError.invalid("double-bond stereo topology")
                    }
                    let axis0 = try delta(b, a)
                    guard axis0.norm > 1e-10 else { throw VivoArtifactValidationError.invalid("zero double bond") }
                    let axis = axis0 / axis0.norm
                    let dx = try delta(x,a), dy = try delta(y,b)
                    let px = dx - axis * dx.dot(axis), py = dy - axis * dy.dot(axis)
                    guard px.norm > 1e-10, py.norm > 1e-10 else { throw VivoArtifactValidationError.invalid("linear double-bond substituents") }
                    let cosine = px.dot(py) / (px.norm * py.norm)
                    guard cosine.isFinite, abs(cosine) > 0.5, (cosine > 0) == sameSide else {
                        throw VivoArtifactValidationError.incompatible("double-bond stereochemistry is twisted or inverted")
                    }
                }
            }
        }
        let compiled: VivoCompiledForceField?
        if let library = request.forceField {
            let assignment = try VivoResidueTemplateAssigner.assignPrepared(structure: output, library: library)
            let options = request.compilationOptions ?? .init()
            guard options.requireAllBondParameters, options.requireAllAngleParameters, options.requireAllProperTorsions else {
                throw VivoArtifactValidationError.invalid("prepared native systems cannot omit required bonded parameters")
            }
            compiled = try VivoForceFieldCompiler.compile(structure: output, library: library, assignment: assignment, options: options)
        } else {
            guard request.compilationOptions == nil else { throw VivoArtifactValidationError.invalid("compilation options require a force-field library") }
            compiled = nil
        }
        return .init(requestFingerprint: try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(request)),
                     sourceFingerprint: try VivoStructureCodec.fingerprint(source), structure: output,
                     sourceToPrepared: rewrite.oldToNew, preparedToSource: reverse,
                     hydrogenIndices: hydrogenIndices, compiledForceField: compiled,
                     requiresCoordinateRelaxation: !request.addHydrogens.isEmpty || !removed.isEmpty || !request.chargeEdits.isEmpty || !request.residueEdits.isEmpty,
                     interpretation: interpretation)
    }
}

public extension VivoResidueTemplateAssigner {
    /// Strict preparation path over the EXISTING library/assigner/compiler.
    /// Incomplete templates and missing bond orders cannot become parameters.
    static func assignPrepared(structure: VivoMolecularStructure,
                               library: VivoForceFieldLibrary,
                               chargeToleranceE: Double = 1e-4) throws -> VivoForceFieldAssignment {
        _ = try VivoStructureValidator.validate(structure); try library.validate()
        guard chargeToleranceE.isFinite, chargeToleranceE >= 0, chargeToleranceE <= 0.01 else {
            throw VivoArtifactValidationError.invalid("native assignment charge tolerance")
        }
        func edge(_ a: String, _ b: String, _ order: VivoBondOrder) -> String {
            // Length-prefixed names prevent atom-name delimiters from colliding.
            let names = [a,b].sorted()
            return "\(names[0].utf8.count):\(names[0])\(names[1].utf8.count):\(names[1]):\(order.rawValue)"
        }
        for residue in structure.residues {
            guard let template = library.residueTemplates.first(where: { $0.residueNames.contains(where: { $0.uppercased() == residue.name.uppercased() }) }) else {
                throw VivoArtifactValidationError.unresolved("no native template for prepared residue \(residue.index) \(residue.name)")
            }
            let names = residue.atomIndices.map { structure.atoms[Int($0)].name }
            guard Set(names).count == names.count, Set(names) == Set(template.atoms.map(\.name)) else {
                throw VivoArtifactValidationError.incompatible("prepared residue \(residue.index) is incomplete or differs from its protonation template")
            }
            let atoms = Set(residue.atomIndices)
            let local = structure.bonds.filter { atoms.contains($0.atomA) && atoms.contains($0.atomB) }
            let actual = Set(local.map { edge(structure.atoms[Int($0.atomA)].name, structure.atoms[Int($0.atomB)].name, $0.order) })
            let expected = Set(template.bonds.map { edge($0.atomA, $0.atomB, $0.order) })
            guard expected.count == template.bonds.count, actual == expected else {
                throw VivoArtifactValidationError.incompatible("prepared residue \(residue.index) connectivity/bond orders differ from its parameter template")
            }
        }
        let assignment = try assign(structure: structure, library: library)
        guard assignment.unresolvedAtoms.isEmpty else { throw VivoArtifactValidationError.unresolved("native assignment remains incomplete") }
        let partial = assignment.chargeByAtomE.reduce(0,+)
        let formal = structure.atoms.reduce(0.0) { $0 + Double($1.formalCharge) }
        guard partial.isFinite, abs(partial - formal) <= chargeToleranceE else {
            throw VivoArtifactValidationError.incompatible("native partial-charge sum does not match the explicit formal-charge microstate")
        }
        return assignment
    }
}
