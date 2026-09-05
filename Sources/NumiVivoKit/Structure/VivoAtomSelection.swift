import Foundation

public indirect enum VivoAtomSelection: Codable, Sendable, Equatable {
    case all
    case none
    case atomIndices([UInt32])
    case elements([UInt16])
    case atomNames([String])
    case residueNames([String])
    case chainIdentifiers([String])
    case residueSequenceRange(Int32, Int32)
    case hetero(Bool)
    case hydrogen(Bool)
    case bonded(to: VivoAtomSelection)
    case within(distanceNM: Double, of: VivoAtomSelection, includeSource: Bool)
    case and([VivoAtomSelection])
    case or([VivoAtomSelection])
    case not(VivoAtomSelection)
}

public struct VivoSelectionResult: Codable, Sendable, Equatable {
    public var atomIndices: [UInt32]
    public init(atomIndices: [UInt32]) { self.atomIndices = atomIndices }
}

public enum VivoSelectionEvaluator {
    public static func evaluate(_ selection: VivoAtomSelection,
                                in structure: VivoMolecularStructure,
                                conformerIndex: Int = 0,
                                periodic: Bool = true) throws -> VivoSelectionResult {
        _ = try VivoStructureValidator.validate(structure)
        let index = try VivoStructureIndex(validated: structure)
        let result = try evaluate(selection, structure: structure, index: index,
                                  conformerIndex: conformerIndex, periodic: periodic, depth: 0)
        return .init(atomIndices: result.sorted())
    }

    private static func evaluate(_ selection: VivoAtomSelection,
                                 structure: VivoMolecularStructure,
                                 index: VivoStructureIndex,
                                 conformerIndex: Int,
                                 periodic: Bool,
                                 depth: Int) throws -> Set<UInt32> {
        guard depth < 64 else { throw VivoArtifactValidationError.invalid("atom selection nesting exceeds 64 levels") }
        let universe = Set(structure.atoms.map(\.index))
        switch selection {
        case .all: return universe
        case .none: return []
        case .atomIndices(let values):
            guard values.allSatisfy({ Int($0) < structure.atoms.count }) else {
                throw VivoArtifactValidationError.invalid("atom selection contains an out-of-range atom index")
            }
            return Set(values)
        case .elements(let values):
            let set = Set(values)
            return Set(structure.atoms.lazy.filter { set.contains($0.element.atomicNumber) }.map(\.index))
        case .atomNames(let values):
            let set = Set(values)
            return Set(structure.atoms.lazy.filter { set.contains($0.name) }.map(\.index))
        case .residueNames(let values):
            let set = Set(values)
            return Set(structure.residues.lazy.filter { set.contains($0.name) }.flatMap(\.atomIndices))
        case .chainIdentifiers(let values):
            let set = Set(values)
            return Set(structure.chains.lazy.filter { set.contains($0.identifier) }
                .flatMap { $0.residueIndices }.flatMap { structure.residues[Int($0)].atomIndices })
        case .residueSequenceRange(let lower, let upper):
            guard lower <= upper else { throw VivoArtifactValidationError.invalid("residue sequence range is reversed") }
            return Set(structure.residues.lazy.filter {
                guard let value = $0.sequenceNumber else { return false }
                return value >= lower && value <= upper
            }.flatMap(\.atomIndices))
        case .hetero(let desired):
            return Set(structure.atoms.lazy.filter { $0.isHetero == desired }.map(\.index))
        case .hydrogen(let desired):
            return Set(structure.atoms.lazy.filter { ($0.element.atomicNumber == 1) == desired }.map(\.index))
        case .bonded(let inner):
            let source = try evaluate(inner, structure: structure, index: index,
                                      conformerIndex: conformerIndex, periodic: periodic, depth: depth + 1)
            var result = Set<UInt32>()
            for atom in source { result.formUnion(index.bondedAtoms[Int(atom)]) }
            return result
        case .within(let distanceNM, let inner, let includeSource):
            guard distanceNM.isFinite, distanceNM >= 0 else {
                throw VivoArtifactValidationError.invalid("within selection requires a finite nonnegative distance")
            }
            guard structure.conformers.indices.contains(conformerIndex) else {
                throw VivoArtifactValidationError.unresolved("within selection requires conformer \(conformerIndex)")
            }
            let source = try evaluate(inner, structure: structure, index: index,
                                      conformerIndex: conformerIndex, periodic: periodic, depth: depth + 1)
            let coordinates = structure.conformers[conformerIndex].positionsNM
            let cutoff2 = distanceNM * distanceNM
            var result = Set<UInt32>()
            for atom in structure.atoms {
                if source.contains(atom.index) {
                    if includeSource { result.insert(atom.index) }
                    continue
                }
                let p = coordinates[Int(atom.index)]
                for origin in source {
                    var delta = p - coordinates[Int(origin)]
                    if periodic, let cell = structure.periodicCell { delta = try cell.minimumImage(delta) }
                    if delta.squaredNorm <= cutoff2 { result.insert(atom.index); break }
                }
            }
            return result
        case .and(let clauses):
            guard let first = clauses.first else { return universe }
            var result = try evaluate(first, structure: structure, index: index,
                                      conformerIndex: conformerIndex, periodic: periodic, depth: depth + 1)
            for clause in clauses.dropFirst() {
                result.formIntersection(try evaluate(clause, structure: structure, index: index,
                                                     conformerIndex: conformerIndex, periodic: periodic, depth: depth + 1))
            }
            return result
        case .or(let clauses):
            var result = Set<UInt32>()
            for clause in clauses {
                result.formUnion(try evaluate(clause, structure: structure, index: index,
                                              conformerIndex: conformerIndex, periodic: periodic, depth: depth + 1))
            }
            return result
        case .not(let clause):
            return universe.subtracting(try evaluate(clause, structure: structure, index: index,
                                                      conformerIndex: conformerIndex, periodic: periodic, depth: depth + 1))
        }
    }
}
