import Foundation

public struct VivoAtomSemanticKey: Codable, Sendable, Equatable, Hashable, CustomStringConvertible {
    public var chainIdentifier: String?
    public var residueSequenceNumber: Int32?
    public var residueInsertionCode: String?
    public var residueName: String?
    public var atomName: String
    public var elementAtomicNumber: UInt16
    public var alternateLocation: String?

    public var description: String {
        [chainIdentifier, residueSequenceNumber.map(String.init), residueInsertionCode,
         residueName, atomName, String(elementAtomicNumber), alternateLocation]
            .compactMap { $0 }.joined(separator: ":")
    }
}

public struct VivoAtomMapping: Codable, Sendable, Equatable {
    public var sourceAtomCount: UInt32
    public var targetAtomCount: UInt32
    /// Dense source index -> target index. Missing mappings are never represented
    /// implicitly; callers that need partial mapping must request it explicitly.
    public var sourceToTarget: [UInt32]

    public init(sourceAtomCount: UInt32, targetAtomCount: UInt32, sourceToTarget: [UInt32]) {
        self.sourceAtomCount = sourceAtomCount
        self.targetAtomCount = targetAtomCount
        self.sourceToTarget = sourceToTarget
    }

    public func validate(requireBijection: Bool = true) throws {
        guard sourceToTarget.count == Int(sourceAtomCount),
              sourceToTarget.allSatisfy({ $0 < targetAtomCount }) else {
            throw VivoArtifactValidationError.invalid("atom mapping shape or target index is invalid")
        }
        if requireBijection {
            guard sourceAtomCount == targetAtomCount,
                  Set(sourceToTarget).count == sourceToTarget.count else {
                throw VivoArtifactValidationError.invalid("bijective atom mapping contains duplicates or unequal atom counts")
            }
        }
    }

    public func remap<T>(_ sourceOrdered: [T]) throws -> [T] {
        try validate(requireBijection: true)
        guard sourceOrdered.count == Int(sourceAtomCount) else {
            throw VivoArtifactValidationError.invalid("source array does not match atom mapping")
        }
        var result = sourceOrdered
        for source in sourceOrdered.indices {
            result[Int(sourceToTarget[source])] = sourceOrdered[source]
        }
        return result
    }
}

public enum VivoAtomMapper {
    /// Builds an exact semantic mapping from chain/residue/atom identity. Ambiguous
    /// keys are a hard error: QM/MM and force-field state must never be assigned by
    /// whichever duplicate happens to occur first in a file.
    public static func semanticBijection(from source: VivoMolecularStructure,
                                         to target: VivoMolecularStructure) throws -> VivoAtomMapping {
        _ = try VivoStructureValidator.validate(source)
        _ = try VivoStructureValidator.validate(target)
        guard source.atoms.count == target.atoms.count else {
            throw VivoArtifactValidationError.incompatible("semantic atom mapping requires equal atom counts")
        }
        let sourceKeys = try semanticKeys(source)
        let targetKeys = try semanticKeys(target)
        let targetIndex = try uniqueIndex(targetKeys, label: "target")
        _ = try uniqueIndex(sourceKeys, label: "source")

        var mapping: [UInt32] = []
        mapping.reserveCapacity(sourceKeys.count)
        for key in sourceKeys {
            guard let index = targetIndex[key] else {
                throw VivoArtifactValidationError.unresolved("target is missing semantic atom key \(key)")
            }
            mapping.append(index)
        }
        let result = VivoAtomMapping(sourceAtomCount: UInt32(source.atoms.count),
                                     targetAtomCount: UInt32(target.atoms.count),
                                     sourceToTarget: mapping)
        try result.validate(requireBijection: true)
        return result
    }

    /// Source serials are useful for reconnecting frames from one file lineage,
    /// but are not considered chemistry identity. This helper is therefore explicit.
    public static func sourceSerialBijection(from source: VivoMolecularStructure,
                                             to target: VivoMolecularStructure) throws -> VivoAtomMapping {
        _ = try VivoStructureValidator.validate(source)
        _ = try VivoStructureValidator.validate(target)
        guard source.atoms.count == target.atoms.count else {
            throw VivoArtifactValidationError.incompatible("serial atom mapping requires equal atom counts")
        }
        func serialIndex(_ structure: VivoMolecularStructure, label: String) throws -> [Int32: UInt32] {
            var result: [Int32: UInt32] = [:]
            for atom in structure.atoms {
                guard let serial = atom.sourceSerial else {
                    throw VivoArtifactValidationError.unresolved("\(label) atom \(atom.index) has no source serial")
                }
                guard result.updateValue(atom.index, forKey: serial) == nil else {
                    throw VivoArtifactValidationError.invalid("\(label) contains duplicate source serial \(serial)")
                }
            }
            return result
        }
        let targetSerials = try serialIndex(target, label: "target")
        _ = try serialIndex(source, label: "source")
        let mapping = try source.atoms.map { atom -> UInt32 in
            guard let serial = atom.sourceSerial, let targetAtom = targetSerials[serial] else {
                throw VivoArtifactValidationError.unresolved("source serial has no target atom")
            }
            return targetAtom
        }
        let result = VivoAtomMapping(sourceAtomCount: UInt32(source.atoms.count),
                                     targetAtomCount: UInt32(target.atoms.count),
                                     sourceToTarget: mapping)
        try result.validate(requireBijection: true)
        return result
    }

    public static func semanticKeys(_ structure: VivoMolecularStructure) throws -> [VivoAtomSemanticKey] {
        _ = try VivoStructureValidator.validate(structure)
        return structure.atoms.map { atom in
            let residue = atom.residueIndex.map { structure.residues[Int($0)] }
            let chain = residue?.chainIndex.map { structure.chains[Int($0)] }
            return VivoAtomSemanticKey(chainIdentifier: chain?.identifier,
                                       residueSequenceNumber: residue?.sequenceNumber,
                                       residueInsertionCode: residue?.insertionCode,
                                       residueName: residue?.name,
                                       atomName: atom.name,
                                       elementAtomicNumber: atom.element.atomicNumber,
                                       alternateLocation: atom.alternateLocation)
        }
    }

    private static func uniqueIndex(_ keys: [VivoAtomSemanticKey], label: String) throws -> [VivoAtomSemanticKey: UInt32] {
        var result: [VivoAtomSemanticKey: UInt32] = [:]
        for (index, key) in keys.enumerated() {
            guard result.updateValue(UInt32(index), forKey: key) == nil else {
                throw VivoArtifactValidationError.invalid("\(label) structure has ambiguous semantic atom key \(key)")
            }
        }
        return result
    }
}
