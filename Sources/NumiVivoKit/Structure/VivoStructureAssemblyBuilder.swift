import Foundation

/// Shared assembly path for text structure formats. It centralizes dense indices,
/// residue/chain membership, and prevents each parser from inventing identity rules.
final class VivoStructureAssemblyBuilder {
    struct ResidueKey: Hashable {
        var chain: String
        var sequence: Int32?
        var insertion: String?
        var name: String
    }

    private(set) var atoms: [VivoMolecularAtom] = []
    private(set) var positions: [VivoVector3D] = []
    private(set) var residues: [VivoMolecularResidue] = []
    private(set) var chains: [VivoMolecularChain] = []
    private var chainMap: [String: UInt32] = [:]
    private var residueMap: [ResidueKey: UInt32] = [:]

    @discardableResult
    func appendAtom(name: String, element: VivoElement, positionNM: VivoVector3D,
                    residueName: String? = nil, chainIdentifier: String? = nil,
                    residueSequence: Int32? = nil, insertionCode: String? = nil,
                    sourceSerial: Int32? = nil, alternateLocation: String? = nil,
                    occupancy: Double? = nil, bFactor: Double? = nil,
                    formalCharge: Int16 = 0, isHetero: Bool = false) throws -> UInt32 {
        guard atoms.count < Int(UInt32.max), positionNM.isFinite else {
            throw VivoArtifactValidationError.invalid("structure atom capacity exceeded or coordinate is nonfinite")
        }
        var residueIndex: UInt32?
        if let residueName {
            let chainName = chainIdentifier ?? ""
            let chainIndex: UInt32
            if let existing = chainMap[chainName] { chainIndex = existing }
            else {
                guard chains.count < Int(UInt32.max) else { throw VivoArtifactValidationError.invalid("chain capacity exceeded") }
                chainIndex = UInt32(chains.count)
                chainMap[chainName] = chainIndex
                chains.append(.init(index: chainIndex, identifier: chainName))
            }
            let key = ResidueKey(chain: chainName, sequence: residueSequence,
                                 insertion: normalized(insertionCode), name: residueName)
            if let existing = residueMap[key] { residueIndex = existing }
            else {
                guard residues.count < Int(UInt32.max) else { throw VivoArtifactValidationError.invalid("residue capacity exceeded") }
                let value = UInt32(residues.count)
                residueMap[key] = value
                residueIndex = value
                residues.append(.init(index: value, name: residueName, chainIndex: chainIndex,
                                      sequenceNumber: residueSequence, insertionCode: normalized(insertionCode)))
                chains[Int(chainIndex)].residueIndices.append(value)
            }
        }
        let index = UInt32(atoms.count)
        atoms.append(.init(index: index, name: name, element: element,
                           formalCharge: formalCharge, residueIndex: residueIndex,
                           sourceSerial: sourceSerial, alternateLocation: normalized(alternateLocation),
                           occupancy: occupancy, bFactor: bFactor, isHetero: isHetero))
        positions.append(positionNM)
        if let residueIndex { residues[Int(residueIndex)].atomIndices.append(index) }
        return index
    }

    func makeStructure(identifier: String, bonds: [VivoMolecularBond] = [],
                       conformers: [VivoMolecularConformer]? = nil,
                       periodicCell: VivoPeriodicCell? = nil,
                       metadata: [String: String] = [:]) throws -> VivoMolecularStructure {
        let resolvedConformers = conformers ?? [.init(positionsNM: positions)]
        let structure = VivoMolecularStructure(identifier: identifier, atoms: atoms, bonds: bonds,
                                               residues: residues, chains: chains,
                                               conformers: resolvedConformers,
                                               periodicCell: periodicCell, metadata: metadata)
        _ = try VivoStructureValidator.validate(structure)
        return structure
    }

    private func normalized(_ value: String?) -> String? {
        guard let value else { return nil }
        let trimmed = value.trimmingCharacters(in: .whitespacesAndNewlines)
        return trimmed.isEmpty || trimmed == "." || trimmed == "?" ? nil : trimmed
    }
}

extension String {
    func vivoField(_ range: Range<Int>) -> String {
        guard range.lowerBound < utf8.count else { return "" }
        let bytes = Array(utf8)
        let upper = min(range.upperBound, bytes.count)
        return String(decoding: bytes[range.lowerBound..<upper], as: UTF8.self)
    }
}
