import Foundation

public enum VivoStructureFormat: String, Codable, Sendable, CaseIterable {
    case vivoJSON
    case pdb
    case mmcif
    case sdf
    case mol2
    case smiles

    public static func infer(path: String) -> VivoStructureFormat? {
        let lower = path.lowercased()
        if lower.hasSuffix(".pdb") || lower.hasSuffix(".ent") { return .pdb }
        if lower.hasSuffix(".cif") || lower.hasSuffix(".mmcif") { return .mmcif }
        if lower.hasSuffix(".sdf") || lower.hasSuffix(".mol") { return .sdf }
        if lower.hasSuffix(".mol2") { return .mol2 }
        if lower.hasSuffix(".smi") || lower.hasSuffix(".smiles") { return .smiles }
        if lower.hasSuffix(".json") || lower.hasSuffix(".vivo-structure") { return .vivoJSON }
        return nil
    }
}

public struct VivoStructureImportResult: Sendable {
    public let structure: VivoMolecularStructure
    public let sourceFingerprint: VivoFingerprint
    public let structureFingerprint: VivoFingerprint
    public let format: VivoStructureFormat
}

public enum VivoStructureCodec {
    public static func decode(_ data: Data, format: VivoStructureFormat,
                              identifier: String = "molecule") throws -> VivoStructureImportResult {
        let sourceFingerprint = try VivoCanonicalJSON.fingerprint(data)
        let structure: VivoMolecularStructure
        switch format {
        case .vivoJSON:
            structure = try VivoCanonicalJSON.decode(VivoMolecularStructure.self, from: data)
        case .pdb:
            structure = try VivoPDB.read(try utf8(data, label: "PDB"), options: .init(identifier: identifier))
        case .mmcif:
            structure = try VivoMMCIF.read(try utf8(data, label: "mmCIF"), identifier: identifier)
        case .sdf:
            structure = try VivoSDF.readFirst(try utf8(data, label: "SDF"), identifier: identifier)
        case .mol2:
            structure = try VivoMOL2.read(try utf8(data, label: "MOL2"), identifier: identifier)
        case .smiles:
            let text = try utf8(data, label: "SMILES")
            guard let first = text.split(whereSeparator: \.isNewline).first,
                  let token = first.split(whereSeparator: \.isWhitespace).first else {
                throw VivoArtifactValidationError.invalid("SMILES input is empty")
            }
            structure = try VivoSMILES.read(String(token), identifier: identifier)
        }
        _ = try VivoStructureValidator.validate(structure)
        let fingerprint = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(structure))
        return .init(structure: structure, sourceFingerprint: sourceFingerprint,
                     structureFingerprint: fingerprint, format: format)
    }

    public static func encode(_ structure: VivoMolecularStructure,
                              format: VivoStructureFormat,
                              conformerIndex: Int = 0) throws -> Data {
        _ = try VivoStructureValidator.validate(structure)
        switch format {
        case .vivoJSON:
            return try VivoCanonicalJSON.encode(structure)
        case .pdb:
            return Data(try VivoPDB.write(structure, conformerIndex: conformerIndex).utf8)
        case .sdf:
            return Data(try VivoSDF.write(structure, conformerIndex: conformerIndex).utf8)
        case .mol2:
            return Data(try VivoMOL2.write(structure, conformerIndex: conformerIndex).utf8)
        case .mmcif:
            throw VivoArtifactValidationError.incompatible("native mmCIF writing is not implemented yet; lossy PDB fallback is never automatic")
        case .smiles:
            throw VivoArtifactValidationError.incompatible("canonical SMILES generation requires authoritative aromaticity/stereo perception and is not implemented yet")
        }
    }

    public static func fingerprint(_ structure: VivoMolecularStructure) throws -> VivoFingerprint {
        _ = try VivoStructureValidator.validate(structure)
        return try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(structure))
    }

    private static func utf8(_ data: Data, label: String) throws -> String {
        guard let text = String(data: data, encoding: .utf8) else {
            throw VivoArtifactValidationError.invalid("\(label) input is not UTF-8")
        }
        return text
    }
}
