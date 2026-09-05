import Foundation

public struct VivoMolecularStructureDocument: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/molecular-structure/v1"
    public var schema: String
    public var sourceFingerprint: VivoFingerprint?
    public var sourceFormat: VivoStructureFormat?
    public var structureFingerprint: VivoFingerprint
    public var structure: VivoMolecularStructure

    public init(structure: VivoMolecularStructure,
                sourceFingerprint: VivoFingerprint? = nil,
                sourceFormat: VivoStructureFormat? = nil) throws {
        self.schema = Self.schema
        self.sourceFingerprint = sourceFingerprint
        self.sourceFormat = sourceFormat
        self.structure = structure
        self.structureFingerprint = try VivoStructureCodec.fingerprint(structure)
        try validate()
    }

    public func validate() throws {
        guard schema == Self.schema else {
            throw VivoArtifactValidationError.incompatible("unsupported molecular structure document schema \(schema)")
        }
        _ = try VivoStructureValidator.validate(structure)
        guard try VivoStructureCodec.fingerprint(structure) == structureFingerprint else {
            throw VivoArtifactValidationError.invalid("molecular structure fingerprint does not match its canonical payload")
        }
    }

    public func canonicalData() throws -> Data {
        try validate()
        return try VivoCanonicalJSON.encode(self)
    }

    public func fingerprint() throws -> VivoFingerprint {
        try VivoCanonicalJSON.fingerprint(canonicalData())
    }

    public static func decode(_ data: Data) throws -> VivoMolecularStructureDocument {
        if let document = try? VivoCanonicalJSON.decode(Self.self, from: data) {
            try document.validate()
            return document
        }
        let structure = try VivoCanonicalJSON.decode(VivoMolecularStructure.self, from: data)
        return try .init(structure: structure)
    }
}
