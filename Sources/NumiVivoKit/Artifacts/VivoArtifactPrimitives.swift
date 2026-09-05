import Foundation

// Shared artifact identity and validation are independent of program compilation
// and host biology. Valid fingerprint encoding is unchanged. Invalid standalone
// fingerprint construction now reports an artifact error, not a pack-header error.
public enum VivoArtifactValidationError: Error, Sendable, CustomStringConvertible {
    case invalid(String)
    case unresolved(String)
    case incompatible(String)

    public var description: String {
        switch self {
        case .invalid(let message): "Invalid artifact: \(message)"
        case .unresolved(let message): "Unresolved artifact reference: \(message)"
        case .incompatible(let message): "Incompatible artifact: \(message)"
        }
    }
}

public struct VivoFingerprint: Hashable, Codable, Sendable, CustomStringConvertible {
    public static let byteCount = 32
    public let bytes: [UInt8]
    public init(bytes: [UInt8]) throws {
        guard bytes.count == Self.byteCount else {
            throw VivoArtifactValidationError.invalid("SHA-256 fingerprints require exactly 32 bytes")
        }
        self.bytes = bytes
    }
    private enum CodingKeys: String, CodingKey { case bytes }
    public init(from decoder: Decoder) throws {
        let container = try decoder.container(keyedBy: CodingKeys.self)
        var array = try container.nestedUnkeyedContainer(forKey: .bytes)
        if let count = array.count, count != Self.byteCount {
            throw DecodingError.dataCorruptedError(forKey: .bytes, in: container,
                                                   debugDescription: "SHA-256 fingerprints require exactly 32 bytes")
        }
        var value: [UInt8] = []
        value.reserveCapacity(Self.byteCount)
        for _ in 0..<Self.byteCount {
            guard !array.isAtEnd else {
                throw DecodingError.dataCorruptedError(forKey: .bytes, in: container, debugDescription: "Truncated fingerprint")
            }
            value.append(try array.decode(UInt8.self))
        }
        guard array.isAtEnd else {
            throw DecodingError.dataCorruptedError(forKey: .bytes, in: container, debugDescription: "Oversized fingerprint")
        }
        bytes = value
    }
    public func encode(to encoder: Encoder) throws {
        var container = encoder.container(keyedBy: CodingKeys.self)
        try container.encode(bytes, forKey: .bytes)
    }
    public var hex: String { bytes.map { String(format: "%02x", $0) }.joined() }
    public var description: String { hex }
}

