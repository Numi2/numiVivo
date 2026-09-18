import Foundation

// Harness-only JSON/hash implementation. Molecular validation and artifact identity
// use exact production sources, not substitute validators or fingerprint types.
#if canImport(CryptoKit)
import CryptoKit
#else
@_silgen_name("SHA256")
private func opensslSHA256(_ input: UnsafeRawPointer?, _ size: Int, _ digest: UnsafeMutablePointer<UInt8>) -> UnsafeMutablePointer<UInt8>?
#endif
public enum VivoCanonicalJSON {
    public static func encode<T: Encodable>(_ value: T) throws -> Data {
        let encoder = JSONEncoder()
        encoder.outputFormatting = [.sortedKeys, .withoutEscapingSlashes]
        encoder.dateEncodingStrategy = .iso8601; encoder.dataEncodingStrategy = .base64
        return try encoder.encode(value)
    }
    public static func decode<T: Decodable>(_ type: T.Type, from data: Data) throws -> T {
        let decoder = JSONDecoder()
        decoder.dateDecodingStrategy = .iso8601; decoder.dataDecodingStrategy = .base64
        return try decoder.decode(type, from: data)
    }
    public static func fingerprint(_ data: Data) throws -> VivoFingerprint {
        #if canImport(CryptoKit)
        return try .init(bytes: Array(SHA256.hash(data: data)))
        #else
        var output = [UInt8](repeating: 0, count: 32)
        let ok = data.withUnsafeBytes { raw in opensslSHA256(raw.baseAddress, raw.count, &output) != nil }
        guard ok else { throw VivoArtifactValidationError.invalid("OpenSSL SHA-256 failed") }
        return try .init(bytes: output)
        #endif
    }
}
