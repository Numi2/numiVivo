import Foundation

// Harness-only interfaces for testing production numerical source on Linux.
// This is NOT the production artifact store, integrity validator, or CryptoKit.
// SHA-256 delegates to system OpenSSL; JSON settings match VivoCanonicalJSON.
#if canImport(CryptoKit)
import CryptoKit
#else
@_silgen_name("SHA256")
private func opensslSHA256(_ input: UnsafeRawPointer?, _ size: Int, _ digest: UnsafeMutablePointer<UInt8>) -> UnsafeMutablePointer<UInt8>?
#endif
public struct VivoFingerprint: Codable, Hashable, Sendable {
    public let bytes: [UInt8]
    public var hex: String { bytes.map { String(format: "%02x", $0) }.joined() }
}
public enum VivoCanonicalJSON {
    public static func encode<T: Encodable>(_ value: T) throws -> Data {
        let e = JSONEncoder()
        e.outputFormatting = [.sortedKeys, .withoutEscapingSlashes]
        e.dateEncodingStrategy = .iso8601; e.dataEncodingStrategy = .base64
        return try e.encode(value)
    }
    public static func decode<T: Decodable>(_ type: T.Type, from data: Data) throws -> T {
        let d = JSONDecoder(); d.dateDecodingStrategy = .iso8601; d.dataDecodingStrategy = .base64
        return try d.decode(type, from: data)
    }
    public static func fingerprint(_ data: Data) throws -> VivoFingerprint {
        #if canImport(CryptoKit)
        return .init(bytes: Array(SHA256.hash(data: data)))
        #else
        var output = [UInt8](repeating: 0, count: 32)
        _ = data.withUnsafeBytes { raw in opensslSHA256(raw.baseAddress, raw.count, &output) }
        return .init(bytes: output)
        #endif
    }
}
