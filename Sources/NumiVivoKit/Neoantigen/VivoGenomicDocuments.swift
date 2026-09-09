import Foundation

/// Bound nesting and reject duplicate object keys before Codable can recurse or
/// different JSON implementations can select different values. Syntax and number
/// representability are then checked by the shared canonical decoder.
public enum VivoGenomicDocuments {
    public static func decode<T: Decodable>(_ type: T.Type, from data: Data) throws -> T {
        try VivoAtlasEvidence.require(data.count <= VivoAtlasEvidence.maximumDocumentBytes, "Genomic JSON exceeds 32 MiB.")
        var scanner = Scanner(bytes: Array(data))
        try scanner.value(depth: 0)
        scanner.whitespace()
        try VivoAtlasEvidence.require(scanner.index == scanner.bytes.count, "Trailing JSON content.")
        return try VivoCanonicalJSON.decode(type, from: data)
    }
    private struct Scanner {
        let bytes: [UInt8]
        var index = 0
        var tokens = 0
        mutating func whitespace() {
            while index < bytes.count && [UInt8(9), 10, 13, 32].contains(bytes[index]) { index += 1 }
        }
        mutating func consume(_ byte: UInt8) throws {
            whitespace()
            try VivoAtlasEvidence.require(index < bytes.count && bytes[index] == byte, "Malformed genomic JSON.")
            index += 1
        }
        mutating func string() throws -> String {
            whitespace(); let start = index
            try consume(34)
            var escaped = false
            while index < bytes.count {
                let byte = bytes[index]; index += 1
                try VivoAtlasEvidence.require(index - start <= 65_536, "JSON string exceeds encoded length limit.")
                if escaped { escaped = false; continue }
                if byte == 92 { escaped = true; continue }
                if byte == 34 {
                    return try JSONDecoder().decode(String.self, from: Data(bytes[start..<index]))
                }
            }
            throw VivoGenomicEvidenceError.invalid("Unterminated JSON string.")
        }
        mutating func value(depth: Int) throws {
            whitespace(); tokens += 1
            try VivoAtlasEvidence.require(depth <= 24 && tokens <= 1_000_000 && index < bytes.count, "JSON nesting, token budget or input boundary exceeded.")
            switch bytes[index] {
            case 123:
                index += 1; whitespace()
                if index < bytes.count && bytes[index] == 125 { index += 1; return }
                var keys = Set<String>()
                while true {
                    let key = try string()
                    try VivoAtlasEvidence.require(keys.insert(key).inserted, "Duplicate JSON object key.")
                    try consume(58); try value(depth: depth + 1); whitespace()
                    if index < bytes.count && bytes[index] == 125 { index += 1; return }
                    try consume(44)
                }
            case 91:
                index += 1; whitespace()
                if index < bytes.count && bytes[index] == 93 { index += 1; return }
                while true {
                    try value(depth: depth + 1); whitespace()
                    if index < bytes.count && bytes[index] == 93 { index += 1; return }
                    try consume(44)
                }
            case 34: _ = try string()
            default:
                let start = index
                while index < bytes.count && ![UInt8(9), 10, 13, 32, 44, 93, 125].contains(bytes[index]) { index += 1 }
                try VivoAtlasEvidence.require(index > start && index - start <= 1024, "Invalid or oversized JSON scalar.")
                // Full JSON syntax is checked by Codable after this bounded pass.
            }
        }
    }
}
