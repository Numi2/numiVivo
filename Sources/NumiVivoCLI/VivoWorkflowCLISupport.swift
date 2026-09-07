import Foundation
import NumiVivoKit
#if canImport(Darwin)
import Darwin
#elseif canImport(Glibc)
import Glibc
#endif

/// Shared CLI-only process termination. Domain code never exits the process.
enum VivoNativeExecution {
    static func exit(_ status: Int32) -> Never {
        #if canImport(Darwin)
        Darwin.exit(status)
        #elseif canImport(Glibc)
        Glibc.exit(status)
        #else
        fatalError("unsupported process-exit platform")
        #endif
    }
}

extension VivoFingerprint {
    /// Strict lowercase/uppercase hexadecimal parser for CLI artifact identifiers.
    init(hex: String) throws {
        guard hex.utf8.count == Self.byteCount*2 else {
            throw VivoArtifactValidationError.invalid("SHA-256 fingerprint hex requires exactly 64 characters")
        }
        let scalars = Array(hex.utf8)
        func nibble(_ byte: UInt8) -> UInt8? {
            switch byte {
            case 48...57: return byte-48
            case 65...70: return byte-65+10
            case 97...102: return byte-97+10
            default: return nil
            }
        }
        var bytes=[UInt8](); bytes.reserveCapacity(Self.byteCount)
        for i in stride(from:0,to:scalars.count,by:2) {
            guard let high=nibble(scalars[i]),let low=nibble(scalars[i+1]) else {
                throw VivoArtifactValidationError.invalid("SHA-256 fingerprint contains non-hexadecimal characters")
            }
            bytes.append((high<<4)|low)
        }
        try self.init(bytes:bytes)
    }
}

extension VivoWorkflowCLICommands {
    static func diagnostics(_ message: String) {
        try? FileHandle.standardError.write(contentsOf: Data((message+"\n").utf8))
    }
}
