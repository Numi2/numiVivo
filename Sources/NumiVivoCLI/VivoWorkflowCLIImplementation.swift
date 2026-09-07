import Foundation
import CryptoKit
import NumiVivoKit

/// One unchanged platform-executable identity schema for workflow and sampling
/// export receipts. Prepared-molecular workflow identities remain independent.
enum VivoWorkflowCLIImplementation {
    private struct Implementation: Codable { let schema: String; let executableSHA256: String; let operatingSystem: String; let executionSemantics: String }
    static func fingerprint() throws -> VivoFingerprint {
        let invocation = CommandLine.arguments[0]
        let executable: URL
        if invocation.contains("/") { executable = URL(fileURLWithPath: invocation).standardizedFileURL.resolvingSymlinksInPath() }
        else {
            guard let path = (ProcessInfo.processInfo.environment["PATH"] ?? "").split(separator: ":").map(String.init)
                .map({ URL(fileURLWithPath: $0.isEmpty ? "." : $0).appendingPathComponent(invocation) })
                .first(where: { FileManager.default.isExecutableFile(atPath: $0.path) }) else {
                throw VivoChemistryError.invalid("cannot resolve executing workflow binary")
            }
            executable = path.standardizedFileURL.resolvingSymlinksInPath()
        }
        let handle = try FileHandle(forReadingFrom: executable); defer { try? handle.close() }
        var digest = SHA256()
        while let bytes = try handle.read(upToCount: 4*1024*1024), !bytes.isEmpty { digest.update(data: bytes) }
        let hash = digest.finalize().map { String(format: "%02x", $0) }.joined()
        return try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(Implementation(schema: "numivivo.org/platform-executable/v1",
            executableSHA256: hash, operatingSystem: ProcessInfo.processInfo.operatingSystemVersionString,
            executionSemantics: "static-native-operations;fp64-electronics;metal-transactional-md;single-owner-artifact-scheduler")))
    }
}
