import Foundation
import NumiVivoKit

/// Repeated calls to the retained production owner, solely for stack sampling.
/// Each output must equal the prior full-cohort receipt; timings are not promoted.
@main struct CountProfile {
    static func main() throws {
        let args = CommandLine.arguments
        guard args.count == 5, let backend = VivoCountStoreNormalizationBackend(rawValue: args[2]) else {
            throw VivoOmicsError.invalid("profile <store> <backend> <reference-normalized> <new-output-root>")
        }
        let store = URL(fileURLWithPath: args[1]), reference = URL(fileURLWithPath: args[3]), root = URL(fileURLWithPath: args[4])
        guard !FileManager.default.fileExists(atPath: root.path) else { throw VivoOmicsError.invalid("profile output already exists") }
        try FileManager.default.createDirectory(at: root, withIntermediateDirectories: false, attributes: [.posixPermissions: 0o700])
        let expected = try VivoCanonicalJSON.decode(VivoCountStoreNormalizationReceipt.self,
            from: Data(contentsOf: reference.appendingPathComponent("receipt.json")))
        let implementation = try VivoFingerprint(bytes: Array(repeating: 0, count: 32))
        var results: [VivoCountStoreNormalizationReceipt] = []
        for i in 0..<8 {
            let result = try VivoH5ADCountStore.normalize(store, target: 10_000, backend: backend,
                implementation: implementation, to: root.appendingPathComponent("run-\(i)"))
            guard result == expected else { throw VivoOmicsError.invalid("profile output differs from qualified receipt") }
            results.append(result)
        }
        try VivoCanonicalJSON.encode(results).write(to: root.appendingPathComponent("checks.json"), options: .withoutOverwriting)
        print("passed: all 8 full-cohort receipts exact")
    }
}
