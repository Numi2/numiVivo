import Foundation
import NumiVivoKit

struct VivoProteinStressCLICommands {
    static func handles(_ command: String?) -> Bool {
        ["protein-stress-example", "protein-stress-validate", "protein-stress-run", "protein-stress-verify",
         "protein-stress-campaign", "protein-stress-help"].contains(command ?? "")
    }
    func run(arguments: [String]) async -> Int32 {
        do {
            guard let command = arguments.first else { throw VivoProteinStressError.invalid("missing command") }
            if command == "protein-stress-help" { print(Self.help); return 0 }
            var options: [String: String] = [:], positions: [String] = [], index = 1
            while index < arguments.count {
                let key = arguments[index]
                if key.hasPrefix("--") {
                    guard ["--output", "--store", "--resume", "--journal"].contains(key),
                          options[key] == nil, index + 1 < arguments.count, !arguments[index + 1].hasPrefix("--") else {
                        throw VivoProteinStressError.invalid("unknown, repeated, or valueless option \(key)")
                    }
                    options[key] = arguments[index + 1]; index += 2
                } else { positions.append(key); index += 1 }
            }
            func required(_ key: String) throws -> String {
                guard let value = options[key], !value.isEmpty else { throw VivoProteinStressError.invalid("missing \(key)") }
                return value
            }
            func rejectUnused(_ allowed: Set<String>) throws {
                guard Set(options.keys).isSubset(of: allowed) else { throw VivoProteinStressError.invalid("option not used by \(command)") }
            }
            func output<T: Encodable>(_ value: T) throws {
                let url = URL(fileURLWithPath: try required("--output"))
                try VivoMDAtomicFileExport.write(VivoCanonicalJSON.encode(value), to: url, overwrite: false)
            }
            if command == "protein-stress-example" {
                try rejectUnused(["--output"])
                guard positions.isEmpty else { throw VivoProteinStressError.invalid("example accepts no input file") }
                try output(VivoProteinStressExample.harmonicFixture()); return 0
            }
            guard positions.count == 1 else { throw VivoProteinStressError.invalid("exactly one request JSON path is required") }
            // Refuse an existing output before beginning expensive simulation.
            let destination = try required("--output")
            guard !FileManager.default.fileExists(atPath: destination) else { throw VivoProteinStressError.invalid("output already exists") }
            let data = try VivoMDAtomicFileExport.read(URL(fileURLWithPath: positions[0]), maximumBytes: 128 << 20)
            if command == "protein-stress-campaign" {
                try rejectUnused(["--store", "--output"])
                let request = try VivoCanonicalJSON.decode(VivoProteinStressCampaignRequest.self, from: data)
                let store = try VivoArtifactStore(rootURL: URL(fileURLWithPath: required("--store")))
                let report = try await VivoProteinStressCampaign.run(request, store: store)
                try output(report); return report.candidates.allSatisfy(\.allReplicasCompleted) ? 0 : 75
            }
            let request = try VivoCanonicalJSON.decode(VivoProteinStressRequest.self, from: data)
            if command == "protein-stress-validate" {
                try rejectUnused(["--output"])
                let compiled = try VivoProteinStressCompilation(request)
                struct Validation: Encodable { let request: VivoFingerprint; let maximumJournalEntries: Int; let status: String }
                try output(Validation(request: request.fingerprint(), maximumJournalEntries: compiled.maximumJournalEntries,
                    status: "input/preflight contract only; no GPU execution or physical qualification")); return 0
            }
            let store = try VivoArtifactStore(rootURL: URL(fileURLWithPath: required("--store")))
            if command == "protein-stress-verify" {
                try rejectUnused(["--output", "--store", "--journal"])
                let tail = try Self.fingerprint(required("--journal"))
                let entries = try await VivoProteinStressRunner.verify(request, store: store, journalTail: tail)
                struct Verification: Encodable { let journalTail: VivoFingerprint; let verifiedEntries: Int; let status: String }
                try output(Verification(journalTail: tail, verifiedEntries: entries.count,
                    status: "identity, clock, transition, and metric reconstruction verified; physical dynamics not independently replayed")); return 0
            }
            try rejectUnused(["--output", "--store", "--resume"])
            let resume = try options["--resume"].map(Self.fingerprint)
            let receipt = try await VivoProteinStressRunner.run(request, store: store, resumeFrom: resume)
            if let tail = receipt.journalTail { FileHandle.standardError.write(Data("Durable journal: \(tail.hex)\n".utf8)) }
            try output(receipt); return receipt.disposition == .completed ? 0 : 75
        } catch {
            FileHandle.standardError.write(Data("\(error)\n".utf8)); return 1
        }
    }
    private static func fingerprint(_ text: String) throws -> VivoFingerprint {
        let input = Array(text.utf8)
        guard input.count == 64, input.allSatisfy({ (48...57).contains($0) || (97...102).contains($0) }) else {
            throw VivoProteinStressError.invalid("expected lowercase SHA-256")
        }
        func nibble(_ value: UInt8) -> UInt8 { value <= 57 ? value - 48 : value - 97 + 10 }
        return try VivoFingerprint(bytes: stride(from: 0, to: 64, by: 2).map { nibble(input[$0]) * 16 + nibble(input[$0 + 1]) })
    }
    static let help = """
    protein-stress-example --output fixture.json
    protein-stress-validate request.json --output validation.json
    protein-stress-run request.json --store artifacts --output receipt.json [--resume JOURNAL_SHA256]
    protein-stress-verify request.json --store artifacts --journal JOURNAL_SHA256 --output verification.json
    protein-stress-campaign campaign.json --store artifacts --output report.json

    Exact prepared checkpoints, fixed-cell NVT, explicit physical-particle selections,
    piecewise-constant harmonic pulling/temperature stages, and bounded storage.
    Existing outputs are never overwritten. Failed/rejected replicas remain visible.
    The example is a harmonic numerical fixture, not a protein benchmark.
    """
}
