import Foundation
import NumiVivoKit

struct VivoMDProtocolCLICommands {
    private static let commands: Set<String> = [
        "md-protocol-template", "md-protocol-validate", "md-protocol-run",
        "md-protocol-resume", "md-protocol-inspect", "md-trajectory-inspect", "md-protocol-help"
    ]
    static func handles(_ command: String?) -> Bool { command.map(commands.contains) ?? false }

    func run(arguments: [String]) async -> Int32 {
        do {
            guard let command = arguments.first else { throw ProtocolCLIError.usage("missing command") }
            if command == "md-protocol-help" {
                FileHandle.standardOutput.write(Data(Self.help.utf8)); return 0
            }
            let args = try ProtocolArguments(Array(arguments.dropFirst()))
            switch command {
            case "md-protocol-template":
                try args.allow(["nvt-steps", "npt-steps", "production-steps"], positionals: 1)
                let system: VivoClassicalSystem = try load(args.positionals[0])
                let plan = VivoMDProtocolPlan.preparationTemplate(systemFingerprint: try system.fingerprint(),
                    nvtSteps: try args.positive("nvt-steps", fallback: 10_000),
                    nptSteps: try args.positive("npt-steps", fallback: 10_000),
                    productionSteps: try args.positive("production-steps", fallback: 100_000))
                try plan.validate(); try printJSON(plan); return 0
            case "md-protocol-validate":
                try args.allow(["system", "state"], positionals: 1)
                let plan: VivoMDProtocolPlan = try load(args.positionals[0])
                let system: VivoClassicalSystem = try load(args.required("system"))
                let state: VivoClassicalInitialState = try load(args.required("state"))
                let reports = try plan.validate(system: system, initialState: state)
                let report = ValidationReport(planFingerprint: try plan.fingerprint(),
                    executable: reports.allSatisfy(\.executable), stages: reports)
                try printJSON(report); return report.executable ? 0 : 65
            case "md-protocol-run", "md-protocol-resume":
                let allowed: Set<String> = command == "md-protocol-run"
                    ? ["system", "state", "store"] : ["system", "store", "checkpoint", "reference"]
                try args.allow(allowed, positionals: 1)
                let plan: VivoMDProtocolPlan = try load(args.positionals[0])
                let system: VivoClassicalSystem = try load(args.required("system"))
                let store = try VivoArtifactStore(rootURL: URL(fileURLWithPath: args.required("store")))
                let runner: VivoMDProtocolRunner
                if command == "md-protocol-run" {
                    let state: VivoClassicalInitialState = try load(args.required("state"))
                    runner = try await .start(system: system, initialState: state, plan: plan, store: store)
                } else {
                    let hash = try await checkpointHash(args, store: store)
                    runner = try await .resume(system: system, plan: plan, store: store, checkpoint: hash)
                }
                diagnostic("checkpoint reference: \(runner.checkpointReferenceName)")
                let receipt = try await runner.run()
                // Scientific outputs stay in the rooted artifact store. stdout is
                // a small receipt that can be redirected by the caller.
                let artifact = try await store.put(data: VivoCanonicalJSON.encode(receipt), kind: "md-protocol-receipt",
                                                     mediaType: "application/vnd.numivivo.md-protocol-receipt+json")
                try printJSON(ExecutionReport(receipt: receipt, receiptArtifact: artifact))
                return receipt.disposition == .completed ? 0 : 75
            case "md-protocol-inspect":
                try args.allow(["store", "checkpoint", "reference"], positionals: 0)
                let store = try VivoArtifactStore(rootURL: URL(fileURLWithPath: args.required("store")), createIfNeeded: false)
                let hash = try await checkpointHash(args, store: store)
                let descriptor = try await store.descriptor(for: hash)
                guard descriptor.kind == "md-protocol-checkpoint", descriptor.byteCount <= 2 * 1024 * 1024 else {
                    throw ProtocolCLIError.usage("reference is not a bounded MD protocol checkpoint")
                }
                let checkpoint = try VivoCanonicalJSON.decode(VivoMDProtocolCheckpoint.self, from: await store.data(for: hash))
                let planDescriptor = try await store.descriptor(for: checkpoint.planFingerprint)
                guard planDescriptor.kind == "md-protocol", planDescriptor.byteCount <= 2 * 1024 * 1024 else {
                    throw ProtocolCLIError.usage("checkpoint's plan is not a bounded MD protocol")
                }
                let plan = try VivoCanonicalJSON.decode(VivoMDProtocolPlan.self, from: await store.data(for: checkpoint.planFingerprint))
                try checkpoint.validate(plan: plan)
                try printJSON(checkpoint); return 0
            case "md-trajectory-inspect":
                try args.allow(["store", "manifest", "maximum-chunks", "verify"], positionals: 0)
                let store = try VivoArtifactStore(rootURL: URL(fileURLWithPath: args.required("store")), createIfNeeded: false)
                let hash = try fingerprint(args.required("manifest"))
                let reader = try await VivoMDTrajectoryArchiveReader.open(store: store, manifest: hash)
                let maximum = try args.positive("maximum-chunks", fallback: 100_000)
                let links = try await reader.index(maximumChunks: maximum)
                if args.verify { try await reader.verify(maximumChunks: maximum) }
                try printJSON(TrajectoryReport(manifest: reader.manifest, indexedChunks: links.count,
                                                allPayloadsVerified: args.verify))
                return 0
            default: throw ProtocolCLIError.usage("unknown MD protocol command")
            }
        } catch let error as ProtocolCLIError {
            diagnostic(error.description); return 64
        } catch let error as VivoMDRuntimeError {
            diagnostic(error.description); return 75
        } catch is CancellationError {
            diagnostic("cancelled before a protocol receipt was available"); return 75
        } catch {
            diagnostic(String(describing: error)); return 65
        }
    }

    private func checkpointHash(_ args: ProtocolArguments, store: VivoArtifactStore) async throws -> VivoFingerprint {
        guard (args.options["checkpoint"] != nil) != (args.options["reference"] != nil) else {
            throw ProtocolCLIError.usage("provide exactly one of --checkpoint or --reference")
        }
        if let raw = args.options["checkpoint"] { return try fingerprint(raw) }
        return try await store.reference(args.required("reference")).artifact.fingerprint
    }
    private func load<T: Decodable & Sendable>(_ path: String) throws -> T {
        let limits = VivoArtifactLoadLimits(maximumBytes: 256 * 1024 * 1024,
                                            maximumNodes: 8_000_000, maximumArrayElements: 2_000_000)
        return try VivoValidatedArtifactLoader.decode(T.self, at: URL(fileURLWithPath: path), limits: limits).value
    }
    private func fingerprint(_ raw: String) throws -> VivoFingerprint {
        let bytes = Array(raw.utf8)
        guard bytes.count == 64, bytes.allSatisfy({ (48...57).contains($0) || (97...102).contains($0) }) else {
            throw ProtocolCLIError.usage("fingerprint must be 64 lowercase hexadecimal characters")
        }
        func digit(_ value: UInt8) -> UInt8 { value <= 57 ? value - 48 : value - 97 + 10 }
        return try .init(bytes: stride(from: 0, to: bytes.count, by: 2).map { digit(bytes[$0]) * 16 + digit(bytes[$0 + 1]) })
    }
    private func printJSON<T: Encodable>(_ value: T) throws {
        FileHandle.standardOutput.write(try VivoCanonicalJSON.encode(value))
        FileHandle.standardOutput.write(Data("\n".utf8))
    }
    private func diagnostic(_ value: String) {
        FileHandle.standardError.write(Data(("numivivo protocol: " + value + "\n").utf8))
    }
    private struct ValidationReport: Encodable {
        let planFingerprint: VivoFingerprint
        let executable: Bool
        let stages: [VivoMDCapabilityReport]
    }
    private struct ExecutionReport: Encodable {
        let receipt: VivoMDProtocolRunReceipt
        let receiptArtifact: VivoStoredArtifact
    }
    private struct TrajectoryReport: Encodable {
        let manifest: VivoMDTrajectoryManifest
        let indexedChunks: Int
        let allPayloadsVerified: Bool
    }

    static let help = """
    NumiVivo staged MD protocols and bounded trajectory archives

      numivivo md-protocol-template system.json > protocol.json
      numivivo md-protocol-validate protocol.json --system system.json --state initial.json
      numivivo md-protocol-run protocol.json --system system.json --state initial.json \
          --store ./md-artifacts > receipt.json
      numivivo md-protocol-resume protocol.json --system system.json --store ./md-artifacts \
          --reference <checkpoint-reference-from-receipt> > resumed-receipt.json
      numivivo md-protocol-inspect --store ./md-artifacts --reference <checkpoint-reference>
      numivivo md-trajectory-inspect --store ./md-artifacts --manifest <sha256> [--verify]

    --checkpoint <sha256> replaces --reference for resume/inspection. A resume forks
    a new run reference; the original prefix is never overwritten. Stages preserve
    accepted time, step, velocities and cell except for declared initialization.
    The template explicitly includes minimization, NVT, NPT and NPT production.
    Its durations are illustrative, not a claim of equilibration or accuracy.

    Position samples are binary chunks in VivoArtifactStore; observation-only ticks
    are persisted even when no positions are sampled. Checkpoint intervals flush
    the archive prefix before publishing the restart cursor. SHA-256 verifies bytes,
    not physical validity. --verify reads every bounded trajectory payload.

    An unsuccessful required minimization or rejected MD candidate blocks the run.
    No automatic timestep adaptation or stochastic retry is performed. Cancellation
    preserves the last exportable accepted boundary. A device failure falls back to
    the previously durable cursor. Interrupted minimization restarts its optimizer
    from the saved accepted geometry, not from serialized line-search history.

    Apple package/Metal compilation and numerical qualification remain required.
    """ + "\n"
}

private enum ProtocolCLIError: Error, CustomStringConvertible {
    case usage(String)
    var description: String { switch self { case .usage(let value): return value } }
}
private struct ProtocolArguments {
    let positionals: [String]
    let options: [String: String]
    let verify: Bool
    init(_ raw: [String]) throws {
        var positional: [String] = [], options: [String: String] = [:], verify = false, i = 0
        while i < raw.count {
            let token = raw[i]
            if token == "--verify" {
                guard !verify else { throw ProtocolCLIError.usage("duplicate --verify") }
                verify = true; i += 1
            } else if token.hasPrefix("--") {
                let key = String(token.dropFirst(2))
                guard options[key] == nil, i + 1 < raw.count, !raw[i + 1].hasPrefix("--") else {
                    throw ProtocolCLIError.usage("missing value or duplicate option \(token)")
                }
                options[key] = raw[i + 1]; i += 2
            } else { positional.append(token); i += 1 }
        }
        self.positionals = positional; self.options = options; self.verify = verify
    }
    func allow(_ keys: Set<String>, positionals count: Int) throws {
        guard positionals.count == count, options.keys.allSatisfy(keys.contains), !verify || keys.contains("verify") else {
            throw ProtocolCLIError.usage("unknown options or wrong positional count; see md-protocol-help")
        }
    }
    func required(_ key: String) throws -> String {
        guard let value = options[key], !value.isEmpty else { throw ProtocolCLIError.usage("missing --\(key)") }
        return value
    }
    func positive(_ key: String, fallback: UInt64) throws -> UInt64 {
        guard let raw = options[key] else { return fallback }
        guard let value = UInt64(raw), value > 0 else { throw ProtocolCLIError.usage("--\(key) requires positive UInt64") }
        return value
    }
}
