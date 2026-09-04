import Foundation
import NumiVivoKit

struct VivoHybridCLICommands {
    static func handles(_ name: String?) -> Bool {
        name == "plan-hybrid" || name == "run-hybrid" || name == "hybrid-help"
    }

    func run(arguments: [String]) async -> Int32 {
        do {
            if arguments.first == "hybrid-help" {
                FileHandle.standardOutput.write(Data(Self.help.utf8))
                return 0
            }
            let args = try HybridArguments(arguments)
            let source = try VivoValidatedArtifactLoader.data(at: URL(fileURLWithPath: args.model))
            let model = try VivoGeneratedFingerprintSourceCompiler.exactSSAModel(from: source).value
            if args.command == "plan-hybrid" {
                try args.allow(["authorities", "maximum-step", "output", "force"])
                let authorities: [String: VivoHybridExecutionMode]
                if let path = args.options["authorities"] {
                    authorities = try load([String: VivoHybridExecutionMode].self, path: path)
                } else { authorities = [:] }
                let plan = try VivoHybridExplicitPlanBuilder.make(
                    model: model, reactionAuthorities: authorities,
                    maximumStep: try args.number("maximum-step", default: 0.01)
                )
                try write(try VivoCanonicalJSON.encode(plan), to: args.options["output"] ?? "-", force: args.force)
                return 0
            }
            try args.allow(["plan", "initial-counts", "restore", "lanes", "steps", "dt", "seed",
                            "exact-events", "exact-dispatches", "publications", "output", "checkpoint", "force"])
            let plan = try load(VivoHybridStochasticPlan.self, path: args.required("plan"))
            let lanes = try args.integer("lanes", default: 1)
            let steps = try args.integer("steps", default: 1)
            guard steps <= 100_000 else { throw HybridCLIError.usage("--steps cannot exceed 100000 per invocation") }
            let dtDouble = try args.number("dt", default: 0.01)
            guard dtDouble <= Double(Float.greatestFiniteMagnitude), Float(dtDouble) > 0 else {
                throw HybridCLIError.usage("--dt must be a positive FP32 value")
            }
            let dt = Float(dtDouble)
            let output = args.options["output"] ?? "-"
            let checkpointPath = args.options["checkpoint"]
            try preflight(output, force: args.force)
            if let checkpointPath {
                guard checkpointPath != "-",
                      URL(fileURLWithPath: checkpointPath).standardizedFileURL != URL(fileURLWithPath: output).standardizedFileURL else {
                    throw HybridCLIError.usage("checkpoint and primary output must have different destinations")
                }
                try preflight(checkpointPath, force: args.force)
            }
            let restore: VivoHybridCheckpoint?
            let counts: [UInt32]
            if let path = args.options["restore"] {
                guard args.options["initial-counts"] == nil else {
                    throw HybridCLIError.usage("--restore and --initial-counts are mutually exclusive")
                }
                let value = try load(VivoHybridCheckpoint.self, path: path)
                guard value.laneCount == lanes, value.speciesCount == UInt32(model.species.count) else {
                    throw HybridCLIError.usage("checkpoint shape does not match --lanes and model species")
                }
                restore = value
                let total = UInt64(lanes) * UInt64(model.species.count)
                guard total <= 16_777_216 else { throw HybridCLIError.usage("CLI initial-state limit is 16777216 elements") }
                counts = [UInt32](repeating: 0, count: Int(total))
            } else {
                restore = nil
                counts = try load([UInt32].self, path: args.required("initial-counts"))
            }
            let publications: [VivoHybridPublicationRequest]
            if let path = args.options["publications"] { publications = try load([VivoHybridPublicationRequest].self, path: path) }
            else { publications = [] }
            let seed: UInt64
            if let raw = args.options["seed"] {
                guard let parsed = UInt64(raw) else { throw HybridCLIError.usage("--seed requires a decimal UInt64") }
                seed = parsed
            } else { seed = restore?.seed ?? 0x4e554d495649564f }
            let configuration = VivoHybridRuntimeConfiguration(
                laneCount: lanes, timeStep: dt,
                minimumTimeStep: min(dt, 1e-7), maximumTimeStep: max(dt, 60), seed: seed,
                exactEventsPerDispatch: try args.integer("exact-events", default: 256),
                maximumExactDispatches: try args.integer("exact-dispatches", default: 64)
            )
            let runtime = try await VivoHybridReactionRuntime.make(
                model: model, plan: plan, configuration: configuration, initialCounts: counts
            )
            if let restore { try await runtime.restore(restore) }
            var chain = VivoDigestChain(genesis: try VivoCanonicalJSON.fingerprint(source))
            var committed: UInt32 = 0
            var last: VivoHybridStepCertificate?
            for _ in 0..<steps {
                try Task.checkCancellation()
                let result = try await runtime.step(deltaTime: dt, publications: publications)
                _ = try chain.append(result)
                last = result
                guard result.disposition == .committed else { break }
                committed += 1
            }
            let snapshot = try await runtime.snapshot()
            var checkpointFingerprint: VivoFingerprint?
            if let checkpointPath {
                let checkpoint = try await runtime.checkpoint()
                checkpointFingerprint = try checkpoint.fingerprint()
                try write(try VivoCanonicalJSON.encode(checkpoint), to: checkpointPath, force: args.force)
            }
            let report = HybridRunReport(schemaVersion: 1, deviceName: runtime.deviceName,
                                          requestedSteps: steps, committedSteps: committed,
                                          modelFingerprint: runtime.execution.modelFingerprint,
                                          executablePlanFingerprint: runtime.execution.planFingerprint,
                                          previousCheckpointFingerprint: try restore?.fingerprint(),
                                          checkpointFingerprint: checkpointFingerprint,
                                          ledger: chain, lastStep: last, snapshot: snapshot)
            try write(try VivoCanonicalJSON.encode(report), to: output, force: args.force)
            return committed == steps ? 0 : 75
        } catch let error as HybridCLIError {
            diagnostic(String(describing: error))
            return 64
        } catch let error as VivoHybridExecutionError {
            diagnostic(error.localizedDescription)
            switch error {
            case .resourceLimit: return 75
            case .metal: return 70
            default: return 65
            }
        } catch is CancellationError {
            diagnostic("cancelled; no in-flight candidate was published")
            return 75
        } catch {
            diagnostic(error.localizedDescription)
            return 65
        }
    }

    private func load<T: Decodable & Sendable>(_ type: T.Type, path: String) throws -> T {
        let bytes = try VivoValidatedArtifactLoader.data(at: URL(fileURLWithPath: path))
        return try VivoValidatedArtifactLoader.decode(type, from: bytes).value
    }
    private func preflight(_ path: String, force: Bool) throws {
        guard path != "-" else { return }
        guard force || !FileManager.default.fileExists(atPath: path) else {
            throw HybridCLIError.usage("output exists: \(path); use --force to replace it")
        }
    }
    private func write(_ data: Data, to path: String, force: Bool) throws {
        if path == "-" {
            FileHandle.standardOutput.write(data)
            FileHandle.standardOutput.write(Data("\n".utf8))
            return
        }
        try preflight(path, force: force)
        let url = URL(fileURLWithPath: path)
        try FileManager.default.createDirectory(at: url.deletingLastPathComponent(), withIntermediateDirectories: true)
        try data.write(to: url, options: .atomic)
    }
    private func diagnostic(_ message: String) {
        FileHandle.standardError.write(Data(("numivivo hybrid: " + message + "\n").utf8))
    }
    static let help = """
    NumiVivo executable hybrid reaction runtime

      numivivo plan-hybrid <model.json> [--authorities <mapping.json>]
          [--maximum-step 0.01] [--output <plan.json>]
      numivivo run-hybrid <model.json> --plan <plan.json>
          --initial-counts <counts.json> --lanes 1 --steps 100 --dt 0.01
          [--publications <requests.json>] [--output <result.json>]
          [--checkpoint <checkpoint.json>] [--seed <decimal UInt64>]
          [--exact-events 256] [--exact-dispatches 64] [--force]
      numivivo run-hybrid <model.json> --plan <plan.json>
          --restore <checkpoint.json> --lanes 1 --steps 100 --dt 0.01

    Unspecified connected components use exactSSA. Authority mapping values are
    exactSSA, tauLeap or deterministicRK2. All reactions connected through any
    species dependency must share one authority. Counts use species-major order.
    Work-budget exhaustion returns exit 75 without advancing the rejected step.
    This path supports the sequence-free VivoExactSSAModel kinetic source format;
    it does not silently execute full VivoProgram rules, delays or spatial terms.
    """ + "\n"
}

private struct HybridRunReport: Encodable {
    let schemaVersion: UInt32
    let deviceName: String
    let requestedSteps: UInt32
    let committedSteps: UInt32
    let modelFingerprint: String
    let executablePlanFingerprint: String
    let previousCheckpointFingerprint: VivoFingerprint?
    let checkpointFingerprint: VivoFingerprint?
    let ledger: VivoDigestChain
    let lastStep: VivoHybridStepCertificate?
    let snapshot: VivoHybridStateSnapshot
}
private enum HybridCLIError: Error, CustomStringConvertible {
    case usage(String)
    var description: String { switch self { case .usage(let reason): return reason } }
}
private struct HybridArguments {
    let command: String
    let model: String
    let options: [String: String]
    let force: Bool
    init(_ raw: [String]) throws {
        guard raw.count >= 2 else { throw HybridCLIError.usage("expected a model path; see numivivo hybrid-help") }
        command = raw[0]
        model = raw[1]
        var values: [String: String] = [:]
        var force = false
        var i = 2
        while i < raw.count {
            let token = raw[i]
            guard token.hasPrefix("--") else { throw HybridCLIError.usage("unexpected positional argument \(token)") }
            let name = String(token.dropFirst(2))
            if name == "force" {
                guard !force else { throw HybridCLIError.usage("duplicate --force") }
                force = true
            } else {
                guard values[name] == nil, i + 1 < raw.count, !raw[i + 1].hasPrefix("--") else {
                    throw HybridCLIError.usage("duplicate option or missing value: \(token)")
                }
                values[name] = raw[i + 1]
                i += 1
            }
            i += 1
        }
        options = values
        self.force = force
    }
    func allow(_ allowed: Set<String>) throws {
        guard options.keys.allSatisfy({ allowed.contains($0) }) else { throw HybridCLIError.usage("unknown option") }
    }
    func required(_ name: String) throws -> String {
        guard let value = options[name], !value.isEmpty else { throw HybridCLIError.usage("missing --\(name)") }
        return value
    }
    func integer(_ name: String, default fallback: UInt32) throws -> UInt32 {
        guard let raw = options[name] else { return fallback }
        guard let value = UInt32(raw), value > 0 else { throw HybridCLIError.usage("--\(name) requires a positive UInt32") }
        return value
    }
    func number(_ name: String, default fallback: Double) throws -> Double {
        guard let raw = options[name] else { return fallback }
        guard let value = Double(raw), value.isFinite, value > 0 else { throw HybridCLIError.usage("--\(name) requires a positive finite number") }
        return value
    }
}
