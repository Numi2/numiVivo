import Foundation
import NumiVivoKit

struct VivoPosteriorCLICommands {
    static func handles(_ name: String?) -> Bool {
        ["engagement-fit", "engagement-predict", "engagement-sensitivity", "posterior-help"].contains(name ?? "")
    }
    func run(arguments: [String]) async -> Int32 {
        do {
            guard let command = arguments.first else { throw CLIError.invalid("missing command") }
            if command == "posterior-help" {
                guard arguments.count == 1 else { throw CLIError.invalid("posterior-help takes no arguments") }
                FileHandle.standardOutput.write(Data(Self.help.utf8)); return 0
            }
            let a = try Arguments(Array(arguments.dropFirst()))
            let store = try a.options["store"].map { try VivoArtifactStore(rootURL: URL(fileURLWithPath: $0)) }
            switch command {
            case "engagement-fit":
                try a.allow(["resume", "checkpoint"])
                let problem: VivoTargetPosteriorProblem = try load(a.input)
                let checkpoint: VivoPosteriorCheckpoint? = try a.options["resume"].map(load)
                let prepared = try VivoPreparedTargetPosterior(problem)
                try checkpoint?.validate(for: prepared.plan)
                let writer = PosteriorCheckpointWriter(store: store,
                    url: a.options["checkpoint"].map { URL(fileURLWithPath: $0) },
                    overwriteInitial: a.force || a.options["resume"].map { a.samePath($0, a.options["checkpoint"] ?? "") } == true)
                let result = try await VivoTargetPosteriorFitter.run(problem, checkpoint: checkpoint,
                    progress: { try await writer.publish($0) })
                try result.validate()
                try await emit(VivoCanonicalJSON.encode(result), arguments: a, store: store, kind: "target-posterior-result")
                return result.posterior.completed ? 0 : 75
            case "engagement-predict":
                try a.allow(["policy"])
                let record: VivoTargetPosteriorRecord = try load(a.input)
                let policy: VivoPosteriorPredictivePolicy = try a.options["policy"].map(load) ?? .init()
                try record.validate(requireComplete: true)
                if let store {
                    _ = try await store.put(data: VivoCanonicalJSON.encode(record), kind: "target-posterior-result", mediaType: "application/json")
                }
                let report = try await VivoTargetPosteriorPredictor.predict(record, policy: policy)
                try await emit(VivoCanonicalJSON.encode(report), arguments: a, store: store, kind: "target-posterior-predictive")
                return report.cases.allSatisfy { $0.failures.isEmpty } ? 0 : 75
            case "engagement-sensitivity":
                try a.allow(["policy"])
                let record: VivoTargetPosteriorRecord = try load(a.input)
                let policy: VivoPosteriorSensitivityPolicy = try a.options["policy"].map(load) ?? .init()
                try record.validate(requireComplete: true)
                if let store {
                    _ = try await store.put(data: VivoCanonicalJSON.encode(record), kind: "target-posterior-result", mediaType: "application/json")
                }
                let report = try VivoTargetPosteriorSensitivity.evaluate(record, policy: policy)
                try await emit(VivoCanonicalJSON.encode(report), arguments: a, store: store, kind: "target-posterior-sensitivity")
                return report.derivativeChecksPassed ? 0 : 75
            default: throw CLIError.invalid("unsupported posterior command")
            }
        } catch is CancellationError {
            error("cancelled; only previously committed checkpoints were published"); return 130
        } catch {
            self.error(error.localizedDescription); return 65
        }
    }

    private func load<T: Decodable>(_ path: String) throws -> T {
        guard path != "-" else { throw CLIError.invalid("a bounded regular input file is required") }
        return try VivoKineticsDocumentIO.read(T.self, from: URL(fileURLWithPath: path))
    }
    private func emit(_ data: Data, arguments: Arguments, store: VivoArtifactStore?, kind: String) async throws {
        if let store { _ = try await store.put(data: data, kind: kind, mediaType: "application/json") }
        if let path = arguments.options["output"], path != "-" {
            try VivoKineticsDocumentIO.write(data, to: URL(fileURLWithPath: path), overwrite: arguments.force)
        } else {
            FileHandle.standardOutput.write(data); FileHandle.standardOutput.write(Data("\n".utf8))
        }
    }
    private func error(_ value: String) { FileHandle.standardError.write(Data("numivivo posterior: \(value)\n".utf8)) }

    private enum CLIError: Error, LocalizedError {
        case invalid(String)
        var errorDescription: String? { switch self { case .invalid(let message): return message } }
    }
    private struct Arguments {
        let input: String
        let options: [String: String]
        let force: Bool
        init(_ raw: [String]) throws {
            var paths: [String] = [], options: [String: String] = [:], force = false, i = 0
            while i < raw.count {
                let token = raw[i]
                if token == "--force" {
                    guard !force else { throw CLIError.invalid("duplicate --force") }
                    force = true; i += 1
                } else if token.hasPrefix("--") {
                    let name = String(token.dropFirst(2))
                    guard !name.isEmpty, options[name] == nil, i + 1 < raw.count,
                          !raw[i + 1].isEmpty, !raw[i + 1].hasPrefix("--") else {
                        throw CLIError.invalid("missing value or duplicate option \(token)")
                    }
                    options[name] = raw[i + 1]; i += 2
                } else { paths.append(token); i += 1 }
            }
            guard paths.count == 1 else { throw CLIError.invalid("one input file required; see posterior-help") }
            self.input = paths[0]; self.options = options; self.force = force
        }
        func samePath(_ a: String, _ b: String) -> Bool {
            URL(fileURLWithPath: a).standardizedFileURL.resolvingSymlinksInPath()
                == URL(fileURLWithPath: b).standardizedFileURL.resolvingSymlinksInPath()
        }
        func allow(_ extra: Set<String>) throws {
            guard options.keys.allSatisfy(extra.union(["output", "store"]).contains) else { throw CLIError.invalid("unsupported option") }
            let inputs = [input] + ["resume", "policy"].compactMap { options[$0] }
            if let output = options["output"], output != "-" {
                guard inputs.allSatisfy({ !samePath($0, output) }), options["checkpoint"].map({ !samePath($0, output) }) ?? true else {
                    throw CLIError.invalid("result output cannot overwrite an input or checkpoint")
                }
            }
            if let checkpoint = options["checkpoint"] {
                guard checkpoint != "-", !samePath(checkpoint, input) else { throw CLIError.invalid("checkpoint must be a separate file") }
            }
        }
    }

    static let help = """
    NumiVivo conditional kinetic posterior inference

      engagement-fit <inference-problem.json> [--resume checkpoint.json]
          [--checkpoint checkpoint.json] [--output posterior.json]
      engagement-predict <posterior.json> [--policy predictive-policy.json]
          [--output predictive.json]
      engagement-sensitivity <posterior.json> [--policy sensitivity-policy.json]
          [--output sensitivity.json]

    Shared options: --store <artifact-root>, --force.
    Fitting uses only study cases labelled calibration. Priors are explicit
    bounded uniform-physical or uniform-log-physical distributions. Correlated
    SMC mutation and posterior prediction preserve JOINT parameter particles.
    The current kinetic likelihood is deterministic native Swift FP64; it does
    not silently replace this likelihood with an approximate Metal evaluation.
    Checkpoints bind training data, priors, sampler and forward numerical policy.
    An interrupted or failed tempered population is not a completed posterior.
    Numerical failures remain failures, not zero likelihood or finite penalties.
    Unknown assay SD cannot enter fitting as zero. Predictive intervals remain
    conditional on the model; no clinical efficacy/safety conclusion is inferred.
    """ + "\n"
}

private actor PosteriorCheckpointWriter {
    private let store: VivoArtifactStore?
    private let url: URL?
    private let overwriteInitial: Bool
    private var published = false
    init(store: VivoArtifactStore?, url: URL?, overwriteInitial: Bool) {
        self.store = store; self.url = url; self.overwriteInitial = overwriteInitial
    }
    func publish(_ checkpoint: VivoPosteriorCheckpoint) async throws {
        guard store != nil || url != nil else { return }
        let data = try VivoCanonicalJSON.encode(checkpoint)
        if let store { _ = try await store.put(data: data, kind: "posterior-checkpoint", mediaType: "application/json") }
        if let url { try VivoKineticsDocumentIO.write(data, to: url, overwrite: overwriteInitial || published) }
        published = true
    }
}
