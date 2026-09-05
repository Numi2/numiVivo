import Foundation
import NumiVivoKit

struct VivoScreenedInferenceCLICommands {
    static func handles(_ name: String?) -> Bool {
        ["engagement-fit-screened", "engagement-screen-check", "screening-help"].contains(name ?? "")
    }
    func run(arguments: [String]) async -> Int32 {
        do {
            guard let command = arguments.first else { throw Usage.invalid("missing command") }
            if command == "screening-help" {
                guard arguments.count == 1 else { throw Usage.invalid("screening-help takes no arguments") }
                FileHandle.standardOutput.write(Data(Self.help.utf8)); return 0
            }
            let a = try Arguments(Array(arguments.dropFirst()))
            try a.validate(fitting: command == "engagement-fit-screened")
            let problem: VivoTargetPosteriorProblem = try load(a.input)
            let policy: VivoMetalTargetScreenConfiguration = try a.options["screen-policy"].map(load) ?? .init()
            let checkPolicy: VivoMetalScreenCheckPolicy = try a.options["check-policy"].map(load) ?? .init()
            let checkpoint: VivoPosteriorCheckpoint? = try a.options["resume"].map(load)
            let gpu = try await VivoMetalTargetLikelihoodScreen.make(problem: problem, configuration: policy)
            let prepared = try VivoPreparedTargetPosterior(problem)
            try checkpoint?.validate(for: prepared.plan.withScreening(gpu.screeningPolicy))
            let store = try a.options["store"].map { try VivoArtifactStore(rootURL: URL(fileURLWithPath: $0)) }
            if let store {
                let descriptor = try await store.put(data: VivoCanonicalJSON.encode(gpu.description),
                    kind: "metal-target-screen-description", mediaType: "application/json")
                guard descriptor.fingerprint == gpu.screeningPolicy.fingerprint else {
                    throw VivoPosteriorError.invalid("stored screen description identity mismatch")
                }
            }
            // Probes exercise permutation, single-candidate padding and repeated
            // batches before the CLI permits using this workspace for inference.
            let checks = try await VivoMetalTargetScreenChecks.run(problem: problem, screen: gpu, policy: checkPolicy)
            if let store {
                _ = try await store.put(data: VivoCanonicalJSON.encode(checks), kind: "metal-screen-probe-report", mediaType: "application/json")
            }
            if command == "engagement-screen-check" {
                try emit(VivoCanonicalJSON.encode(checks), to: a.options["output"], overwrite: a.force)
                return checks.deterministicProbeChecksPassed ? 0 : 75
            }
            guard checks.deterministicProbeChecksPassed else {
                throw VivoPosteriorError.numerical("screen batch-independence probes failed; inference was not started")
            }
            guard let store else { throw Usage.invalid("screened fitting requires --store to retain screen and probe provenance") }
            let checkpointURL = a.options["checkpoint"].map { URL(fileURLWithPath: $0) }
            let resumeSame = a.options["resume"].map { a.samePath($0, a.options["checkpoint"] ?? "") } ?? false
            let sink = VivoScreenedCheckpointSink(store: store, url: checkpointURL, overwriteInitial: a.force || resumeSame)
            let result = try await VivoTargetPosteriorFitter.run(problem, checkpoint: checkpoint,
                progress: { try await sink.publish($0) }, screen: gpu.screen())
            try result.validate()
            let data = try VivoCanonicalJSON.encode(result)
            _ = try await store.put(data: data, kind: "target-posterior-result", mediaType: "application/json")
            try emit(data, to: a.options["output"], overwrite: a.force)
            return result.posterior.completed ? 0 : 75
        } catch is CancellationError {
            stderr("cancelled; completed-stage checkpoints remain available"); return 130
        } catch let error as Usage {
            stderr(error.localizedDescription); return 64
        } catch let error as VivoPosteriorError {
            stderr(error.localizedDescription)
            switch error { case .numerical, .budget: return 75; case .invalid, .busy: return 65 }
        } catch {
            stderr(error.localizedDescription); return 65
        }
    }
    private func load<T: Decodable>(_ path: String) throws -> T {
        guard path != "-" else { throw Usage.invalid("provide a bounded regular input file") }
        return try VivoKineticsDocumentIO.read(T.self, from: URL(fileURLWithPath: path))
    }
    private func emit(_ data: Data, to path: String?, overwrite: Bool) throws {
        if let path, path != "-" {
            try VivoKineticsDocumentIO.write(data, to: URL(fileURLWithPath: path), overwrite: overwrite)
        } else {
            FileHandle.standardOutput.write(data); FileHandle.standardOutput.write(Data("\n".utf8))
        }
    }
    private func stderr(_ message: String) {
        FileHandle.standardError.write(Data("numivivo screening: \(message)\n".utf8))
    }
    private enum Usage: Error, LocalizedError {
        case invalid(String)
        var errorDescription: String? { switch self { case .invalid(let message): return message } }
    }
    private struct Arguments {
        let input: String
        let options: [String: String]
        let force: Bool
        init(_ raw: [String]) throws {
            var inputs: [String] = [], options: [String: String] = [:], force = false, i = 0
            while i < raw.count {
                let token = raw[i]
                if token == "--force" {
                    guard !force else { throw Usage.invalid("duplicate --force") }
                    force = true; i += 1
                } else if token.hasPrefix("--") {
                    let name = String(token.dropFirst(2))
                    guard !name.isEmpty, options[name] == nil, i + 1 < raw.count,
                          !raw[i + 1].hasPrefix("--"), !raw[i + 1].isEmpty else {
                        throw Usage.invalid("missing value or duplicate option: \(token)")
                    }
                    options[name] = raw[i + 1]; i += 2
                } else { inputs.append(token); i += 1 }
            }
            guard inputs.count == 1 else { throw Usage.invalid("one inference-problem file is required") }
            self.input = inputs[0]; self.options = options; self.force = force
        }
        func samePath(_ a: String, _ b: String) -> Bool {
            guard !a.isEmpty, !b.isEmpty else { return false }
            return URL(fileURLWithPath: a).standardizedFileURL.resolvingSymlinksInPath()
                == URL(fileURLWithPath: b).standardizedFileURL.resolvingSymlinksInPath()
        }
        func validate(fitting: Bool) throws {
            var allowed: Set<String> = ["screen-policy", "check-policy", "output", "store"]
            if fitting { allowed.formUnion(["resume", "checkpoint"]) }
            guard options.keys.allSatisfy(allowed.contains), !fitting || options["store"] != nil else {
                throw Usage.invalid("unsupported option, or missing --store for screened fitting")
            }
            let sourcePaths = [input] + ["screen-policy", "check-policy", "resume"].compactMap { options[$0] }
            if let output = options["output"], output != "-" {
                guard sourcePaths.allSatisfy({ !samePath(output, $0) }),
                      options["checkpoint"].map({ !samePath(output, $0) }) ?? true else {
                    throw Usage.invalid("output must differ from all inputs and checkpoint")
                }
            }
            if let checkpoint = options["checkpoint"] {
                guard checkpoint != "-", !samePath(input, checkpoint),
                      ["screen-policy", "check-policy"].compactMap({ options[$0] }).allSatisfy({ !samePath(checkpoint, $0) }) else {
                    throw Usage.invalid("checkpoint would replace an input configuration")
                }
            }
        }
    }
    static let help = """
    NumiVivo Metal screening with authoritative FP64 delayed acceptance

      engagement-screen-check <inference-problem.json>
          [--screen-policy policy.json] [--check-policy probes.json]
          [--output checks.json] [--store artifact-root]

      engagement-fit-screened <inference-problem.json> --store artifact-root
          [--screen-policy policy.json] [--check-policy probes.json]
          [--resume checkpoint.json] [--checkpoint checkpoint.json]
          [--output posterior.json] [--force]

    Fitting keeps the FP64 likelihood, priors and assay model authoritative.
    Metal only screens proposals; survivors receive an FP64 acceptance correction.
    The returned record works with engagement-predict and engagement-design.
    The GPU screen supports exact-valued observations. Within-case correlation
    is omitted only by the screen and retained in the FP64 correction.
    Censored-data screening is rejected. There is no silent backend fallback.
    Probes run before fitting; finite probes are not universal qualification.
    Device, executable, shaders, fixed prior-bounded schedule and memory policy
    are fingerprinted. Changed screening identity cannot resume an old run.
    No speedup or biological validity follows merely from completing this command.
    """ + "\n"
}

private actor VivoScreenedCheckpointSink {
    let store: VivoArtifactStore
    let url: URL?
    let overwriteInitial: Bool
    var published = false
    init(store: VivoArtifactStore, url: URL?, overwriteInitial: Bool) {
        self.store = store; self.url = url; self.overwriteInitial = overwriteInitial
    }
    func publish(_ checkpoint: VivoPosteriorCheckpoint) async throws {
        let data = try VivoCanonicalJSON.encode(checkpoint)
        _ = try await store.put(data: data, kind: "posterior-checkpoint", mediaType: "application/json")
        if let url { try VivoKineticsDocumentIO.write(data, to: url, overwrite: overwriteInitial || published) }
        published = true
    }
}
