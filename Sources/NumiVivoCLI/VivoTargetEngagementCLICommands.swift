import Foundation
import NumiVivoKit

struct VivoTargetEngagementCLICommands {
    static let commands: Set<String> = ["engagement-help", "engagement-validate", "engagement-source",
        "engagement-compile", "engagement-run", "engagement-batch", "engagement-study", "engagement-rate", "engagement-apply-rate"]
    static func handles(_ name: String?) -> Bool { commands.contains(name ?? "") }

    func run(arguments: [String]) async -> Int32 {
        do {
            guard let command = arguments.first else { throw Usage.invalid("missing command") }
            if command == "engagement-help" {
                guard arguments.count == 1 else { throw Usage.invalid("engagement-help accepts no arguments") }
                FileHandle.standardOutput.write(Data(Self.help.utf8)); return 0
            }
            let a = try Arguments(Array(arguments.dropFirst()))
            switch command {
            case "engagement-validate":
                try a.allow([])
                let experiment: VivoTargetEngagementExperiment = try load(a.input)
                try experiment.validate()
                let report = ValidationReport(schemaVersion: 1,
                    experimentFingerprint: try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(experiment)),
                    containsAssumptions: experiment.kinetics.containsAssumptions || experiment.exposure.origin == .assumed,
                    unknownParameterUncertainty: experiment.kinetics.hasUnknownParameterUncertainty,
                    status: "valid-input-contract; execution and biological validity not established")
                try await emit(VivoCanonicalJSON.encode(report), arguments: a, kind: "target-engagement-validation")
            case "engagement-source":
                try a.allow([])
                let experiment: VivoTargetEngagementExperiment = try load(a.input)
                let e = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(experiment))
                let k = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(experiment.kinetics))
                let source = try VivoTargetEngagementProgramSource.make(experiment, kineticsFingerprint: k.hex, experimentFingerprint: e.hex)
                try await emit(source, arguments: a, kind: "target-engagement-source")
            case "engagement-compile":
                try a.allow([])
                guard a.options["output"] != nil, a.options["output"] != "-" else {
                    throw Usage.invalid("binary ProgramPack requires --output <file>")
                }
                let experiment: VivoTargetEngagementExperiment = try load(a.input)
                let compiled = try VivoTargetEngagementCompiler.compile(experiment)
                try await emit(compiled.pack.data, arguments: a, kind: "target-engagement-program",
                               mediaType: "application/octet-stream")
            case "engagement-run":
                try a.allow(["backend", "policy"])
                let experiment: VivoTargetEngagementExperiment = try load(a.input)
                switch a.options["backend"] ?? "reference" {
                case "reference":
                    let policy: VivoTargetEngagementNumerics = try a.options["policy"].map(load) ?? .init()
                    let record = try VivoTargetEngagementRunRecord.reference(experiment, numerics: policy)
                    try record.validate()
                    try await emit(VivoCanonicalJSON.encode(record), arguments: a, kind: "target-engagement-reference-result")
                case "metal":
                    let policy: VivoTargetEngagementMetalPolicy = try a.options["policy"].map(load) ?? .init()
                    let result = try await VivoTargetEngagementMetalRunner.run([experiment], policy: policy)
                    try await emit(VivoCanonicalJSON.encode(result), arguments: a, kind: "target-engagement-metal-result")
                default: throw Usage.invalid("--backend must be reference or metal; no silent fallback")
                }
            case "engagement-batch":
                try a.allow(["policy"])
                let experiments: [VivoTargetEngagementExperiment] = try load(a.input)
                let policy: VivoTargetEngagementMetalPolicy = try a.options["policy"].map(load) ?? .init()
                let result = try await VivoTargetEngagementMetalRunner.run(experiments, policy: policy)
                try await emit(VivoCanonicalJSON.encode(result), arguments: a, kind: "target-engagement-metal-result")
            case "engagement-study":
                try a.allow(["policy"])
                let study: VivoTargetEngagementStudy = try load(a.input)
                let policy: VivoTargetEngagementNumerics = try a.options["policy"].map(load) ?? .init()
                let result = try VivoTargetEngagementStudyEvaluator.evaluate(study, numerics: policy)
                try await emit(VivoCanonicalJSON.encode(result), arguments: a, kind: "target-engagement-study-result")
                return result.cases.contains(where: { $0.failure != nil }) ? 75 : 0
            case "engagement-rate":
                try a.allow([])
                let request: VivoTransitionStateRateRequest = try load(a.input)
                let derivation = try VivoTransitionStateDerivation.calculate(request)
                try await emit(derivation.evidenceData, arguments: a, kind: "conditional-transition-state-derivation")
            case "engagement-apply-rate":
                try a.allow(["rate"])
                guard let ratePath = a.options["rate"], let storePath = a.options["store"] else {
                    throw Usage.invalid("applying a rate requires --rate <request.json> and --store <artifact-root>")
                }
                let experiment: VivoTargetEngagementExperiment = try load(a.input)
                let request: VivoTransitionStateRateRequest = try load(ratePath)
                let derivation = try VivoTransitionStateDerivation.calculate(request)
                let model = try derivation.applyingInactivation(to: experiment.kinetics)
                let updated = VivoTargetEngagementExperiment(kinetics: model, exposure: experiment.exposure,
                    initial: experiment.initial, sampleTimesSeconds: experiment.sampleTimesSeconds)
                try updated.validate()
                let store = try VivoArtifactStore(rootURL: URL(fileURLWithPath: storePath))
                let descriptor = try await store.put(data: derivation.evidenceData,
                    kind: "conditional-transition-state-derivation", mediaType: "application/json")
                guard descriptor.fingerprint == derivation.evidenceFingerprint else {
                    throw VivoKineticsError.invalid("stored derivation identity mismatch")
                }
                try await emit(VivoCanonicalJSON.encode(updated), arguments: a, kind: "target-engagement-experiment")
            default: throw Usage.invalid("unknown engagement command")
            }
            return 0
        } catch let error as Usage {
            stderr(error.localizedDescription); return 64
        } catch let error as VivoKineticsError {
            stderr(error.localizedDescription)
            switch error { case .numerical, .capacity: return 75; case .invalid, .unsupported: return 65 }
        } catch is CancellationError {
            stderr("cancelled; no incomplete result was published"); return 130
        } catch {
            stderr(error.localizedDescription); return 65
        }
    }

    private struct ValidationReport: Encodable {
        let schemaVersion: UInt32
        let experimentFingerprint: VivoFingerprint
        let containsAssumptions: Bool
        let unknownParameterUncertainty: Bool
        let status: String
    }
    private func load<T: Decodable>(_ path: String) throws -> T {
        guard path != "-" else { throw Usage.invalid("provide a bounded regular input file, not stdin") }
        return try VivoKineticsDocumentIO.read(T.self, from: URL(fileURLWithPath: path))
    }
    private func emit(_ data: Data, arguments: Arguments, kind: String,
                      mediaType: String = "application/json") async throws {
        guard data.count <= 536_870_912 else { throw VivoKineticsError.capacity("output document") }
        // Immutable store publication and file output are individually atomic,
        // not a cross-filesystem transaction. A failed file output may leave a
        // valid content-addressed object in the explicitly requested store.
        if let root = arguments.options["store"] {
            let store = try VivoArtifactStore(rootURL: URL(fileURLWithPath: root))
            _ = try await store.put(data: data, kind: kind, mediaType: mediaType)
        }
        if let destination = arguments.options["output"], destination != "-" {
            try VivoKineticsDocumentIO.write(data, to: URL(fileURLWithPath: destination), overwrite: arguments.force)
        } else {
            FileHandle.standardOutput.write(data)
            FileHandle.standardOutput.write(Data("\n".utf8))
        }
    }
    private func stderr(_ message: String) { FileHandle.standardError.write(Data("numivivo engagement: \(message)\n".utf8)) }

    private enum Usage: Error, LocalizedError {
        case invalid(String)
        var errorDescription: String? { switch self { case .invalid(let message): return message } }
    }
    private struct Arguments {
        let input: String
        let options: [String: String]
        let force: Bool
        init(_ raw: [String]) throws {
            var paths: [String] = [], options: [String: String] = [:], force = false, index = 0
            while index < raw.count {
                let token = raw[index]
                if token == "--force" {
                    guard !force else { throw Usage.invalid("duplicate --force") }
                    force = true; index += 1
                } else if token.hasPrefix("--") {
                    let name = String(token.dropFirst(2))
                    guard !name.isEmpty, options[name] == nil, index + 1 < raw.count,
                          !raw[index + 1].hasPrefix("--"), !raw[index + 1].isEmpty else {
                        throw Usage.invalid("missing value or duplicate option: \(token)")
                    }
                    options[name] = raw[index + 1]; index += 2
                } else { paths.append(token); index += 1 }
            }
            guard paths.count == 1 else { throw Usage.invalid("one input file is required; see engagement-help") }
            self.input = paths[0]; self.options = options; self.force = force
        }
        func allow(_ extra: Set<String>) throws {
            let allowed = extra.union(["output", "store"])
            guard options.keys.allSatisfy(allowed.contains) else { throw Usage.invalid("unsupported option for command") }
            if let output = options["output"], output != "-" {
                let outputURL = URL(fileURLWithPath: output).standardizedFileURL.resolvingSymlinksInPath()
                for path in [input] + ["policy", "rate"].compactMap({ options[$0] }) {
                    let inputURL = URL(fileURLWithPath: path).standardizedFileURL.resolvingSymlinksInPath()
                    guard inputURL != outputURL else { throw Usage.invalid("input and output paths must differ") }
                }
            }
        }
    }

    static let help = """
    NumiVivo target engagement and conditional kinetics

      engagement-validate <experiment.json>
      engagement-source <experiment.json> --output program.json
      engagement-compile <experiment.json> --output program.nvpack
      engagement-run <experiment.json> [--backend reference|metal] [--policy policy.json]
      engagement-batch <experiments.json> [--policy metal-policy.json]
      engagement-study <study.json> [--policy reference-policy.json]
      engagement-rate <activation-free-energy-request.json>
      engagement-apply-rate <experiment.json> --rate request.json --store artifact-root

    Output options: --output <path|->, --store <artifact-root>, --force.
    Reference is the default; Metal never silently falls back. Batch uses common
    kinetics/initial states/observation times with different unbound exposures.
    Inputs use seconds, molar unbound concentrations and explicitly sourced rates.
    Source/compile commands reuse the existing typed F1 ProgramPack compiler.
    No command generates quantum chemistry, predicts clinical efficacy, or gives
    a patient-specific treatment instruction. Missing uncertainty stays unknown.
    """ + "\n"
}
