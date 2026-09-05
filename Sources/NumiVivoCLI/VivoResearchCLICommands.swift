import Foundation
import NumiVivoKit

struct VivoResearchCLICommands {
    static func handles(_ name: String?) -> Bool {
        ["finite-drug-run", "engagement-design", "occupancy-benchmark-inspect", "occupancy-benchmark-compare", "research-help"].contains(name ?? "")
    }
    func run(arguments: [String]) async -> Int32 {
        do {
            guard let command = arguments.first else { throw Usage.invalid("missing command") }
            if command == "research-help" {
                guard arguments.count == 1 else { throw Usage.invalid("research-help takes no arguments") }
                FileHandle.standardOutput.write(Data(Self.help.utf8)); return 0
            }
            let a = try Arguments(Array(arguments.dropFirst()))
            let output: Data, kind: String
            switch command {
            case "finite-drug-run":
                try a.allow([])
                let experiment: VivoFiniteDrugExperiment = try load(a.input)
                output = try VivoCanonicalJSON.encode(VivoFiniteDrugRunRecord.run(experiment)); kind = "finite-drug-result"
            case "engagement-design":
                try a.allow(["candidates", "policy"])
                guard let path = a.options["candidates"] else { throw Usage.invalid("--candidates <file> is required") }
                let posterior: VivoTargetPosteriorRecord = try load(a.input)
                let candidates: [VivoTargetDesignCandidate] = try load(path)
                let policy: VivoTargetDesignPolicy = try a.options["policy"].map(load) ?? .init()
                output = try VivoCanonicalJSON.encode(VivoTargetExperimentDesigner.rank(record: posterior, candidates: candidates, policy: policy))
                kind = "target-measurement-design"
            case "occupancy-benchmark-inspect":
                try a.allow([])
                let benchmark: VivoAggregateOccupancyBenchmark = try load(a.input)
                try benchmark.validate()
                output = try VivoCanonicalJSON.encode(benchmark); kind = "aggregate-occupancy-benchmark"
            case "occupancy-benchmark-compare":
                try a.allow([])
                let request: VivoAggregateBenchmarkRequest = try load(a.input)
                output = try VivoCanonicalJSON.encode(VivoAggregateOccupancyEvaluator.compare(request)); kind = "aggregate-occupancy-comparison"
            default: throw Usage.invalid("unsupported command")
            }
            if let path = a.options["store"] {
                let store = try VivoArtifactStore(rootURL: URL(fileURLWithPath: path))
                _ = try await store.put(data: output, kind: kind, mediaType: "application/json")
            }
            if let path = a.options["output"], path != "-" {
                try VivoKineticsDocumentIO.write(output, to: URL(fileURLWithPath: path), overwrite: a.force)
            } else { FileHandle.standardOutput.write(output); FileHandle.standardOutput.write(Data("\n".utf8)) }
            return 0
        } catch is CancellationError { error("cancelled; no partial numerical result published"); return 130 }
        catch { self.error(error.localizedDescription); return 65 }
    }
    private func load<T: Decodable>(_ path: String) throws -> T {
        guard path != "-" else { throw Usage.invalid("bounded input file required") }
        return try VivoKineticsDocumentIO.read(T.self, from: URL(fileURLWithPath: path))
    }
    private func error(_ value: String) { FileHandle.standardError.write(Data("numivivo research: \(value)\n".utf8)) }
    private enum Usage: Error, LocalizedError {
        case invalid(String)
        var errorDescription: String? { switch self { case .invalid(let value): return value } }
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
                    guard !force else { throw Usage.invalid("duplicate --force") }
                    force = true; i += 1
                } else if token.hasPrefix("--") {
                    let key = String(token.dropFirst(2))
                    guard !key.isEmpty, options[key] == nil, i + 1 < raw.count,
                          !raw[i + 1].isEmpty, !raw[i + 1].hasPrefix("--") else { throw Usage.invalid("option missing value or duplicated") }
                    options[key] = raw[i + 1]; i += 2
                } else { paths.append(token); i += 1 }
            }
            guard paths.count == 1 else { throw Usage.invalid("one input file required; see research-help") }
            self.input = paths[0]; self.options = options; self.force = force
        }
        func allow(_ additional: Set<String>) throws {
            guard options.keys.allSatisfy(additional.union(["output", "store"]).contains) else { throw Usage.invalid("unknown option") }
            if let output = options["output"], output != "-" {
                let target = URL(fileURLWithPath: output).standardizedFileURL.resolvingSymlinksInPath()
                for path in [input] + ["policy", "candidates"].compactMap({ options[$0] }) {
                    guard URL(fileURLWithPath: path).standardizedFileURL.resolvingSymlinksInPath() != target else {
                        throw Usage.invalid("output cannot overwrite an input")
                    }
                }
            }
        }
    }
    static let help = """
    NumiVivo portable research operations
      finite-drug-run <finite-experiment.json>
      engagement-design <posterior.json> --candidates <candidate-array.json> [--policy policy.json]
      occupancy-benchmark-inspect <aggregate-observations.json>
      occupancy-benchmark-compare <benchmark-comparison-request.json>
    Output options: --output <path|->, --store <artifact-root>, --force.
    Finite drug is a local reaction operator, not a circulation or treatment model.
    Design ranks proposed single measurements; it does not execute experiments.
    Aggregate cohorts are not individual measurements or fitted kinetic data.
    """ + "\n"
}
