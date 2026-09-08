import Foundation
import NumiVivoKit

struct VivoMolecularSamplingExportCLICommands {
    static func handles(_ name: String?) -> Bool {
        ["molecule-sampling-export", "molecule-sampling-export-verify"].contains(name ?? "")
    }

    func run(arguments raw: [String]) async -> Int32 {
        do {
            try Task.checkCancellation()
            let arguments = try Arguments(raw)
            let limits = try VivoMolecularSamplingCLILimits.load(arguments.options["--read-limits"])
            let root = try VivoWorkflowCLIDocumentPaths.canonicalURL(URL(fileURLWithPath: arguments.options["--store"]!))
            try Task.checkCancellation()
            let outputs = try OutputPlan(arguments: arguments, storeRoot: root)
            try Task.checkCancellation()
            // Reading and selection require an existing store. No typo in its
            // path creates a new empty archive during export or verification.
            let store = try VivoArtifactStore(rootURL: root, createIfNeeded: false)
            let implementation = try VivoWorkflowCLIImplementation.fingerprint()
            if arguments.command == "molecule-sampling-export" {
                let checkpoint: VivoFingerprint
                if let hash = arguments.options["--checkpoint"] { checkpoint = try .init(hex: hash) }
                else {
                    // Resolve once through this actor; never follow it again
                    // during validation, publication or optional file writes.
                    checkpoint = try await store.reference(arguments.options["--reference"]!,
                        maximumObjectBytes: limits.maximumCursorBytes).artifact.fingerprint
                }
                let reader = try await VivoMolecularSamplingArchiveReader.open(store: store, checkpoint: checkpoint, limits: limits)
                let receipt = try arguments.options["--receipt"].map { try VivoFingerprint(hex: $0) }
                let selection = try await reader.selectReplica(index: Int(arguments.options["--replica"]!)!,
                    validation: arguments.switches.contains("--verify-all-payloads") ? .allPayloads : .restart, receipt: receipt)
                if arguments.switches.contains("--require-converged"), !selection.inspection.declaredObservableCriteriaSatisfied {
                    throw VivoChemistryError.invalid("selected prefix does not satisfy the declared consecutive observable criteria")
                }
                try Task.checkCancellation()
                let budget = VivoChemistryBudget(maximumBytes: arguments.options["--maximum-export-bytes"].flatMap(Int.init)
                    ?? VivoChemistryBudget().maximumBytes)
                let publication = try await VivoMolecularSamplingExporter.publish(selection, implementationFingerprint: implementation, budget: budget)
                let workflow = VivoChemistryWorkflow(store: store)
                for name in ["checkpoint", "snapshot", "mapping"] where outputs.destinations[name] != nil {
                    try Task.checkCancellation()
                    guard let port = publication.receipt.outputs.first(where: { $0.name == name }) else {
                        throw VivoChemistryError.invalid("sampling export omitted the requested \(name) output")
                    }
                    let payload = try await workflow.payload(artifact: port.artifact, expectedKind: port.kind)
                    guard publication.receipt.payloadFingerprints[name] == (try VivoCanonicalJSON.fingerprint(payload)) else {
                        throw VivoChemistryError.invalid("sampling export payload differs from its retained identity")
                    }
                    try outputs.write(payload, name: name)
                }
                // A successful final output never precedes optional file writes.
                try outputs.write(VivoCanonicalJSON.encode(publication.receipt), name: "receipt")
            } else {
                let receipt = try VivoKineticsDocumentIO.read(VivoMolecularSamplingExportReceipt.self,
                    from: URL(fileURLWithPath: arguments.positionals[0]),
                    maximumBytes: min(limits.maximumDiagnosticBytes, 512 * 1024 * 1024))
                let verification = try await VivoMolecularSamplingExporter.verify(receipt, store: store,
                    implementationFingerprint: implementation, limits: limits,
                    minimumValidation: arguments.switches.contains("--verify-all-payloads") ? .allPayloads : nil)
                if arguments.switches.contains("--require-converged"), !verification.declaredObservableCriteriaSatisfied {
                    throw VivoChemistryError.invalid("verified prefix does not satisfy the declared consecutive observable criteria")
                }
                try outputs.write(VivoCanonicalJSON.encode(verification), name: "receipt")
            }
            return 0
        } catch is CancellationError {
            try? FileHandle.standardError.write(contentsOf: Data("Sampling export cancelled; no successful final output was published.\n".utf8))
            return 130
        } catch {
            try? FileHandle.standardError.write(contentsOf: Data("Sampling export failed: \(error)\n".utf8))
            return 65
        }
    }

    private struct Arguments {
        let command: String
        var positionals: [String] = []
        var options: [String: String] = [:]
        var switches: Set<String> = []

        init(_ raw: [String]) throws {
            guard let command = raw.first, VivoMolecularSamplingExportCLICommands.handles(command) else {
                throw VivoChemistryError.invalid("unknown sampling export command")
            }
            self.command = command
            let values: Set<String> = ["--store", "--checkpoint", "--reference", "--replica", "--receipt",
                "--read-limits", "--maximum-export-bytes", "--output", "--checkpoint-output", "--snapshot-output", "--mapping-output"]
            let flags: Set<String> = ["--require-converged", "--verify-all-payloads"]
            var cursor = 1
            while cursor < raw.count {
                let argument = raw[cursor]
                if flags.contains(argument) {
                    guard switches.insert(argument).inserted else { throw VivoChemistryError.invalid("duplicate sampling export switch") }
                    cursor += 1
                } else if argument.hasPrefix("--") {
                    guard values.contains(argument), options[argument] == nil,
                          cursor + 1 < raw.count, !raw[cursor + 1].hasPrefix("--") else {
                        throw VivoChemistryError.invalid("unknown, duplicate or incomplete sampling export option")
                    }
                    options[argument] = raw[cursor + 1]; cursor += 2
                } else { positionals.append(argument); cursor += 1 }
            }
            guard let store = options["--store"], !store.isEmpty else { throw VivoChemistryError.invalid("sampling export requires --store") }
            if command == "molecule-sampling-export" {
                guard positionals.isEmpty, let replica = options["--replica"], let index = Int(replica), index >= 0,
                      (options["--checkpoint"] != nil) != (options["--reference"] != nil) else {
                    throw VivoChemistryError.invalid("sampling export requires --replica INDEX and exactly one of --checkpoint HASH or --reference NAME")
                }
                for key in ["--checkpoint", "--receipt"] {
                    if let hash = options[key] { _ = try VivoFingerprint(hex: hash) }
                }
                if let bytes = options["--maximum-export-bytes"] {
                    guard let count = Int(bytes), count > 0 else { throw VivoChemistryError.invalid("--maximum-export-bytes must be a positive integer") }
                }
                if let reference = options["--reference"], reference.isEmpty { throw VivoChemistryError.invalid("empty sampling checkpoint reference") }
            } else {
                let allowed: Set<String> = ["--store", "--read-limits", "--output"]
                guard positionals.count == 1, Set(options.keys).isSubset(of: allowed) else {
                    throw VivoChemistryError.invalid("sampling export verification requires one export JSON file and optional --output/--read-limits")
                }
            }
            guard options.values.allSatisfy({ !$0.isEmpty }) else { throw VivoChemistryError.invalid("empty sampling export option") }
        }
    }

    private struct OutputPlan {
        let destinations: [String: URL]
        private let prepared: [String: VivoKineticsDocumentIO.PreparedOutput]
        init(arguments: Arguments, storeRoot: URL) throws {
            let options = [("receipt", "--output"), ("checkpoint", "--checkpoint-output"),
                           ("snapshot", "--snapshot-output"), ("mapping", "--mapping-output")]
            var destinations: [String: URL] = [:]
            let inputs = (arguments.positionals + [arguments.options["--read-limits"]].compactMap { $0 })
                .map { URL(fileURLWithPath: $0) }
            let root = try VivoWorkflowCLIDocumentPaths.canonicalURL(storeRoot)
            let prefix = root.path == "/" ? "/" : root.path + "/"
            for (name, option) in options {
                guard let path = arguments.options[option] else { continue }
                let target = try VivoWorkflowCLIDocumentPaths.canonicalURL(URL(fileURLWithPath: path))
                guard target != root, !target.path.hasPrefix(prefix) else {
                    throw VivoChemistryError.invalid("sampling export output aliases artifact storage")
                }
                for input in inputs {
                    if try VivoWorkflowCLIDocumentPaths.aliases(target, input) {
                        throw VivoChemistryError.invalid("sampling export output aliases an input")
                    }
                }
                for other in destinations.values {
                    if try VivoWorkflowCLIDocumentPaths.aliases(target, other) {
                        throw VivoChemistryError.invalid("sampling export destinations alias each other")
                    }
                    guard !target.path.hasPrefix(other.path + "/"), !other.path.hasPrefix(target.path + "/") else {
                        throw VivoChemistryError.invalid("sampling export file destinations cannot contain each other")
                    }
                }
                guard !FileManager.default.fileExists(atPath: target.path) else {
                    throw VivoChemistryError.invalid("sampling export output exists; outputs are no-clobber")
                }
                destinations[name] = target
            }
            self.destinations = destinations
            // Complete collective preflight before capturing any directory
            // authority. Keep these handles across asynchronous publication so
            // a later ancestor replacement cannot redirect a local output.
            self.prepared = try destinations.mapValues {
                try VivoKineticsDocumentIO.prepareNoClobberOutput(to: $0)
            }
        }

        func write(_ data: Data, name: String) throws {
            try Task.checkCancellation()
            if let destination = prepared[name] {
                try destination.write(data)
            } else if name == "receipt" {
                try FileHandle.standardOutput.write(contentsOf: data)
                try FileHandle.standardOutput.write(contentsOf: Data("\n".utf8))
            }
        }
    }

    static let help = """
    Select one exact accepted replica at the last complete cross-replica block:
      numivivo molecule-sampling-export --store STORE --checkpoint SHA256 --replica 0 --output export.json
      numivivo molecule-sampling-export --store STORE --reference NAME --replica 1 [--receipt SHA256] \
          [--require-converged] [--verify-all-payloads] [--read-limits limits.json] [--maximum-export-bytes BYTES] [--output export.json] \
          [--checkpoint-output checkpoint.json] [--snapshot-output snapshot.json] [--mapping-output mapping.json]
      numivivo molecule-sampling-export-verify export.json --store STORE [--output verification.json] \
          [--require-converged] [--verify-all-payloads] [--read-limits limits.json]
    Export never advances, thermalizes or minimizes MD. Replica indices are explicit
    and zero-based; a block-zero cursor has no accepted production state to export.
    A reference is resolved once. Source termination is not inferred when no bound
    receipt is supplied. Its preliminary cursor integrity read uses maximumCursorBytes
    and is separate from the reader's cumulative budget. Partial prefixes export (0),
    with their original termination and declared-observable criteria retained.
    --require-converged rejects unmet criteria (65). Cooperative task cancellation
    returns 130; process signal handling follows the existing CLI behavior.
    All file outputs are no-clobber and outside source inputs/artifact storage.
    The export receipt is written last. Optional files are individually published;
    a later failure can leave completed files and immutable intermediate artifacts.
    Snapshot output contains all particles; mapping retains the structure-atom map.
    Workflow envelopes expose the receipt's named typed outputs without renaming
    checkpoint object kinds. Downstream electronic recipes still require explicit
    basis/electron settings; periodic geometry requires an explicit QM/MM policy.
    --maximum-export-bytes controls the separate host publication budget (default
    268435456). It changes allocation admission, never simulation parameters.
    Verification freshly validates the archive and export; --verify-all-payloads
    strengthens payload verification without changing the recorded export scope.
    This checks an accepted prefix and declared scalar criteria, not ensemble
    representativeness, all earlier stopping opportunities, rates or kinetics.

    \(VivoMolecularSamplingCLILimits.help)
    """
}
