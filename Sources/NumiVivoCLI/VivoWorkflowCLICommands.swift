import Foundation
import NumiVivoKit
#if canImport(Metal)
import Metal
#endif

struct VivoWorkflowCLICommands {
    static func handles(_ command: String?) -> Bool {
        ["workflow-catalog", "workflow-template", "workflow-plan", "workflow-run", "workflow-import", "workflow-export", "workflow-verify", "workflow-help"].contains(command ?? "")
    }
    private struct Identity: Codable {
        let binary: VivoFingerprint
        let system: String
        let architecture: String
        let metalDevice: String?
        let metalRegistryID: UInt64?
    }
    private struct Receipt: Codable {
        let schema: String
        let reportArtifact: VivoFingerprint
        let recipeFingerprint: VivoFingerprint
        let allTasksSucceeded: Bool
    }
    private struct Options {
        let positional: [String]
        let values: [String: String]
        let force: Bool
        init(_ args: [String]) throws {
            var p: [String] = [], v: [String: String] = [:], force = false, i = 0
            while i < args.count {
                let token = args[i]
                if token == "--force" {
                    guard !force else { throw VivoChemistryError.invalid("duplicate --force") }; force = true; i += 1
                } else if token.hasPrefix("--") {
                    guard ["--output", "--store", "--kind"].contains(token), v[token] == nil,
                          i + 1 < args.count, !args[i + 1].isEmpty, !args[i + 1].hasPrefix("--") else {
                        throw VivoChemistryError.invalid("unknown, duplicated or valueless workflow option")
                    }
                    v[token] = args[i + 1]; i += 2
                } else { p.append(token); i += 1 }
            }
            positional = p; values = v; self.force = force
        }
        func require(_ count: Int, _ allowed: Set<String>) throws {
            guard positional.count == count, Set(values.keys).isSubset(of: allowed), !force || values["--output"] != nil else {
                throw VivoChemistryError.invalid("workflow argument/option mismatch; see workflow-help")
            }
        }
    }
    private func canonicalPath(_ url: URL) -> URL {
        var root = url.standardizedFileURL, suffix: [String] = []
        while !FileManager.default.fileExists(atPath: root.path), root.path != "/" {
            suffix.append(root.lastPathComponent); root.deleteLastPathComponent()
        }
        var result = root.resolvingSymlinksInPath().standardizedFileURL
        for part in suffix.reversed() { result.appendPathComponent(part) }
        return result.standardizedFileURL
    }
    private func checkDestination(_ url: URL, inputs: [URL], store: URL?) throws {
        let target = canonicalPath(url)
        if let store {
            let root = canonicalPath(store), prefix = root.path == "/" ? "/" : root.path + "/"
            guard target != root, !target.path.hasPrefix(prefix) else { throw VivoChemistryError.invalid("workflow export cannot overwrite the artifact store") }
        }
        for input in inputs {
            guard target != canonicalPath(input) else { throw VivoChemistryError.invalid("workflow output/receipt aliases an input") }
            if let a = try? FileManager.default.attributesOfItem(atPath: target.path),
               let b = try? FileManager.default.attributesOfItem(atPath: input.path),
               let ai = a[.systemFileNumber] as? NSNumber, let bi = b[.systemFileNumber] as? NSNumber,
               let ad = a[.systemNumber] as? NSNumber, let bd = b[.systemNumber] as? NSNumber, ai == bi, ad == bd {
                throw VivoChemistryError.invalid("workflow output/receipt is hard-linked to an input")
            }
        }
    }
    private func output(_ data: Data, destination: String?, force: Bool) throws {
        guard let destination, destination != "-" else {
            FileHandle.standardOutput.write(data); FileHandle.standardOutput.write(Data("\n".utf8)); return
        }
        // Reuse the existing no-clobber/pinned-root publication primitive.
        let url = URL(fileURLWithPath: destination)
        let parent = canonicalPath(url.deletingLastPathComponent())
        try VivoKineticsDocumentIO.write(data, to: parent.appendingPathComponent(url.lastPathComponent), overwrite: force)
    }
    private func fingerprint(_ hex: String) throws -> VivoFingerprint {
        let encoded = Array(hex.utf8)
        guard encoded.count == 64 else { throw VivoChemistryError.invalid("a 64-character SHA-256 is required") }
        func nibble(_ byte: UInt8) throws -> UInt8 {
            switch byte {
            case 48...57: return byte - 48
            case 65...70: return byte - 65 + 10
            case 97...102: return byte - 97 + 10
            default: throw VivoChemistryError.invalid("invalid ASCII SHA-256 hex")
            }
        }
        let bytes = try stride(from: 0, to: 64, by: 2).map { i in
            try (nibble(encoded[i]) << 4) | nibble(encoded[i + 1])
        }
        return try .init(bytes: bytes)
    }
    private func registry() throws -> VivoWorkflowRegistry {
        guard let binary = Bundle.main.executableURL else { throw VivoChemistryError.invalid("workflow executable identity unavailable") }
        #if arch(arm64)
        let architecture = "arm64"
        #elseif arch(x86_64)
        let architecture = "x86_64"
        #else
        let architecture = "other"
        #endif
        #if canImport(Metal)
        let device = MTLCreateSystemDefaultDevice(), name = device?.name, registryID = device?.registryID
        #else
        let name: String? = nil, registryID: UInt64? = nil
        #endif
        let identity = try Identity(binary: VivoCanonicalJSON.fingerprint(VivoMDAtomicFileExport.read(binary, maximumBytes: 1 << 30)),
            system: ProcessInfo.processInfo.operatingSystemVersionString, architecture: architecture, metalDevice: name, metalRegistryID: registryID)
        return try VivoPlatformOperations.registry(implementationFingerprint: VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(identity)))
    }
    func run(arguments: [String]) async -> Int32 {
        do {
            guard let command = arguments.first else { throw VivoChemistryError.invalid("missing workflow command") }
            let args = try Options(Array(arguments.dropFirst())), destination = args.values["--output"]
            if command == "workflow-help" {
                try args.require(0, []); FileHandle.standardOutput.write(Data(Self.help.utf8)); return 0
            }
            let usesStore = ["workflow-run", "workflow-import", "workflow-export", "workflow-verify"].contains(command)
            let root = usesStore ? canonicalPath(URL(fileURLWithPath: args.values["--store"] ?? ".numivivo/chemistry-artifacts")) : nil
            let inputFiles = ["workflow-plan", "workflow-run", "workflow-import"].contains(command) ? args.positional.map { URL(fileURLWithPath: $0) } : []
            if let destination, destination != "-" {
                for name in [destination, destination + ".receipt.json"] {
                    let url = URL(fileURLWithPath: name)
                    try checkDestination(url, inputs: inputFiles, store: root)
                    guard args.force || !FileManager.default.fileExists(atPath: url.path) else { throw VivoChemistryError.invalid("workflow output already exists; use --force explicitly") }
                }
            }
            switch command {
            case "workflow-catalog":
                try args.require(0, ["--output"])
                try output(VivoCanonicalJSON.encode(registry().catalog), destination: destination, force: args.force)
            case "workflow-template":
                try args.require(1, ["--output"])
                try output(VivoCanonicalJSON.encode(VivoWorkflowTemplates.make(args.positional[0])), destination: destination, force: args.force)
            case "workflow-plan", "workflow-run":
                try args.require(1, command == "workflow-plan" ? ["--output"] : ["--output", "--store"])
                let data = try VivoMDAtomicFileExport.read(URL(fileURLWithPath: args.positional[0]), maximumBytes: 128 << 20)
                let recipe = try VivoCanonicalJSON.decode(VivoWorkflowRecipe.self, from: data), registry = try registry()
                let plan = try VivoWorkflowPlanner.compile(recipe, registry: registry)
                if command == "workflow-plan" { try output(VivoCanonicalJSON.encode(plan), destination: destination, force: args.force); return 0 }
                let store = try VivoArtifactStore(rootURL: root!)
                let execution = try await VivoWorkflowExecutor(store: store, registry: registry).run(recipe)
                try publish(execution, destination: destination, force: args.force)
                return execution.report.allTasksSucceeded ? 0 : 1
            case "workflow-import":
                try args.require(1, ["--output", "--store", "--kind"])
                guard let kind = args.values["--kind"], !kind.isEmpty else { throw VivoChemistryError.invalid("workflow-import requires --kind") }
                let payload = try VivoCanonicalJSON.decode(VivoJSONValue.self, from: VivoMDAtomicFileExport.read(URL(fileURLWithPath: args.positional[0]), maximumBytes: 128 << 20))
                let store = try VivoArtifactStore(rootURL: root!)
                let artifact = try await store.put(data: VivoCanonicalJSON.encode(payload), kind: kind, mediaType: "application/json")
                try output(VivoCanonicalJSON.encode(artifact), destination: destination, force: args.force)
            case "workflow-export":
                try args.require(1, ["--output", "--store", "--kind"])
                guard let kind = args.values["--kind"] else { throw VivoChemistryError.invalid("workflow-export requires the expected --kind") }
                let store = try VivoArtifactStore(rootURL: root!, createIfNeeded: false), workflow = VivoChemistryWorkflow(store: store)
                let payload = try await workflow.payload(artifact: fingerprint(args.positional[0]), expectedKind: kind)
                try output(payload, destination: destination, force: args.force)
            case "workflow-verify":
                try args.require(1, ["--output", "--store"])
                let store = try VivoArtifactStore(rootURL: root!, createIfNeeded: false)
                let execution = try await VivoWorkflowExecutor(store: store, registry: registry()).verify(fingerprint(args.positional[0]))
                try publish(execution, destination: destination, force: args.force)
                return execution.report.allTasksSucceeded ? 0 : 1
            default: throw VivoChemistryError.invalid("unknown workflow command")
            }
            return 0
        } catch is CancellationError {
            FileHandle.standardError.write(Data("Workflow cancelled. Accepted per-task artifacts remain resumable; no completed run report was published.\n".utf8)); return 130
        } catch {
            FileHandle.standardError.write(Data("Workflow failed: \(error)\n".utf8)); return 65
        }
    }
    private func publish(_ execution: VivoWorkflowExecution, destination: String?, force: Bool) throws {
        try output(VivoCanonicalJSON.encode(execution.report), destination: destination, force: force)
        if let destination, destination != "-" {
            let receipt = Receipt(schema: "numivivo.org/workflow-cli-receipt/v1", reportArtifact: execution.artifact.fingerprint,
                recipeFingerprint: execution.report.recipeFingerprint, allTasksSucceeded: execution.report.allTasksSucceeded)
            try output(VivoCanonicalJSON.encode(receipt), destination: destination + ".receipt.json", force: force)
        }
        let reused = execution.report.nodes.filter { $0.outcome.reused }.count
        FileHandle.standardError.write(Data("Workflow report \(execution.artifact.fingerprint.hex): \(execution.report.nodes.count) nodes, \(reused) reused, allTasksSucceeded=\(execution.report.allTasksSucceeded).\n".utf8))
    }
    static let help = """
    General native workflow recipes
      numivivo workflow-catalog
      numivivo workflow-template molecular-analysis --output recipe.json
      numivivo workflow-template md-segments --output md.json
      numivivo workflow-plan recipe.json
      numivivo workflow-run recipe.json --store ./artifacts --output report.json
      numivivo workflow-import input.json --kind vivo.electronic-system --store ./artifacts
      numivivo workflow-export SHA256 --kind vivo.many-body-result --store ./artifacts --output result.json
      numivivo workflow-verify REPORT_SHA256 --store ./artifacts --output verified.json

    The JSON recipe chooses registered operations, typed inputs, dependencies,
    numerical settings, resource limits and named outputs. No shell commands or
    automatic method fallback. Independent branches use bounded concurrency;
    failed branches remain visible and their descendants are marked blocked.
    Rerun the same recipe/store to resume from verified accepted task outputs.
    A run with failures returns exit 1 and retains its partial report. Invalid
    recipes return 65 before execution. Exports never overwrite inputs or store
    objects, and existing files require --force. Plan does not execute kernels.
    Cancellation is cooperative at task boundaries; already accepted task outputs
    may remain cached. MD replay validation is explicitly distinct from integrity.
    This interface does not upgrade numerical results into qualified rates or
    biological predictions. The existing stage engines retain their own limits.

    """
}
