import Foundation
import CryptoKit
import NumiVivoKit

struct VivoWorkflowCLICommands {
    static func handles(_ command: String?) -> Bool {
        ["workflow-help", "workflow-catalog", "workflow-template", "workflow-plan", "workflow-run", "workflow-verify", "workflow-export", "artifact-put", "artifact-show", "campaign-template", "campaign-plan", "campaign-run", "campaign-resume"].contains(command ?? "")
    }
    private struct Implementation: Codable { let schema: String; let executableSHA256: String; let operatingSystem: String; let executionSemantics: String }
    private static func implementation() throws -> VivoFingerprint {
        let invocation = CommandLine.arguments[0]
        let executable: URL
        if invocation.contains("/") { executable = URL(fileURLWithPath: invocation).standardizedFileURL.resolvingSymlinksInPath() }
        else {
            guard let path = (ProcessInfo.processInfo.environment["PATH"] ?? "").split(separator: ":").map(String.init)
                .map({ URL(fileURLWithPath: $0.isEmpty ? "." : $0).appendingPathComponent(invocation) })
                .first(where: { FileManager.default.isExecutableFile(atPath: $0.path) }) else {
                throw VivoChemistryError.invalid("cannot resolve executing workflow binary")
            }
            executable = path.standardizedFileURL.resolvingSymlinksInPath()
        }
        let handle = try FileHandle(forReadingFrom: executable); defer { try? handle.close() }
        var digest = SHA256()
        while let bytes = try handle.read(upToCount: 4*1024*1024), !bytes.isEmpty { digest.update(data: bytes) }
        let hash = digest.finalize().map { String(format: "%02x", $0) }.joined()
        return try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(Implementation(schema: "numivivo.org/platform-executable/v1",
            executableSHA256: hash, operatingSystem: ProcessInfo.processInfo.operatingSystemVersionString,
            executionSemantics: "static-native-operations;fp64-electronics;metal-transactional-md;single-owner-artifact-scheduler")))
    }
    private struct Arguments {
        let command: String
        var positionals: [String] = []
        var options: [String: String] = [:]
        var switches: Set<String> = []
        init(_ raw: [String]) throws {
            guard let command = raw.first else { throw VivoChemistryError.invalid("missing workflow command") }
            self.command = command; var cursor = 1
            let flags: Set<String> = ["--store", "--output", "--name", "--kind", "--media-type"]
            while cursor < raw.count {
                let argument = raw[cursor]
                if argument == "--force" { guard switches.insert(argument).inserted else { throw VivoChemistryError.invalid("duplicate --force") }; cursor += 1 }
                else if argument.hasPrefix("--") {
                    guard flags.contains(argument), options[argument] == nil, cursor+1 < raw.count, !raw[cursor+1].hasPrefix("--") else {
                        throw VivoChemistryError.invalid("unknown, duplicate or incomplete workflow option \(argument)")
                    }
                    options[argument] = raw[cursor+1]; cursor += 2
                } else { positionals.append(argument); cursor += 1 }
            }
        }
        func require(_ count: Int, allowed: Set<String>) throws {
            guard positionals.count == count, Set(options.keys).isSubset(of: allowed),
                  switches.isEmpty || (allowed.contains("--force") && switches == ["--force"]) else {
                throw VivoChemistryError.invalid("invalid arguments for \(command); use workflow-help")
            }
        }
    }
    private static func file(_ path: String) throws -> Data {
        let url = URL(fileURLWithPath: path)
        let values = try url.resourceValues(forKeys: [.isRegularFileKey, .fileSizeKey])
        guard values.isRegularFile == true, let size = values.fileSize, size >= 0, size <= 256*1024*1024 else {
            throw VivoChemistryError.invalid("workflow input must be a bounded regular file")
        }
        let handle = try FileHandle(forReadingFrom: url); defer { try? handle.close() }
        let data = try handle.read(upToCount: 256*1024*1024+1) ?? Data()
        guard data.count <= 256*1024*1024 else { throw VivoChemistryError.resourceLimit("workflow file grew beyond its read limit") }
        return data
    }
    private static func output<T: Encodable>(_ value: T, arguments: Arguments, source: String? = nil) throws {
        if let path = arguments.options["--output"] {
            let destination = URL(fileURLWithPath: path).standardizedFileURL.resolvingSymlinksInPath()
            if let source, destination == URL(fileURLWithPath: source).standardizedFileURL.resolvingSymlinksInPath() {
                throw VivoChemistryError.invalid("refusing to overwrite workflow input")
            }
            if FileManager.default.fileExists(atPath: destination.path), !arguments.switches.contains("--force") {
                throw VivoChemistryError.invalid("output exists; use --force explicitly")
            }
            try VivoCanonicalJSON.encode(value).write(to: destination, options: .atomic)
        } else { print(String(decoding: try VivoCanonicalJSON.encode(value), as: UTF8.self)) }
    }
    private static func publish(_ execution: VivoWorkflowExecution, arguments: Arguments) throws {
        try output(execution.report, arguments: arguments)
        if let path = arguments.options["--output"] {
            let receipt = URL(fileURLWithPath: path).appendingPathExtension("receipt.json")
            if FileManager.default.fileExists(atPath: receipt.path), !arguments.switches.contains("--force") {
                throw VivoChemistryError.invalid("receipt output exists; use --force explicitly")
            }
            struct Receipt: Codable { let schema: String; let reportArtifact: VivoFingerprint; let recipeFingerprint: VivoFingerprint; let allTasksSucceeded: Bool }
            try VivoCanonicalJSON.encode(Receipt(schema: "numivivo.org/workflow-export-receipt/v1", reportArtifact: execution.artifact.fingerprint,
                recipeFingerprint: execution.report.recipeFingerprint, allTasksSucceeded: execution.report.allTasksSucceeded)).write(to: receipt, options: .atomic)
        }
        if !execution.report.allTasksSucceeded { throw WorkflowFailure.tasks }
    }
    private static func publishCampaign(_ execution: VivoAdaptiveCampaignExecution,arguments: Arguments) throws {
        try output(execution.report,arguments: arguments)
        if let path = arguments.options["--output"] {
            let receipt = URL(fileURLWithPath: path).appendingPathExtension("receipt.json")
            if FileManager.default.fileExists(atPath: receipt.path), !arguments.switches.contains("--force") {
                throw VivoChemistryError.invalid("campaign receipt exists; use --force explicitly")
            }
            struct Receipt: Codable { let schema: String; let reportArtifact: VivoFingerprint; let checkpointArtifact: VivoFingerprint; let criteriaSatisfied: Bool }
            try VivoCanonicalJSON.encode(Receipt(schema: "numivivo.org/campaign-export-receipt/v1",reportArtifact: execution.artifact.fingerprint,
                checkpointArtifact: execution.report.checkpointArtifact,criteriaSatisfied: execution.report.finalDecision.allCriteriaSatisfied)).write(to: receipt,options: .atomic)
        }
        if execution.report.termination != .criteriaSatisfied { throw WorkflowFailure.tasks }
    }
    private enum WorkflowFailure: Error { case tasks }
    private static func diagnostics(_ error: Error) {
        try? FileHandle.standardError.write(contentsOf: Data("\(error)\n".utf8))
    }
    static func run(_ raw: [String]) async {
        do { try await execute(Arguments(raw)) }
        catch WorkflowFailure.tasks { diagnostics("workflow or campaign is incomplete; report and checkpoint retain failures and unexecuted work"); VivoNativeExecution.exit(2) }
        catch { diagnostics(error); VivoNativeExecution.exit(1) }
    }
    private static func execute(_ arguments: Arguments) async throws {
        if arguments.command == "workflow-help" { try arguments.require(0, allowed: []); print(help); return }
        let id = try implementation(), registry = try VivoPlatformOperations.registry(implementationFingerprint: id)
        let rootStore = URL(fileURLWithPath: arguments.options["--store"] ?? ".numivivo/workflow-artifacts").standardizedFileURL.resolvingSymlinksInPath()
        if ["workflow-run", "workflow-verify", "campaign-run", "campaign-resume"].contains(arguments.command), let path = arguments.options["--output"] {
            let target = URL(fileURLWithPath: path).standardizedFileURL.resolvingSymlinksInPath()
            let receipt = target.appendingPathExtension("receipt.json")
            if !arguments.switches.contains("--force"), FileManager.default.fileExists(atPath: receipt.path) {
                throw VivoChemistryError.invalid("workflow receipt exists; use --force explicitly")
            }
        }
        if let path = arguments.options["--output"] {
            let target = URL(fileURLWithPath: path).standardizedFileURL.resolvingSymlinksInPath()
            if arguments.options["--store"] != nil || ["workflow-run","workflow-verify","workflow-export","artifact-put","artifact-show","campaign-run","campaign-resume"].contains(arguments.command) {
                guard target.path != rootStore.path, !target.path.hasPrefix(rootStore.path + "/") else {
                    throw VivoChemistryError.invalid("workflow exports must be outside the artifact store")
                }
            }
            if ["workflow-plan", "workflow-run", "artifact-put", "campaign-plan", "campaign-run"].contains(arguments.command), let input = arguments.positionals.first {
                guard target != URL(fileURLWithPath: input).standardizedFileURL.resolvingSymlinksInPath() else {
                    throw VivoChemistryError.invalid("refusing to overwrite workflow source input")
                }
            }
            if FileManager.default.fileExists(atPath: target.path), !arguments.switches.contains("--force") {
                throw VivoChemistryError.invalid("output exists; use --force explicitly")
            }
        }
        switch arguments.command {
        case "campaign-template":
            try arguments.require(1,allowed: ["--output","--force"])
            try output(VivoAdaptiveCampaignExamples.make(arguments.positionals[0]),arguments: arguments)
        case "campaign-plan", "campaign-run":
            try arguments.require(1,allowed: arguments.command == "campaign-run" ? ["--store","--output","--force"] : ["--output","--force"])
            let request = try VivoCanonicalJSON.decode(VivoAdaptiveCampaignRequest.self,from: file(arguments.positionals[0]))
            if arguments.command == "campaign-plan" {
                try output(VivoAdaptiveCampaignPlanner.compile(request,registry: registry),arguments: arguments,source: arguments.positionals[0])
            } else {
                let store = try VivoArtifactStore(rootURL: rootStore)
                let executor = VivoAdaptiveCampaignExecutor(store: store,registry: registry)
                try await publishCampaign(executor.run(request),arguments: arguments)
            }
        case "campaign-resume":
            try arguments.require(1,allowed: ["--store","--output","--force"])
            let store = try VivoArtifactStore(rootURL: rootStore)
            let executor = VivoAdaptiveCampaignExecutor(store: store,registry: registry)
            try await publishCampaign(executor.resume(VivoFingerprint(hex: arguments.positionals[0])),arguments: arguments)
        case "workflow-catalog":
            try arguments.require(0, allowed: ["--output", "--force"]); try output(registry.catalog, arguments: arguments)
        case "workflow-template":
            try arguments.require(1, allowed: ["--output", "--force"])
            try output(VivoWorkflowTemplates.make(arguments.positionals[0]), arguments: arguments)
        case "workflow-plan":
            try arguments.require(1, allowed: ["--output", "--force"])
            let recipe = try VivoCanonicalJSON.decode(VivoWorkflowRecipe.self, from: file(arguments.positionals[0]))
            try output(VivoWorkflowPlanner.compile(recipe, registry: registry), arguments: arguments, source: arguments.positionals[0])
        case "workflow-run":
            try arguments.require(1, allowed: ["--store", "--output", "--force"])
            let recipe = try VivoCanonicalJSON.decode(VivoWorkflowRecipe.self, from: file(arguments.positionals[0]))
            let store = try VivoArtifactStore(rootURL: rootStore)
            let executor = VivoWorkflowExecutor(store: store, registry: registry)
            try await publish(executor.run(recipe), arguments: arguments)
        case "workflow-verify":
            try arguments.require(1, allowed: ["--store", "--output", "--force"])
            let store = try VivoArtifactStore(rootURL: rootStore)
            let executor = VivoWorkflowExecutor(store: store, registry: registry)
            try await publish(executor.verify(VivoFingerprint(hex: arguments.positionals[0])), arguments: arguments)
        case "workflow-export":
            try arguments.require(1, allowed: ["--store", "--output", "--name", "--force"])
            guard let name = arguments.options["--name"], let path = arguments.options["--output"] else {
                throw VivoChemistryError.invalid("workflow-export requires --name and --output")
            }
            let store = try VivoArtifactStore(rootURL: rootStore)
            let runID = try VivoFingerprint(hex: arguments.positionals[0])
            let executor = VivoWorkflowExecutor(store: store, registry: registry)
            let verified = try await executor.verify(runID)
            guard let port = verified.report.exports.first(where: { $0.name == name })?.artifact else {
                throw VivoChemistryError.invalid("requested workflow export is unavailable")
            }
            let reader = VivoChemistryWorkflow(store: store)
            let bytes = try await reader.payload(artifact: port.artifact, expectedKind: port.kind)
            let destination = URL(fileURLWithPath: path)
            if FileManager.default.fileExists(atPath: destination.path), !arguments.switches.contains("--force") {
                throw VivoChemistryError.invalid("export output exists; use --force explicitly")
            }
            let receipt = destination.appendingPathExtension("receipt.json")
            if FileManager.default.fileExists(atPath: receipt.path), !arguments.switches.contains("--force") { throw VivoChemistryError.invalid("export receipt exists") }
            try bytes.write(to: destination, options: .atomic)
            struct Receipt: Codable { let schema: String; let verifiedRun: VivoFingerprint; let output: VivoChemistryOutputReceipt }
            try VivoCanonicalJSON.encode(Receipt(schema: "numivivo.org/workflow-payload-receipt/v1", verifiedRun: runID, output: port)).write(to: receipt, options: .atomic)
        case "artifact-put":
            try arguments.require(1, allowed: ["--store", "--output", "--kind", "--media-type", "--force"])
            guard let kind = arguments.options["--kind"], !kind.isEmpty else { throw VivoChemistryError.invalid("artifact-put requires --kind") }
            let store = try VivoArtifactStore(rootURL: rootStore)
            let descriptor = try await store.put(data: file(arguments.positionals[0]), kind: kind,
                mediaType: arguments.options["--media-type"] ?? "application/octet-stream")
            try output(descriptor, arguments: arguments, source: arguments.positionals[0])
        case "artifact-show":
            try arguments.require(1, allowed: ["--store", "--output", "--force"])
            let store = try VivoArtifactStore(rootURL: rootStore)
            let fingerprint = try VivoFingerprint(hex: arguments.positionals[0])
            let descriptor = try await store.descriptor(for: fingerprint)
            _ = try await store.data(for: fingerprint, verify: true)
            try output(descriptor, arguments: arguments)
        default: throw VivoChemistryError.invalid("unknown workflow command")
        }
    }
    static let help = """
    NumiVivo reusable native workflow DAG

    workflow-catalog [--output catalog.json]
    workflow-template molecular-analysis|md-segments|native-electronic|qmmm-force|qmmm-minimize|qmmm-pmf|qmmm-surface|qmmm-transmission|qmmm-transmission-apply|qmmm-pmf-rate|qmmm-replicated-rate|qmmm-rate-apply|qmmm-pmf-rate-apply|property-refinement|kinetics [--output recipe.json]
    workflow-plan recipe.json [--output plan.json]
    workflow-run recipe.json [--store DIRECTORY] [--output report.json]
    workflow-verify REPORT_ARTIFACT_SHA256 [--store DIRECTORY] [--output report.json]
    workflow-export REPORT_ARTIFACT_SHA256 --name EXPORT --output payload.json [--store DIRECTORY]
    artifact-put FILE --kind KIND [--media-type TYPE] [--store DIRECTORY] [--output descriptor.json]
    artifact-show ARTIFACT_SHA256 [--store DIRECTORY] [--output descriptor.json]
    campaign-template equilibrium-correction|electronic-crosscheck [--output campaign.json]
    campaign-plan campaign.json [--output campaign.plan.json]
    campaign-run campaign.json [--store DIRECTORY] [--output campaign.report.json]
    campaign-resume CHECKPOINT_ARTIFACT_SHA256 [--store DIRECTORY] [--output resumed.report.json]

    Existing output files require explicit --force. Failed or blocked nodes remain
    in the run report; unrelated branches can finish. Reuse validates artifacts and
    existing native numerical contracts. Campaigns select complete predeclared
    recipes, retain native observable-bound metrics, and stop on failed confirmation.
    Acceptance covers declared criteria, not unexamined chemistry or biological validity.
    Availability does not imply chemical, kinetic, biological, or Apple performance
    qualification. Templates are small examples, not a specific paper's inputs.
    """
}
