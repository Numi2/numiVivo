import Foundation
import NumiVivoKit

struct VivoGenomicEvidenceCLICommands {
    static func handles(_ name: String?) -> Bool {
        ["neoantigen-atlas-request", "neoantigen-splice-job", "neoantigen-evidence-import", "neoantigen-evidence-verify", "neoantigen-evidence-review", "neoantigen-evidence-example", "neoantigen-evidence-help"].contains(name ?? "")
    }
    private struct Completion: Codable { var schema: String; var files: [String: String] }
    private func read(_ path: String, maximum: Int = 128 * 1024) throws -> Data {
        try VivoSingleCellCampaignIO.readDocument(URL(fileURLWithPath: path), maximumBytes: maximum)
    }
    private func load<T: Decodable>(_ type: T.Type, _ path: String, maximum: Int = 128 * 1024) throws -> T {
        try VivoGenomicDocuments.decode(type, from: read(path, maximum: maximum))
    }
    private func printJSON<T: Encodable>(_ value: T) throws {
        FileHandle.standardOutput.write(try VivoCanonicalJSON.encode(value)); FileHandle.standardOutput.write(Data("\n".utf8))
    }
    private func options(_ args: ArraySlice<String>) throws -> [String: String] {
        var values: [String: String] = [:], index = args.startIndex
        while index < args.endIndex {
            let key = args[index]
            guard ["--store", "--output", "--selectors", "--atlas-bundle", "--splice-bundle", "--receipt"].contains(key), values[key] == nil,
                  index + 1 < args.endIndex, !args[index + 1].isEmpty, args[index + 1] != "-", !args[index + 1].hasPrefix("--") else {
                throw VivoGenomicEvidenceError.invalid("Unknown, duplicate or incomplete option.\n" + Self.help)
            }
            values[key] = args[index + 1]; index += 2
        }
        return values
    }
    private func check(_ options: [String: String], required: Set<String>, optional: Set<String> = []) throws {
        guard required.isSubset(of: Set(options.keys)), Set(options.keys).isSubset(of: required.union(optional)) else { throw VivoGenomicEvidenceError.invalid(Self.help) }
    }
    private func destination(_ path: String, store: String) throws -> URL {
        let url = URL(fileURLWithPath: path).standardizedFileURL.resolvingSymlinksInPath()
        let root = URL(fileURLWithPath: store).standardizedFileURL.resolvingSymlinksInPath().path
        guard !FileManager.default.fileExists(atPath: url.path), url.path != root, !url.path.hasPrefix(root == "/" ? "/" : root + "/") else {
            throw VivoGenomicEvidenceError.invalid("Output must be new and outside the artifact store.")
        }
        return url
    }
    private func bundle(_ path: String, atlas: Bool) throws -> [String: Data] {
        let root = URL(fileURLWithPath: path)
        let completion = try load(Completion.self, root.appendingPathComponent("complete.json").path)
        guard completion.schema == "numivivo.org/external-capture-complete/v1", !completion.files.isEmpty,
              Set(completion.files.keys).isSubset(of: VivoGenomicEvidenceRevisionBuilder.inputNames) else {
            throw VivoGenomicEvidenceError.invalid("Incomplete or unsupported evidence bundle. The completion record must be written last.")
        }
        var result: [String: Data] = [:]
        for name in completion.files.keys.sorted() {
            let bytes = try read(root.appendingPathComponent(name).path, maximum: VivoAtlasEvidence.maximumDocumentBytes)
            guard try VivoAtlasEvidence.digest(bytes) == completion.files[name] else { throw VivoGenomicEvidenceError.invalid("Evidence bundle digest mismatch: " + name) }
            if name.hasPrefix("atlas-") == atlas { result[name] = bytes }
        }
        guard !result.isEmpty else { throw VivoGenomicEvidenceError.invalid("The bundle has no inputs for the selected evidence type.") }
        return result
    }
    private func sealed(_ inputs: [String: Data]) throws -> [String: Data] {
        var files = inputs
        files["complete.json"] = try VivoCanonicalJSON.encode(Completion(schema: "numivivo.org/external-capture-complete/v1", files: inputs.mapValues(VivoAtlasEvidence.digest)))
        return files
    }
    private func partial(_ revision: VivoGenomicEvidenceRevision) -> Bool {
        revision.atlasCapture?.outcomes.contains { $0.status == .failed } ?? false
    }
    func run(arguments: [String]) async -> Int32 {
        do {
            guard let command = arguments.first else { throw VivoGenomicEvidenceError.invalid(Self.help) }
            if command == "neoantigen-evidence-help" {
                guard arguments.count == 1 else { throw VivoGenomicEvidenceError.invalid(Self.help) }
                FileHandle.standardOutput.write(Data(Self.help.utf8)); return 0
            }
            let example = command == "neoantigen-evidence-example"
            guard example || (arguments.count >= 2 && !arguments[1].hasPrefix("--") && arguments[1] != "-") else { throw VivoGenomicEvidenceError.invalid(Self.help) }
            let opts = try options(arguments.dropFirst(example ? 1 : 2))
            let implementation = try VivoWorkflowCLIImplementation.fingerprint()
            if example {
                try check(opts, required: ["--store", "--output"])
                let output = try destination(opts["--output"]!, store: opts["--store"]!)
                let store = try VivoArtifactStore(rootURL: URL(fileURLWithPath: opts["--store"]!))
                let demo = try VivoNeoantigenExample.inputs()
                let (parentReceipt, parent) = try await VivoNeoantigenArtifacts.publish(manifestData: VivoCanonicalJSON.encode(demo.manifest), tsv: demo.tsv, implementation: implementation, store: store)
                let inputs = try VivoGenomicEvidenceExample.inputs(binding: VivoGenomicEvidenceArtifacts.binding(for: parent), parentInputs: VivoGenomicEvidenceArtifacts.parentInputs(parent), implementationSHA256: implementation.hex)
                let (receipt, revision, _) = try await VivoGenomicEvidenceArtifacts.publish(parentReceipt: parentReceipt, files: inputs, implementation: implementation, store: store)
                var files = try sealed(inputs)
                files["case.json"] = try VivoCanonicalJSON.encode(demo.manifest); files["all_epitopes.tsv"] = demo.tsv
                files["parent-receipt.json"] = try VivoCanonicalJSON.encode(parentReceipt); files["receipt.json"] = try VivoCanonicalJSON.encode(receipt)
                files["revision.json"] = try VivoCanonicalJSON.encode(revision)
                files["report.html"] = try VivoGenomicEvidenceHTML.render(revision, parentData: VivoCanonicalJSON.encode(parent))
                try VivoOmicsDirectoryExport.write(files, to: output)
                try printJSON(["status": "synthetic-example-only", "open": output.appendingPathComponent("report.html").path]); return 0
            }
            switch command {
            case "neoantigen-atlas-request", "neoantigen-splice-job":
                let atlas = command == "neoantigen-atlas-request"
                try check(opts, required: atlas ? ["--store", "--output", "--selectors"] : ["--store", "--output"], optional: atlas ? ["--splice-bundle"] : [])
                let output = try destination(opts["--output"]!, store: opts["--store"]!)
                let store = try VivoArtifactStore(rootURL: URL(fileURLWithPath: opts["--store"]!), createIfNeeded: false)
                let receipt = try load(VivoNeoantigenReceipt.self, arguments[1])
                let parent = try await VivoNeoantigenArtifacts.verify(receipt, implementation: implementation, store: store)
                if atlas {
                    let selectors = try load(VivoAtlasSelectors.self, opts["--selectors"]!)
                    let spliceFiles = try opts["--splice-bundle"].map { try bundle($0, atlas: false) } ?? [:]
                    let request = try VivoGenomicEvidenceArtifacts.request(parent: parent, selectors: selectors, spliceFiles: spliceFiles, implementation: implementation)
                    try VivoOmicsDirectoryExport.write(sealed(["atlas-request.json": try VivoCanonicalJSON.encode(request)]), to: output)
                    try printJSON(["status": "prepared-not-queried", "variants": String(request.variants.count), "excludedCandidates": String(request.excludedCandidates.count)])
                } else {
                    let template = try VivoGenomicEvidenceArtifacts.spliceJobTemplate(parent: parent)
                    try VivoOmicsDirectoryExport.write(["splice-job.json": template], to: output)
                    try printJSON(["status": "template-only-fill-explicit-inputs-and-resources", "directory": output.path])
                }
                return 0
            case "neoantigen-evidence-import":
                try check(opts, required: ["--store", "--output"], optional: ["--atlas-bundle", "--splice-bundle"])
                guard opts["--atlas-bundle"] != nil || opts["--splice-bundle"] != nil else { throw VivoGenomicEvidenceError.invalid("At least one completed evidence bundle is required.") }
                let output = try destination(opts["--output"]!, store: opts["--store"]!)
                let store = try VivoArtifactStore(rootURL: URL(fileURLWithPath: opts["--store"]!), createIfNeeded: false)
                let parentReceipt = try load(VivoNeoantigenReceipt.self, arguments[1])
                var inputs = try opts["--splice-bundle"].map { try bundle($0, atlas: false) } ?? [:]
                let atlasInputs: [String: Data] = try opts["--atlas-bundle"].map { try bundle($0, atlas: true) } ?? [:]
                for (name, bytes) in atlasInputs { inputs[name] = bytes }
                let (receipt, revision, parent) = try await VivoGenomicEvidenceArtifacts.publish(parentReceipt: parentReceipt, files: inputs, implementation: implementation, store: store)
                var files = try sealed(inputs); files["receipt.json"] = try VivoCanonicalJSON.encode(receipt)
                files["revision.json"] = try VivoCanonicalJSON.encode(revision)
                files["report.html"] = try VivoGenomicEvidenceHTML.render(revision, parentData: VivoCanonicalJSON.encode(parent))
                try VivoOmicsDirectoryExport.write(files, to: output)
                try printJSON(["status": partial(revision) ? "imported-with-failed-atlas-queries" : "research-evidence-imported", "revisionSHA256": receipt.revision.hex, "open": output.appendingPathComponent("report.html").path])
                return partial(revision) ? 2 : 0
            case "neoantigen-evidence-verify":
                try check(opts, required: ["--store"])
                let store = try VivoArtifactStore(rootURL: URL(fileURLWithPath: opts["--store"]!), createIfNeeded: false)
                let receipt = try load(VivoGenomicEvidenceReceipt.self, arguments[1])
                let (revision, _) = try await VivoGenomicEvidenceArtifacts.verify(receipt, implementation: implementation, store: store)
                try printJSON(["integrity": "verified-and-reconstructed", "clinicalQualification": "none", "atlasQueries": partial(revision) ? "contains-failures" : "no-failures-recorded"])
                return partial(revision) ? 2 : 0
            case "neoantigen-evidence-review":
                try check(opts, required: ["--store", "--receipt"])
                let store = try VivoArtifactStore(rootURL: URL(fileURLWithPath: opts["--store"]!), createIfNeeded: false)
                let receipt = try load(VivoGenomicEvidenceReceipt.self, opts["--receipt"]!)
                let review = try load(VivoGenomicEvidenceReview.self, arguments[1], maximum: 1024 * 1024)
                let stored = try await VivoGenomicEvidenceArtifacts.record(review, receipt: receipt, implementation: implementation, store: store)
                try printJSON(["status": "research-review-recorded", "reviewSHA256": stored.fingerprint.hex, "reviewerAuthentication": "self-declared"]); return 0
            default: throw VivoGenomicEvidenceError.invalid(Self.help)
            }
        } catch {
            FileHandle.standardError.write(Data("genomic-evidence: \(error)\n".utf8)); return 1
        }
    }
    static let help = """
    Genomic evidence extensions (public-reference research, not patient care):
      neoantigen-evidence-example --store <directory> --output <new-directory>
      neoantigen-atlas-request <parent-receipt.json> --store <directory> --selectors <selectors.json> --output <new-directory> [--splice-bundle <directory>]
      neoantigen-splice-job <parent-receipt.json> --store <directory> --output <new-directory>
      neoantigen-evidence-import <parent-receipt.json> --store <directory> --output <new-directory> [--atlas-bundle <directory>] [--splice-bundle <directory>]
      neoantigen-evidence-verify <receipt.json> --store <directory>
      neoantigen-evidence-review <review.json> --receipt <receipt.json> --store <directory>
      neoantigen-evidence-help

    Selectors contain requestedScorers and ontologyTerms, using exact captured SDK metadata names.
    External retrieval/execution adapters are in ReferenceAdapters/AlphaGenomeAtlas.
    No network requests are made by these native commands. Bundles require complete.json.
    Use originating executable identities; reimport parent inputs after rebuilding.
    Exit 0: completed; 1: invalid input/I/O failure; 2: recorded Atlas query failures.
    The demonstration intentionally contains a failed query and invented predictions.

    """
}
