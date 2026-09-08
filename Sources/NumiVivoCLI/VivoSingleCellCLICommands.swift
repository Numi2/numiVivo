import Foundation
import NumiVivoKit

struct VivoSingleCellCLICommands {
    static func handles(_ name: String?) -> Bool {
        ["singlecell-run", "singlecell-verify", "singlecell-export", "singlecell-help"].contains(name ?? "")
    }
    private struct Verification: Encodable {
        let schemaVersion: Int
        let status: String
        let evidence: VivoOmicsEvidence
        let cells: Int
        let features: Int
        let pseudobulkGroups: Int
        let numericalProfile: String
    }
    private func canonicalURL(_ url: URL) throws -> URL {
        let manager = FileManager.default
        var ancestor = url.absoluteURL, suffix: [String] = []
        guard ancestor.isFileURL, ancestor.path.utf8.count <= 8192 else {
            throw VivoOmicsError.invalid("bounded local output path required")
        }
        while !manager.fileExists(atPath: ancestor.path) {
            if let attributes = try? manager.attributesOfItem(atPath: ancestor.path),
               attributes[.type] as? FileAttributeType == .typeSymbolicLink {
                throw VivoOmicsError.invalid("dangling symbolic link in output path")
            }
            let parent = ancestor.deletingLastPathComponent()
            guard parent.path != ancestor.path, suffix.count < 4096 else {
                throw VivoOmicsError.invalid("output has no resolvable ancestor")
            }
            suffix.append(ancestor.lastPathComponent); ancestor = parent
        }
        var resolved = ancestor.resolvingSymlinksInPath().standardizedFileURL
        for part in suffix.reversed() { resolved.appendPathComponent(part) }
        return resolved.standardizedFileURL
    }
    func run(arguments: [String]) async -> Int32 {
        do {
            guard let command = arguments.first else { throw VivoOmicsError.invalid("missing command") }
            if command == "singlecell-help" {
                guard arguments.count == 1 else { throw VivoOmicsError.invalid("singlecell-help takes no arguments") }
                FileHandle.standardOutput.write(Data(Self.help.utf8)); return 0
            }
            guard arguments.count >= 4, !arguments[1].hasPrefix("--"), arguments[1] != "-" else {
                throw VivoOmicsError.invalid("one manifest/receipt file and --store <directory> are required")
            }
            let inputURL = URL(fileURLWithPath: arguments[1]).standardizedFileURL
            var options: [String: String] = [:], index = 2
            while index < arguments.count {
                let key = arguments[index]
                guard ["--store", "--output"].contains(key), options[key] == nil,
                      index + 1 < arguments.count, !arguments[index + 1].isEmpty,
                      !arguments[index + 1].hasPrefix("--") else { throw VivoOmicsError.invalid("unknown, duplicate or incomplete option") }
                options[key] = arguments[index + 1]; index += 2
            }
            guard let storePath = options["--store"], storePath != "-" else { throw VivoOmicsError.invalid("--store <directory> is required") }
            let storeURL = URL(fileURLWithPath: storePath).standardizedFileURL
            if let output = options["--output"], output != "-" {
                let url = try canonicalURL(URL(fileURLWithPath: output))
                let root = try canonicalURL(storeURL).path
                guard url != inputURL.resolvingSymlinksInPath(), url.path != root,
                      !url.path.hasPrefix(root == "/" ? "/" : root + "/"), !FileManager.default.fileExists(atPath: url.path) else {
                    throw VivoOmicsError.invalid("output must be new, outside the artifact store, and not an input")
                }
            }
            let implementation = try VivoWorkflowCLIImplementation.fingerprint()
            let output: Data
            switch command {
            case "singlecell-run":
                let input = try VivoSingleCellCampaignIO.snapshot(manifestURL: inputURL)
                let store = try VivoArtifactStore(rootURL: storeURL)
                let receipt = try await VivoSingleCellArtifacts.publish(input: input, implementation: implementation, store: store)
                output = try VivoCanonicalJSON.encode(receipt)
            case "singlecell-verify", "singlecell-export":
                let receipt = try VivoKineticsDocumentIO.read(VivoSingleCellRunReceipt.self, from: inputURL)
                let store = try VivoArtifactStore(rootURL: storeURL, createIfNeeded: false)
                let report = try await VivoSingleCellArtifacts.verify(receipt: receipt, implementation: implementation, store: store)
                if command == "singlecell-export" { output = try VivoCanonicalJSON.encode(report) }
                else { output = try VivoCanonicalJSON.encode(Verification(schemaVersion: 1,
                    status: "verified-native-reconstruction-not-biological-validation", evidence: report.dataset.evidence,
                    cells: report.dataset.cells.count, features: report.dataset.features.count,
                    pseudobulkGroups: report.pseudobulk.groups.count, numericalProfile: report.numericalProfile)) }
            default: throw VivoOmicsError.invalid("unsupported single-cell command")
            }
            try Task.checkCancellation()
            if let path = options["--output"], path != "-" {
                try VivoKineticsDocumentIO.write(output, to: URL(fileURLWithPath: path), overwrite: false)
            } else {
                FileHandle.standardOutput.write(output); FileHandle.standardOutput.write(Data("\n".utf8))
            }
            return 0
        } catch is CancellationError {
            FileHandle.standardError.write(Data("numivivo singlecell: cancelled; no success receipt published\n".utf8)); return 130
        } catch {
            FileHandle.standardError.write(Data("numivivo singlecell: \(error)\n".utf8)); return 65
        }
    }
    static let help = """
    NumiVivo single-cell count workflows
      singlecell-run <manifest.json> --store <artifact-directory> [--output <new-receipt.json|->]
      singlecell-verify <receipt.json> --store <artifact-directory> [--output <new-summary.json|->]
      singlecell-export <receipt.json> --store <artifact-directory> [--output <new-report.json|->]
    Source paths are relative to the manifest directory; symlinks and ../ are rejected.
    Inputs: uncompressed 10x integer/general Matrix Market, three-column Gene Expression
    features.tsv and barcodes.tsv. Use explicit sample, replicate, evidence and count units.
    Raw counts remain exact UInt64. Optional log normalization is a separate FP64 view.
    No filtering, inferred cell types, batch correction or differential-expression test is performed.
    Verification needs the recorded executable/OS and reconstructs results from immutable source bytes.
    This is bounded native CPU processing, not a GPU/atlas-scale or biological-validity claim.
    """ + "\n"
}
