import Foundation
import NumiVivoKit

struct VivoNeoantigenCLICommands {
    static func handles(_ name: String?) -> Bool {
        ["neoantigen-example", "neoantigen-import", "neoantigen-verify", "neoantigen-review", "neoantigen-help"].contains(name ?? "")
    }
    private func read(_ path: String, maximum: Int) throws -> Data {
        try VivoSingleCellCampaignIO.readDocument(URL(fileURLWithPath: path), maximumBytes: maximum)
    }
    private func printJSON<T: Encodable>(_ value: T) throws {
        FileHandle.standardOutput.write(try VivoCanonicalJSON.encode(value))
        FileHandle.standardOutput.write(Data("\n".utf8))
    }
    private func destination(_ path: String, outside store: String? = nil) throws -> URL {
        guard path != "-", !path.isEmpty else { throw VivoNeoantigenError.invalid("A new local output directory is required.") }
        let output = URL(fileURLWithPath: path).standardizedFileURL.resolvingSymlinksInPath()
        guard !FileManager.default.fileExists(atPath: output.path) else { throw VivoNeoantigenError.invalid("Output already exists; source data and prior reports are never overwritten.") }
        if let store {
            let root = URL(fileURLWithPath: store).standardizedFileURL.resolvingSymlinksInPath().path
            guard output.path != root, !output.path.hasPrefix(root == "/" ? "/" : root + "/") else {
                throw VivoNeoantigenError.invalid("Output must be outside the artifact store.")
            }
        }
        return output
    }
    func run(arguments: [String]) async -> Int32 {
        do {
            guard let command = arguments.first else { throw VivoNeoantigenError.invalid(Self.help) }
            if command == "neoantigen-help" {
                guard arguments.count == 1 else { throw VivoNeoantigenError.invalid(Self.help) }
                FileHandle.standardOutput.write(Data(Self.help.utf8)); return 0
            }
            let implementation = try VivoWorkflowCLIImplementation.fingerprint()
            switch command {
            case "neoantigen-example":
                guard arguments.count == 3, arguments[1] == "--output" else { throw VivoNeoantigenError.invalid(Self.help) }
                let output = try destination(arguments[2]), example = try VivoNeoantigenExample.inputs()
                let report = try VivoNeoantigenWorkbench.analyze(manifest: example.manifest, tsv: example.tsv, implementationSHA256: implementation.hex)
                try VivoOmicsDirectoryExport.write([
                    "case.json": try VivoCanonicalJSON.encode(example.manifest), "all_epitopes.tsv": example.tsv,
                    "report.json": try VivoCanonicalJSON.encode(report), "report.html": try VivoNeoantigenHTML.render(report)
                ], to: output)
                try printJSON(["status": "synthetic-demonstration-only", "open": output.appendingPathComponent("report.html").path])
                return 0
            case "neoantigen-import":
                guard arguments.count == 8, arguments[2] == "--tsv", arguments[4] == "--store", arguments[6] == "--output" else { throw VivoNeoantigenError.invalid(Self.help) }
                let output = try destination(arguments[7], outside: arguments[5])
                let manifest = try read(arguments[1], maximum: 128 * 1024)
                let tsv = try read(arguments[3], maximum: VivoNeoantigenWorkbench.maximumTSVBytes)
                let store = try VivoArtifactStore(rootURL: URL(fileURLWithPath: arguments[5]))
                let result = try await VivoNeoantigenArtifacts.publish(manifestData: manifest, tsv: tsv, implementation: implementation, store: store)
                try VivoOmicsDirectoryExport.write([
                    "case.json": manifest, "all_epitopes.tsv": tsv,
                    "report.json": try VivoCanonicalJSON.encode(result.report), "report.html": try VivoNeoantigenHTML.render(result.report),
                    "receipt.json": try VivoCanonicalJSON.encode(result.receipt)
                ], to: output)
                try printJSON(["status": result.report.state.rawValue, "reportSHA256": result.receipt.report.hex,
                               "open": output.appendingPathComponent("report.html").path])
                return result.report.state == .blocked ? 2 : 0
            case "neoantigen-verify":
                guard arguments.count == 4, arguments[2] == "--store" else { throw VivoNeoantigenError.invalid(Self.help) }
                let receipt = try VivoCanonicalJSON.decode(VivoNeoantigenReceipt.self, from: read(arguments[1], maximum: 128 * 1024))
                let store = try VivoArtifactStore(rootURL: URL(fileURLWithPath: arguments[3]), createIfNeeded: false)
                let report = try await VivoNeoantigenArtifacts.verify(receipt, implementation: implementation, store: store)
                try printJSON(["integrity": "verified-and-reconstructed", "researchState": report.state.rawValue, "clinicalQualification": "none"])
                return report.state == .blocked ? 2 : 0
            case "neoantigen-review":
                guard arguments.count == 6, arguments[2] == "--receipt", arguments[4] == "--store" else { throw VivoNeoantigenError.invalid(Self.help) }
                let review = try VivoCanonicalJSON.decode(VivoNeoantigenReview.self, from: read(arguments[1], maximum: 1024 * 1024))
                let receipt = try VivoCanonicalJSON.decode(VivoNeoantigenReceipt.self, from: read(arguments[3], maximum: 128 * 1024))
                let store = try VivoArtifactStore(rootURL: URL(fileURLWithPath: arguments[5]), createIfNeeded: false)
                let saved = try await VivoNeoantigenArtifacts.record(review, receipt: receipt, implementation: implementation, store: store)
                try printJSON(["status": "research-review-recorded", "reviewSHA256": saved.fingerprint.hex, "reviewerAuthentication": "self-declared"])
                return 0
            default: throw VivoNeoantigenError.invalid(Self.help)
            }
        } catch {
            FileHandle.standardError.write(Data("neoantigen: \(error)\n".utf8)); return 1
        }
    }
    static let help = """
    Reference-data neoantigen workbench (research only; no patient treatment):
      neoantigen-example --output <new-directory>
      neoantigen-import <case.json> --tsv <all_epitopes.tsv> --store <directory> --output <new-directory>
      neoantigen-verify <receipt.json> --store <directory>
      neoantigen-review <review.json> --receipt <receipt.json> --store <directory>
      neoantigen-help

    Open report.html locally to inspect evidence and export research decisions.
    Imports are class-I, unaggregated pVACseq reports, not raw-read analysis.
    The example is invented; no pVACseq run or biological validation is claimed.
    Exit 0: completed, 1: invalid input/I/O failure, 2: report has blocking findings.
    Verification checks integrity and reconstruction, never clinical suitability.

    """
}
