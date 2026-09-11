// Baseline runner linked against the frozen pre-refactor native library.
// Calls the actual existing statistical owner; it is not a replacement model.
import Foundation
@testable import NumiVivoKit

@main struct NativeReference {
    static func main() throws {
        let arguments = Array(CommandLine.arguments.dropFirst())
        guard arguments.count == 3 else { fatalError("NativeReference <resident-count-bundle> <contrast.json> <new-report.json>") }
        let source = URL(fileURLWithPath: arguments[0]), output = URL(fileURLWithPath: arguments[2]), toStdout = arguments[2] == "-"
        guard toStdout || !FileManager.default.fileExists(atPath: output.path) else { throw VivoOmicsError.invalid("baseline output exists") }
        let snapshot = FileManager.default.temporaryDirectory.appendingPathComponent("numivivo-expression-baseline-" + UUID().uuidString)
        try FileManager.default.createDirectory(at: snapshot, withIntermediateDirectories: false)
        defer { try? FileManager.default.removeItem(at: snapshot) }
        func copy(_ name: String, _ limit: Int) throws -> VivoFingerprint {
            try VivoOmicsFileSnapshot.fingerprint(source.appendingPathComponent(name), copyTo: snapshot.appendingPathComponent(name), maximumBytes: limit)
        }
        _ = try copy("receipt.json", 65_536)
        let receipt = try VivoCanonicalJSON.decode(VivoCountStreamReceipt.self, from: Data(contentsOf: snapshot.appendingPathComponent("receipt.json")))
        guard try copy("plan.json", 268_435_456) == receipt.plan, try copy("report.json", 536_870_912) == receipt.report else {
            throw VivoOmicsError.invalid("baseline count artifact fingerprint")
        }
        let plan = try VivoCanonicalJSON.decode(VivoCountStreamPlan.self, from: Data(contentsOf: snapshot.appendingPathComponent("plan.json")))
        let report = try VivoCanonicalJSON.decode(VivoCountStreamReport.self, from: Data(contentsOf: snapshot.appendingPathComponent("report.json")))
        guard try plan.validate() == report.canonicalNonzeros else { throw VivoOmicsError.invalid("baseline source axes") }
        let contrast = try VivoCanonicalJSON.decode(VivoOmicsExpressionContrast.self, from: Data(contentsOf: URL(fileURLWithPath: arguments[1])))
        let result = try VivoPseudobulkDifferentialExpression.evaluate(metadata: plan.metadata, bulk: report.pseudobulk, contrast: contrast)
        let raw = try VivoCanonicalJSON.encode(result)
        if toStdout { FileHandle.standardOutput.write(raw) }
        else { try raw.write(to: output, options: .withoutOverwriting) }
        FileHandle.standardError.write(Data("completed-native-resident-expression \(result.features.count) \(result.testedFeatures)\n".utf8))
    }
}
