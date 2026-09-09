import Foundation

@main struct NeoantigenPortableCheck {
    static func main() throws {
        let implementation = try VivoCanonicalJSON.fingerprint(Data(contentsOf: URL(fileURLWithPath: CommandLine.arguments[0]))).hex
        let fixture = try VivoNeoantigenExample.inputs()
        var checks = 0
        func check(_ condition: Bool, _ name: String) throws {
            guard condition else { throw VivoNeoantigenError.invalid("FAIL: " + name) }; checks += 1
        }
        func expectFailure(_ name: String, _ operation: () throws -> Void) throws {
            do { try operation() } catch { checks += 1; return }
            throw VivoNeoantigenError.invalid("FAIL: expected rejection: " + name)
        }
        func analyze(_ manifest: VivoNeoantigenCase, _ tsv: Data) throws -> VivoNeoantigenReport {
            try VivoNeoantigenWorkbench.analyze(manifest: manifest, tsv: tsv, implementationSHA256: implementation)
        }
        func revised(_ column: String, _ value: String, row: Int = 1) throws -> (VivoNeoantigenCase, Data) {
            var lines = String(decoding: fixture.tsv, as: UTF8.self).split(separator: "\n").map { $0.components(separatedBy: "\t") }
            guard let index = lines[0].firstIndex(of: column) else { throw VivoNeoantigenError.invalid("test column missing") }
            lines[row][index] = value
            let bytes = Data((lines.map { $0.joined(separator: "\t") }.joined(separator: "\n") + "\n").utf8)
            var manifest = fixture.manifest; manifest.externalRun.reportSHA256 = try VivoNeoantigenWorkbench.digest(bytes)
            return (manifest, bytes)
        }
        let report = try analyze(fixture.manifest, fixture.tsv)
        try check(report.state == .researchReview && report.candidates.count == 3, "reference-fixture import")
        try check(report.candidates[0].sourceLine == 2 && report.candidates[0].start == 100, "source row and zero-based coordinates")
        try check(report.candidates[1].geneExpression == nil && report.candidates[1].tumorRNAFraction == nil, "missing is not zero")
        try check(report.candidates[2].geneExpression == 0 && report.candidates[2].tumorRNAFraction == 0, "zero retained")
        try check(report.candidates[2].evidenceGaps.contains { $0.contains("unchanged from wildtype") }, "unchanged peptide flagged")
        try check(report.candidates[0].sourceFields["Evaluation"] == "Accept" && report.candidates[0].evidenceGaps.contains { $0.contains("not been experimentally") }, "upstream Accept is not experimental evidence")
        try check(report == analyze(fixture.manifest, fixture.tsv), "deterministic reconstruction")
        try check(Set(report.candidates.map(\.id)).count == 3, "distinct row identities")
        var altered = fixture.manifest; altered.samples[0].assembly = "GRCh37"
        try check(analyze(altered, fixture.tsv).state == .blocked, "mixed assembly blocked")
        altered = fixture.manifest; altered.samples[0].subjectID = "other-subject"
        try check(analyze(altered, fixture.tsv).state == .blocked, "declared subject mismatch blocked")
        altered = fixture.manifest; altered.samples.removeLast()
        try check(analyze(altered, fixture.tsv).state == .blocked, "missing RNA record blocked")
        altered = fixture.manifest; altered.externalRun.version = "latest"
        try check(analyze(altered, fixture.tsv).state == .blocked, "unpinned external version blocked")
        altered = fixture.manifest; altered.externalRun.sampleInputs["demo-tumorDNA"] = String(repeating: "0", count: 64)
        try check(analyze(altered, fixture.tsv).state == .blocked, "input mapping mismatch blocked")
        try check(analyze(fixture.manifest, fixture.tsv + Data("\n".utf8)).state == .blocked, "TSV byte tampering blocked")
        for (column, bad) in [("Median MT IC50 Score", "NaN"), ("Median MT IC50 Score", "Infinity"), ("Median MT IC50 Score", "-1"), ("Median MT IC50 Score", "0"), ("Tumor DNA VAF", "1.1"), ("Tumor RNA Depth", "0.5"), ("Tumor RNA Depth", "0"), ("HLA Allele", "HLA-B*07:02"), ("Peptide Length", "8"), ("Start", "-1"), ("Stop", "99"), ("MT Epitope Seq", "ACDXFGHIK")] {
            let (m, t) = try revised(column, bad)
            try expectFailure(column + "=" + bad) { _ = try analyze(m, t) }
        }
        let header = fixture.tsv.split(separator: 10)[0]
        let empty = Data(header) + Data([10]); altered = fixture.manifest
        altered.externalRun.reportSHA256 = try VivoNeoantigenWorkbench.digest(empty)
        try check(analyze(altered, empty).state == .noCandidates, "header-only report is explicitly empty")
        let duplicateHeader = Data(("Chromosome\t" + String(decoding: fixture.tsv, as: UTF8.self)).utf8)
        altered = fixture.manifest; altered.externalRun.reportSHA256 = try VivoNeoantigenWorkbench.digest(duplicateHeader)
        try expectFailure("duplicate header") { _ = try analyze(altered, duplicateHeader) }
        let duplicateRows = fixture.tsv + Data(fixture.tsv.split(separator: 10)[1]) + Data([10])
        altered = fixture.manifest; altered.externalRun.reportSHA256 = try VivoNeoantigenWorkbench.digest(duplicateRows)
        try expectFailure("duplicate row") { _ = try analyze(altered, duplicateRows) }
        let crlf = Data(String(decoding: fixture.tsv, as: UTF8.self).replacingOccurrences(of: "\n", with: "\r\n").utf8)
        altered = fixture.manifest; altered.externalRun.reportSHA256 = try VivoNeoantigenWorkbench.digest(crlf)
        try check(analyze(altered, crlf).candidates.count == 3, "CRLF supported")
        let invalidUTF8 = Data([255, 254]); altered = fixture.manifest
        altered.externalRun.reportSHA256 = try VivoNeoantigenWorkbench.digest(invalidUTF8)
        try expectFailure("invalid UTF-8") { _ = try analyze(altered, invalidUTF8) }
        let overLimit = Data(repeating: 65, count: VivoNeoantigenWorkbench.maximumTSVBytes + 1)
        try expectFailure("byte bound") { _ = try analyze(fixture.manifest, overLimit) }
        var review = VivoNeoantigenReview(schema: "numivivo.org/neoantigen-review/v1", reportSHA256: try VivoNeoantigenWorkbench.fingerprint(report),
            reviewerID: "test-reviewer", reviewedAt: Date(timeIntervalSince1970: 1_700_000_000),
            decisions: [.init(candidateID: report.candidates[0].id, disposition: .deferReview, rationale: "Experimental recognition evidence is missing.")])
        try VivoNeoantigenWorkbench.validate(review, against: report); checks += 1
        review.reportSHA256 = String(repeating: "0", count: 64)
        try expectFailure("stale report review") { try VivoNeoantigenWorkbench.validate(review, against: report) }
        review.reportSHA256 = try VivoNeoantigenWorkbench.fingerprint(report)
        review.decisions[0].candidateID = String(repeating: "0", count: 64)
        try expectFailure("foreign candidate") { try VivoNeoantigenWorkbench.validate(review, against: report) }
        review.decisions[0].candidateID = report.candidates[0].id; review.decisions[0].rationale = "  "
        try expectFailure("empty rationale") { try VivoNeoantigenWorkbench.validate(review, against: report) }
        altered = fixture.manifest; altered.sourceCitation = "</script><script>globalThis.__neoInjected=true</script>"
        let injectionHTML = try VivoNeoantigenHTML.render(analyze(altered, fixture.tsv))
        try check(!String(decoding: injectionHTML, as: UTF8.self).contains("globalThis.__neoInjected=true"), "HTML source injection encoded")
        if CommandLine.arguments.count == 3, CommandLine.arguments[1] == "--write-demo" {
            let directory = URL(fileURLWithPath: CommandLine.arguments[2])
            try FileManager.default.createDirectory(at: directory, withIntermediateDirectories: false)
            try VivoCanonicalJSON.encode(fixture.manifest).write(to: directory.appendingPathComponent("case.json"))
            try fixture.tsv.write(to: directory.appendingPathComponent("all_epitopes.tsv"))
            try VivoCanonicalJSON.encode(report).write(to: directory.appendingPathComponent("report.json"))
            try VivoNeoantigenHTML.render(report).write(to: directory.appendingPathComponent("report.html"))
            try injectionHTML.write(to: directory.appendingPathComponent("injection-check.html"))
        } else if CommandLine.arguments.count == 3, CommandLine.arguments[1] == "--check-review" {
            let bytes = try Data(contentsOf: URL(fileURLWithPath: CommandLine.arguments[2]))
            let imported = try VivoCanonicalJSON.decode(VivoNeoantigenReview.self, from: bytes)
            try VivoNeoantigenWorkbench.validate(imported, against: report); checks += 1
        } else if CommandLine.arguments.count != 1 {
            throw VivoNeoantigenError.invalid("Usage: neoantigen-check [--write-demo <new-directory> | --check-review <review.json>]")
        }
        print("PASS: \(checks) portable checks; no Apple package, artifact-store, external predictor or clinical validation claimed.")
    }
}
