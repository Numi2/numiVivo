import Foundation
@testable import NumiVivoKit

/// Execute the production numerical owner against complete, previously qualified
/// PCA inputs. Full artifact lifecycle is exercised separately through the CLI.
@main struct MNNGaussianMain {
    static func main() throws {
        guard CommandLine.arguments.count == 4 else { throw VivoOmicsError.invalid("mnn-gaussian-check COHORT_ROOT OUTPUT MODE(tiled|scan)") }
        let root = URL(fileURLWithPath: CommandLine.arguments[1]), out = URL(fileURLWithPath: CommandLine.arguments[2])
        let mode = CommandLine.arguments[3]
        guard mode == "tiled" || mode == "scan" else { throw VivoOmicsError.invalid("matching mode") }
        try VivoH5ADCountStore.requireNew(out)
        let original = root.appendingPathComponent("pca"), prior = root.appendingPathComponent("mnn")
        let metadata = try VivoPCAIntegration.read(VivoSingleCellCountMetadata.self, original, "metadata.json", maximum: 536_870_912)
        let old = try VivoPCAIntegration.read(VivoMNNIntegrationReport.self, prior, "report.json", maximum: 16_777_216)
        let plan = try VivoPCAIntegration.read(VivoPCAIntegrationPlan.self, prior, "plan.json", maximum: 65_536)
        guard var options = plan.mnn, options.kernel == nil else { throw VivoOmicsError.invalid("original exhaustive MNN expected") }
        let oldEncoding = try VivoCanonicalJSON.encode(options)
        guard !String(decoding: oldEncoding, as: UTF8.self).contains("kernel") else { throw VivoOmicsError.invalid("legacy option encoding") }
        options.kernel = mode == "tiled" ? .tiledGaussian : nil
        let n = metadata.cells.count, d = old.components
        guard n == old.cells else { throw VivoOmicsError.invalid("cohort axis") }
        let reader = try VivoPCAScoreReader(original.appendingPathComponent("scores.bin"), rows: n, columns: d)
        var x: [Double] = []; x.reserveCapacity(n*d)
        for start in stride(from: 0, to: n, by: 2_048) { x.append(contentsOf: try reader.readRows(start..<min(n,start+2_048))) }
        let cells = metadata.cells.map { VivoOmicsCellIdentity(sampleID: $0.sampleID, barcode: $0.barcode) }
        let start = Date()
        let result = try VivoMNNIntegration.run(cells: cells, scores: x, dimensions: d, samples: metadata.samples, options: options)
        let seconds = Date().timeIntervalSince(start)
        guard result.report.alignments == old.alignments, result.report.assemblyOrder == old.assemblyOrder,
              result.report.steps.count == old.steps.count, result.report.panoramas == old.panoramas,
              result.report.medianRowNorm == old.medianRowNorm, result.report.cellLevels == old.cellLevels,
              result.report.zeroDirectionCells == old.zeroDirectionCells else { throw VivoOmicsError.invalid("original assembly differs") }
        var maximumStepError = 0.0
        for (a,b) in zip(result.report.steps,old.steps) {
            guard a.alignment == b.alignment, a.targetLevels == b.targetLevels,
                  a.referenceLevels == b.referenceLevels, a.anchors == b.anchors,
                  a.correctedCells == b.correctedCells, a.zeroWeightCells == b.zeroWeightCells,
                  a.status == b.status else { throw VivoOmicsError.invalid("step structure differs") }
            for (x,y) in [(a.minimumWeight,b.minimumWeight),(a.maximumWeight,b.maximumWeight),(a.correctionRMS,b.correctionRMS)] {
                let error = abs(x-y)/(1+abs(y)); maximumStepError = max(maximumStepError,error)
                guard error <= 1e-10 else { throw VivoOmicsError.invalid("step numerical tolerance") }
            }
        }
        let oldReader = try VivoPCAScoreReader(prior.appendingPathComponent("scores.bin"), rows: n, columns: d)
        var maximumScoreError = 0.0, maximumOriginalScore = 0.0, exact = true
        for start in stride(from: 0, to: n, by: 2_048) {
            let values = try oldReader.readRows(start..<min(n,start+2_048))
            for (j,value) in values.enumerated() {
                let actual = result.scores[start*d+j]
                maximumScoreError = max(maximumScoreError,abs(actual-value))
                maximumOriginalScore = max(maximumOriginalScore,abs(value))
                exact = exact && actual.bitPattern == value.bitPattern
            }
        }
        guard maximumScoreError <= 1e-10*(1+maximumOriginalScore) else { throw VivoOmicsError.invalid("complete coordinate tolerance") }
        try FileManager.default.createDirectory(at: out, withIntermediateDirectories: false)
        let scores = try VivoCountRecordWriter(out.appendingPathComponent("scores.bin"))
        for i in 0..<n { for j in 0..<d { try scores.append(row: i, feature: j, bits: result.scores[i*d+j].bitPattern) } }
        _ = try scores.finish()
        let anchors = try VivoCountRecordWriter(out.appendingPathComponent("anchors.bin"))
        for (i,a) in result.anchors.enumerated() { try anchors.append(row: a.first, feature: a.second, bits: result.anchorDistances[i].bitPattern) }
        _ = try anchors.finish()
        for name in (mode == "scan" ? ["scores.bin","anchors.bin"] : ["anchors.bin"]) {
            guard try VivoOmicsFileSnapshot.fingerprint(prior.appendingPathComponent(name), maximumBytes: 1_000_000_000) == VivoOmicsFileSnapshot.fingerprint(out.appendingPathComponent(name), maximumBytes: 1_000_000_000) else { throw VivoOmicsError.invalid("complete original bytes differ: \(name)") }
        }
        if mode == "scan" { guard result.report == old else { throw VivoOmicsError.invalid("legacy report changed") } }
        _ = try VivoPCAIntegration.write(result.report,out,"report.json",maximum:16_777_216)
        let record: [String:Any] = ["status":"passed","mode":mode,"cells":n,"components":d,"anchors":result.anchors.count,"seconds":seconds,"allOriginalAnchorBytesExact":true,"allOriginalScoreBytesExact":exact,"maximumScoreError":maximumScoreError,"maximumOriginalScore":maximumOriginalScore,"maximumRelativeStepError":maximumStepError,"allAssemblyStructuresExact":true,"legacyOptionEncodingPreserved":true,"scope":"Production numerical owner on complete original PCA; artifact lifecycle and biological evidence remain separate."]
        try JSONSerialization.data(withJSONObject:record,options:[.sortedKeys,.prettyPrinted]).write(to:out.appendingPathComponent("checks.json"))
        print("\(mode) \(n) cells: all anchors/assembly exact; score max error \(maximumScoreError); \(seconds) seconds")
    }
}
