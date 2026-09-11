import Foundation
import NumiVivoKit

struct Case: Codable {
    let id: String
    let counts: [UInt64]
    let libraryCounts: [UInt64]
    let cellDispersion: Double
    let gammaPriorShape: Double
    let gammaPriorRatePerCPM: Double
    let plannedLibraryCounts: [UInt64]
}
struct Result: Codable {
    let id: String
    let posterior: VivoCountObservationPosterior?
    let prediction: VivoCountObservationPredictiveMoments?
    let error: String?
}
@main struct Main {
    static func main() throws {
        let args=CommandLine.arguments
        guard args.count==3 else { fatalError("count-observation-check INPUT.json OUTPUT.jsonl") }
        let cases=try JSONDecoder().decode([Case].self,from: Data(contentsOf: URL(fileURLWithPath: args[1])))
        guard !FileManager.default.fileExists(atPath: args[2]) else { fatalError("output exists") }
        _=FileManager.default.createFile(atPath: args[2],contents: nil)
        let output=try FileHandle(forWritingTo: URL(fileURLWithPath: args[2]))
        defer { try? output.close() }
        let encoder=JSONEncoder();encoder.outputFormatting=[.sortedKeys]
        var failed=0
        for c in cases {
            let r: Result
            do {
                let p=try VivoCountObservation.posterior(counts: c.counts,libraryCounts: c.libraryCounts,cellDispersion: c.cellDispersion,gammaPriorShape: c.gammaPriorShape,gammaPriorRatePerCPM: c.gammaPriorRatePerCPM)
                r = .init(id: c.id,posterior: p,prediction: try VivoCountObservation.predictiveMoments(p,plannedLibraryCounts: c.plannedLibraryCounts),error: nil)
            } catch { failed += 1;r = .init(id: c.id,posterior: nil,prediction: nil,error: String(describing: error)) }
            var data=try encoder.encode(r);data.append(10);try output.write(contentsOf: data)
        }
        print("cases=\(cases.count) failed=\(failed)")
        if failed>0 { exit(1) }
    }
}
