import Foundation
import CryptoKit
import NumiVivoKit

struct Input: Codable {
    let origin: String
    let featureIDs: [String]
    let groups: [VivoPairedCountGroup]
    let controlConditionID: String
    let treatedConditionID: String
}
struct Result: Codable {
    let origin: String
    let excludedDonorID: String?
    let report: VivoPairedCountReport
}
@main struct Main {
    static func main() throws {
        let bytes=try FileHandle.standardInput.readToEnd()!
        let input=try JSONDecoder().decode(Input.self,from: bytes)
        let source=try VivoFingerprint(bytes: Array(SHA256.hash(data: bytes)))
        let encoder=JSONEncoder();encoder.outputFormatting=[.sortedKeys]
        let donors=Set(input.groups.map(\.donorID)).sorted()
        // Every paired donor is omitted once for sensitivity, not prediction.
        for excluded in [nil]+donors.map({ Optional($0) }) {
            try autoreleasepool {
                let groups=input.groups.filter { $0.donorID != excluded }
                let report=try VivoPairedCountMoments.analyze(groups: groups,featureIDs: input.featureIDs,
                    controlConditionID: input.controlConditionID,treatedConditionID: input.treatedConditionID,trainingSource: source)
                var line=try encoder.encode(Result(origin: input.origin,excludedDonorID: excluded,report: report));line.append(10)
                try FileHandle.standardOutput.write(contentsOf: line)
            }
        }
    }
}
