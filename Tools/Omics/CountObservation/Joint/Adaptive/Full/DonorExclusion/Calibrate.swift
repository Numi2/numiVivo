import Foundation
import CryptoKit
import NumiVivoKit

struct Input: Decodable {
    let origin: String
    let excludedDonorID: String
    let featureIDs: [String]
    let groups: [VivoPairedCountGroup]
    let controlConditionID: String
    let treatedConditionID: String
}
struct Output: Encodable {
    let origin: String
    let excludedDonorID: String
    let inputSHA256: String
    let control: VivoCountObservationCalibrationModel
    let treated: VivoCountObservationCalibrationModel
}
@main struct Main {
    static func main() throws {
        let bytes=try FileHandle.standardInput.readToEnd()!
        let input=try JSONDecoder().decode(Input.self,from: bytes)
        let ids=Set(input.groups.map(\.donorID)).sorted()
        guard !ids.contains(input.excludedDonorID),input.controlConditionID != input.treatedConditionID,
              input.groups.count==ids.count*2 else { throw VivoOmicsError.invalid("excluded donor or paired condition contract") }
        let source=try VivoFingerprint(bytes: Array(SHA256.hash(data: bytes)))
        func fit(_ condition: String) throws -> VivoCountObservationCalibrationModel {
            let groups=try ids.map { id -> VivoCellCountMoments in
                let selected=input.groups.filter { $0.donorID==id && $0.conditionID==condition }
                guard selected.count==1 else { throw VivoOmicsError.invalid("missing or duplicate training donor condition") }
                return selected[0].moments
            }
            return try VivoCountObservationCalibration.fit(groups,featureIDs: input.featureIDs,donorIDs: ids,conditionID: condition,trainingSource: source)
        }
        let output=try Output(origin: input.origin,excludedDonorID: input.excludedDonorID,inputSHA256: source.bytes.map { String(format:"%02x",$0) }.joined(),control: fit(input.controlConditionID),treated: fit(input.treatedConditionID))
        let encoder=JSONEncoder();encoder.outputFormatting=[.sortedKeys]
        try FileHandle.standardOutput.write(contentsOf: encoder.encode(output))
    }
}
