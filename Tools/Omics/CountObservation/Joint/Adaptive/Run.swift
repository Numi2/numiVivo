import Foundation
import CryptoKit
import NumiVivoKit
struct Gene: Decodable {
    let featureID: String
    let controlCellDispersion: Double?
    let treatedCellDispersion: Double?
    let controlCalibrationStatus: String
    let treatedCalibrationStatus: String
    let pairs: [VivoJointCountPair]
}
struct Query: Decodable {
    let id: String
    let featureID: String
    let donorID: String
    let control: VivoCountDepthStratum
    let plannedLibraryCounts: [UInt64]
}
struct Input: Decodable { let origin: String;let genes: [Gene];let queries: [Query];let gridSizes: [Int] }
struct Row: Encodable {
    let kind: String;let modelID: String;let queryID: String?;let status: String
    let model: VivoAdaptiveJointCountModel?;let prediction: VivoJointCountPrediction?;let detail: String?
}
@main struct Main {
    static func main() throws {
        let bytes=try FileHandle.standardInput.readToEnd()!,input=try JSONDecoder().decode(Input.self,from: bytes)
        let source=try VivoFingerprint(bytes: Array(SHA256.hash(data: bytes)))
        let encoder=JSONEncoder();encoder.outputFormatting=[.sortedKeys]
        func emit(_ row: Row) throws { var d=try encoder.encode(row);d.append(10);try FileHandle.standardOutput.write(contentsOf: d) }
        var failures=0
        for g in input.genes { for grid in [9] {
            try autoreleasepool {
                let id=input.origin+":"+g.featureID+":grid="+String(grid)
                var model: VivoAdaptiveJointCountModel?,status="unavailableCellDispersion"
                var detail: String?="control="+g.controlCalibrationStatus+"; treated="+g.treatedCalibrationStatus
                if let cp=g.controlCellDispersion,let tp=g.treatedCellDispersion {
                    do {
                        model=try VivoAdaptiveJointCountResponse.fit(pairs: g.pairs,featureID: g.featureID,controlConditionID: "control",treatedConditionID: "IFNB",
                            controlCellDispersion: cp,treatedCellDispersion: tp,trainingSource: source,plan: .init(initialGridPointsPerAxis: grid))
                        status=model!.status;detail=nil
                    } catch { failures+=1;status="error";detail=String(describing: error) }
                }
                try emit(.init(kind: "fit",modelID: id,queryID: nil,status: status,model: model,prediction: nil,detail: detail))
                for q in input.queries where q.featureID==g.featureID {
                    if let model,status=="boundedContinuousLikelihood" {
                        do {
                            let p=try VivoAdaptiveJointCountResponse.predict(control: q.control,featureID: q.featureID,queryDonorID: q.donorID,controlConditionID: "control",querySource: source,model: model,plannedTreatedLibraryCounts: q.plannedLibraryCounts)
                            try emit(.init(kind: "query",modelID: id,queryID: q.id,status: "conditionalPrediction",model: nil,prediction: p,detail: nil))
                        } catch { failures+=1;try emit(.init(kind: "query",modelID: id,queryID: q.id,status: "error",model: nil,prediction: nil,detail: String(describing: error))) }
                    } else { try emit(.init(kind: "query",modelID: id,queryID: q.id,status: status,model: nil,prediction: nil,detail: detail)) }
                }
            }
        } }
        if failures>0 { exit(1) }
    }
}
