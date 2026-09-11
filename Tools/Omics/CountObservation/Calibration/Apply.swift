import Foundation
import NumiVivoKit
struct Query: Codable { let id: String; let featureID: String; let donorID: String; let counts: [UInt64]; let libraryCounts: [UInt64]; let plannedLibraryCounts: [UInt64] }
struct Input: Codable { let models: [VivoCountObservationCalibrationModel]; let queries: [Query] }
struct Row: Codable {
    let model: Int
    let id: String
    let status: String
    let posterior: VivoCountObservationPosterior?
    let prediction: VivoCountObservationPredictiveMoments?
    let failure: String?
}
@main struct Main {
    static func main() throws {
        let input=try JSONDecoder().decode(Input.self,from: Data(contentsOf: URL(fileURLWithPath: CommandLine.arguments[1])))
        let encoder=JSONEncoder();encoder.outputFormatting=[.sortedKeys]
        var errors=0
        for (i,model) in input.models.enumerated() {
            for q in input.queries {
                let result: Row
                if let feature=model.features.first(where: { $0.featureID==q.featureID }),feature.status != "availableConditionalMomentCalibration" {
                    result = .init(model: i,id: q.id,status: "unavailableCalibration",posterior: nil,prediction: nil,failure: feature.status)
                } else {
                    do {
                        let p=try VivoCountObservationCalibration.posterior(counts: q.counts,libraryCounts: q.libraryCounts,featureID: q.featureID,queryDonorID: q.donorID,conditionID: "control",model: model)
                        result = .init(model: i,id: q.id,status: "availableConditionalPosterior",posterior: p,prediction: try VivoCountObservation.predictiveMoments(p,plannedLibraryCounts: q.plannedLibraryCounts),failure: nil)
                    } catch { errors+=1;result = .init(model: i,id: q.id,status: "error",posterior: nil,prediction: nil,failure: String(describing: error)) }
                }
                var data=try encoder.encode(result);data.append(10);try FileHandle.standardOutput.write(contentsOf: data)
            }
        }
        if errors>0 { exit(1) }
    }
}
