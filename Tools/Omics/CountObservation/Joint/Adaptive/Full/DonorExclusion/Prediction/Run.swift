import Foundation
import CryptoKit
import NumiVivoKit
struct Header: Decodable {
    let origin: String
    let queryDonorID: String
    let trainingDonorIDs: [String]
    let firstFeatureIndex: Int
    let featureIDs: [String]
    let libraryCounts: [UInt64]
    let cellsPerLibrary: [Int]
}
struct Gene: Decodable {
    let featureIndex: Int
    let featureID: String
    let status: String
    let model: VivoAdaptiveJointCountModel?
    let bins: [Int]
    let counts: [UInt64]
}
struct Row: Encodable {
    let featureIndex: Int
    let featureID: String
    let status: String
    let prediction: VivoJointCountPrediction?
}
@main struct Main {
    static func main() throws {
        let decoder=JSONDecoder(),encoder=JSONEncoder();encoder.outputFormatting=[.sortedKeys]
        guard let line=readLine() else { throw VivoOmicsError.invalid("missing prediction header") }
        let headerBytes=Data(line.utf8),h=try decoder.decode(Header.self,from: headerBytes)
        guard !h.featureIDs.isEmpty,h.featureIDs.count<=64,Set(h.featureIDs).count==h.featureIDs.count,
              h.firstFeatureIndex>=0,h.trainingDonorIDs.count>=3,
              Set(h.trainingDonorIDs).count==h.trainingDonorIDs.count,
              !h.trainingDonorIDs.contains(h.queryDonorID),
              h.libraryCounts.count==h.cellsPerLibrary.count else {
            throw VivoOmicsError.invalid("prediction donor or feature header")
        }
        var processed=0
        while let line=readLine() {
            try autoreleasepool {
                let bytes=Data(line.utf8),g=try decoder.decode(Gene.self,from: bytes)
                guard processed<h.featureIDs.count,g.featureIndex==h.firstFeatureIndex+processed,
                      g.featureID==h.featureIDs[processed],g.bins.count==g.counts.count else {
                    throw VivoOmicsError.invalid("prediction feature sequence")
                }
                var prediction: VivoJointCountPrediction?
                if g.status=="boundedContinuousLikelihood" {
                    guard let m=g.model,m.status==g.status,
                          m.model.trainingDonorIDs.sorted()==h.trainingDonorIDs.sorted() else {
                        throw VivoOmicsError.invalid("prediction model training identity")
                    }
                    var counts=Array(repeating: UInt64(0),count: h.libraryCounts.count),previous = -1
                    for i in g.bins.indices {
                        let b=g.bins[i]
                        guard b>previous,b<counts.count,g.counts[i]>0 else { throw VivoOmicsError.invalid("prediction count bins") }
                        counts[b]=g.counts[i];previous=b
                    }
                    var hash=SHA256();hash.update(data: headerBytes);hash.update(data: bytes)
                    prediction=try VivoAdaptiveJointCountResponse.predict(
                        control: .init(libraryCounts: h.libraryCounts,cellsPerLibrary: h.cellsPerLibrary,geneCountsPerLibrary: counts),
                        featureID: g.featureID,queryDonorID: h.queryDonorID,controlConditionID: "control",
                        querySource: try .init(bytes: Array(hash.finalize())),model: m,
                        plannedTreatedLibraryCounts: [10_000])
                } else {
                    guard g.bins.isEmpty,g.counts.isEmpty else { throw VivoOmicsError.invalid("unavailable prediction has counts") }
                }
                var output=try encoder.encode(Row(featureIndex: g.featureIndex,featureID: g.featureID,
                    status: prediction == nil ? g.status : "conditionalPrediction",prediction: prediction))
                output.append(10);try FileHandle.standardOutput.write(contentsOf: output);processed+=1
            }
        }
        guard processed==h.featureIDs.count else { throw VivoOmicsError.invalid("truncated prediction stream") }
    }
}
