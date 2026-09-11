import Foundation
import CryptoKit
import NumiVivoKit
struct Group: Decodable {
    let donorID: String
    let conditionID: String
    let libraryCounts: [UInt64]
    let cellsPerLibrary: [Int]
}
struct Header: Decodable {
    let origin: String
    let firstFeatureIndex: Int
    let featureIDs: [String]
    let groups: [Group]
    let cacheSHA256: String
    let manifestSHA256: String
}
struct SparseCounts: Decodable { let bins: [Int];let counts: [UInt64] }
struct Gene: Decodable {
    let featureIndex: Int
    let featureID: String
    let controlCellDispersion: Double?
    let treatedCellDispersion: Double?
    let controlCalibrationStatus: String
    let treatedCalibrationStatus: String
    let groups: [SparseCounts]
}
struct Row: Encodable {
    let origin: String
    let featureIndex: Int
    let featureID: String
    let status: String
    let model: VivoAdaptiveJointCountModel?
    let detail: String?
}
@main struct Main {
    static func main() throws {
        guard let first=readLine(),let headerData=first.data(using: .utf8) else { throw VivoOmicsError.invalid("missing full-gene stream header") }
        let decoder=JSONDecoder(),encoder=JSONEncoder();encoder.outputFormatting=[.sortedKeys]
        let header=try decoder.decode(Header.self,from: headerData)
        guard !header.featureIDs.isEmpty,header.featureIDs.count<=256,Set(header.featureIDs).count==header.featureIDs.count,
              header.firstFeatureIndex>=0,header.groups.count>=6,header.groups.count<=128,
              header.cacheSHA256.count==64,header.manifestSHA256.count==64 else { throw VivoOmicsError.invalid("full-gene header identities or budget") }
        let donors=Set(header.groups.map(\.donorID)).sorted()
        var pairIndices=[(String,Int,Int)]()
        for donor in donors {
            let c=header.groups.indices.filter { header.groups[$0].donorID==donor && header.groups[$0].conditionID=="control" }
            let t=header.groups.indices.filter { header.groups[$0].donorID==donor && header.groups[$0].conditionID=="IFNB" }
            guard c.count==1,t.count==1 else { throw VivoOmicsError.invalid("full-gene donor conditions") }
            pairIndices.append((donor,c[0],t[0]))
        }
        guard pairIndices.count*2==header.groups.count else { throw VivoOmicsError.invalid("full-gene extra condition") }
        var processed=0
        while let line=readLine() {
            try autoreleasepool {
                let bytes=Data(line.utf8),gene=try decoder.decode(Gene.self,from: bytes)
                guard processed<header.featureIDs.count,gene.featureIndex==header.firstFeatureIndex+processed,
                      gene.featureID==header.featureIDs[processed] else { throw VivoOmicsError.invalid("full-gene ordering or identity") }
                var digest=SHA256();digest.update(data: headerData);digest.update(data: bytes)
                let source=try VivoFingerprint(bytes: Array(digest.finalize()))
                var model: VivoAdaptiveJointCountModel?,status="unavailableCellDispersion",detail: String?="control="+gene.controlCalibrationStatus+"; treated="+gene.treatedCalibrationStatus
                if let cp=gene.controlCellDispersion,let tp=gene.treatedCellDispersion {
                    guard gene.groups.count==header.groups.count else { throw VivoOmicsError.invalid("full-gene group count") }
                    var strata=[VivoCountDepthStratum]()
                    for i in header.groups.indices {
                        let g=header.groups[i],s=gene.groups[i]
                        guard s.bins.count==s.counts.count else { throw VivoOmicsError.invalid("full-gene sparse count dimensions") }
                        var previous = -1,counts=Array(repeating: UInt64(0),count: g.libraryCounts.count)
                        for j in s.bins.indices {
                            let bin=s.bins[j];guard bin>previous,bin<counts.count,s.counts[j]>0 else { throw VivoOmicsError.invalid("full-gene sparse count order or range") }
                            counts[bin]=s.counts[j];previous=bin
                        }
                        strata.append(.init(libraryCounts: g.libraryCounts,cellsPerLibrary: g.cellsPerLibrary,geneCountsPerLibrary: counts))
                    }
                    let pairs=pairIndices.map { VivoJointCountPair(donorID: $0.0,control: strata[$0.1],treated: strata[$0.2]) }
                    do {
                        model=try VivoAdaptiveJointCountResponse.fit(pairs: pairs,featureID: gene.featureID,controlConditionID: "control",treatedConditionID: "IFNB",controlCellDispersion: cp,treatedCellDispersion: tp,trainingSource: source)
                        status=model!.status;detail=nil
                    } catch { status="fitError";detail=String(describing: error) }
                } else { guard gene.groups.isEmpty else { throw VivoOmicsError.invalid("unavailable gene unexpectedly carries fit counts") } }
                var output=try encoder.encode(Row(origin: header.origin,featureIndex: gene.featureIndex,featureID: gene.featureID,status: status,model: model,detail: detail));output.append(10);try FileHandle.standardOutput.write(contentsOf: output);processed+=1
            }
        }
        guard processed==header.featureIDs.count else { throw VivoOmicsError.invalid("truncated full-gene stream") }
    }
}
