import Foundation
import CryptoKit
import NumiVivoKit
struct Group: Codable { let donorID: String; let conditionID: String; let sourceRows: [Int] }
struct Metadata: Codable { let origin: String; let featureIDs: [String]; let groups: [Group] }
struct Report: Codable {
    let origin: String
    let metadataSHA256: String
    let streamSHA256: String
    let bindingSHA256: String
    let cells: Int
    let nonzeros: Int
    let moments: [VivoCellCountMoments]
    let models: [VivoCountObservationCalibrationModel]
}
@main struct Main {
    static func main() throws {
        guard CommandLine.arguments.count==2 else { fatalError("calibrate-cells METADATA.json < framed-cells.bin") }
        let bytes=try Data(contentsOf: URL(fileURLWithPath: CommandLine.arguments[1]))
        let meta=try JSONDecoder().decode(Metadata.self,from: bytes)
        guard !meta.groups.isEmpty,meta.groups.count<=128,Set(meta.groups.flatMap(\.sourceRows)).count==meta.groups.reduce(0,{$0+$1.sourceRows.count}) else { fatalError("overlapping source cells") }
        var remaining=meta.groups.map { Set($0.sourceRows) }
        let accumulators=try meta.groups.map { _ in try VivoCellCountMomentAccumulator(featureCount: meta.featureIDs.count) }
        let stdin=FileHandle.standardInput
        var stream=SHA256(),binding=SHA256();binding.update(data: bytes)
        func read(_ n: Int,allowEOF: Bool = false) throws -> Data {
            var data=Data()
            while data.count<n {
                guard let part=try stdin.read(upToCount: n-data.count),!part.isEmpty else {
                    if allowEOF && data.isEmpty { return data };throw NSError(domain: "truncated cell stream",code: 1)
                }
                data.append(part)
            }
            stream.update(data: data);binding.update(data: data);return data
        }
        func u32(_ d: Data,_ offset: Int) -> Int { d.withUnsafeBytes { Int(UInt32(littleEndian: $0.loadUnaligned(fromByteOffset: offset,as: UInt32.self))) } }
        var cells=0,nonzeros=0
        while try autoreleasepool(invoking: {
            let header=try read(12,allowEOF: true);if header.isEmpty { return false }
            let group=u32(header,0),row=u32(header,4),n=u32(header,8)
            guard remaining.indices.contains(group),remaining[group].remove(row) != nil,n>0,n<=meta.featureIDs.count else { fatalError("cell stream identity or sparse shape") }
            let data=try read(n*8)
            var ids: [Int]=[],counts: [UInt64]=[];ids.reserveCapacity(n);counts.reserveCapacity(n)
            for k in 0..<n { ids.append(u32(data,k*8));counts.append(UInt64(u32(data,k*8+4))) }
            try accumulators[group].appendCell(featureIndices: ids,counts: counts)
            cells+=1;nonzeros+=n
            return true
        }) {}
        guard remaining.allSatisfy(\.isEmpty) else { fatalError("missing selected source cells") }
        let moments=try accumulators.map { try $0.finish() },digest=binding.finalize()
        let fingerprint=try VivoFingerprint(bytes: Array(digest))
        var models: [VivoCountObservationCalibrationModel]=[]
        for condition in Set(meta.groups.map(\.conditionID)).sorted() {
            let indices=meta.groups.indices.filter { meta.groups[$0].conditionID==condition }
            models.append(try VivoCountObservationCalibration.fit(indices.map { moments[$0] },featureIDs: meta.featureIDs,
                donorIDs: indices.map { meta.groups[$0].donorID },conditionID: condition,trainingSource: fingerprint))
        }
        func hex<D: Sequence>(_ d: D) -> String where D.Element == UInt8 { d.map { String(format: "%02x",$0) }.joined() }
        let report=Report(origin: meta.origin,metadataSHA256: hex(SHA256.hash(data: bytes)),streamSHA256: hex(stream.finalize()),bindingSHA256: hex(digest),cells: cells,nonzeros: nonzeros,moments: moments,models: models)
        let encoder=JSONEncoder();encoder.outputFormatting=[.sortedKeys]
        try FileHandle.standardOutput.write(contentsOf: encoder.encode(report))
    }
}
