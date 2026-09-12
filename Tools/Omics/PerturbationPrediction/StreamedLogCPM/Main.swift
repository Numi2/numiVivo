import Foundation
struct Plan: Decodable { let featureIDs: [String]; let groupIDs: [String]; let rowGroups: [Int]; let rowTotals: [UInt64] }
@main struct Main {
 static func main() throws {
  let root=URL(fileURLWithPath:CommandLine.arguments[1])
  let p=try JSONDecoder().decode(Plan.self,from:Data(contentsOf:root.appendingPathComponent("plan.json")))
  var a=try VivoStreamedLogCPM(featureIDs:p.featureIDs,groupIDs:p.groupIDs,rowGroups:p.rowGroups,rowTotals:p.rowTotals)
  let h=try FileHandle(forReadingFrom:root.appendingPathComponent("counts.bin"))
  defer { try? h.close() }
  while let data=try h.read(upToCount:1048576), !data.isEmpty {
   guard data.count % 16 == 0 else { fatalError("truncated") }
   try data.withUnsafeBytes { b in
    for offset in stride(from:0,to:data.count,by:16) {
     let row=Int(UInt32(littleEndian:b.loadUnaligned(fromByteOffset:offset,as:UInt32.self)))
     let feature=Int(UInt32(littleEndian:b.loadUnaligned(fromByteOffset:offset+4,as:UInt32.self)))
     let count=UInt64(littleEndian:b.loadUnaligned(fromByteOffset:offset+8,as:UInt64.self))
     try a.add(row:row,feature:feature,count:count)
    }
   }
  }
  let result=try a.finish()
  try JSONEncoder().encode(result).write(to:root.appendingPathComponent("native.json"),options:.withoutOverwriting)
  func rejected(_ f: () throws -> Void) { do { try f(); fatalError("accepted invalid input") } catch {} }
  rejected { _ = try a.finish() }
  rejected { try a.add(row:0,feature:0,count:1) }
  var z=try VivoStreamedLogCPM(featureIDs:["a"],groupIDs:["g"],rowGroups:[0,0],rowTotals:[0,2])
  try z.add(row:1,feature:0,count:2)
  let zr=try z.finish(); precondition(zr.cellCounts == [2] && zr.zeroCellCounts == [1] && abs(zr.means[0][0]-log1p(1e6)/2)<1e-12)
  var empty=try VivoStreamedLogCPM(featureIDs:["a"],groupIDs:["g"],rowGroups:[0],rowTotals:[0])
  let er = try empty.finish(); precondition(er.means == [[0]])
  for records in [[(0,0,UInt64(1)),(0,0,1)],[(0,0,3)],[(1,0,1)],[(0,1,1)],[(0,0,0)]] {
   var bad=try VivoStreamedLogCPM(featureIDs:["a"],groupIDs:["g"],rowGroups:[0],rowTotals:[2])
   rejected { for r in records { try bad.add(row:r.0,feature:r.1,count:r.2) }; _ = try bad.finish() }
   rejected { _ = try bad.finish() }
  }
  var missing=try VivoStreamedLogCPM(featureIDs:["a"],groupIDs:["g"],rowGroups:[0],rowTotals:[2])
  rejected { _ = try missing.finish() }
  print("PASS native real counts and malformed/zero-cell checks")
 }
}
