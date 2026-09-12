import Foundation
import NumiVivoKit
struct Plan: Decodable { let featureIDs: [String]; let groupIDs: [String]; let rowGroups: [Int]; let rowTotals: [UInt64] }
@main struct Main {
 static func main() throws {
  let p=try JSONDecoder().decode(Plan.self,from:Data(contentsOf:URL(fileURLWithPath:CommandLine.arguments[1])))
  var a=try VivoStreamedLogCPM(featureIDs:p.featureIDs,groupIDs:p.groupIDs,rowGroups:p.rowGroups,rowTotals:p.rowTotals)
  var buffer=Data()
  while true {
   let chunk=try FileHandle.standardInput.read(upToCount:1048576) ?? Data()
   if chunk.isEmpty { break }
   buffer.append(chunk)
   let complete=buffer.count / 16 * 16
   try buffer.withUnsafeBytes { b in
    for offset in stride(from:0,to:complete,by:16) {
     try a.add(row:Int(UInt32(littleEndian:b.loadUnaligned(fromByteOffset:offset,as:UInt32.self))),
       feature:Int(UInt32(littleEndian:b.loadUnaligned(fromByteOffset:offset+4,as:UInt32.self))),
       count:UInt64(littleEndian:b.loadUnaligned(fromByteOffset:offset+8,as:UInt64.self)))
    }
   }
   buffer=Data(buffer.suffix(buffer.count-complete))
  }
  guard buffer.isEmpty else { throw VivoStreamedLogCPM.Failure.invalidRecord }
  let result=try a.finish()
  try JSONEncoder().encode(result).write(to:URL(fileURLWithPath:CommandLine.arguments[2]),options:.withoutOverwriting)
 }
}
