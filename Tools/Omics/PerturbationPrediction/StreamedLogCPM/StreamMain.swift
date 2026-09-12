import Foundation
import NumiVivoKit
struct Plan: Decodable { let featureIDs: [String]; let groupIDs: [String]; let rowGroups: [Int]; let rowTotals: [UInt64] }
@main struct Main {
 static func main() throws {
  let p=try JSONDecoder().decode(Plan.self,from:Data(contentsOf:URL(fileURLWithPath:CommandLine.arguments[1])))
  let result=try VivoStreamedLogCPM.consume(featureIDs:p.featureIDs,groupIDs:p.groupIDs,rowGroups:p.rowGroups,rowTotals:p.rowTotals) {
   try FileHandle.standardInput.read(upToCount:1048576) ?? Data()
  }
  try JSONEncoder().encode(result).write(to:URL(fileURLWithPath:CommandLine.arguments[2]),options:.withoutOverwriting)
 }
}
