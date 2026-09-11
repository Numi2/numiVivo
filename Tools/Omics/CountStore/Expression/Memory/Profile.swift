import Foundation
import Darwin
import NumiVivoKit
@main struct Profile {
 static func mark(_ stage: String) {
  var usage = rusage(); _ = getrusage(RUSAGE_SELF, &usage)
  FileHandle.standardOutput.write(Data("{\"stage\":\"\(stage)\",\"peakRSS\":\(usage.ru_maxrss),\"unix\":\(Date().timeIntervalSince1970)}\n".utf8))
 }
 static func main() throws {
  let args = CommandLine.arguments
  let plan = try VivoCanonicalJSON.decode(VivoFileExpressionPlan.self, from: Data(contentsOf: URL(fileURLWithPath: args[2])))
  mark("before-evaluate")
  let result = try VivoFileExpression.run(source: URL(fileURLWithPath: args[1]), plan: plan)
  mark("after-evaluate")
  let raw = try VivoCanonicalJSON.encode(result)
  mark("after-whole-report-encode")
  let digest = try VivoCanonicalJSON.fingerprint(raw).bytes.map { String(format: "%02x", $0) }.joined()
  guard digest == "4919f0c96dbe114707dd58b5939d708b6a3e9e7a712cb91380edd2424d92ed11" else { fatalError("Report changed") }
  print("verified-report \(raw.count) \(digest)")
 }
}
