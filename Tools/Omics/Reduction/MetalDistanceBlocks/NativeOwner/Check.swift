import Foundation
@testable import NumiVivoKit
@main struct Check {
 static func main() async throws {
  let root=URL(fileURLWithPath:CommandLine.arguments[1]);try FileManager.default.createDirectory(at:root,withIntermediateDirectories:false)
  let old=try JSONSerialization.jsonObject(with:JSONEncoder().encode(VivoPCANeighborExecution())) as! [String:Any]
  precondition(Set(old.keys) == Set(["workers","queryBlockRows","candidateBlockRows"]))
  var execution=VivoPCANeighborExecution();execution.backend = .metalFP32
  do {try execution.validate();fatalError("workers admitted")} catch {}
  execution.workers=1;execution.queryBlockRows=128;execution.candidateBlockRows=8192
  try execution.validate()
  var options=VivoSingleCellNeighborOptions();options.neighbors=20;options.maximumDistancePairs=100_000_000
  let input=URL(fileURLWithPath:"/Users/n/numivivo-native-mnn-20260910/hagai/pca/scores.bin")
  let begin=Date();var ids=[UInt32](),distances=[Double]()
  try VivoMetalPCANeighbors.stream(source:input,rows:13863,dimensions:20,options:options,execution:execution) {row,i,d in
   precondition(i[0]==row && d[0]==0 && i.count==20);ids += i.map(UInt32.init);distances += d
  }
  try ids.withUnsafeBytes {try Data($0).write(to:root.appendingPathComponent("indices.bin"))}
  try distances.withUnsafeBytes {try Data($0).write(to:root.appendingPathComponent("distances.bin"))}
  print("complete-real-owner-seconds",Date().timeIntervalSince(begin))
  var budget=options;budget.maximumDistancePairs=1
  do {try VivoMetalPCANeighbors.stream(source:input,rows:13863,dimensions:20,options:budget,execution:execution){_,_,_ in fatalError("budget sink")};fatalError("budget admitted")} catch {}
  let cancelOptions=options,cancelExecution=execution
  let task=Task.detached {try VivoMetalPCANeighbors.stream(source:input,rows:13863,dimensions:20,options:cancelOptions,execution:cancelExecution){_,_,_ in}}
  try await Task.sleep(nanoseconds:100_000_000);task.cancel()
  do {try await task.value;fatalError("cancellation ignored")} catch is CancellationError {print("cancellation PASS")}
  print("default-encoding invalid-workers work-budget PASS")
 }
}
