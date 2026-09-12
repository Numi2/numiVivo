import Foundation
@testable import NumiVivoKit
@main struct Controls {
 static func main() async throws {
  let url=URL(fileURLWithPath:CommandLine.arguments[1])
  let n=VivoWindowedCountRecords.windowBytes/16+3
  var data=Data()
  for i in 0..<n {
   var row=UInt32(i%31).littleEndian,feature=UInt32(i%17).littleEndian,bits=UInt64(i).littleEndian
   withUnsafeBytes(of:&row){data.append(contentsOf:$0)}
   withUnsafeBytes(of:&feature){data.append(contentsOf:$0)}
   withUnsafeBytes(of:&bits){data.append(contentsOf:$0)}
  }
  try data.write(to:url,options:.withoutOverwriting)
  let reader=try VivoWindowedCountRecords(url,entries:n)
  _ = try reader.record(n-1)
  for _ in 0..<2 {
   var seen=0
   try reader.forEachRecord { row,feature,bits in
    precondition(row==seen%31 && feature==seen%17 && bits==UInt64(seen));seen+=1
   }
   precondition(seen==n)
  }
  let task=Task {
   withUnsafeCurrentTask { $0?.cancel() }
   let cancelledReader=try VivoWindowedCountRecords(url,entries:n)
   do {try cancelledReader.forEachRecord{_,_,_ in fatalError("cancelled traversal visited record")};fatalError("cancellation missing")}
   catch is CancellationError {}
  }
  try await task.value
  print("PASS: two complete traversals, partial final window, prior random access, cancellation before traversal")
 }
}
