import Foundation
@testable import NumiVivoKit
@main struct Controls {
 static func main() async throws {
  let url=URL(fileURLWithPath:CommandLine.arguments[1])
  let n=VivoWindowedCountRecords.windowBytes/16+3
  var data=Data()
  for i in 0..<n {
   var row=UInt32(i%31).littleEndian,feature=UInt32(i%17).littleEndian,bits=Double(i+1).bitPattern.littleEndian
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
    precondition(row==seen%31 && feature==seen%17 && bits==Double(seen+1).bitPattern);seen+=1
   }
   precondition(seen==n)
  }
  let entries=try VivoMappedReductionEntries(url:url,expectedEntries:n,rows:31,columns:17)
  let v=(0..<17).map { Double($0)-8.5 }, u=(0..<31).map { Double($0)/8-2 }
  var expectedProject=[Double](repeating:-0.25,count:31),expectedTranspose=[Double](repeating:0.75,count:17)
  for i in 0..<n {
   expectedProject[i%31] += Double(i+1)*v[i%17]
   expectedTranspose[i%17] += Double(i+1)*u[i%31]
  }
  let project=try entries.project(v,shift:0.25,rows:31)
  let transpose=try entries.transpose(u,initial:[Double](repeating:0.75,count:17))
  precondition(project==expectedProject && transpose==expectedTranspose && entries.visits==2*n)
  print("PASS: full-window unsorted sparse project and transpose match sequential FP64 reference")
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
