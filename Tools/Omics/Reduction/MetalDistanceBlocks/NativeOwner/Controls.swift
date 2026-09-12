import Foundation
@testable import NumiVivoKit
@main struct Controls {
 static func main() throws {
  let root=URL(fileURLWithPath:CommandLine.arguments[1]);try FileManager.default.createDirectory(at:root,withIntermediateDirectories:false)
  func write(_ values:[Double],_ name:String) throws -> URL {
   var data=Data()
   for (i,v) in values.enumerated() {
    var row=UInt32(i),col:UInt32=0,x=v
    withUnsafeBytes(of:&row){data.append(contentsOf:$0)};withUnsafeBytes(of:&col){data.append(contentsOf:$0)};withUnsafeBytes(of:&x){data.append(contentsOf:$0)}
   }
   let url=root.appendingPathComponent(name);try data.write(to:url);return url
  }
  var execution=VivoPCANeighborExecution();execution.backend = .metalFP32;execution.workers=1;execution.queryBlockRows=3;execution.candidateBlockRows=2
  var options=VivoSingleCellNeighborOptions();options.neighbors=3
  let duplicate=try write([0,0,0,1],"duplicates.bin");var rows=[[Int]]()
  try VivoMetalPCANeighbors.stream(source:duplicate,rows:4,dimensions:1,options:options,execution:execution){_,ids,_ in rows.append(ids)}
  precondition(rows == [[0,1,2],[1,0,2],[2,0,1],[3,0,1]])
  options.neighbors=2
  for (name,values) in [("conversion-overflow",[Double.greatestFiniteMagnitude,0]),("distance-overflow",[1e30,0]),("nonfinite",[Double.nan,0])] {
   let file=try write(values,name)
   do {try VivoMetalPCANeighbors.stream(source:file,rows:2,dimensions:1,options:options,execution:execution){_,_,_ in fatalError("invalid sink")};fatalError("invalid admitted")} catch {print(name,"rejected",error)}
  }
  print("duplicate ties and partial query/candidate tiles PASS")
 }
}
