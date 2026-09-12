import Foundation
@main struct Main {
 static func main() throws {
  func make(_ ordering:VivoStreamedLogCPM.Ordering) throws -> VivoStreamedLogCPM {
   try .init(featureIDs:["a","b"],groupIDs:["x"],rowGroups:[0,0,0],rowTotals:[3,4,0],ordering:ordering)
  }
  var row=try make(.rowMajor);for (r,f,c) in [(0,0,1),(0,1,2),(1,0,4)] {try row.add(row:r,feature:f,count:UInt64(c))}
  var col=try make(.featureMajor);for (r,f,c) in [(0,0,1),(1,0,4),(0,1,2)] {try col.add(row:r,feature:f,count:UInt64(c))}
  let expected=try row.finish(),actual=try col.finish();precondition(expected==actual)
  for records in [[(0,0,1),(0,0,1)],[(1,0,4),(0,0,1)],[(0,1,2),(1,0,4)],[(0,0,4)],[(2,0,1)],[(0,0,1)]] {
   var a=try make(.featureMajor);var rejected=false
   do {for (r,f,c) in records {try a.add(row:r,feature:f,count:UInt64(c))};_ = try a.finish()}
   catch {rejected=true}
   precondition(rejected)
   do {_ = try a.finish();fatalError("failed accumulator reused")} catch VivoStreamedLogCPM.Failure.closed {}
  }
  print("PASS layout equality, zero cells, duplicate, ordering, excess, missing and failed reuse")
 }
}
