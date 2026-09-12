import Foundation
@main struct Main {
 static func main() async throws {
  func consume(_ read: () throws -> Data) throws -> VivoStreamedLogCPM.Result {
   try VivoStreamedLogCPM.consume(featureIDs:["g"],groupIDs:["c"],rowGroups:[0],rowTotals:[0],read:read)
  }
  do { _ = try consume { Data(repeating:0,count:1_048_577) }; fatalError("oversize accepted") }
  catch VivoStreamedLogCPM.Failure.invalidRecord {}
  do { var n=0; _ = try consume { n+=1; return n==1 ? Data([0]) : Data() }; fatalError("partial accepted") }
  catch VivoStreamedLogCPM.Failure.invalidRecord {}
  enum Injected: Error {case read}
  do {_ = try consume { throw Injected.read };fatalError("read error swallowed")}
  catch Injected.read {}
  let task = Task { () throws -> Void in
   withUnsafeCurrentTask { $0?.cancel() }
   _ = try consume { fatalError("read after cancellation") }
  }
  do {try await task.value;fatalError("cancellation ignored")}
  catch is CancellationError {}
  print("PASS oversize partial read-error and cancellation")
 }
}
