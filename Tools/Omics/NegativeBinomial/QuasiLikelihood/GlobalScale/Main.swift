import Foundation
struct Input: Decodable {
 let counts: [[UInt64]], design: [[Double]], offsets: [Double], contrast: [Double]
 let trendDispersions: [Double], abundanceCovariates: [Double]
}
struct Output: Encodable { let fit: VivoOmicsNBQLGlobalFit?; let error: String? }
@main struct Main {
 static func main() throws {
  let i = try JSONDecoder().decode(Input.self,from: FileHandle.standardInput.readDataToEndOfFile())
  let out: Output
  do { out = .init(fit: try VivoOmicsNBQLGlobalScale.fit(counts: i.counts,design: i.design,offsets: i.offsets,contrast: i.contrast,trendDispersions: i.trendDispersions,abundanceCovariates: i.abundanceCovariates),error: nil) }
  catch { out = .init(fit: nil,error: error.localizedDescription) }
  let encoder = JSONEncoder();encoder.outputFormatting = [.sortedKeys]
  FileHandle.standardOutput.write(try encoder.encode(out))
 }
}
