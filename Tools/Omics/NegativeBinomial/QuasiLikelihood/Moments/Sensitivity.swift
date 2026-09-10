import Foundation
struct Input: Decodable { let counts: [UInt64]; let means: [Double]; let design: [[Double]]; let dispersion: Double; let averageQuasiDispersion: Double }
struct Result: Encodable { let value: VivoOmicsNBAdjustedResiduals?; let error: String? }
@main struct Main {
 static func main() throws {
  let input = try JSONDecoder().decode([Input].self,from: FileHandle.standardInput.readDataToEndOfFile())
  let result = input.map { v -> Result in
   do { return .init(value: try VivoOmicsNBResidualAdjustment.adjustedResiduals(counts: v.counts,means: v.means,design: v.design,dispersion: v.dispersion,averageQuasiDispersion: v.averageQuasiDispersion,relativeTolerance: 1e-10,maximumTerms: 10_000_000),error: nil) }
   catch { return .init(value: nil,error: error.localizedDescription) }
  }
  let encoder = JSONEncoder();encoder.outputFormatting = [.sortedKeys];FileHandle.standardOutput.write(try encoder.encode(result))
 }
}
