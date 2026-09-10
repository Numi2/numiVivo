import Foundation
struct MomentInput: Decodable { let mean: Double; let dispersion: Double; let relativeTolerance: Double?; let maximumTerms: Int? }
struct ResidualInput: Decodable { let counts: [UInt64]; let means: [Double]; let design: [[Double]]; let dispersion: Double; let averageQuasiDispersion: Double }
struct Input: Decodable { let moments: [MomentInput]?; let residuals: [ResidualInput]? }
struct Result<T: Encodable>: Encodable { let value: T?; let error: String? }
struct Output: Encodable { let moments: [Result<VivoOmicsNBDevianceMoments>]; let residuals: [Result<VivoOmicsNBAdjustedResiduals>] }
@main struct Main {
 static func main() throws {
  let input = try JSONDecoder().decode(Input.self,from: FileHandle.standardInput.readDataToEndOfFile())
  let moments = (input.moments ?? []).map { v -> Result<VivoOmicsNBDevianceMoments> in
   do { return .init(value: try VivoOmicsNBResidualAdjustment.moments(mean: v.mean,dispersion: v.dispersion,relativeTolerance: v.relativeTolerance ?? 1e-10,maximumTerms: v.maximumTerms ?? 1_000_000,method: .adaptive),error: nil) }
   catch { return .init(value: nil,error: error.localizedDescription) }
  }
  let residuals = (input.residuals ?? []).map { v -> Result<VivoOmicsNBAdjustedResiduals> in
   do { return .init(value: try VivoOmicsNBResidualAdjustment.adjustedResiduals(counts: v.counts,means: v.means,design: v.design,dispersion: v.dispersion,averageQuasiDispersion: v.averageQuasiDispersion,method: .adaptive),error: nil) }
   catch { return .init(value: nil,error: error.localizedDescription) }
  }
  let encoder = JSONEncoder();encoder.outputFormatting = [.sortedKeys]
  FileHandle.standardOutput.write(try encoder.encode(Output(moments: moments,residuals: residuals)))
 }
}
