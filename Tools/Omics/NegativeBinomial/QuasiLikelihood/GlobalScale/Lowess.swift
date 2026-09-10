import Foundation
struct Input: Decodable { let x: [Double]; let y: [Double]; let span: Double?; let robustnessIterations: Int?; let deltaFraction: Double? }
struct Output: Encodable { let fit: VivoOmicsLowessFit?; let error: String? }
@main struct Main {
 static func main() throws {
  let cases = try JSONDecoder().decode([Input].self,from: FileHandle.standardInput.readDataToEndOfFile())
  let output = cases.map { i -> Output in
   do { return .init(fit: try VivoOmicsRobustLowess.fit(x: i.x,y: i.y,span: i.span ?? 0.5,robustnessIterations: i.robustnessIterations ?? 3,deltaFraction: i.deltaFraction ?? 0.01),error: nil) }
   catch { return .init(fit: nil,error: error.localizedDescription) }
  }
  let encoder = JSONEncoder();encoder.outputFormatting = [.sortedKeys]
  FileHandle.standardOutput.write(try encoder.encode(output))
 }
}
