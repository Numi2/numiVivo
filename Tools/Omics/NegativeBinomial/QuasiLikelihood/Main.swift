import Foundation
struct Row: Decodable { let counts: [UInt64]; let means: [Double]; let dispersion: Double }
struct Result: Encodable { let values: [Double]?; let error: String? }
@main struct Main {
    static func main() throws {
        let rows = try JSONDecoder().decode([Row].self,from: FileHandle.standardInput.readDataToEndOfFile())
        let results = rows.map { row -> Result in
            do {
                guard row.counts.count == row.means.count else { throw VivoOmicsStatisticsError.invalid("count/mean dimensions") }
                return .init(values: try zip(row.counts,row.means).map {
                    try VivoOmicsNegativeBinomial.unitDeviance(count: $0.0,mean: $0.1,dispersion: row.dispersion)
                },error: nil)
            } catch { return .init(values: nil,error: error.localizedDescription) }
        }
        let encoder = JSONEncoder();encoder.outputFormatting = [.sortedKeys]
        FileHandle.standardOutput.write(try encoder.encode(results))
    }
}
