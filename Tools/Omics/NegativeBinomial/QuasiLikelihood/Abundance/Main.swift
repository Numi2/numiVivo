import Foundation
struct Row: Decodable {
    let counts: [UInt64], offsets: [Double], dispersion: Double
    let priorCount: Double?
    let maximumScoreEvaluations: Int?
}
struct Input: Decodable {
    let rows: [Row]?
    let counts: [[UInt64]]?, design: [[Double]]?, offsets: [Double]?, contrast: [Double]?, trendDispersions: [Double]?
}
struct RowOutput: Encodable { let fit: VivoOmicsNBAbundanceFit?; let error: String? }
struct Output: Encodable { let rows: [RowOutput]?; let fit: VivoOmicsNBQLNativeAbundanceFit?; let error: String? }
@main struct Main {
    static func main() throws {
        let input = try JSONDecoder().decode(Input.self,from: FileHandle.standardInput.readDataToEndOfFile())
        let output: Output
        if let rows = input.rows {
            output = .init(rows: rows.map { row in
                do { return RowOutput(fit: try VivoOmicsNBAbundance.fit(counts: row.counts,offsets: row.offsets,
                    dispersion: row.dispersion,priorCount: row.priorCount ?? 2,
                    maximumScoreEvaluations: row.maximumScoreEvaluations ?? 128),error: nil) }
                catch { return RowOutput(fit: nil,error: error.localizedDescription) }
            },fit: nil,error: nil)
        } else {
            do {
                guard let counts = input.counts, let design = input.design, let offsets = input.offsets,
                      let contrast = input.contrast, let dispersions = input.trendDispersions else {
                    throw VivoOmicsStatisticsError.invalid("missing native-abundance family input")
                }
                output = .init(rows: nil,fit: try VivoOmicsNBQLGlobalScale.fitEstimatingAbundance(
                    counts: counts,design: design,offsets: offsets,contrast: contrast,trendDispersions: dispersions),error: nil)
            } catch { output = .init(rows: nil,fit: nil,error: error.localizedDescription) }
        }
        let encoder = JSONEncoder(); encoder.outputFormatting = [.sortedKeys]
        FileHandle.standardOutput.write(try encoder.encode(output))
    }
}
