import Foundation
struct SmoothInput: Decodable { let x: [Double], y: [Double], weights: [Double], span: Double }
struct TailInput: Decodable { let logStatistic: Double, numeratorDF: Double, denominatorDF: Double }
struct Input: Decodable { let smoothers: [SmoothInput], shapes: [Double], tails: [TailInput] }
struct SmoothOutput: Encodable { let fit: VivoOmicsPrecisionLowessFit?; let error: String? }
struct TailOutput: Encodable { let fit: VivoOmicsFLogTails?; let error: String? }
struct ShapeOutput: Encodable { let logMinusDigamma: Double, trigamma: Double }
struct Output: Encodable { let smoothers: [SmoothOutput], shapes: [ShapeOutput], tails: [TailOutput] }
@main struct Main {
    static func main() throws {
        let input = try JSONDecoder().decode(Input.self,from: FileHandle.standardInput.readDataToEndOfFile())
        let smoothers = input.smoothers.map { row -> SmoothOutput in
            do { return .init(fit: try VivoOmicsPrecisionLowess.fit(x: row.x,y: row.y,weights: row.weights,span: row.span),error: nil) }
            catch { return .init(fit: nil,error: error.localizedDescription) }
        }
        let shapes = try input.shapes.map { ShapeOutput(logMinusDigamma: try VivoOmicsQLSpecialFunctions.logMinusDigamma($0),
                                                        trigamma: try VivoOmicsQLSpecialFunctions.trigamma($0)) }
        let tails = input.tails.map { row -> TailOutput in
            do { return .init(fit: try VivoOmicsQLSpecialFunctions.logFTails(logStatistic: row.logStatistic,
                numeratorDF: row.numeratorDF,denominatorDF: row.denominatorDF),error: nil) }
            catch { return .init(fit: nil,error: error.localizedDescription) }
        }
        let encoder = JSONEncoder(); encoder.outputFormatting = [.sortedKeys]
        FileHandle.standardOutput.write(try encoder.encode(Output(smoothers: smoothers,shapes: shapes,tails: tails)))
    }
}
