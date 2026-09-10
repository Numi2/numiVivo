import Foundation
struct Input: Decodable {
    let featureIndices: [Int]
    let counts: [[UInt64]], design: [[Double]], offsets: [Double], contrast: [Double], trendDispersions: [Double]
    let expectedRefits: [VivoOmicsNBFit]?
    let expectedModeration: VivoOmicsQLModeratedVariances?
    let expectedResidualDF: [Double]?
    let expectedAbundance: [Double]?
    let maximumNullIterations: Int?
}
struct Output: Encodable {
    let featureIndices: [Int]
    let completed: Bool
    let inference: VivoOmicsNBQLTestFamily?
    let failures: [VivoOmicsNBQLFailure]
    let refitsExactlyEqual: Bool?
    let residualDFExactlyEqual: Bool?
    let abundanceExactlyEqual: Bool?
    let moderationExactlyEqual: Bool?
    let averageQLDispersion: Double?
    let seconds: Double
}
@main struct Main {
    static func main() throws {
        let input = try JSONDecoder().decode(Input.self,from: FileHandle.standardInput.readDataToEndOfFile())
        guard input.featureIndices.count == input.counts.count, Set(input.featureIndices).count == input.counts.count else {
            throw VivoOmicsStatisticsError.invalid("QL test feature identity")
        }
        let start = Date()
        let fit = try VivoOmicsNBQLInference.fit(counts: input.counts,design: input.design,offsets: input.offsets,
            contrast: input.contrast,trendDispersions: input.trendDispersions,
            maximumNullIterations: input.maximumNullIterations ?? 100)
        let global = fit.native.globalFit
        let output = Output(featureIndices: input.featureIndices,completed: fit.completed,inference: fit.inference,failures: fit.failures,
            refitsExactlyEqual: input.expectedRefits.map { global?.refittedFits == $0.map(Optional.some) },
            residualDFExactlyEqual: input.expectedResidualDF.map { global?.adjustedResiduals.map { $0?.degreesOfFreedom } == $0.map(Optional.some) },
            abundanceExactlyEqual: input.expectedAbundance.map { fit.native.abundanceFits.map { $0?.log2CountsPerMillion } == $0.map(Optional.some) },
            moderationExactlyEqual: input.expectedModeration.map { fit.moderation == $0 },
            averageQLDispersion: global?.averageQuasiDispersion,seconds: Date().timeIntervalSince(start))
        let encoder = JSONEncoder(); encoder.outputFormatting = [.sortedKeys]
        FileHandle.standardOutput.write(try encoder.encode(output))
    }
}
