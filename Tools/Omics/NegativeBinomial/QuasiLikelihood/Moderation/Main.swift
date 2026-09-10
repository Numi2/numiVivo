import Foundation

struct Input: Decodable {
    let featureIndices: [Int]?
    let variances: [Double]
    let degreesOfFreedom: [Double]
    let abundance: [Double]?
    let robust: Bool?
    let maximumProfileEvaluations: Int?
    let maximumFeatureEvaluations: Int?
    let maximumNeighborhoodVisits: Int?
}
struct Output: Encodable {
    let featureIndices: [Int]
    let fit: VivoOmicsQLModeratedVariances?
    let error: String?
    let seconds: Double
}
@main struct Main {
    static func main() throws {
        let input = try JSONDecoder().decode(Input.self, from: FileHandle.standardInput.readDataToEndOfFile())
        let indices = input.featureIndices ?? Array(input.variances.indices)
        guard indices.count == input.variances.count, Set(indices).count == indices.count else {
            throw VivoOmicsStatisticsError.invalid("moderation feature identity length or duplicates")
        }
        let start = Date(), output: Output
        do {
            let fit = try VivoOmicsQLModeration.fit(variances: input.variances,
                degreesOfFreedom: input.degreesOfFreedom, abundance: input.abundance,
                robust: input.robust ?? true,
                maximumProfileEvaluations: input.maximumProfileEvaluations ?? 128,
                maximumFeatureEvaluations: input.maximumFeatureEvaluations ?? 100_000_000,
                maximumNeighborhoodVisits: input.maximumNeighborhoodVisits ?? 100_000_000)
            output = .init(featureIndices: indices, fit: fit, error: nil, seconds: Date().timeIntervalSince(start))
        } catch {
            output = .init(featureIndices: indices, fit: nil, error: String(describing: error), seconds: Date().timeIntervalSince(start))
        }
        let encoder = JSONEncoder(); encoder.outputFormatting = [.sortedKeys]
        FileHandle.standardOutput.write(try encoder.encode(output))
    }
}
