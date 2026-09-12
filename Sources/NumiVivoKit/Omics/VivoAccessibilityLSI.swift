import Foundation

public struct VivoAccessibilityLSIResult: Codable, Sendable {
    public let method: String
    public let cells: [VivoOmicsCellIdentity]
    public let featureIDs: [String]
    public let singularValues: [Double]
    public let leftSingularVectors: [[Double]]
    public let standardizedEmbeddings: [[Double]]
    public let featureLoadings: [[Double]]
    public let relativeResiduals: [Double]
    public let maximumLoadingOrthogonalityError: Double
    public let retainedEnergyFraction: Double
    public let zeroRows: [Int]
    public let basisSize: Int
    public let operatorScans: Int
}

/// Uncentered sparse SVD of a frozen TF-IDF matrix. The caller owns its immutable
/// snapshot and identity binding; every scan must replay the same sparse values.
/// Reuses the existing residual-qualified Krylov solver with zero centers.
public enum VivoAccessibilityLSI {
    public static func fit(cells: [VivoOmicsCellIdentity], featureIDs: [String],
                           components: Int = 30, maximumBasis: Int = 256,
                           relativeResidualTolerance: Double = 1e-5, seed: UInt64 = 7,
                           scan: @escaping (_ accept: (Int, Int, Double) throws -> Void) throws -> Void) throws -> VivoAccessibilityLSIResult {
        let n = cells.count, width = featureIDs.count
        guard n > components, n <= 10_000_000, width >= components, width <= 200_000,
              Set(cells).count == n else { throw VivoOmicsError.invalid("LSI axes") }
        var options = VivoSingleCellReductionOptions()
        options.components = components; options.maximumBasis = maximumBasis
        options.relativeResidualTolerance = relativeResidualTolerance; options.seed = seed
        options.featurePanel = featureIDs
        try options.validate()
        var energy = 0.0, records = 0
        var rowEnergy = [Double](repeating: 0, count: n)
        var previous = [Int](repeating: -1, count: n)
        try scan { row, feature, value in
            guard row >= 0, row < n, feature >= 0, feature < width,
                  feature > previous[row], value.isFinite, value > 0 else {
                throw VivoOmicsError.invalid("LSI canonical sparse value")
            }
            previous[row] = feature; records += 1
            energy += value * value; rowEnergy[row] += value * value
        }
        guard energy.isFinite, energy > 0 else { throw VivoOmicsError.invalid("LSI zero or nonfinite energy") }
        var scans = 0
        let result = try VivoSingleCellReduction.fit(cells: cells, statistics: [], selected: Array(0..<width),
            centers: [Double](repeating: 0, count: width), totalVariance: energy / Double(n - 1), options: options,
            project: { vector, offset in
                guard offset == 0 else { throw VivoOmicsError.invalid("LSI must remain uncentered") }
                var output = [Double](repeating: 0, count: n), seen = 0
                try scan { row, feature, value in
                    guard row >= 0, row < n, feature >= 0, feature < width,
                          value.isFinite, value > 0, seen < records else { throw VivoOmicsError.invalid("LSI replay") }
                    output[row] += value * vector[feature]; seen += 1
                }
                guard seen == records else { throw VivoOmicsError.invalid("LSI replay cardinality") }
                scans += 1; return output
            }, transpose: { vector, initial in
                var output = initial, seen = 0
                try scan { row, feature, value in
                    guard row >= 0, row < n, feature >= 0, feature < width,
                          value.isFinite, value > 0, seen < records else { throw VivoOmicsError.invalid("LSI replay") }
                    output[feature] += value * vector[row]; seen += 1
                }
                guard seen == records else { throw VivoOmicsError.invalid("LSI replay cardinality") }
                scans += 1; return output
            })
        let singular = result.explainedVariance.map { sqrt($0 * Double(n - 1)) }
        var vectors = result.scores
        for row in 0..<n { for k in 0..<components { vectors[row][k] /= singular[k] } }
        var standardized = vectors
        for k in 0..<components {
            // Welford preserves exactly constant columns: summing then dividing
            // can round their mean away from the repeated value and amplify
            // that rounding residue into a spurious standardized component.
            var mean = 0.0, m2 = 0.0
            for row in 0..<n {
                let delta = vectors[row][k] - mean
                mean += delta / Double(row + 1)
                m2 += delta * (vectors[row][k] - mean)
            }
            let variance = m2 / Double(n - 1)
            guard variance.isFinite, variance > 0 else { throw VivoOmicsError.invalid("LSI constant embedding") }
            let deviation = sqrt(variance)
            for row in 0..<n { standardized[row][k] = (vectors[row][k] - mean) / deviation }
        }
        return .init(method: "uncentered-tfidf-krylov-SVD-v1", cells: cells, featureIDs: featureIDs,
            singularValues: singular, leftSingularVectors: vectors, standardizedEmbeddings: standardized,
            featureLoadings: result.loadings, relativeResiduals: result.relativeResiduals,
            maximumLoadingOrthogonalityError: result.maximumLoadingOrthogonalityError,
            retainedEnergyFraction: result.explainedVarianceRatio.reduce(0, +),
            zeroRows: rowEnergy.indices.filter { rowEnergy[$0] == 0 }, basisSize: result.basisSize, operatorScans: scans)
    }
}
