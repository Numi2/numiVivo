import Foundation

/// Imported, frozen gene-space reference. Scores are uncalibrated annotations,
/// never authoritative cell identities. No fitting or majority voting occurs here.
public struct VivoCellTypistReferenceModel: Codable, Sendable {
    public let sourceSHA256: String
    public let classes: [String]
    public let features: [String]
    public let coefficients: [[Double]]
    public let intercepts: [Double]
    public let means: [Double]
    public let scales: [Double]
    public let withMean: Bool
    public let withStd: Bool
    public let multiClass: String
}

public struct VivoCellTypistReferenceResult: Sendable {
    public let label: String
    public let decisions: [Double]
    /// Independent sigmoid scores; these do not sum to one.
    public let probabilities: [Double]
}

public enum VivoCellTypistReferenceError: Error {
    case invalidModel, invalidAxis, missingFeatures, invalidCounts, nonfiniteScore
}

/// Bounded per-row inference: O(source features + model parameters) retained state,
/// O(model features + sparse row) work before the class-by-feature dot products.
/// Full source-row totals determine normalization, including non-model genes.
public struct VivoCellTypistReference: Sendable {
    private let model: VivoCellTypistReferenceModel
    private let sourceToModel: [Int]
    private let zeros: [Double]

    public init(model: VivoCellTypistReferenceModel, sourceFeatures: [String]) throws {
        let g = model.features.count, c = model.classes.count
        guard (1...100_000).contains(g), (2...256).contains(c),
              g * c <= 10_000_000,
              model.multiClass == "ovr", model.withStd,
              model.sourceSHA256.count == 64,
              model.sourceSHA256.allSatisfy({ "0123456789abcdef".contains($0) }),
              Set(model.features).count == g, Set(model.classes).count == c,
              model.features.allSatisfy({ !$0.isEmpty }), model.classes.allSatisfy({ !$0.isEmpty }),
              model.coefficients.count == c, model.intercepts.count == c,
              model.means.count == g, model.scales.count == g,
              model.coefficients.allSatisfy({ $0.count == g && $0.allSatisfy(\.isFinite) }),
              model.intercepts.allSatisfy(\.isFinite), model.means.allSatisfy(\.isFinite),
              model.scales.allSatisfy({ $0.isFinite && $0 > 0 }) else {
            throw VivoCellTypistReferenceError.invalidModel
        }
        guard !sourceFeatures.isEmpty, sourceFeatures.count <= 1_000_000,
              Set(sourceFeatures).count == sourceFeatures.count,
              sourceFeatures.allSatisfy({ !$0.isEmpty }) else {
            throw VivoCellTypistReferenceError.invalidAxis
        }
        let indices = Dictionary(uniqueKeysWithValues: model.features.enumerated().map { ($1, $0) })
        let mapping = sourceFeatures.map { indices[$0] ?? -1 }
        // This full-panel path deliberately rejects missing genes. It does not
        // silently emulate CellTypist's optional intersection-panel behavior.
        guard mapping.filter({ $0 >= 0 }).count == g else {
            throw VivoCellTypistReferenceError.missingFeatures
        }
        self.model = model
        self.sourceToModel = mapping
        self.zeros = (0..<g).map { min(10, (model.withMean ? -model.means[$0] : 0) / model.scales[$0]) }
        guard zeros.allSatisfy(\.isFinite) else { throw VivoCellTypistReferenceError.invalidModel }
    }

    /// Canonical sparse row: strictly increasing source coordinates, positive UMI counts.
    /// An empty row is retained and normalized to all-zero log expression.
    public func predict(indices: [Int], counts: [UInt64]) throws -> VivoCellTypistReferenceResult {
        guard indices.count == counts.count else { throw VivoCellTypistReferenceError.invalidCounts }
        var total: UInt64 = 0, previous = -1
        for (index, count) in zip(indices, counts) {
            let sum = total.addingReportingOverflow(count)
            guard index > previous, index < sourceToModel.count, count > 0, !sum.overflow else {
                throw VivoCellTypistReferenceError.invalidCounts
            }
            previous = index; total = sum.partialValue
        }
        var values = zeros
        if total > 0 {
            for (index, count) in zip(indices, counts) {
                let j = sourceToModel[index]
                if j >= 0 {
                    let expression = log1p(Double(count) / Double(total) * 10_000)
                    values[j] = min(10, (expression - (model.withMean ? model.means[j] : 0)) / model.scales[j])
                }
            }
        }
        var decisions = model.intercepts
        for c in decisions.indices {
            for j in values.indices { decisions[c] += values[j] * model.coefficients[c][j] }
        }
        guard decisions.allSatisfy(\.isFinite) else { throw VivoCellTypistReferenceError.nonfiniteScore }
        let probabilities = decisions.map { x in
            x >= 0 ? 1 / (1 + exp(-x)) : exp(x) / (1 + exp(x))
        }
        var best = 0
        for c in decisions.indices where decisions[c] > decisions[best] { best = c }
        return VivoCellTypistReferenceResult(label: model.classes[best], decisions: decisions, probabilities: probabilities)
    }
}
