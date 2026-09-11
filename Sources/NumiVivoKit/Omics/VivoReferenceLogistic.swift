import Foundation

public struct VivoReferenceLogisticOptions: Codable, Sendable, Equatable {
    public var penalty: Double = 1
    public var gradientTolerance: Double = 1e-7
    public var maximumIterations: Int = 2000
    public var maximumWork: Int = 20_000_000_000
    public init() {}
    private enum CodingKeys: String, CodingKey { case penalty, gradientTolerance, maximumIterations, maximumWork }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["penalty", "gradientTolerance", "maximumIterations", "maximumWork"])
        let c = try decoder.container(keyedBy: CodingKeys.self)
        penalty = try c.decodeIfPresent(Double.self, forKey: .penalty) ?? 1
        gradientTolerance = try c.decodeIfPresent(Double.self, forKey: .gradientTolerance) ?? 1e-7
        maximumIterations = try c.decodeIfPresent(Int.self, forKey: .maximumIterations) ?? 2000
        maximumWork = try c.decodeIfPresent(Int.self, forKey: .maximumWork) ?? 20_000_000_000
    }
    func validate() throws {
        guard penalty.isFinite, (1e-6...1e6).contains(penalty), gradientTolerance.isFinite,
              (1e-10...1e-4).contains(gradientTolerance), (1...10_000).contains(maximumIterations),
              (1...100_000_000_000).contains(maximumWork) else { throw VivoOmicsError.invalid("reference logistic options") }
    }
}

public struct VivoReferenceLogisticModel: Codable, Sendable, Equatable {
    public let method: String
    public let means: [Double]
    public let scales: [Double]
    public let classCounts: [Int]
    public let classWeights: [Double]
    /// Class-major standardized coefficients; last entry of each row is intercept.
    public let parameters: [[Double]]
    public let objectives: [Double]
    public let gradientMaximum: Double
    public let iterations: Int
    public let evaluations: Int
    public let chargedWork: Int
}

/// Convex weighted multinomial likelihood in frozen training PCA coordinates.
/// Dense state is cells-by-components, never cells-by-genes or cells-by-cells.
public enum VivoReferenceLogistic {
    static func dot(_ a: [Double], _ b: [Double]) -> Double {
        var value = 0.0
        for j in a.indices { value += a[j] * b[j] }
        return value
    }
    static func softmax(_ logits: [Double]) throws -> (values: [Double], logSum: Double) {
        guard let maximum = logits.max(), maximum.isFinite else { throw VivoOmicsError.invalid("logistic logits") }
        let values = logits.map { exp($0 - maximum) }, sum = values.reduce(0, +)
        guard sum.isFinite, sum > 0 else { throw VivoOmicsError.invalid("logistic normalization") }
        return (values.map { $0 / sum }, maximum + log(sum))
    }
    /// Weighted mean negative log likelihood plus penalty/(2n) ||coefficients||².
    /// No penalty on intercepts. Used by the optimizer and independent derivative tests.
    static func objective(_ parameters: [Double], x: [[Double]], labels: [Int], classes: Int,
                          weights: [Double], penalty: Double) throws -> (value: Double, gradient: [Double]) {
        let n = x.count, d = x[0].count, width = d + 1
        var loss = 0.0, gradient = [Double](repeating: 0, count: parameters.count)
        for i in 0..<n {
            if i % 1024 == 0 { try Task.checkCancellation() }
            var logits = [Double](repeating: 0, count: classes)
            for c in 0..<classes {
                let base = c * width
                logits[c] = parameters[base + d]
                for j in 0..<d { logits[c] += parameters[base + j] * x[i][j] }
            }
            let probability = try softmax(logits), weight = weights[labels[i]] / Double(n)
            loss += weight * (probability.logSum - logits[labels[i]])
            for c in 0..<classes {
                let error = weight * (probability.values[c] - (labels[i] == c ? 1 : 0)), base = c * width
                for j in 0..<d { gradient[base + j] += error * x[i][j] }
                gradient[base + d] += error
            }
        }
        for c in 0..<classes { for j in 0..<d {
            let index = c * width + j, scale = penalty / Double(n)
            loss += 0.5 * scale * parameters[index] * parameters[index]
            gradient[index] += scale * parameters[index]
        } }
        guard loss.isFinite, gradient.allSatisfy(\.isFinite) else { throw VivoOmicsError.invalid("nonfinite logistic objective") }
        return (loss, gradient)
    }
    public static func fit(scores: [[Double]], labels: [Int], classes: Int, options: VivoReferenceLogisticOptions) throws -> VivoReferenceLogisticModel {
        try options.validate()
        let n = scores.count, d = scores.first?.count ?? 0
        guard (2...100_000).contains(n), (1...64).contains(d), (2...256).contains(classes), labels.count == n,
              scores.allSatisfy({ $0.count == d && $0.allSatisfy(\.isFinite) }), labels.allSatisfy({ (0..<classes).contains($0) }) else {
            throw VivoOmicsError.invalid("reference logistic training axes")
        }
        var counts = [Int](repeating: 0, count: classes), means = [Double](repeating: 0, count: d)
        for i in 0..<n {
            counts[labels[i]] += 1
            for j in 0..<d { means[j] += scores[i][j] }
        }
        for j in means.indices { means[j] /= Double(n) }
        guard counts.allSatisfy({ $0 > 0 }) else { throw VivoOmicsError.invalid("reference logistic class absent") }
        var variances = [Double](repeating: 0, count: d)
        for row in scores { for j in 0..<d { let delta = row[j] - means[j]; variances[j] += delta * delta / Double(n) } }
        // Treat variance below the floating-point error bound as constant, as
        // StandardScaler does; do not amplify cancellation into huge weights.
        let scales = variances.indices.map { j in
            let meanError = Double(n) * means[j] * Double.ulpOfOne
            let bound = Double(n) * Double.ulpOfOne * variances[j] + meanError * meanError
            return variances[j] <= bound ? 1 : sqrt(variances[j])
        }
        guard means.allSatisfy(\.isFinite), scales.allSatisfy({ $0.isFinite && $0 > 0 }) else { throw VivoOmicsError.invalid("reference logistic scaling") }
        let x = scores.map { row in row.indices.map { (row[$0] - means[$0]) / scales[$0] } }
        let weights = counts.map { Double(n) / (Double(classes) * Double($0)) }, width = d + 1
        let evaluationWork = 2 * n * classes * width
        var evaluations = 0, work = 0
        func evaluate(_ p: [Double]) throws -> (value: Double, gradient: [Double]) {
            guard evaluationWork <= options.maximumWork - work else { throw VivoOmicsError.limit("reference logistic training work") }
            work += evaluationWork; evaluations += 1
            return try objective(p, x: x, labels: labels, classes: classes, weights: weights, penalty: options.penalty)
        }
        var parameters = [Double](repeating: 0, count: classes * width), state = try evaluate(parameters)
        var objectives = [state.value], historyS: [[Double]] = [], historyY: [[Double]] = [], rho: [Double] = []
        var iteration = 0
        while state.gradient.map(abs).max()! > options.gradientTolerance {
            try Task.checkCancellation()
            guard iteration < options.maximumIterations else { throw VivoOmicsError.invalid("reference logistic did not converge before iteration limit") }
            var direction = state.gradient, alpha = [Double](repeating: 0, count: rho.count)
            for i in rho.indices.reversed() {
                alpha[i] = rho[i] * dot(historyS[i], direction)
                for j in direction.indices { direction[j] -= alpha[i] * historyY[i][j] }
            }
            if let y = historyY.last, let s = historyS.last {
                let scale = dot(s, y) / dot(y, y)
                for j in direction.indices { direction[j] *= scale }
            }
            for i in rho.indices {
                let beta = rho[i] * dot(historyY[i], direction)
                for j in direction.indices { direction[j] += historyS[i][j] * (alpha[i] - beta) }
            }
            direction = direction.map { -$0 }
            var slope = dot(state.gradient, direction)
            if !slope.isFinite || slope >= 0 {
                historyS.removeAll(); historyY.removeAll(); rho.removeAll()
                direction = state.gradient.map { -$0 }; slope = -dot(state.gradient, state.gradient)
            }
            var step = 1.0, accepted: (parameters: [Double], value: Double, gradient: [Double])?
            for _ in 0..<50 {
                let candidate = parameters.indices.map { parameters[$0] + step * direction[$0] }
                let proposed = try evaluate(candidate)
                if proposed.value <= state.value + 1e-4 * step * slope {
                    accepted = (candidate, proposed.value, proposed.gradient); break
                }
                step *= 0.5
            }
            guard let accepted else { throw VivoOmicsError.invalid("reference logistic line search did not converge") }
            let s = parameters.indices.map { accepted.parameters[$0] - parameters[$0] }
            let y = parameters.indices.map { accepted.gradient[$0] - state.gradient[$0] }
            let curvature = dot(s, y)
            if curvature > 1e-14 * sqrt(dot(s, s) * dot(y, y)), curvature.isFinite {
                if rho.count == 10 { historyS.removeFirst(); historyY.removeFirst(); rho.removeFirst() }
                historyS.append(s); historyY.append(y); rho.append(1 / curvature)
            }
            parameters = accepted.parameters; state = (accepted.value, accepted.gradient)
            objectives.append(state.value); iteration += 1
        }
        return .init(method: "balanced-multinomial-logistic-L2-training-standardization-LBFGS-v1", means: means, scales: scales,
            classCounts: counts, classWeights: weights, parameters: (0..<classes).map { Array(parameters[($0 * width)..<(($0 + 1) * width)]) },
            objectives: objectives, gradientMaximum: state.gradient.map(abs).max()!, iterations: iteration, evaluations: evaluations, chargedWork: work)
    }
    public static func probabilities(_ score: [Double], model: VivoReferenceLogisticModel) throws -> [Double] {
        let d = model.means.count
        guard score.count == d, model.scales.count == d, score.allSatisfy(\.isFinite), model.means.allSatisfy(\.isFinite),
              model.scales.allSatisfy({ $0.isFinite && $0 > 0 }), (2...256).contains(model.parameters.count),
              model.parameters.allSatisfy({ $0.count == d + 1 && $0.allSatisfy(\.isFinite) }) else { throw VivoOmicsError.invalid("reference logistic query axes or values") }
        let x = score.indices.map { (score[$0] - model.means[$0]) / model.scales[$0] }
        return try softmax(model.parameters.map { row in
            var value = row[d]
            for j in 0..<d { value += row[j] * x[j] }
            return value
        }).values
    }
}
