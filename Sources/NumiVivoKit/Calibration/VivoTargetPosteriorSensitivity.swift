import Foundation

public struct VivoPosteriorSensitivityPolicy: Codable, Equatable, Sendable {
    public let priorCoordinateStep: Double
    public let relativeRankThreshold: Double
    public let maximumDerivativeChange: Double
    public init(priorCoordinateStep: Double = 1e-3, relativeRankThreshold: Double = 1e-8,
                maximumDerivativeChange: Double = 0.05) {
        self.priorCoordinateStep = priorCoordinateStep; self.relativeRankThreshold = relativeRankThreshold
        self.maximumDerivativeChange = maximumDerivativeChange
    }
    public func validate() throws {
        guard priorCoordinateStep.isFinite, (1e-6...0.05).contains(priorCoordinateStep),
              relativeRankThreshold.isFinite, (1e-12...1e-2).contains(relativeRankThreshold),
              maximumDerivativeChange.isFinite, (1e-6...0.5).contains(maximumDerivativeChange) else {
            throw VivoPosteriorError.invalid("local sensitivity policy")
        }
    }
}

public struct VivoLocalWeakParameterDirection: Codable, Equatable, Sendable {
    public let relativeInformation: Double
    public let priorCoordinateCoefficients: [Double]
}

public struct VivoTargetPosteriorSensitivityReport: Codable, Equatable, Sendable {
    public let schemaVersion: UInt32
    public let posteriorFingerprint: VivoFingerprint
    public let policy: VivoPosteriorSensitivityPolicy
    public let parameterOrder: [String]
    public let evaluationPriorCoordinates: [Double]
    public let calibrationObservationCount: Int
    public let whitenedColumnNorms: [Double]
    public let relativeDerivativeChanges: [Double]
    public let fisherInformationRowMajor: [Double]
    public let eigenvaluesDescending: [Double]
    public let numericalRank: Int?
    public let weakDirections: [VivoLocalWeakParameterDirection]
    public let derivativeChecksPassed: Bool
    public let limitations: [String]
}

/// Local likelihood information J'J in prior coordinates, NOT a global or
/// structural-identifiability proof and not a posterior-convergence diagnostic.
/// Uses h and h/2 perturbations; refuses to report numerical rank when local
/// derivative estimates fail the declared consistency threshold.
public enum VivoTargetPosteriorSensitivity {
    public static func evaluate(_ record: VivoTargetPosteriorRecord,
                                policy: VivoPosteriorSensitivityPolicy = .init()) throws -> VivoTargetPosteriorSensitivityReport {
        try record.validate(requireComplete: true); try policy.validate()
        guard !record.problem.usesExtendedTrainingAssays else {
            throw VivoPosteriorError.invalid("mean-only J'J sensitivity does not represent fitted/correlated/censored assay information; covariance-aware information is required")
        }
        let prepared = try VivoPreparedTargetPosterior(record.problem)
        guard let checkpoint = record.posterior.checkpoint else { throw VivoPosteriorError.invalid("missing posterior") }
        let d = prepared.plan.parameters.count, n = Double(checkpoint.particles.count)
        var center = [Double](repeating: 0, count: d)
        for particle in checkpoint.particles { for j in 0..<d { center[j] += particle.coordinates[j] / n } }
        let observations = prepared.calibrationCases.flatMap(\.observations)
        let scales = try observations.map { observation -> Double in
            guard let sd = observation.standardDeviation, sd > 0 else { throw VivoPosteriorError.invalid("measurement SD") }
            return sd
        }
        let m = scales.count
        var columns: [[Double]] = [], norms: [Double] = [], changes: [Double] = []
        for parameter in 0..<d {
            try Task.checkCancellation()
            var estimates: [[Double]] = []
            for h in [policy.priorCoordinateStep, policy.priorCoordinateStep * 0.5] {
                var lower = center, upper = center
                lower[parameter] = max(0, center[parameter] - h)
                upper[parameter] = min(1, center[parameter] + h)
                let width = upper[parameter] - lower[parameter]
                guard width > 0 else { throw VivoPosteriorError.numerical("unresolvable derivative interval") }
                let low = try predictions(lower, prepared: prepared), high = try predictions(upper, prepared: prepared)
                guard low.count == m, high.count == m else { throw VivoPosteriorError.numerical("derivative observation shape") }
                let estimate = (0..<m).map { (high[$0] - low[$0]) / width / scales[$0] }
                guard estimate.allSatisfy({ $0.isFinite && abs($0) <= 1e140 }) else {
                    throw VivoPosteriorError.numerical("whitened derivative outside supported numerical range")
                }
                estimates.append(estimate)
            }
            let norm = sqrt(estimates[1].reduce(0) { $0 + $1 * $1 })
            let otherNorm = sqrt(estimates[0].reduce(0) { $0 + $1 * $1 })
            let difference = sqrt(zip(estimates[0], estimates[1]).reduce(0) {
                $0 + ($1.0 - $1.1) * ($1.0 - $1.1)
            })
            let denominator = max(norm, otherNorm)
            changes.append(denominator > 0 ? difference / denominator : 0)
            norms.append(norm); columns.append(estimates[1])
        }
        var fisher = [Double](repeating: 0, count: d * d)
        for i in 0..<d { for j in 0...i {
            let product = zip(columns[i], columns[j]).reduce(0) { $0 + $1.0 * $1.1 }
            guard product.isFinite else { throw VivoPosteriorError.numerical("information matrix overflow") }
            fisher[i * d + j] = product; fisher[j * d + i] = product
        } }
        let (eigenvalues, vectors) = try symmetricSpectrum(fisher, dimension: d)
        let maximum = eigenvalues.first ?? 0
        let stable = changes.allSatisfy { $0 <= policy.maximumDerivativeChange }
        let threshold = maximum * policy.relativeRankThreshold
        let rank = stable ? eigenvalues.filter { $0 > threshold }.count : nil
        let weak = stable ? eigenvalues.indices.filter { eigenvalues[$0] <= threshold }.map { index in
            VivoLocalWeakParameterDirection(relativeInformation: maximum > 0 ? eigenvalues[index] / maximum : 0,
                                           priorCoordinateCoefficients: vectors[index])
        } : []
        return try .init(schemaVersion: 1,
            posteriorFingerprint: VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(record)), policy: policy,
            parameterOrder: prepared.plan.parameters.map(\.identifier), evaluationPriorCoordinates: center,
            calibrationObservationCount: m, whitenedColumnNorms: norms, relativeDerivativeChanges: changes,
            fisherInformationRowMajor: fisher, eigenvaluesDescending: eigenvalues, numericalRank: rank,
            weakDirections: weak, derivativeChecksPassed: stable, limitations: [
                "Local Gaussian-likelihood information at the posterior mean in prior coordinates; other regions may differ.",
                "A full local rank does not establish global identifiability, adequate posterior exploration or biological correctness.",
                "Weak directions combine parameter changes in the declared prior coordinates; logarithmic and linear priors have different scales.",
                "Finite-difference consistency at h and h/2 is a numerical diagnostic, not a bound on forward-model or sampling error.",
                "Training observations only. Held-out outcomes do not select parameters or improve this information matrix.",
                "Baseline target abundance is excluded from occupancy-only fitting because this reservoir fraction model does not identify it."
            ])
    }

    private static func predictions(_ coordinates: [Double], prepared: VivoPreparedTargetPosterior) throws -> [Double] {
        let values = try zip(prepared.plan.parameters, coordinates).map { try $0.value(at: $1) }
        var output: [Double] = []
        for item in prepared.calibrationCases {
            try Task.checkCancellation()
            let experiment = try prepared.experiment(for: item, values: values)
            let result = try VivoTargetEngagementReference.run(experiment, numerics: prepared.problem.numerics)
            let samples = Dictionary(uniqueKeysWithValues: result.samples.map { ($0.timeSeconds, $0) })
            for observation in item.observations {
                guard let sample = samples[observation.timeSeconds] else { throw VivoPosteriorError.numerical("missing sensitivity sample") }
                output.append(observation.observable.value(sample))
            }
        }
        return output
    }

    private static func symmetricSpectrum(_ matrix: [Double], dimension d: Int) throws -> ([Double], [[Double]]) {
        let scale = matrix.map(abs).max() ?? 0
        var vectors = [Double](repeating: 0, count: d * d)
        for i in 0..<d { vectors[i * d + i] = 1 }
        if scale == 0 { return ([Double](repeating: 0, count: d), (0..<d).map { j in (0..<d).map { vectors[$0 * d + j] } }) }
        var a = matrix.map { $0 / scale }, converged = d == 1
        for _ in 0..<(100 * d * d) {
            var p = 0, q = 0, largest = 0.0
            for i in 0..<d { for j in 0..<i {
                if abs(a[i * d + j]) > largest { largest = abs(a[i * d + j]); p = j; q = i }
            } }
            if largest <= 1e-13 { converged = true; break }
            let app = a[p * d + p], aqq = a[q * d + q], apq = a[p * d + q]
            let angle = 0.5 * atan2(2 * apq, aqq - app), c = cos(angle), s = sin(angle)
            for k in 0..<d where k != p && k != q {
                let akp = a[k * d + p], akq = a[k * d + q]
                a[k * d + p] = c * akp - s * akq; a[p * d + k] = a[k * d + p]
                a[k * d + q] = s * akp + c * akq; a[q * d + k] = a[k * d + q]
            }
            a[p * d + p] = c * c * app - 2 * s * c * apq + s * s * aqq
            a[q * d + q] = s * s * app + 2 * s * c * apq + c * c * aqq
            a[p * d + q] = 0; a[q * d + p] = 0
            for k in 0..<d {
                let vkp = vectors[k * d + p], vkq = vectors[k * d + q]
                vectors[k * d + p] = c * vkp - s * vkq
                vectors[k * d + q] = s * vkp + c * vkq
            }
        }
        guard converged else { throw VivoPosteriorError.numerical("local information eigensolver did not converge") }
        let order = (0..<d).sorted { a[$0 * d + $0] > a[$1 * d + $1] }
        let values = try order.map { index -> Double in
            let value = a[index * d + index]
            guard value.isFinite, value >= -1e-10 else { throw VivoPosteriorError.numerical("negative information eigenvalue") }
            // J'J is positive semidefinite by construction; only tiny eigensolver
            // roundoff below zero is replaced with zero, not biological state.
            return max(0, value) * scale
        }
        return (values, order.map { j in (0..<d).map { vectors[$0 * d + j] } })
    }
}
