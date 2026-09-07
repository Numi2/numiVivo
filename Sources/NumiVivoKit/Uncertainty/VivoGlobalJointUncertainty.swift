import Foundation

public struct VivoJointProbabilityModel: Codable, Sendable, Equatable {
    public let identifier: String
    /// Explicit mixture probability, not proportional to the number of draws.
    public let probability: Double
    public let evidenceFingerprint: VivoFingerprint
    public let interpretation: String
    public init(identifier: String, probability: Double, evidenceFingerprint: VivoFingerprint, interpretation: String) {
        self.identifier = identifier; self.probability = probability
        self.evidenceFingerprint = evidenceFingerprint; self.interpretation = interpretation
    }
}
public struct VivoJointParameterDraw: Codable, Sendable, Equatable {
    public let identifier: String
    public let modelIdentifier: String
    /// Declared resampling/dependence unit. Repeated states must not be relabeled
    /// independent merely to increase an effective-sample diagnostic.
    public let dependenceBlockIdentifier: String
    public let logWeight: Double
    public let values: [Double]
    public init(identifier: String, modelIdentifier: String, dependenceBlockIdentifier: String,
                logWeight: Double = 0, values: [Double]) {
        self.identifier = identifier; self.modelIdentifier = modelIdentifier
        self.dependenceBlockIdentifier = dependenceBlockIdentifier; self.logWeight = logWeight; self.values = values
    }
}
public struct VivoGlobalJointUncertaintyRequest: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/global-joint-uncertainty-request/v1"
    public var schema: String = Self.schema
    public var contextFingerprint: VivoFingerprint
    public var parameterIdentifiers: [String]
    public var models: [VivoJointProbabilityModel]
    public var draws: [VivoJointParameterDraw]
    /// These remain unresolved in the RESULT, even if sampled variance is zero.
    public var unresolvedComponents: [String]
    public var assumptions: [String]
    public var maximumPrimitiveWork: Int
    public init(contextFingerprint: VivoFingerprint, parameterIdentifiers: [String],
                models: [VivoJointProbabilityModel], draws: [VivoJointParameterDraw],
                unresolvedComponents: [String], assumptions: [String], maximumPrimitiveWork: Int = 100_000_000) {
        self.contextFingerprint = contextFingerprint; self.parameterIdentifiers = parameterIdentifiers
        self.models = models; self.draws = draws; self.unresolvedComponents = unresolvedComponents
        self.assumptions = assumptions; self.maximumPrimitiveWork = maximumPrimitiveWork
    }
    public func validate() throws {
        let ids = models.map(\.identifier), n = parameterIdentifiers.count
        func validName(_ s: String) -> Bool { !s.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty && s.utf8.count <= 4096 }
        guard schema == Self.schema, (1...4096).contains(n), (1...128).contains(models.count),
              (2...65_536).contains(draws.count), maximumPrimitiveWork > 0,
              Set(parameterIdentifiers).count == n, parameterIdentifiers.allSatisfy(validName),
              Set(ids).count == ids.count, models.allSatisfy({ validName($0.identifier) && validName($0.interpretation) && $0.probability.isFinite && $0.probability > 0 }),
              abs(models.reduce(0) { $0+$1.probability }-1) <= 1e-12,
              Set(draws.map(\.identifier)).count == draws.count,
              draws.allSatisfy({ validName($0.identifier) && validName($0.dependenceBlockIdentifier)
                && ids.contains($0.modelIdentifier) && $0.logWeight.isFinite && $0.values.count == n && $0.values.allSatisfy(\.isFinite) }),
              Set(draws.map(\.modelIdentifier)) == Set(ids),
              unresolvedComponents.allSatisfy(validName), assumptions.allSatisfy(validName),
              Double(n)*Double(draws.count) <= Double(maximumPrimitiveWork) else {
            throw VivoChemistryError.invalid("global joint uncertainty dimensions, mixture weights, complete draws or context")
        }
    }
}
public struct VivoJointObservableSummary: Codable, Sendable, Equatable {
    public let identifier: String
    public let definedProbability: Double
    public let mean: Double?
    /// Posterior/predictive distribution variance, NOT the variance of its mean.
    public let variance: Double?
    public let withinModelVariance: Double?
    public let betweenModelVariance: Double?
    public let quantile025: Double?
    public let median: Double?
    public let quantile975: Double?
}
public struct VivoJointObservableDraw: Codable, Sendable, Equatable {
    public let parameterDrawIdentifier: String
    public let modelIdentifier: String
    public let probability: Double
    public let values: [Double?]
}
public struct VivoGlobalJointUncertaintyResult: Codable, Sendable, Equatable {
    public let requestFingerprint: VivoFingerprint
    public let contextFingerprint: VivoFingerprint
    public let summaries: [VivoJointObservableSummary]
    public let draws: [VivoJointObservableDraw]
    /// Unconditional covariance is absent if any observable is undefined in any
    /// positive-weight draw. Pairwise deletion would silently change the measure.
    public let covariance: VivoQMMatrix?
    public let weightEffectiveDraws: Double
    public let declaredBlockEffectiveCount: Double
    public let unresolvedComponents: [String]
    public let assumptions: [String]
    public let chargedPrimitiveWork: Int
    public let evidenceFingerprint: VivoFingerprint
    public let interpretation: String
}
public struct VivoJointForwardEvaluation: Sendable {
    public let values: [Double?]
    public let chargedPrimitiveWork: Int
    public init(values: [Double?], chargedPrimitiveWork: Int) { self.values = values; self.chargedPrimitiveWork = chargedPrimitiveWork }
}

/// Global nonlinear propagation of supplied JOINT samples. No independent
/// marginal resampling, linearization, likelihood fitting or posterior invention.
public enum VivoGlobalJointUncertainty {
    public static let interpretation = "Full nonlinear push-forward of evidence-bound joint draws with explicit model-mixture probabilities. Within/between-model variance uses total variance, and shared draw/block identities preserve supplied dependencies. Weight/block effective counts are diagnostics, not autocorrelation measurements, Monte Carlo error bars, confidence coverage or proof of chemical/model completeness. Undefined observables and unresolved uncertainty sources remain explicit."
    public static func propagate(_ request: VivoGlobalJointUncertaintyRequest, observableIdentifiers: [String],
        evaluate: ([Double], String, Int) throws -> VivoJointForwardEvaluation) throws -> VivoGlobalJointUncertaintyResult {
        try request.validate()
        let m = observableIdentifiers.count, count = request.draws.count
        guard (1...512).contains(m), Set(observableIdentifiers).count == m,
              observableIdentifiers.allSatisfy({ !$0.isEmpty && $0.utf8.count <= 4096 }) else {
            throw VivoChemistryError.invalid("joint observable identifiers")
        }
        let rawReserve = Double(count)*Double(m*m+request.parameterIdentifiers.count+4*m)
        guard rawReserve < Double(request.maximumPrimitiveWork) else { throw VivoChemistryError.resourceLimit("joint uncertainty storage/statistics budget") }
        var work = Int(rawReserve)
        var probabilities = [Double](repeating: 0,count: count)
        for model in request.models {
            let indices = request.draws.indices.filter { request.draws[$0].modelIdentifier == model.identifier }
            let maximum = indices.map { request.draws[$0].logWeight }.max()!
            let weights = indices.map { exp(request.draws[$0].logWeight-maximum) }, total = weights.reduce(0,+)
            guard total.isFinite, total > 0 else { throw VivoChemistryError.invalid("joint within-model weights") }
            for (slot,i) in indices.enumerated() { probabilities[i] = model.probability*weights[slot]/total }
        }
        // Tiny numerical residual from normalized model weights is removed once.
        let totalProbability = probabilities.reduce(0,+)
        probabilities = probabilities.map { $0/totalProbability }
        var outputs: [VivoJointObservableDraw] = []
        for (i, draw) in request.draws.enumerated() {
            let remaining = request.maximumPrimitiveWork-work
            guard remaining > 0 else { throw VivoChemistryError.resourceLimit("joint forward-model aggregate budget") }
            let output = try evaluate(draw.values,draw.modelIdentifier,remaining)
            guard output.chargedPrimitiveWork >= 0, output.chargedPrimitiveWork <= remaining,
                  output.values.count == m, output.values.allSatisfy({ $0?.isFinite != false }) else {
                throw VivoChemistryError.invalid("joint forward model returned nonfinite/incomplete output or exceeded budget")
            }
            work += output.chargedPrimitiveWork
            outputs.append(.init(parameterDrawIdentifier: draw.identifier, modelIdentifier: draw.modelIdentifier,
                probability: probabilities[i], values: output.values))
        }
        var summaries: [VivoJointObservableSummary] = []
        for column in 0..<m {
            let support = outputs.indices.filter { outputs[$0].values[column] != nil && probabilities[$0] > 0 }
            let mass = support.reduce(0) { $0+probabilities[$1] }
            if mass == 0 {
                summaries.append(.init(identifier: observableIdentifiers[column], definedProbability: 0, mean: nil, variance: nil,
                    withinModelVariance: nil,betweenModelVariance: nil,quantile025: nil,median: nil,quantile975: nil)); continue
            }
            let anchor = outputs[support[0]].values[column]!
            let mean = anchor+support.reduce(0) { $0+probabilities[$1]/mass*(outputs[$1].values[column]!-anchor) }
            let variance = support.reduce(0) { $0+probabilities[$1]/mass*pow(outputs[$1].values[column]!-mean,2) }
            var within = 0.0, between = 0.0
            for model in request.models {
                let selected = support.filter { outputs[$0].modelIdentifier == model.identifier }
                let weight = selected.reduce(0) { $0+probabilities[$1] }
                if weight == 0 { continue }
                let localAnchor = outputs[selected[0]].values[column]!
                let mu = localAnchor+selected.reduce(0) { $0+probabilities[$1]/weight*(outputs[$1].values[column]!-localAnchor) }
                within += selected.reduce(0) { $0+probabilities[$1]/mass*pow(outputs[$1].values[column]!-mu,2) }
                between += weight/mass*pow(mu-mean,2)
            }
            guard [mean,variance,within,between].allSatisfy(\.isFinite),
                  abs(variance-within-between) <= 1e-9*max(1,variance) else {
                throw VivoChemistryError.convergence("joint total-variance accounting overflow")
            }
            let sorted = support.sorted { outputs[$0].values[column]! < outputs[$1].values[column]! }
            func quantile(_ probability: Double) -> Double {
                var cumulative = 0.0
                for i in sorted { cumulative += probabilities[i]/mass; if cumulative >= probability { return outputs[i].values[column]! } }
                return outputs[sorted.last!].values[column]!
            }
            summaries.append(.init(identifier: observableIdentifiers[column], definedProbability: mass,
                mean: mean,variance: variance,withinModelVariance: within,betweenModelVariance: between,
                quantile025: quantile(0.025),median: quantile(0.5),quantile975: quantile(0.975)))
        }
        var covariance: VivoQMMatrix?
        if outputs.allSatisfy({ $0.values.allSatisfy({ $0 != nil }) }) {
            var matrix = VivoQMMatrix(m,m)
            for i in 0..<m { for j in 0...i {
                let value = outputs.indices.reduce(0.0) {
                    $0+probabilities[$1]*(outputs[$1].values[i]!-summaries[i].mean!)*(outputs[$1].values[j]!-summaries[j].mean!)
                }
                guard value.isFinite else { throw VivoChemistryError.convergence("joint observable covariance overflow") }
                matrix[i,j] = value; matrix[j,i] = value
            } }
            covariance = matrix
        }
        var blocks: [String:Double] = [:]
        for i in request.draws.indices { blocks[request.draws[i].dependenceBlockIdentifier,default: 0] += probabilities[i] }
        let weightESS = 1/probabilities.reduce(0) { $0+$1*$1 }
        let blockESS = 1/blocks.values.reduce(0) { $0+$1*$1 }
        let requestID = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(request))
        struct Evidence: Encodable { let request: VivoFingerprint; let summaries: [VivoJointObservableSummary]; let draws: [VivoJointObservableDraw]; let covariance: VivoQMMatrix? }
        let evidence = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(Evidence(request: requestID,summaries: summaries,draws: outputs,covariance: covariance)))
        return .init(requestFingerprint: requestID,contextFingerprint: request.contextFingerprint,summaries: summaries,
            draws: outputs,covariance: covariance,weightEffectiveDraws: weightESS,declaredBlockEffectiveCount: blockESS,
            unresolvedComponents: request.unresolvedComponents,assumptions: request.assumptions,chargedPrimitiveWork: work,
            evidenceFingerprint: evidence,interpretation: interpretation)
    }
}

public enum VivoJointGaussianDraws {
    /// Optional explicit distribution model. A supplied covariance is NOT
    /// inferred from marginal SDs. Singular PSD matrices are supported; genuinely
    /// indefinite covariance is rejected, not repaired by jitter.
    public static func generate(mean: [Double], covariance: VivoQMMatrix, count: Int, seed: UInt64,
                                modelIdentifier: String, maximumPrimitiveWork: Int = 100_000_000) throws -> [VivoJointParameterDraw] {
        let n = mean.count
        guard (1...256).contains(n), (2...65_536).contains(count), covariance.rows == n, covariance.columns == n,
              mean.allSatisfy(\.isFinite), !modelIdentifier.isEmpty,
              Double(n*n)*Double(n+count) <= Double(maximumPrimitiveWork) else {
            throw VivoChemistryError.invalid("joint Gaussian distribution dimensions/budget")
        }
        let scale = covariance.values.map(abs).max() ?? 0
        guard scale.isFinite else { throw VivoChemistryError.invalid("nonfinite joint covariance") }
        let normalized = try VivoQMMatrix(rows: n,columns: n,values: covariance.values.map { scale == 0 ? 0 : $0/scale })
        let eig = try VivoQMDenseAlgebra.symmetricEigen(normalized)
        guard eig.values.allSatisfy({ $0 >= -1e-12 }) else { throw VivoChemistryError.invalid("joint covariance is not positive semidefinite") }
        let roots = eig.values.map { sqrt(max(0,$0))*sqrt(scale) }
        var rng = VivoSplitMix64(state: seed), result: [VivoJointParameterDraw] = []
        for i in 0..<count {
            let z = (0..<n).map { _ in rng.normal() }
            let values = (0..<n).map { row in mean[row]+(0..<n).reduce(0) { $0+eig.vectors[row,$1]*roots[$1]*z[$1] } }
            guard values.allSatisfy(\.isFinite) else { throw VivoChemistryError.convergence("joint Gaussian sample overflow") }
            let id = "\(modelIdentifier):seed-\(seed):draw-\(i)"
            result.append(.init(identifier: id,modelIdentifier: modelIdentifier,dependenceBlockIdentifier: id,values: values))
        }
        return result
    }
}
