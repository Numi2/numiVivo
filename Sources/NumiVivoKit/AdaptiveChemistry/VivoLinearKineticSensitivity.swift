import Foundation

public struct VivoKineticGeneratorDerivative: Codable, Sendable, Equatable {
    public let identifier: String
    public let parameterization: String
    public let generatorDerivativePerSecond: VivoQMMatrix
    public let initialProbabilityDerivative: [Double]
    public init(identifier: String, parameterization: String = "natural-log-rate",
                generatorDerivativePerSecond: VivoQMMatrix, initialProbabilityDerivative: [Double] = []) {
        self.identifier = identifier; self.parameterization = parameterization
        self.generatorDerivativePerSecond = generatorDerivativePerSecond
        self.initialProbabilityDerivative = initialProbabilityDerivative
    }
}
public struct VivoKineticSensitivityConfiguration: Codable, Sendable, Equatable {
    public var absoluteProbabilityTolerance: Double
    public var absoluteSensitivityTolerance: Double
    public var maximumPrimitiveWork: Int
    public var maximumPoissonTerms: Int
    public init(absoluteProbabilityTolerance: Double = 1e-12, absoluteSensitivityTolerance: Double = 1e-10,
                maximumPrimitiveWork: Int = 100_000_000, maximumPoissonTerms: Int = 512) {
        self.absoluteProbabilityTolerance = absoluteProbabilityTolerance
        self.absoluteSensitivityTolerance = absoluteSensitivityTolerance
        self.maximumPrimitiveWork = maximumPrimitiveWork; self.maximumPoissonTerms = maximumPoissonTerms
    }
}
public struct VivoKineticObservableDerivative: Codable, Sendable, Equatable {
    public let parameterIdentifier: String
    public let probabilityByState: [Double]
    public let survivalProbability: Double
    public let reactedProbability: Double
    public let logSurvivalProbability: Double?
    public let reactiveFluxPerSecond: Double
    public let hazardPerSecond: Double?
}
public struct VivoLinearKineticObservation: Codable, Sendable, Equatable {
    public let timeSeconds: Double
    public let probabilityByState: [Double]
    public let survivalProbability: Double
    public let reactedProbability: Double
    public let reactiveFluxPerSecond: Double
    public let hazardPerSecond: Double?
    public let derivatives: [VivoKineticObservableDerivative]
}
public struct VivoLinearKineticSensitivityResult: Codable, Sendable, Equatable {
    public let observations: [VivoLinearKineticObservation]
    public let chargedPrimitiveWork: Int
    /// Analytic truncation estimates only; not a bound on FP64 roundoff or the
    /// model's physical/parameter uncertainty.
    public let probabilityPoissonTailBound: Double
    public let sensitivityPoissonTailBounds: [Double]
    public let interpretation: String
}

/// One uniformization authority for transient finite-state dynamics and their
/// forward generator sensitivities. Lambda is held fixed when differentiating;
/// differentiating the argmax used to choose lambda would be erroneous.
public enum VivoLinearKineticSensitivity {
    public static let interpretation = "Transient row-generator evolution and local parameter derivatives through the same uniformization recurrence. Poisson tail estimates exclude floating-point roundoff. Parameter derivatives are local sensitivities, not uncertainty distributions, mechanism completeness, or experimental validation."
    public static func calculate(generator q: VivoQMMatrix, initialProbability: [Double],
                                 observationTimesSeconds times: [Double],
                                 derivatives: [VivoKineticGeneratorDerivative] = [],
                                 configuration cfg: VivoKineticSensitivityConfiguration = .init()) throws -> VivoLinearKineticSensitivityResult {
        let n = q.rows, m = derivatives.count
        guard (1...128).contains(n), q.columns == n, initialProbability.count == n,
              q.values.allSatisfy(\.isFinite), initialProbability.allSatisfy({ $0.isFinite && $0 >= 0 }),
              initialProbability.reduce(0,+) <= 1+1e-12,
              (1...16384).contains(times.count), times.allSatisfy({ $0.isFinite && $0 >= 0 }),
              zip(times,times.dropFirst()).allSatisfy({ $0 < $1 }), m <= 256,
              Set(derivatives.map(\.identifier)).count == m,
              cfg.absoluteProbabilityTolerance.isFinite, cfg.absoluteProbabilityTolerance > 0,
              cfg.absoluteSensitivityTolerance.isFinite, cfg.absoluteSensitivityTolerance > 0,
              cfg.maximumPrimitiveWork > 0, (16...4096).contains(cfg.maximumPoissonTerms) else {
            throw VivoChemistryError.invalid("linear kinetic sensitivity dimensions, tolerances or initial state")
        }
        guard n*n*max(1,m) <= cfg.maximumPrimitiveWork,
              times.count*n*max(1,m) <= cfg.maximumPrimitiveWork else {
            throw VivoChemistryError.resourceLimit("linear kinetic sensitivity storage")
        }
        var loss = [Double](repeating: 0,count: n), lambda = 0.0
        for i in 0..<n {
            let rowSum = (0..<n).reduce(0.0) { $0+q[i,$1] }
            let rowScale = (0..<n).reduce(0.0) { $0+abs(q[i,$1]) }
            guard q[i,i] <= 0, rowSum <= 64*Double.ulpOfOne*rowScale else {
                throw VivoChemistryError.invalid("transient generator must have nonpositive row sums")
            }
            for j in 0..<n where i != j && q[i,j] < 0 { throw VivoChemistryError.invalid("negative kinetic transition") }
            loss[i] = max(0,-rowSum); lambda = max(lambda,-q[i,i])
        }
        var sensitivities: [[Double]] = [], derivativeLoss: [[Double]] = [], derivativeNorms: [Double] = []
        for derivative in derivatives {
            let b = derivative.generatorDerivativePerSecond
            guard !derivative.identifier.isEmpty, derivative.identifier.utf8.count <= 512,
                  ["natural-log-rate","dimensionless-parameter"].contains(derivative.parameterization),
                  b.rows == n, b.columns == n, b.values.allSatisfy(\.isFinite),
                  derivative.initialProbabilityDerivative.isEmpty ||
                    (derivative.initialProbabilityDerivative.count == n && derivative.initialProbabilityDerivative.allSatisfy(\.isFinite)) else {
                throw VivoChemistryError.invalid("kinetic generator derivative")
            }
            sensitivities.append(derivative.initialProbabilityDerivative.isEmpty ? [Double](repeating: 0,count: n) : derivative.initialProbabilityDerivative)
            derivativeLoss.append((0..<n).map { i in -(0..<n).reduce(0.0) { $0+b[i,$1] } })
            derivativeNorms.append((0..<n).map { i in (0..<n).reduce(0.0) { $0+abs(b[i,$1]) } }.max()!)
        }
        let finalTime = times.last!, maxB = derivativeNorms.max() ?? 0
        let rawChunks = ceil(lambda*finalTime/16)+Double(times.count)
        guard rawChunks.isFinite, rawChunks < Double(Int.max), finalTime*maxB < Double.greatestFiniteMagnitude/4 else {
            throw VivoChemistryError.resourceLimit("kinetic sensitivity time/work scale")
        }
        let chunkBound = max(1,Int(rawChunks))
        let pTolerance = min(cfg.absoluteProbabilityTolerance/(2*Double(chunkBound)),
                             cfg.absoluteSensitivityTolerance/(2*Double(chunkBound)*(1+finalTime*maxB)))
        let sTolerance = cfg.absoluteSensitivityTolerance/(2*Double(chunkBound))
        guard pTolerance > 0, sTolerance > 0 else { throw VivoChemistryError.resourceLimit("unrepresentable kinetic truncation tolerance") }
        var work = 0
        func charge(_ count: Int) throws {
            guard count <= cfg.maximumPrimitiveWork-work else { throw VivoChemistryError.resourceLimit("kinetic aggregate primitive work") }
            work += count
        }
        func product(_ row: [Double],_ matrix: VivoQMMatrix) throws -> [Double] {
            try charge(n*n)
            var value = [Double](repeating: 0,count: n)
            for i in 0..<n where row[i] != 0 { for j in 0..<n { value[j] += row[i]*matrix[i,j] } }
            return value
        }
        var transition = VivoQMMatrix(n,n), scaledDerivatives: [VivoQMMatrix] = []
        if lambda > 0 {
            for i in 0..<n { for j in 0..<n {
                transition[i,j] = (i == j ? 1.0 : 0.0)+q[i,j]/lambda
                guard transition[i,j] >= 0 else { throw VivoChemistryError.convergence("negative uniformized transition") }
            } }
            scaledDerivatives = derivatives.map { $0.generatorDerivativePerSecond.scaled(1/lambda) }
        }
        var probability = initialProbability, currentTime = 0.0, pError = 0.0
        var sErrors = [Double](repeating: 0,count: m), observations: [VivoLinearKineticObservation] = []
        func step(_ deltaTime: Double) throws {
            if lambda == 0 {
                for j in 0..<m {
                    let change = try product(probability,derivatives[j].generatorDerivativePerSecond)
                    for i in 0..<n { sensitivities[j][i] += deltaTime*change[i] }
                }
                return
            }
            let x = lambda*deltaTime
            var v = probability, s = sensitivities, weight = exp(-x)
            var output = v.map { $0*weight }, derivativeOutput = s.map { $0.map { $0*weight } }
            let initialMass = probability.reduce(0,+), initialNorms = s.map { $0.reduce(0) { $0+abs($1) } }
            var tailP = Double.infinity, tailS = [Double](repeating: Double.infinity,count: m), finished = false
            for term in 1...cfg.maximumPoissonTerms {
                let old = v; v = try product(old,transition)
                for parameter in 0..<m {
                    let propagated = try product(s[parameter],transition)
                    let forcing = try product(old,scaledDerivatives[parameter])
                    s[parameter] = zip(propagated,forcing).map(+)
                }
                weight *= x/Double(term)
                for i in 0..<n { output[i] += weight*v[i] }
                for j in 0..<m { for i in 0..<n { derivativeOutput[j][i] += weight*s[j][i] } }
                if Double(term+2) > x {
                    let next = weight*x/Double(term+1), ratio = x/Double(term+2)
                    let tailMass = next/(1-ratio)
                    let tailFirstMoment = next*(Double(term+1)/(1-ratio)+ratio/((1-ratio)*(1-ratio)))
                    tailP = initialMass*tailMass
                    tailS = (0..<m).map { tailMass*initialNorms[$0]+initialMass*tailFirstMoment*derivativeNorms[$0]/lambda }
                    if tailP <= pTolerance && tailS.allSatisfy({ $0 <= sTolerance }) { finished = true; break }
                }
            }
            guard finished, output.allSatisfy({ $0.isFinite && $0 >= 0 }),
                  derivativeOutput.flatMap({ $0 }).allSatisfy(\.isFinite) else {
                throw VivoChemistryError.convergence("kinetic sensitivity Poisson series failed its joint tolerance")
            }
            for j in 0..<m { sErrors[j] += tailS[j]+deltaTime*derivativeNorms[j]*pError }
            pError += tailP; probability = output; sensitivities = derivativeOutput
        }
        for time in times {
            let duration = time-currentTime
            let chunks = max(1,Int(ceil(lambda*duration/16)))
            if duration > 0 { for _ in 0..<chunks { try step(duration/Double(chunks)) } }
            currentTime = time
            let survival = probability.reduce(0,+), flux = zip(probability,loss).reduce(0) { $0+$1.0*$1.1 }
            guard survival.isFinite, survival >= 0, survival <= 1+1e-10, flux.isFinite, flux >= 0 else {
                throw VivoChemistryError.convergence("kinetic probability or flux invariant")
            }
            let hazard: Double? = survival > 0 ? flux/survival : nil
            var observed: [VivoKineticObservableDerivative] = []
            for j in 0..<m {
                let ds = sensitivities[j].reduce(0,+)
                let df = zip(sensitivities[j],loss).reduce(0) { $0+$1.0*$1.1 }
                       + zip(probability,derivativeLoss[j]).reduce(0) { $0+$1.0*$1.1 }
                let dh: Double? = survival > 0 ? (df-hazard!*ds)/survival : nil
                guard ds.isFinite, df.isFinite, dh == nil || dh!.isFinite else {
                    throw VivoChemistryError.convergence("kinetic derivative overflow")
                }
                observed.append(.init(parameterIdentifier: derivatives[j].identifier,probabilityByState: sensitivities[j],
                    survivalProbability: ds,reactedProbability: -ds,logSurvivalProbability: survival > 0 ? ds/survival : nil,
                    reactiveFluxPerSecond: df,hazardPerSecond: dh))
            }
            observations.append(.init(timeSeconds: time,probabilityByState: probability,survivalProbability: survival,
                reactedProbability: max(0,1-survival),reactiveFluxPerSecond: flux,hazardPerSecond: hazard,derivatives: observed))
        }
        guard pError <= cfg.absoluteProbabilityTolerance, sErrors.allSatisfy({ $0 <= cfg.absoluteSensitivityTolerance }) else {
            throw VivoChemistryError.convergence("kinetic cumulative truncation tolerance")
        }
        return .init(observations: observations,chargedPrimitiveWork: work,probabilityPoissonTailBound: pError,
                     sensitivityPoissonTailBounds: sErrors,interpretation: interpretation)
    }
}
