import Foundation

/// These are one-dimensional barrier models, not classical transmission
/// probabilities or multidimensional quantum dynamical rate calculations.
public enum VivoBarrierTunnellingModel: String, Codable, Sendable {
    case wigner, asymmetricEckart
}
public struct VivoBarrierTunnellingRequest: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/barrier-tunnelling-request/v1"
    public var schema: String = Self.schema
    public var model: VivoBarrierTunnellingModel
    public var temperatureK: Double
    /// Magnitude of the single imaginary mode, in inverse centimetres.
    public var imaginaryWavenumberPerCM: Double
    /// Positive barriers from the two asymptotes, using ONE energy/ZPE convention.
    public var forwardBarrierKJPerMol: Double
    public var reverseBarrierKJPerMol: Double
    public var hamiltonianFingerprint: VivoFingerprint
    public var stationaryPointEvidence: VivoFingerprint
    public var energyConvention: String
    public var relativeTolerance: Double
    public var maximumEvaluations: Int
    public init(model: VivoBarrierTunnellingModel, temperatureK: Double,
                imaginaryWavenumberPerCM: Double, forwardBarrierKJPerMol: Double,
                reverseBarrierKJPerMol: Double, hamiltonianFingerprint: VivoFingerprint,
                stationaryPointEvidence: VivoFingerprint, energyConvention: String,
                relativeTolerance: Double = 1e-8, maximumEvaluations: Int = 262_144) {
        self.model = model; self.temperatureK = temperatureK
        self.imaginaryWavenumberPerCM = imaginaryWavenumberPerCM
        self.forwardBarrierKJPerMol = forwardBarrierKJPerMol
        self.reverseBarrierKJPerMol = reverseBarrierKJPerMol
        self.hamiltonianFingerprint = hamiltonianFingerprint
        self.stationaryPointEvidence = stationaryPointEvidence; self.energyConvention = energyConvention
        self.relativeTolerance = relativeTolerance; self.maximumEvaluations = maximumEvaluations
    }
    public func validate() throws {
        guard schema == Self.schema, temperatureK.isFinite, temperatureK > 0,
              imaginaryWavenumberPerCM.isFinite, imaginaryWavenumberPerCM > 0,
              imaginaryWavenumberPerCM <= 100_000,
              forwardBarrierKJPerMol.isFinite, forwardBarrierKJPerMol > 0,
              reverseBarrierKJPerMol.isFinite, reverseBarrierKJPerMol > 0,
              !energyConvention.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty,
              energyConvention.utf8.count <= 1024,
              relativeTolerance.isFinite, relativeTolerance >= 1e-12, relativeTolerance <= 1e-3,
              (256...4_194_304).contains(maximumEvaluations) else {
            throw VivoChemistryError.invalid("barrier tunnelling identity, positive barriers, units or quadrature limits")
        }
    }
}
public struct VivoBarrierTunnellingResult: Codable, Sendable, Equatable {
    public let request: VivoBarrierTunnellingRequest
    public let logFactor: Double
    public let crossoverTemperatureK: Double
    public let evaluations: Int
    public let quadratureRelativeChange: Double
    public let logAbsoluteTailBound: Double?
    public let warnings: [String]
    public let evidenceFingerprint: VivoFingerprint
    public var factor: Double? {
        let x = exp(logFactor)
        return x.isFinite && x > 0 ? x : nil
    }
}

/// Native FP64 Wigner and analytic Eckart transmission with a dimensionless
/// thermal integral. Numerical convergence does NOT qualify the barrier model.
/// Wigner is refused outside its small-correction regime instead of extrapolated.
public enum VivoBarrierTunnelling {
    static let rtPerK = VivoAtomicUnits.gasConstantJPerMolK / 1000
    static let hnuPerCM = VivoAtomicUnits.planckJS * 299_792_458 * 100 * 6.02214076e23 / 1000
    private static func logSinh(_ x: Double) -> Double {
        x + log(-expm1(-2*x)) - log(2)
    }
    private static func logCosh(_ x: Double) -> Double {
        let a = abs(x); return a + log1p(exp(-2*a)) - log(2)
    }
    private static func logAdd(_ a: Double, _ b: Double) -> Double {
        if a == -.infinity { return b }; if b == -.infinity { return a }
        let m = max(a,b); return m + log(exp(a-m)+exp(b-m))
    }
    /// E is relative to the reactant asymptote. Energies below either open
    /// channel have ZERO transmission; abs(negative channel energy) is invalid.
    public static func transmissionProbability(energyKJPerMol e: Double,
                                               request r: VivoBarrierTunnellingRequest) throws -> Double {
        try r.validate()
        guard r.model == .asymmetricEckart, e.isFinite else {
            throw VivoChemistryError.invalid("Eckart transmission requires finite energy and Eckart model")
        }
        return exp(try logTransmission(e, forward: r.forwardBarrierKJPerMol,
            reverse: r.reverseBarrierKJPerMol, hnu: hnuPerCM*r.imaginaryWavenumberPerCM))
    }
    private static func logTransmission(_ e: Double, forward v1: Double, reverse v2: Double,
                                        hnu: Double) throws -> Double {
        let productEnergy = v1-v2
        if e <= max(0,productEnergy) { return -.infinity }
        let a1 = 2*Double.pi*v1/hnu, a2 = 2*Double.pi*v2/hnu
        let denominator = 1/sqrt(a1)+1/sqrt(a2)
        let a = 2*sqrt(a1*(e/v1))/denominator
        let b = 2*sqrt(a1*((e-productEnergy)/v1))/denominator
        let d2 = a1*a2 - Double.pi*Double.pi/4
        guard [a,b,d2].allSatisfy(\.isFinite), a > 0, b > 0 else {
            throw VivoChemistryError.convergence("Eckart channel scale overflow")
        }
        let s = logCosh(a+b), logDen: Double
        if d2 >= 0 { logDen = logAdd(s,logCosh(2*sqrt(d2))) }
        else {
            // Analytic continuation for a shallow barrier uses cos, not cosh(abs).
            let correction = cos(2*sqrt(-d2))*exp(-s)
            guard correction > -1 else { throw VivoChemistryError.convergence("unresolved shallow Eckart denominator") }
            logDen = s + log1p(correction)
        }
        let value = log(2) + logSinh(a) + logSinh(b) - logDen
        guard value.isFinite, value <= 1e-9 else { throw VivoChemistryError.convergence("Eckart transmission outside [0,1]") }
        return min(0,value)
    }
    public static func calculate(_ r: VivoBarrierTunnellingRequest) throws -> VivoBarrierTunnellingResult {
        try r.validate()
        let rt = rtPerK*r.temperatureK, hnu = hnuPerCM*r.imaginaryWavenumberPerCM
        guard rt.isFinite, rt > 0, hnu.isFinite, hnu > 0 else { throw VivoChemistryError.invalid("nuclear thermal scale") }
        let crossover = hnu/(2*Double.pi*rtPerK)
        var count = 0, change = 0.0, tail: Double?, value: Double
        var warnings = ["One-dimensional separable barrier approximation; not a full quantum-nuclear rate or classical recrossing coefficient."]
        switch r.model {
        case .wigner:
            let x = hnu/rt
            guard x <= 1 else { throw VivoChemistryError.unsupported("Wigner h*nu/(kBT) exceeds 1; use a justified nonperturbative model") }
            value = log1p(x*x/24)
        case .asymmetricEckart:
            // Choose the lower asymptote as zero. Forward/reverse correction
            // factors are identical; differing thresholds must not break balance.
            let hi = max(r.forwardBarrierKJPerMol,r.reverseBarrierKJPerMol)
            let lo = min(r.forwardBarrierKJPerMol,r.reverseBarrierKJPerMol)
            let threshold = hi-lo, barrier = lo/rt
            guard barrier.isFinite, barrier <= 10_000 else { throw VivoChemistryError.resourceLimit("Eckart thermal dynamic range") }
            let xMax = barrier + max(40,-log(r.relativeTolerance)+8)
            let tMax = sqrt(xMax)
            tail = barrier-xMax // since 0 <= T(E) <= 1
            func integral(_ n: Int) throws -> Double {
                guard n+1 <= r.maximumEvaluations-count else { throw VivoChemistryError.resourceLimit("Eckart quadrature budget") }
                count += n+1
                let dt = tMax/Double(n)
                var sum = -Double.infinity
                for i in 1...n {
                    let t = Double(i)*dt, x = t*t
                    let logT = try logTransmission(threshold+rt*x, forward: hi, reverse: lo, hnu: hnu)
                    let w = i == n ? 1.0 : (i % 2 == 0 ? 2.0 : 4.0)
                    sum = logAdd(sum,logT + barrier-x + log(2*t*w))
                }
                return sum + log(dt/3)
            }
            var n = 64, previous = try integral(64), converged = false
            value = previous
            while n <= (r.maximumEvaluations-1)/2 {
                n *= 2; value = try integral(n)
                change = abs(expm1(previous-value))
                if n >= 256 && change <= r.relativeTolerance && tail! - value < log(r.relativeTolerance) {
                    converged = true; break
                }
                previous = value
            }
            guard converged else { throw VivoChemistryError.convergence("Eckart thermal integral did not converge") }
            if r.temperatureK < crossover { warnings.append("Below crossover: a fitted one-dimensional barrier can miss multidimensional tunnelling paths.") }
        }
        guard value.isFinite else { throw VivoChemistryError.convergence("nonfinite logarithmic nuclear correction") }
        struct Evidence: Codable { let request: VivoBarrierTunnellingRequest; let logFactor: Double; let evaluations: Int; let change: Double }
        let fp = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(Evidence(request: r, logFactor: value, evaluations: count, change: change)))
        return .init(request: r, logFactor: value, crossoverTemperatureK: crossover, evaluations: count,
            quadratureRelativeChange: change, logAbsoluteTailBound: tail, warnings: warnings, evidenceFingerprint: fp)
    }
    public static func validate(_ result: VivoBarrierTunnellingResult) throws {
        guard result == (try calculate(result.request)) else { throw VivoChemistryError.invalid("nuclear correction does not reconstruct") }
    }
}

/// A separate, explicit model composition. Never overwrites a PMF's classical
/// kappa, infers separability, or attaches itself to a kinetic pack automatically.
public struct VivoNuclearRateComposition: Codable, Sendable, Equatable {
    public let baselineFingerprint: VivoFingerprint
    public let hamiltonianFingerprint: VivoFingerprint
    public let temperatureK: Double
    public let baselineLogRatePerSecond: Double
    public let correction: VivoBarrierTunnellingResult
    public let separabilityAssumption: String
    public let correctedLogRatePerSecond: Double
    public init(baselineFingerprint: VivoFingerprint, hamiltonianFingerprint: VivoFingerprint,
                temperatureK: Double, baselineLogRatePerSecond: Double,
                baselineAlreadyContainsTunnelling: Bool, correction: VivoBarrierTunnellingResult,
                separabilityAssumption: String) throws {
        try VivoBarrierTunnelling.validate(correction)
        guard !baselineAlreadyContainsTunnelling, baselineLogRatePerSecond.isFinite,
              temperatureK == correction.request.temperatureK,
              hamiltonianFingerprint == correction.request.hamiltonianFingerprint,
              !separabilityAssumption.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty,
              separabilityAssumption.utf8.count <= 4096,
              (baselineLogRatePerSecond+correction.logFactor).isFinite else {
            throw VivoChemistryError.invalid("nuclear rate composition double counting, context or undeclared coupling assumption")
        }
        self.baselineFingerprint = baselineFingerprint; self.hamiltonianFingerprint = hamiltonianFingerprint
        self.temperatureK = temperatureK; self.baselineLogRatePerSecond = baselineLogRatePerSecond
        self.correction = correction; self.separabilityAssumption = separabilityAssumption
        correctedLogRatePerSecond = baselineLogRatePerSecond+correction.logFactor
    }
}
