import Foundation

public enum VivoTransmissionEvidenceKind: String, Codable, Sendable {
    /// kappa=1 or another explicitly assumed value. This is a model assumption,
    /// not independent dynamical evidence.
    case assumed
    /// Coefficient obtained from an explicitly identified dynamical calculation.
    case computedDynamics
    /// Coefficient constrained by an explicitly identified experimental source.
    case experimental
}

public struct VivoTransmissionCoefficientEvidence: Codable, Sendable, Equatable {
    public let coefficient: Double
    public let kind: VivoTransmissionEvidenceKind
    public let sourceIdentifier: String
    public init(coefficient: Double = 1, kind: VivoTransmissionEvidenceKind = .assumed,
                sourceIdentifier: String = "classical no-recrossing TST assumption") {
        self.coefficient = coefficient; self.kind = kind; self.sourceIdentifier = sourceIdentifier
    }
    public func validate() throws {
        guard coefficient.isFinite, coefficient > 0, coefficient <= 1e6,
              !sourceIdentifier.isEmpty, sourceIdentifier.utf8.count <= 4096 else {
            throw VivoChemistryError.invalid("transmission coefficient, evidence kind or source identity")
        }
        if kind == .assumed, sourceIdentifier.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty {
            throw VivoChemistryError.invalid("assumed transmission coefficient requires an explicit model label")
        }
    }
}

public struct VivoTransitionStateTheoryRequest: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/transition-state-theory-rate/v1"
    public let schema: String
    public let connectivity: VivoReactionConnectivityResult
    public let reactantEndpointIdentifier: String
    public let transmission: VivoTransmissionCoefficientEvidence
    public init(connectivity: VivoReactionConnectivityResult, reactantEndpointIdentifier: String,
                transmission: VivoTransmissionCoefficientEvidence = .init()) {
        schema = Self.schema; self.connectivity = connectivity
        self.reactantEndpointIdentifier = reactantEndpointIdentifier; self.transmission = transmission
    }
    public func validate() throws {
        try transmission.validate()
        guard schema == Self.schema, connectivity.converged,
              !reactantEndpointIdentifier.isEmpty,
              connectivity.request.endpoints.contains(where: { $0.identifier == reactantEndpointIdentifier }) else {
            throw VivoChemistryError.invalid("TST rate requires a converged mapped connection and declared reactant endpoint")
        }
    }
}

public struct VivoTransitionStateTheoryResult: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/transition-state-theory-rate-result/v1"
    public let schema: String
    public let request: VivoTransitionStateTheoryRequest
    public let barrier: VivoHarmonicBarrierEstimate
    public let temperatureK: Double
    public let molecularity: Int
    public let transmissionCoefficient: Double
    public let standardStateValue: Double
    public let standardStateUnits: String
    public let rateUnits: String
    public let thermalFrequencyPerSecond: Double
    public let activationExponent: Double
    public let logRateConstant: Double
    public let rateConstant: Double
    public let hasIndependentTransmissionEvidence: Bool
    public let interpretation: String
}

/// Canonical Eyring/TST conversion after three separate prerequisites:
///   1) a reconstructed first-order saddle/RRHO barrier,
///   2) independently mapped two-branch reaction connectivity,
///   3) an explicit transmission-coefficient evidence record.
///
/// The returned numerical rate is NOT promoted to experimental/predictive truth.
/// `hasIndependentTransmissionEvidence` only distinguishes an externally sourced
/// kappa from a classical assumption; it does not validate the electronic model.
public enum VivoTransitionStateTheory {
    public static let planckJouleSecond = 6.62607015e-34
    public static let interpretation = "Eyring transition-state-theory standard-state rate reconstructed from a qualified local saddle, mapped two-branch connectivity, ideal RRHO activation Gibbs energy and an explicit transmission coefficient; no claim of tunneling/recrossing accuracy, solvent dynamical completeness, experimental validation or paper reproduction"

    private static func endpoint(_ request: VivoTransitionStateTheoryRequest) throws -> VivoMappedReactionEndpoint {
        guard let endpoint = request.connectivity.request.endpoints.first(where: { $0.identifier == request.reactantEndpointIdentifier }) else {
            throw VivoChemistryError.invalid("TST reactant endpoint disappeared from connectivity request")
        }
        return endpoint
    }

    public static func estimate(_ request: VivoTransitionStateTheoryRequest) throws -> VivoTransitionStateTheoryResult {
        try request.validate()
        // A stored `converged=true` flag is not sufficient. Reconstruct the
        // displacement/step refinement and mapped endpoint assignment.
        try VivoReactionConnectivity.validate(request.connectivity, request: request.connectivity.request)
        let reactantEndpoint = try endpoint(request)
        let reactants = reactantEndpoint.components.map(\.point)
        let saddle = request.connectivity.request.saddle
        let barrier = try VivoHarmonicBarrier.estimate(saddle: saddle, reactants: reactants)
        let cfg = saddle.thermochemistry.configuration, temperature = cfg.temperatureK
        guard barrier.reactantMolecularity == reactants.count, barrier.activationGibbsHartree.isFinite,
              temperature.isFinite, temperature > 0 else {
            throw VivoChemistryError.invalid("TST barrier, molecularity or temperature")
        }
        let thermal = VivoAtomicUnits.boltzmannJPerK * temperature / planckJouleSecond
        let exponent = -barrier.activationGibbsHartree / (VivoNuclearUnits.kHartree * temperature)
        let standardValue: Double, standardUnits: String, rateUnits: String
        switch cfg.standardState {
        case .idealGas(let pressure):
            standardValue = pressure; standardUnits = "Pa"
            rateUnits = reactants.count == 1 ? "s^-1" : "Pa^\(1-reactants.count) s^-1"
        case .concentration(let concentration, _):
            standardValue = concentration; standardUnits = "mol/L"
            rateUnits = reactants.count == 1 ? "s^-1" : "(mol/L)^\(1-reactants.count) s^-1"
        }
        guard thermal.isFinite, thermal > 0, standardValue.isFinite, standardValue > 0 else {
            throw VivoChemistryError.convergence("TST thermal prefactor or standard-state scale")
        }
        let logRate = log(request.transmission.coefficient) + log(thermal) + exponent
            + Double(1 - reactants.count) * log(standardValue)
        guard logRate.isFinite else { throw VivoChemistryError.convergence("TST logarithmic rate overflow") }
        let minimumLog = log(Double.leastNonzeroMagnitude), maximumLog = log(Double.greatestFiniteMagnitude)
        let rate: Double
        if logRate < minimumLog { rate = 0 }
        else if logRate > maximumLog { throw VivoChemistryError.convergence("TST rate exceeds FP64 range") }
        else { rate = exp(logRate) }
        return .init(schema: VivoTransitionStateTheoryResult.schema, request: request, barrier: barrier,
            temperatureK: temperature, molecularity: reactants.count,
            transmissionCoefficient: request.transmission.coefficient,
            standardStateValue: standardValue, standardStateUnits: standardUnits, rateUnits: rateUnits,
            thermalFrequencyPerSecond: thermal, activationExponent: exponent,
            logRateConstant: logRate, rateConstant: rate,
            hasIndependentTransmissionEvidence: request.transmission.kind != .assumed,
            interpretation: interpretation)
    }

    public static func validate(_ result: VivoTransitionStateTheoryResult,
                                request: VivoTransitionStateTheoryRequest) throws {
        guard result.schema == VivoTransitionStateTheoryResult.schema, result.request == request,
              result.interpretation == interpretation else {
            throw VivoChemistryError.invalid("TST result request, schema or interpretation binding")
        }
        let rebuilt = try estimate(request)
        guard rebuilt == result else {
            throw VivoChemistryError.invalid("TST connectivity, barrier, standard-state or rate reconstruction differs")
        }
    }
}
