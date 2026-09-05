import Foundation

public enum VivoBarrierQuantity: String, Codable, Sendable {
    case electronicEnergyDifference, activationGibbsFreeEnergy
}
public enum VivoBarrierReferenceState: String, Codable, Sendable {
    case preReactiveBoundComplex, separatedReactants
}
public enum VivoMolarEnergyUnit: String, Codable, Sendable {
    case joulesPerMol = "J/mol"
    case kilojoulesPerMol = "kJ/mol"
    case kilocaloriesPerMol = "kcal/mol"
    public var joulesPerMol: Double {
        switch self { case .joulesPerMol: return 1; case .kilojoulesPerMol: return 1000; case .kilocaloriesPerMol: return 4184 }
    }
}

/// Distinguishes a sampled activation free energy from an electronic energy
/// difference. This representation does not perform QM or generate either one.
public struct VivoActivationBarrier: Codable, Equatable, Sendable {
    public let context: VivoKineticContext
    public let quantity: VivoBarrierQuantity
    public let referenceState: VivoBarrierReferenceState
    public let value: Double
    public let unit: VivoMolarEnergyUnit
    public let conditionalStandardDeviation: Double?
    public let method: String
    public let samplingDescription: String
    public let origin: VivoKineticOrigin
    public let evidence: VivoKineticEvidence
    public init(context: VivoKineticContext, quantity: VivoBarrierQuantity,
                referenceState: VivoBarrierReferenceState, value: Double, unit: VivoMolarEnergyUnit,
                conditionalStandardDeviation: Double? = nil, method: String,
                samplingDescription: String, origin: VivoKineticOrigin, evidence: VivoKineticEvidence) {
        self.context = context; self.quantity = quantity; self.referenceState = referenceState
        self.value = value; self.unit = unit; self.conditionalStandardDeviation = conditionalStandardDeviation
        self.method = method; self.samplingDescription = samplingDescription; self.origin = origin; self.evidence = evidence
    }
}

public struct VivoTransitionStateRateRequest: Codable, Equatable, Sendable {
    public let schemaVersion: UInt32
    public let barrier: VivoActivationBarrier
    /// Classical product-formation probability. Quantum tunnelling corrections
    /// are a different model and are not silently overloaded into this number.
    public let transmissionProbability: Double
    public let transmissionOrigin: VivoKineticOrigin
    public let transmissionEvidence: VivoKineticEvidence
    public init(barrier: VivoActivationBarrier, transmissionProbability: Double,
                transmissionOrigin: VivoKineticOrigin, transmissionEvidence: VivoKineticEvidence) {
        schemaVersion = 1; self.barrier = barrier; self.transmissionProbability = transmissionProbability
        self.transmissionOrigin = transmissionOrigin; self.transmissionEvidence = transmissionEvidence
    }
}

public struct VivoTransitionStateRateEstimate: Codable, Equatable, Sendable {
    public let ratePerSecond: Double
    public let naturalLogRatePerSecond: Double
    /// Propagation of the supplied barrier SD ONLY, conditional on T, kappa and
    /// model assumptions. Nil remains unknown; this is not total prediction error.
    public let conditionalLogRateStandardDeviation: Double?
    public let containsAssumptions: Bool
    public let limitations: [String]
}

/// Unimolecular conventional transition-state relation for a pre-reactive bound
/// complex. k = kappa*(kB*T/h)*exp(-deltaG/(R*T)). It is NOT kon, a free-binding
/// energy conversion, or a kinetic interpretation of biased trajectory time.
/// Primary definitions: IUPAC Gold Book T06470 and G02631.
public enum VivoTransitionStateRateEstimator {
    public static func estimate(_ request: VivoTransitionStateRateRequest) throws -> VivoTransitionStateRateEstimate {
        let b = request.barrier
        try b.context.validate(); try b.evidence.validate(origin: b.origin)
        try request.transmissionEvidence.validate(origin: request.transmissionOrigin)
        guard request.schemaVersion == 1, request.transmissionProbability.isFinite,
              request.transmissionProbability > 0, request.transmissionProbability <= 1,
              b.value.isFinite, b.value >= 0, !b.method.isEmpty, !b.samplingDescription.isEmpty,
              b.method.utf8.count <= 4096, b.samplingDescription.utf8.count <= 16_384 else {
            throw VivoKineticsError.invalid("transition-state input, sampling declaration or classical transmission probability")
        }
        guard b.quantity == .activationGibbsFreeEnergy else {
            throw VivoKineticsError.unsupported("electronic energy alone is not an activation Gibbs free energy")
        }
        guard b.referenceState == .preReactiveBoundComplex else {
            throw VivoKineticsError.unsupported("separated-reactant standard states require a different kinetic order and standard-state treatment")
        }
        let boltzmann = 1.380649e-23, planck = 6.62607015e-34, avogadro = 6.02214076e23
        let rt = boltzmann * avogadro * b.context.temperatureK
        let barrier = b.value * b.unit.joulesPerMol
        guard rt.isFinite, rt > 0, barrier.isFinite else { throw VivoKineticsError.numerical("molar energy conversion overflow") }
        let logRate = log(request.transmissionProbability) + log(boltzmann / planck)
            + log(b.context.temperatureK) - barrier / rt
        guard logRate.isFinite, logRate >= log(Double.leastNonzeroMagnitude),
              logRate <= log(Double.greatestFiniteMagnitude) else {
            throw VivoKineticsError.numerical("conditional rate is not representable in FP64")
        }
        let rate = exp(logRate)
        guard rate.isFinite, rate > 0 else { throw VivoKineticsError.numerical("rate overflow/underflow") }
        var logSD: Double?
        if let sd = b.conditionalStandardDeviation {
            guard sd.isFinite, sd >= 0 else { throw VivoKineticsError.invalid("barrier conditional SD") }
            let propagated = sd * b.unit.joulesPerMol / rt
            guard propagated.isFinite else { throw VivoKineticsError.numerical("conditional uncertainty overflow") }
            logSD = propagated
        }
        return .init(ratePerSecond: rate, naturalLogRatePerSecond: logRate,
                     conditionalLogRateStandardDeviation: logSD,
                     containsAssumptions: b.origin == .assumed || request.transmissionOrigin == .assumed,
                     limitations: [
                        "A conditional unimolecular rate from an activation free energy for the specified bound complex.",
                        "No native electronic structure, reaction-path search or sampling is executed by this conversion.",
                        "Barrier SD propagation does not include missing mechanism, sampling bias, transmission uncertainty or context transfer error.",
                        "Multiple reactive conformers or pathways require their populations and kinetic network; averaging barriers is not implemented.",
                        "Classical transmission is restricted to (0,1]; tunnelling-enhanced or barrierless kinetics are outside this contract."
                     ])
    }
}
