import Foundation

public struct VivoQMMMChemicalThermodynamicState: Codable, Sendable, Equatable {
    public var identifier: String
    /// Number of bound protons relative to an arbitrary common reference. Only
    /// differences matter, so adding the same integer to every state is invariant.
    public var boundProtonOffset: Int
    /// Relative semigrand Gibbs free energy at the request's reference pH.
    /// A common additive free-energy offset is immaterial.
    public var relativeSemigrandFreeEnergyKJPerMol: Double
    public var origin: VivoKineticOrigin
    public var evidence: VivoKineticEvidence

    public init(identifier: String, boundProtonOffset: Int,
                relativeSemigrandFreeEnergyKJPerMol: Double,
                origin: VivoKineticOrigin, evidence: VivoKineticEvidence) {
        self.identifier = identifier
        self.boundProtonOffset = boundProtonOffset
        self.relativeSemigrandFreeEnergyKJPerMol = relativeSemigrandFreeEnergyKJPerMol
        self.origin = origin
        self.evidence = evidence
    }

    public func validate() throws {
        guard !identifier.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty,
              identifier.utf8.count <= 512,
              (-128...128).contains(boundProtonOffset),
              relativeSemigrandFreeEnergyKJPerMol.isFinite else {
            throw VivoKineticsError.invalid("QM/MM chemical thermodynamic state")
        }
        try evidence.validate(origin: origin)
    }
}

public struct VivoQMMMChemicalStateThermodynamicsRequest: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/qmmm-chemical-state-thermodynamics/v1"
    public var schema: String
    public var identifier: String
    public var temperatureK: Double
    public var referencePH: Double
    public var targetPH: Double
    public var states: [VivoQMMMChemicalThermodynamicState]

    public init(identifier: String, temperatureK: Double, referencePH: Double,
                targetPH: Double, states: [VivoQMMMChemicalThermodynamicState]) {
        self.schema = Self.schema
        self.identifier = identifier
        self.temperatureK = temperatureK
        self.referencePH = referencePH
        self.targetPH = targetPH
        self.states = states
    }
}

public struct VivoQMMMChemicalStatePopulation: Codable, Sendable, Equatable {
    public let identifier: String
    public let boundProtonOffset: Int
    public let relativeSemigrandFreeEnergyKJPerMolAtTargetPH: Double
    public let logStatisticalWeight: Double
    public let population: Double
}

public struct VivoQMMMChemicalStateThermodynamicsResult: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/qmmm-chemical-state-thermodynamics-result/v1"
    public let schema: String
    public let requestFingerprint: VivoFingerprint
    public let populations: [VivoQMMMChemicalStatePopulation]
    public let mostPopulatedStateIdentifier: String
    public let interpretation: String
    public let evidenceFingerprint: VivoFingerprint

    public func population(identifier: String) -> Double? {
        populations.first(where: { $0.identifier == identifier })?.population
    }

    public var populationEvidence: VivoKineticEvidence {
        .init(source: "NumiVivo semigrand chemical-state population calculation",
              locator: "evidence-bound relative state free energies with explicit proton stoichiometry and pH shift",
              sourceFingerprint: evidenceFingerprint.hex)
    }
}

/// Reweights an explicitly enumerated state set across pH in the semigrand
/// canonical ensemble. Input free energies are relative semigrand Gibbs free
/// energies at `referencePH`. For state i with relative bound-proton count n_i,
/// moving to another pH adds n_i R T ln(10) (pH-referencePH). Therefore
/// log w_i = -G_i(reference)/(RT) - n_i ln(10)(pH-referencePH).
///
/// The operation does not enumerate states, predict pKa values, calculate the
/// supplied relative free energies, or claim constant-pH dynamical sampling.
public enum VivoQMMMChemicalStateThermodynamics {
    private static let gasConstantKJ = 0.00831446261815324
    public static let interpretation = "Semigrand-canonical pH reweighting of an explicit evidence-backed chemical-state set. Relative free energies are defined at a declared reference pH and shifted by each state's relative bound-proton stoichiometry. State enumeration, pKa prediction, conformational completeness and free-energy accuracy remain separate evidence requirements."

    public static func calculate(_ request: VivoQMMMChemicalStateThermodynamicsRequest) throws -> VivoQMMMChemicalStateThermodynamicsResult {
        guard request.schema == VivoQMMMChemicalStateThermodynamicsRequest.schema,
              !request.identifier.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty,
              request.identifier.utf8.count <= 512,
              request.temperatureK.isFinite, request.temperatureK > 0,
              request.referencePH.isFinite, request.targetPH.isFinite,
              (-10...30).contains(request.referencePH), (-10...30).contains(request.targetPH),
              !request.states.isEmpty, request.states.count <= 4096,
              Set(request.states.map(\.identifier)).count == request.states.count else {
            throw VivoKineticsError.invalid("QM/MM chemical-state thermodynamics request")
        }
        for state in request.states { try state.validate() }
        let rt = gasConstantKJ * request.temperatureK
        let deltaPH = request.targetPH - request.referencePH
        let ln10 = log(10.0)
        let targetFreeEnergies = request.states.map { state in
            state.relativeSemigrandFreeEnergyKJPerMol + Double(state.boundProtonOffset) * rt * ln10 * deltaPH
        }
        guard targetFreeEnergies.allSatisfy(\.isFinite) else {
            throw VivoKineticsError.numerical("chemical-state pH-shifted semigrand free energy")
        }
        let rawLogWeights = targetFreeEnergies.map { -$0 / rt }
        guard let maxLog = rawLogWeights.max(), maxLog.isFinite else {
            throw VivoKineticsError.numerical("chemical-state statistical weights")
        }
        let shifted = rawLogWeights.map { exp($0 - maxLog) }
        let sum = shifted.reduce(0, +)
        guard sum.isFinite, sum > 0 else {
            throw VivoKineticsError.numerical("chemical-state partition normalization")
        }
        let probabilities = shifted.map { $0 / sum }
        let probabilitySum = probabilities.reduce(0, +)
        guard probabilities.allSatisfy({ $0.isFinite && $0 >= 0 && $0 <= 1 }),
              abs(probabilitySum - 1) <= 1e-12 else {
            throw VivoKineticsError.numerical("chemical-state normalized populations")
        }
        let minimumFreeEnergy = targetFreeEnergies.min()!
        let populations = request.states.indices.map { index in
            VivoQMMMChemicalStatePopulation(identifier: request.states[index].identifier,
                boundProtonOffset: request.states[index].boundProtonOffset,
                relativeSemigrandFreeEnergyKJPerMolAtTargetPH: targetFreeEnergies[index] - minimumFreeEnergy,
                logStatisticalWeight: rawLogWeights[index] - maxLog,
                population: probabilities[index])
        }
        let maximumIndex = probabilities.indices.max { probabilities[$0] < probabilities[$1] }!
        let requestID = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(request))
        struct Evidence: Codable {
            let schema: String
            let request: VivoQMMMChemicalStateThermodynamicsRequest
            let populations: [VivoQMMMChemicalStatePopulation]
        }
        let evidenceID = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(Evidence(
            schema: "numivivo.org/qmmm-chemical-state-thermodynamics-evidence/v1",
            request: request, populations: populations)))
        return .init(schema: VivoQMMMChemicalStateThermodynamicsResult.schema,
                     requestFingerprint: requestID, populations: populations,
                     mostPopulatedStateIdentifier: request.states[maximumIndex].identifier,
                     interpretation: interpretation, evidenceFingerprint: evidenceID)
    }

    public static func validate(_ result: VivoQMMMChemicalStateThermodynamicsResult,
                                request: VivoQMMMChemicalStateThermodynamicsRequest) throws {
        guard result == (try calculate(request)) else {
            throw VivoKineticsError.invalid("chemical-state thermodynamic populations do not reconstruct")
        }
    }
}
