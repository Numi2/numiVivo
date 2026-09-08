import Foundation

/// Identifies one mapped component of the declared reactant endpoint. Equal
/// molecules in different endpoint slots remain distinct tagged/bath roles.
public struct VivoConditionalEncounterComponent: Codable, Sendable, Equatable {
    public let index: Int
    public let atomIdentifiers: [String]
    public let qualifiedPointFingerprint: VivoFingerprint
    public init(index: Int, atomIdentifiers: [String], qualifiedPointFingerprint: VivoFingerprint) {
        self.index = index; self.atomIdentifiers = atomIdentifiers
        self.qualifiedPointFingerprint = qualifiedPointFingerprint
    }
}

public enum VivoConditionalEncounterReservoirUnit: String, Codable, Sendable {
    case pascal = "Pa"
    case molar = "mol/L"
}

/// A constant externally maintained partial pressure or ideal concentration.
/// It is not a concentration calculated from the sampled MD cell or a finite
/// pool that this first-event model depletes.
public struct VivoConditionalEncounterReservoir: Codable, Sendable, Equatable {
    public let component: VivoConditionalEncounterComponent
    public let value: Double
    public let unit: VivoConditionalEncounterReservoirUnit
    public init(component: VivoConditionalEncounterComponent, value: Double, unit: VivoConditionalEncounterReservoirUnit) {
        self.component = component; self.value = value; self.unit = unit
    }
}

public struct VivoConditionalEncounterRequest: Codable, Sendable, Equatable {
    public static let schemaID = "numivivo.org/conditional-encounter/v1"
    public let schema: String
    public let identifier: String
    public let reactantEndpointIdentifier: String
    public let temperatureK: Double
    public let taggedComponent: VivoConditionalEncounterComponent
    public let reservoirs: [VivoConditionalEncounterReservoir]
    public let observationTimesSeconds: [Double]
    public let includeReservoirSensitivities: Bool
    public let numerics: VivoKineticSensitivityConfiguration

    public init(identifier: String, reactantEndpointIdentifier: String, temperatureK: Double,
                taggedComponent: VivoConditionalEncounterComponent, reservoirs: [VivoConditionalEncounterReservoir],
                observationTimesSeconds: [Double], includeReservoirSensitivities: Bool = false,
                numerics: VivoKineticSensitivityConfiguration = .init()) {
        schema = Self.schemaID; self.identifier = identifier; self.reactantEndpointIdentifier = reactantEndpointIdentifier
        self.temperatureK = temperatureK; self.taggedComponent = taggedComponent; self.reservoirs = reservoirs
        self.observationTimesSeconds = observationTimesSeconds
        self.includeReservoirSensitivities = includeReservoirSensitivities; self.numerics = numerics
    }
}

public struct VivoConditionalEncounterResult: Codable, Sendable, Equatable {
    public static let schemaID = "numivivo.org/conditional-encounter-result/v1"
    public let schema: String
    public let request: VivoConditionalEncounterRequest
    /// Canonical complete TST result identity. The workflow separately retains
    /// exact source wrapper bytes and its input artifact identity.
    public let sourceTransitionStateFingerprint: VivoFingerprint
    public let reactantComponents: [VivoConditionalEncounterComponent]
    public let molecularity: Int
    public let standardState: VivoThermoStandardState
    public let standardStateValue: Double
    public let standardStateUnits: String
    public let rateConstant: Double
    public let logRateConstant: Double
    public let rateUnits: String
    /// Unchanged source declaration; the encounter model supplies no new
    /// independent transmission measurement or uncertainty distribution.
    public let transmission: VivoTransmissionCoefficientEvidence
    public let conditionalHazardPerSecond: Double
    /// Nil exactly when a declared zero reservoir makes the hazard zero.
    public let logConditionalHazardPerSecond: Double?
    public let kinetics: VivoLinearKineticSensitivityResult
    public let interpretation: String
}

public enum VivoConditionalEncounter {
    public static let interpretation = "First-event probability, survival, flux and hazard for one tagged reactant with every other mapped reactant component externally maintained at its declared constant ideal pressure or concentration. The original molecularity, standard state and transmission declaration are retained. No finite-reservoir depletion, reverse events, occupancy, efficacy, equilibrium population, independent transmission evidence, or physical/model uncertainty qualification is supplied. Local log sensitivities are not uncertainty distributions."

    /// Descriptive bindings for constructing a request, not verified source
    /// authority. Every calculation freshly validates the complete TST result.
    public static func componentBindings(in source: VivoTransitionStateTheoryResult) throws -> [VivoConditionalEncounterComponent] {
        try Task.checkCancellation()
        try source.request.validate()
        let connection = source.request.connectivity.request
        try VivoReactionConnectivity.validateRequest(connection)
        guard let endpoint = connection.endpoints.first(where: { $0.identifier == source.request.reactantEndpointIdentifier }),
              (1...256).contains(endpoint.components.count) else {
            throw invalid("reactant component count exceeds the log-sensitivity interface")
        }
        return try endpoint.components.enumerated().map { index, component in
            try Task.checkCancellation()
            return try .init(index: index, atomIdentifiers: component.atomIndices.map { connection.atomIdentifiers[$0] },
                qualifiedPointFingerprint: VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(component.point)))
        }
    }

    /// Structural request construction only; this does not qualify a reaction.
    public static func componentBinding(index: Int, source: VivoTransitionStateTheoryResult) throws -> VivoConditionalEncounterComponent {
        let components = try componentBindings(in: source)
        guard components.indices.contains(index) else { throw invalid("component index outside the reactant endpoint") }
        return components[index]
    }

    /// Request construction/admission only. Runs the bounded probability owner
    /// to check its existing time/work limits, but does not replay or qualify
    /// the supplied reaction. A successful check is not verified source evidence.
    public static func validateConditions(_ request: VivoConditionalEncounterRequest,
                                           source: VivoTransitionStateTheoryResult) throws {
        let prepared = try prepareConditions(request, source: source)
        _ = try prepared.plan.propagate(request)
        try Task.checkCancellation()
    }

    public static func calculate(_ request: VivoConditionalEncounterRequest,
                                 source: VivoTransitionStateTheoryResult) throws -> VivoConditionalEncounterResult {
        let prepared = try prepareConditions(request, source: source)
        // Admit time/work through the existing numerical owner before expensive
        // reaction replay. Its output is discarded unless full validation passes.
        let kinetics = try prepared.plan.propagate(request)
        try Task.checkCancellation()
        // Neither decoded result summaries nor descriptive component bindings
        // grant authority. This reconstructs connectivity, RRHO and TST before
        // any successful output, including a zero-bath result.
        try VivoTransitionStateTheory.validate(source, request: source.request)
        try Task.checkCancellation()
        return .init(schema: VivoConditionalEncounterResult.schemaID, request: request,
            sourceTransitionStateFingerprint: try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(source)),
            reactantComponents: prepared.components, molecularity: source.molecularity, standardState: prepared.standardState,
            standardStateValue: source.standardStateValue, standardStateUnits: source.standardStateUnits,
            rateConstant: source.rateConstant, logRateConstant: source.logRateConstant, rateUnits: source.rateUnits,
            transmission: source.request.transmission, conditionalHazardPerSecond: prepared.plan.hazard,
            logConditionalHazardPerSecond: prepared.plan.logHazard, kinetics: kinetics, interpretation: interpretation)
    }

    private static func prepareConditions(_ request: VivoConditionalEncounterRequest, source: VivoTransitionStateTheoryResult) throws
        -> (components: [VivoConditionalEncounterComponent], standardState: VivoThermoStandardState, plan: VivoConditionalEncounterPlan) {
        try Task.checkCancellation()
        let components = try componentBindings(in: source)
        let standardState = source.request.connectivity.request.saddle.thermochemistry.configuration.standardState
        let plan = try VivoConditionalEncounterPlan.prepare(request, components: components,
            sourceEndpoint: source.request.reactantEndpointIdentifier, temperatureK: source.temperatureK,
            molecularity: source.molecularity, standardState: standardState,
            standardStateValue: source.standardStateValue, standardStateUnits: source.standardStateUnits,
            rateUnits: source.rateUnits, logRateConstant: source.logRateConstant)
        return (components, standardState, plan)
    }

    public static func validate(_ result: VivoConditionalEncounterResult, request: VivoConditionalEncounterRequest,
                                source: VivoTransitionStateTheoryResult) throws {
        guard result == (try calculate(request, source: source)) else {
            throw invalid("result does not reconstruct from its complete source reaction and conditions")
        }
    }

    private static func invalid(_ message: String) -> VivoChemistryError {
        .invalid("conditional encounter: " + message)
    }
}

/// Numerical/admission owner used only around full source validation above.
/// This plan is not a verified reaction token and cannot create a public result.
/// Keeping these cheap checks separate permits bounded host conformance tests.
struct VivoConditionalEncounterPlan {
    let hazard: Double
    let logHazard: Double?
    let reservoirIndices: [Int]

    static func prepare(_ request: VivoConditionalEncounterRequest,
                        components: [VivoConditionalEncounterComponent], sourceEndpoint: String,
                        temperatureK: Double, molecularity: Int, standardState: VivoThermoStandardState,
                        standardStateValue: Double, standardStateUnits: String,
                        rateUnits: String, logRateConstant: Double) throws -> Self {
        func invalid(_ message: String) -> VivoChemistryError { .invalid("conditional encounter: " + message) }
        guard request.schema == VivoConditionalEncounterRequest.schemaID,
              !request.identifier.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty,
              request.identifier.utf8.count <= 4096, request.reactantEndpointIdentifier == sourceEndpoint,
              !sourceEndpoint.isEmpty, sourceEndpoint.utf8.count <= 1024,
              request.temperatureK.isFinite, request.temperatureK > 0, request.temperatureK == temperatureK,
              (1...256).contains(molecularity), components.count == molecularity,
              components.enumerated().allSatisfy({ $0.offset == $0.element.index }),
              components.allSatisfy({ !$0.atomIdentifiers.isEmpty && $0.atomIdentifiers.allSatisfy({ !$0.isEmpty && $0.utf8.count <= 1024 }) }),
              Set(components.flatMap(\.atomIdentifiers)).count == components.reduce(0, { $0 + $1.atomIdentifiers.count }),
              request.reservoirs.count == molecularity - 1, logRateConstant.isFinite else {
            throw invalid("schema, endpoint, temperature or mapped component partition")
        }
        let tagged = request.taggedComponent.index
        guard components.indices.contains(tagged), request.taggedComponent == components[tagged] else {
            throw invalid("tagged component identity differs from the source endpoint")
        }
        let expectedUnit: VivoConditionalEncounterReservoirUnit, expectedValue: Double
        switch standardState {
        case .idealGas(let pressure): expectedUnit = .pascal; expectedValue = pressure
        case .concentration(let concentration, let referencePressure):
            guard referencePressure.isFinite, referencePressure > 0 else { throw invalid("gas reference pressure") }
            expectedUnit = .molar; expectedValue = concentration
        }
        let expectedRateUnits = molecularity == 1 ? "s^-1" :
            (expectedUnit == .pascal ? "Pa^\(1-molecularity) s^-1" : "(mol/L)^\(1-molecularity) s^-1")
        guard expectedValue.isFinite, expectedValue > 0, standardStateValue == expectedValue,
              standardStateUnits == expectedUnit.rawValue, rateUnits == expectedRateUnits else {
            throw invalid("source standard state, molecularity or rate units")
        }
        var seen = Set([tagged]), logarithms = [logRateConstant], zero = false
        for reservoir in request.reservoirs {
            let index = reservoir.component.index
            guard components.indices.contains(index), reservoir.component == components[index],
                  seen.insert(index).inserted, reservoir.unit == expectedUnit,
                  reservoir.value.isFinite, reservoir.value >= 0 else {
                throw invalid("missing, duplicate or mismatched reservoir component/value/units")
            }
            if reservoir.value == 0 {
                guard !request.includeReservoirSensitivities else { throw invalid("a log-reservoir derivative requires a positive reservoir value") }
                zero = true
            } else { logarithms.append(log(reservoir.value)) }
        }
        guard seen.count == components.count else { throw invalid("incomplete component partition") }
        let cfg = request.numerics, times = request.observationTimesSeconds
        guard (1...16384).contains(times.count), times.allSatisfy({ $0.isFinite && $0 >= 0 }),
              zip(times, times.dropFirst()).allSatisfy({ $0 < $1 }),
              cfg.absoluteProbabilityTolerance.isFinite, cfg.absoluteProbabilityTolerance > 0,
              cfg.absoluteSensitivityTolerance.isFinite, cfg.absoluteSensitivityTolerance > 0,
              cfg.maximumPrimitiveWork > 0, (16...4096).contains(cfg.maximumPoissonTerms) else {
            throw invalid("observation times or kinetic numerical budget")
        }
        // Logarithmic products avoid overflow/underflow of an intermediate
        // reservoir product when the final conditional hazard is representable.
        var sum = 0.0, correction = 0.0
        for value in logarithms {
            let next = sum + value
            correction += abs(sum) >= abs(value) ? (sum - next) + value : (value - next) + sum
            sum = next
        }
        let logarithm = sum + correction
        let hazard: Double
        if zero { hazard = 0 }
        else {
            guard logarithm.isFinite, logarithm >= log(Double.leastNonzeroMagnitude),
                  logarithm <= log(Double.greatestFiniteMagnitude) else {
                throw VivoChemistryError.convergence("conditional encounter positive hazard exceeds FP64 range")
            }
            hazard = exp(logarithm)
            guard hazard.isFinite, hazard > 0 else { throw VivoChemistryError.convergence("conditional encounter hazard cannot be represented") }
        }
        return .init(hazard: hazard, logHazard: zero ? nil : logarithm,
                     reservoirIndices: request.reservoirs.map { $0.component.index }.sorted())
    }

    func propagate(_ request: VivoConditionalEncounterRequest) throws -> VivoLinearKineticSensitivityResult {
        try Task.checkCancellation()
        var generator = VivoQMMatrix(1, 1)
        generator[0, 0] = -hazard
        var derivatives = [VivoKineticGeneratorDerivative(identifier: "natural-log-rate", generatorDerivativePerSecond: generator)]
        if request.includeReservoirSensitivities {
            derivatives += reservoirIndices.map {
                .init(identifier: "natural-log-reservoir-\($0)", generatorDerivativePerSecond: generator)
            }
        }
        return try VivoLinearKineticSensitivity.calculate(generator: generator, initialProbability: [1],
            observationTimesSeconds: request.observationTimesSeconds, derivatives: derivatives, configuration: request.numerics)
    }
}
