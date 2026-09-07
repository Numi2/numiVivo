import Foundation

public struct VivoQMMMChemicalExchangeState: Codable, Sendable, Equatable {
    public var identifier: String
    public var initialPopulation: Double
    public var populationOrigin: VivoKineticOrigin
    public var populationEvidence: VivoKineticEvidence
    public var pathways: [VivoQMMMPathwayRate]

    public init(identifier: String, initialPopulation: Double,
                populationOrigin: VivoKineticOrigin,
                populationEvidence: VivoKineticEvidence,
                pathways: [VivoQMMMPathwayRate]) {
        self.identifier = identifier
        self.initialPopulation = initialPopulation
        self.populationOrigin = populationOrigin
        self.populationEvidence = populationEvidence
        self.pathways = pathways
    }

    public func validate() throws {
        guard !identifier.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty,
              identifier.utf8.count <= 512,
              initialPopulation.isFinite,
              initialPopulation >= 0,
              initialPopulation <= 1,
              !pathways.isEmpty,
              Set(pathways.map(\.identifier)).count == pathways.count else {
            throw VivoKineticsError.invalid("QM/MM exchange-state identity, initial population or pathway set")
        }
        try populationEvidence.validate(origin: populationOrigin)
        for pathway in pathways { try pathway.validate() }
    }
}

public struct VivoQMMMChemicalExchangeEdge: Codable, Sendable, Equatable {
    public var fromStateIdentifier: String
    public var toStateIdentifier: String
    public var rate: VivoKineticParameter

    public init(fromStateIdentifier: String, toStateIdentifier: String, rate: VivoKineticParameter) {
        self.fromStateIdentifier = fromStateIdentifier
        self.toStateIdentifier = toStateIdentifier
        self.rate = rate
    }

    public func validate() throws {
        guard !fromStateIdentifier.isEmpty,
              !toStateIdentifier.isEmpty,
              fromStateIdentifier != toStateIdentifier else {
            throw VivoKineticsError.invalid("chemical-state exchange edge identity")
        }
        try rate.validate(unit: .perSecond, label: "chemical-state interconversion", positive: true)
    }
}

public struct VivoQMMMChemicalExchangeNetworkRequest: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/qmmm-chemical-exchange-network/v2"
    public var schema: String
    public var identifier: String
    public var states: [VivoQMMMChemicalExchangeState]
    public var exchangeEdges: [VivoQMMMChemicalExchangeEdge]
    public var observationTimesSeconds: [Double]
    public var maximumUniformizationWork: Int

    public init(identifier: String,
                states: [VivoQMMMChemicalExchangeState],
                exchangeEdges: [VivoQMMMChemicalExchangeEdge],
                observationTimesSeconds: [Double],
                maximumUniformizationWork: Int = 2_000_000) {
        self.schema = Self.schema
        self.identifier = identifier
        self.states = states
        self.exchangeEdges = exchangeEdges
        self.observationTimesSeconds = observationTimesSeconds
        self.maximumUniformizationWork = maximumUniformizationWork
    }
}

public struct VivoQMMMChemicalExchangeObservation: Codable, Sendable, Equatable {
    public let timeSeconds: Double
    public let unreactedProbabilityByState: [Double]
    public let survivalProbability: Double
    public let reactedProbability: Double
    public let instantaneousHazardPerSecond: Double
    public let apparentFirstOrderRatePerSecond: Double?
}

public struct VivoQMMMChemicalExchangeNetworkResult: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/qmmm-chemical-exchange-network-result/v2"
    public let schema: String
    public let requestFingerprint: VivoFingerprint
    public let stateIdentifiers: [String]
    public let stateChemicalRatesPerSecond: [Double]
    public let observations: [VivoQMMMChemicalExchangeObservation]
    public let interpretation: String
    public let evidenceFingerprint: VivoFingerprint
}

/// Explicit transient continuous-time Markov treatment for bound chemical states
/// whose exchange is not fast enough for rapid-pre-equilibrium averaging. Each
/// transient state carries an evidence-bound initial population and independently
/// replicated QM/MM chemical loss rate; directed exchange edges compete with that
/// loss. The calculation returns the full non-single-exponential survival curve
/// and therefore does not manufacture one global rate constant when unsupported.
public enum VivoQMMMChemicalExchangeNetwork {
    public static let interpretation = "Explicit finite-state continuous-time chemical-state exchange with evidence-bound initial populations and absorbing QM/MM chemical conversion. State exchange and reaction compete on their declared timescales; reported survival, hazard and apparent first-order rate are time dependent unless the network actually reduces to a single exponential."

    private static func sameEnvironment(_ a: VivoKineticContext, _ b: VivoKineticContext) -> Bool {
        a.compound == b.compound && a.target == b.target && a.targetVariant == b.targetVariant &&
        a.site == b.site && a.hostContext == b.hostContext && a.temperatureK == b.temperatureK &&
        a.pH == b.pH && a.ionicStrengthM == b.ionicStrengthM
    }

    struct ResolvedNetwork {
        let identifiers: [String]
        let chemicalRates: [Double]
        let generator: VivoQMMatrix
    }
    /// One validation and generator-construction authority for propagation and
    /// local kinetic sensitivity. Derivative adapters may not bypass rate proof.
    static func resolve(_ request: VivoQMMMChemicalExchangeNetworkRequest) throws -> ResolvedNetwork {
        guard request.schema == VivoQMMMChemicalExchangeNetworkRequest.schema,
              !request.identifier.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty,
              request.identifier.utf8.count <= 512,
              (1...64).contains(request.states.count),
              Set(request.states.map(\.identifier)).count == request.states.count,
              request.exchangeEdges.count <= 4096,
              (1...16_384).contains(request.observationTimesSeconds.count),
              request.maximumUniformizationWork > 0 else {
            throw VivoKineticsError.invalid("chemical-state exchange network identity or capacity")
        }
        for state in request.states { try state.validate() }
        for edge in request.exchangeEdges { try edge.validate() }
        let population = request.states.reduce(0.0) { $0 + $1.initialPopulation }
        guard population.isFinite, abs(population - 1) <= 1e-10 else {
            throw VivoKineticsError.invalid("chemical-state initial populations must sum to one")
        }
        var previous = -Double.infinity
        for time in request.observationTimesSeconds {
            guard time.isFinite, time >= 0, time > previous else {
                throw VivoKineticsError.invalid("chemical-state observation times must be finite and strictly increasing")
            }
            previous = time
        }

        let identifiers = request.states.map(\.identifier)
        let index = Dictionary(uniqueKeysWithValues: identifiers.enumerated().map { ($1, $0) })
        var seenEdges = Set<String>()
        for edge in request.exchangeEdges {
            guard index[edge.fromStateIdentifier] != nil, index[edge.toStateIdentifier] != nil else {
                throw VivoKineticsError.invalid("chemical-state exchange edge references an unknown state")
            }
            let key = edge.fromStateIdentifier + "\u{0}" + edge.toStateIdentifier
            guard seenEdges.insert(key).inserted else {
                throw VivoKineticsError.invalid("duplicate directed chemical-state exchange edge")
            }
        }

        guard let firstPath = request.states.first?.pathways.first,
              let firstRateRequest = firstPath.request.replicas.first else {
            throw VivoKineticsError.invalid("chemical-state exchange network has no reference pathway")
        }
        let referenceContext = firstRateRequest.context
        var chemicalRates = [Double](repeating: 0, count: request.states.count)
        for (stateIndex, state) in request.states.enumerated() {
            for pathway in state.pathways {
                guard let rateRequest = pathway.request.replicas.first,
                      sameEnvironment(referenceContext, rateRequest.context),
                      rateRequest.context.chemicalState == state.identifier else {
                    throw VivoKineticsError.invalid("chemical-state exchange pathway context differs from state/environment identity")
                }
                chemicalRates[stateIndex] += pathway.result.geometricMeanRatePerSecond
            }
            guard chemicalRates[stateIndex].isFinite, chemicalRates[stateIndex] > 0 else {
                throw VivoKineticsError.numerical("chemical-state total conditional reaction rate")
            }
        }

        var generator = VivoQMMatrix(request.states.count,request.states.count)
        for edge in request.exchangeEdges {
            let i = index[edge.fromStateIdentifier]!, j = index[edge.toStateIdentifier]!
            generator[i,j] += edge.rate.value
        }
        for i in 0..<generator.rows {
            let exchangeOut = (0..<generator.columns).reduce(0.0) { $0+generator[i,$1] }
            let totalOut = exchangeOut+chemicalRates[i]
            guard totalOut.isFinite, totalOut > 0 else { throw VivoKineticsError.numerical("chemical-state generator overflow") }
            generator[i,i] = -totalOut
        }
        return .init(identifiers: identifiers,chemicalRates: chemicalRates,generator: generator)
    }

    public static func calculate(_ request: VivoQMMMChemicalExchangeNetworkRequest) throws -> VivoQMMMChemicalExchangeNetworkResult {
        let resolved = try resolve(request), identifiers = resolved.identifiers, chemicalRates = resolved.chemicalRates
        let trajectory: VivoLinearKineticSensitivityResult
        do {
            trajectory = try VivoLinearKineticSensitivity.calculate(generator: resolved.generator,
                initialProbability: request.states.map(\.initialPopulation), observationTimesSeconds: request.observationTimesSeconds,
                configuration: .init(maximumPrimitiveWork: request.maximumUniformizationWork))
        } catch VivoChemistryError.resourceLimit(let reason) { throw VivoKineticsError.capacity(reason) }
          catch VivoChemistryError.invalid(let reason) { throw VivoKineticsError.invalid(reason) }
          catch VivoChemistryError.convergence(let reason) { throw VivoKineticsError.numerical(reason) }
        let observations = trajectory.observations.map { point in
            let survival = min(1,max(0,point.survivalProbability))
            return VivoQMMMChemicalExchangeObservation(timeSeconds: point.timeSeconds,
                unreactedProbabilityByState: point.probabilityByState, survivalProbability: survival,
                reactedProbability: max(0,1-survival),instantaneousHazardPerSecond: point.hazardPerSecond ?? 0,
                apparentFirstOrderRatePerSecond: point.timeSeconds > 0 && survival > 0 ? -log(survival)/point.timeSeconds : nil)
        }

        let requestID = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(request))
        struct Evidence: Codable {
            let schema: String
            let request: VivoQMMMChemicalExchangeNetworkRequest
            let stateChemicalRatesPerSecond: [Double]
            let observations: [VivoQMMMChemicalExchangeObservation]
        }
        let evidence = Evidence(schema: "numivivo.org/qmmm-chemical-exchange-network-evidence/v2",
                                request: request, stateChemicalRatesPerSecond: chemicalRates,
                                observations: observations)
        let evidenceID = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(evidence))
        return .init(schema: VivoQMMMChemicalExchangeNetworkResult.schema,
                     requestFingerprint: requestID, stateIdentifiers: identifiers,
                     stateChemicalRatesPerSecond: chemicalRates, observations: observations,
                     interpretation: interpretation, evidenceFingerprint: evidenceID)
    }

    public static func validate(_ result: VivoQMMMChemicalExchangeNetworkResult,
                                request: VivoQMMMChemicalExchangeNetworkRequest) throws {
        guard result == (try calculate(request)) else {
            throw VivoKineticsError.invalid("chemical-state exchange network does not reconstruct")
        }
    }
}
