import Foundation

public struct VivoQMMMChemicalExchangeState: Codable, Sendable, Equatable {
    public var identifier: String
    public var initialPopulation: Double
    public var pathways: [VivoQMMMPathwayRate]

    public init(identifier: String, initialPopulation: Double, pathways: [VivoQMMMPathwayRate]) {
        self.identifier = identifier
        self.initialPopulation = initialPopulation
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
    public static let schema = "numivivo.org/qmmm-chemical-exchange-network/v1"
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
    public static let schema = "numivivo.org/qmmm-chemical-exchange-network-result/v1"
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
/// transient state carries its independently replicated QM/MM chemical loss rate;
/// directed exchange edges compete with that loss. The calculation returns the
/// full non-single-exponential survival curve and therefore does not manufacture
/// one global rate constant when the kinetics do not support one.
public enum VivoQMMMChemicalExchangeNetwork {
    public static let interpretation = "Explicit finite-state continuous-time chemical-state exchange with absorbing QM/MM chemical conversion. State exchange and reaction compete on their declared timescales; reported survival, hazard and apparent first-order rate are time dependent unless the network actually reduces to a single exponential."

    private static func sameEnvironment(_ a: VivoKineticContext, _ b: VivoKineticContext) -> Bool {
        a.compound == b.compound && a.target == b.target && a.targetVariant == b.targetVariant &&
        a.site == b.site && a.hostContext == b.hostContext && a.temperatureK == b.temperatureK &&
        a.pH == b.pH && a.ionicStrengthM == b.ionicStrengthM
    }

    private static func multiply(_ row: [Double], _ matrix: [[Double]]) -> [Double] {
        var out = [Double](repeating: 0, count: row.count)
        for i in row.indices where row[i] != 0 {
            for j in out.indices where matrix[i][j] != 0 {
                out[j] += row[i] * matrix[i][j]
            }
        }
        return out
    }

    private static func uniformizedStep(_ initial: [Double], transition: [[Double]], x: Double) throws -> [Double] {
        if x == 0 { return initial }
        guard x.isFinite, x > 0, x <= 32 else {
            throw VivoKineticsError.numerical("chemical-state uniformization step size")
        }
        var term = initial
        var weight = exp(-x)
        var cumulative = weight
        var output = initial.map { $0 * weight }
        for n in 1...512 {
            term = multiply(term, transition)
            weight *= x / Double(n)
            cumulative += weight
            for i in output.indices { output[i] += weight * term[i] }
            if 1 - cumulative <= 1e-14 && n > Int(x) { break }
        }
        guard output.allSatisfy({ $0.isFinite && $0 >= -1e-13 }), cumulative > 1 - 1e-11 else {
            throw VivoKineticsError.numerical("chemical-state uniformization convergence")
        }
        return output.map { max(0, $0) }
    }

    private static func advance(_ initial: [Double], transition: [[Double]], lambda: Double,
                                deltaTime: Double, maximumWork: Int) throws -> [Double] {
        if deltaTime == 0 { return initial }
        let scaled = lambda * deltaTime
        guard scaled.isFinite, scaled >= 0 else {
            throw VivoKineticsError.numerical("chemical-state exchange time scale")
        }
        let rawChunks = ceil(scaled / 32)
        guard rawChunks.isFinite, rawChunks <= Double(Int.max) else {
            throw VivoKineticsError.capacity("chemical-state exchange interval requires unrepresentable uniformization work")
        }
        let chunks = max(1, Int(rawChunks))
        let work = chunks.multipliedReportingOverflow(by: max(1, initial.count * initial.count))
        guard !work.overflow, work.partialValue <= maximumWork else {
            throw VivoKineticsError.capacity("chemical-state exchange uniformization work exceeds declared bound")
        }
        let x = scaled / Double(chunks)
        var state = initial
        for _ in 0..<chunks { state = try uniformizedStep(state, transition: transition, x: x) }
        return state
    }

    public static func calculate(_ request: VivoQMMMChemicalExchangeNetworkRequest) throws -> VivoQMMMChemicalExchangeNetworkResult {
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

        var generator = [[Double]](repeating: [Double](repeating: 0, count: request.states.count),
                                   count: request.states.count)
        for edge in request.exchangeEdges {
            let i = index[edge.fromStateIdentifier]!, j = index[edge.toStateIdentifier]!
            generator[i][j] += edge.rate.value
        }
        var lambda = 0.0
        for i in generator.indices {
            let exchangeOut = generator[i].reduce(0, +)
            let totalOut = exchangeOut + chemicalRates[i]
            generator[i][i] = -totalOut
            lambda = max(lambda, totalOut)
        }
        guard lambda.isFinite, lambda > 0 else {
            throw VivoKineticsError.numerical("chemical-state exchange generator")
        }
        var transition = [[Double]](repeating: [Double](repeating: 0, count: request.states.count),
                                    count: request.states.count)
        for i in transition.indices {
            for j in transition.indices {
                transition[i][j] = (i == j ? 1.0 : 0.0) + generator[i][j] / lambda
                guard transition[i][j].isFinite, transition[i][j] >= -1e-13 else {
                    throw VivoKineticsError.numerical("chemical-state uniformized transition matrix")
                }
                transition[i][j] = max(0, transition[i][j])
            }
        }

        var probability = request.states.map(\.initialPopulation)
        var currentTime = 0.0
        var observations: [VivoQMMMChemicalExchangeObservation] = []
        observations.reserveCapacity(request.observationTimesSeconds.count)
        for time in request.observationTimesSeconds {
            probability = try advance(probability, transition: transition, lambda: lambda,
                                      deltaTime: time - currentTime, maximumWork: request.maximumUniformizationWork)
            currentTime = time
            let survival = probability.reduce(0, +)
            guard survival.isFinite, survival >= 0, survival <= 1 + 1e-10 else {
                throw VivoKineticsError.numerical("chemical-state survival probability")
            }
            let boundedSurvival = min(1, max(0, survival))
            let reacted = max(0, 1 - boundedSurvival)
            let reactiveFlux = zip(probability, chemicalRates).reduce(0.0) { $0 + $1.0 * $1.1 }
            let hazard = boundedSurvival > 0 ? reactiveFlux / boundedSurvival : 0
            let apparent: Double? = time > 0 && boundedSurvival > 0 ? -log(boundedSurvival) / time : nil
            guard hazard.isFinite, hazard >= 0, apparent == nil || (apparent!.isFinite && apparent! >= 0) else {
                throw VivoKineticsError.numerical("chemical-state hazard or apparent rate")
            }
            observations.append(.init(timeSeconds: time, unreactedProbabilityByState: probability,
                                      survivalProbability: boundedSurvival, reactedProbability: reacted,
                                      instantaneousHazardPerSecond: hazard, apparentFirstOrderRatePerSecond: apparent))
        }

        let requestID = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(request))
        struct Evidence: Codable {
            let schema: String
            let request: VivoQMMMChemicalExchangeNetworkRequest
            let stateChemicalRatesPerSecond: [Double]
            let observations: [VivoQMMMChemicalExchangeObservation]
        }
        let evidence = Evidence(schema: "numivivo.org/qmmm-chemical-exchange-network-evidence/v1",
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
