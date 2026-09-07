import Foundation

public enum VivoJointKineticParameter: Codable, Sendable, Equatable {
    case stateFreeEnergyShiftKJPerMol(state: Int)
    case initialLogPopulationShift(state: Int)
    /// PMF/flux-normalization uncertainty EXCLUDING separately sampled kappa.
    case pathwayLogPMFFactor(state: Int, pathway: Int)
    case pathwayBarrierModelShiftKJPerMol(state: Int, pathway: Int)
    /// log(kappa_draw/kappa_baseline); the resulting classical kappa stays <= 1.
    case pathwayLogTransmissionFactor(state: Int, pathway: Int)
    /// A distinct separability model input, never overwrites classical kappa.
    case pathwayLogNuclearFactor(state: Int, pathway: Int)
    case exchangeLogRateShift(edge: Int)

    public var identifier: String {
        switch self {
        case .stateFreeEnergyShiftKJPerMol(let i): return "state[\(i)].deltaG-kJ/mol"
        case .initialLogPopulationShift(let i): return "state[\(i)].delta-log-initial-population"
        case .pathwayLogPMFFactor(let s,let p): return "state[\(s)].path[\(p)].delta-log-PMF-flux"
        case .pathwayBarrierModelShiftKJPerMol(let s,let p): return "state[\(s)].path[\(p)].delta-model-barrier-kJ/mol"
        case .pathwayLogTransmissionFactor(let s,let p): return "state[\(s)].path[\(p)].delta-log-classical-kappa"
        case .pathwayLogNuclearFactor(let s,let p): return "state[\(s)].path[\(p)].log-nuclear-factor"
        case .exchangeLogRateShift(let i): return "exchange[\(i)].delta-log-rate"
        }
    }
}
public struct VivoJointKineticScenario: Codable, Sendable, Equatable {
    public let modelIdentifier: String
    public let network: VivoQMMMChemicalExchangeNetworkRequest
    public let thermodynamics: VivoQMMMChemicalStateThermodynamicsRequest?
    public init(modelIdentifier: String, network: VivoQMMMChemicalExchangeNetworkRequest,
                thermodynamics: VivoQMMMChemicalStateThermodynamicsRequest? = nil) {
        self.modelIdentifier = modelIdentifier; self.network = network; self.thermodynamics = thermodynamics
    }
}
public struct VivoGlobalKineticUncertaintyRequest: Codable, Sendable, Equatable {
    public let scenarios: [VivoJointKineticScenario]
    public let parameters: [VivoJointKineticParameter]
    public let ensemble: VivoGlobalJointUncertaintyRequest
    public init(scenarios: [VivoJointKineticScenario], parameters: [VivoJointKineticParameter],
                ensemble: VivoGlobalJointUncertaintyRequest) {
        self.scenarios = scenarios; self.parameters = parameters; self.ensemble = ensemble
    }
}

/// Baseline rates must already pass the existing independent-replica and
/// Hamiltonian/context gates. Draws are probabilistic perturbations of those
/// baselines, NOT relabeled qualified trajectories or kinetic-pack parameters.
public enum VivoGlobalKineticUncertainty {
    public static func context(scenarios: [VivoJointKineticScenario], parameters: [VivoJointKineticParameter]) throws -> VivoFingerprint {
        struct Context: Encodable { let schema: String; let scenarios: [VivoJointKineticScenario]; let parameters: [VivoJointKineticParameter] }
        return try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(Context(
            schema: "numivivo.org/global-kinetic-context/v1",scenarios: scenarios,parameters: parameters)))
    }
    public static func calculate(_ request: VivoGlobalKineticUncertaintyRequest) throws -> VivoGlobalJointUncertaintyResult {
        let scenarios = request.scenarios, parameters = request.parameters
        try request.ensemble.validate()
        guard !scenarios.isEmpty, scenarios.count <= 128,
              Set(scenarios.map(\.modelIdentifier)).count == scenarios.count,
              Set(scenarios.map(\.modelIdentifier)) == Set(request.ensemble.models.map(\.identifier)),
              request.ensemble.contextFingerprint == (try context(scenarios: scenarios,parameters: parameters)),
              request.ensemble.parameterIdentifiers == parameters.map(\.identifier),
              Set(parameters.map(\.identifier)).count == parameters.count else {
            throw VivoChemistryError.invalid("joint kinetic scenario/parameter/context binding")
        }
        let times = scenarios[0].network.observationTimesSeconds
        guard times.count <= 128 else { throw VivoChemistryError.resourceLimit("joint kinetic observation covariance capacity") }
        let stateIDs = scenarios[0].network.states.map(\.identifier)
        guard times.count*(4+stateIDs.count) <= 512 else { throw VivoChemistryError.resourceLimit("joint state/time observable covariance capacity") }
        let pathIDs = scenarios[0].network.states.map { $0.pathways.map(\.identifier) }
        let edgeIDs = scenarios[0].network.exchangeEdges.map { $0.fromStateIdentifier+"->"+$0.toStateIdentifier }
        let baselineContext = scenarios[0].network.states.first?.pathways.first?.request.replicas.first?.context
        guard let baselineContext else { throw VivoChemistryError.invalid("missing joint kinetic reaction context") }
        var byID: [String:VivoJointKineticScenario] = [:]
        for scenario in scenarios {
            _ = try VivoQMMMChemicalExchangeNetwork.resolve(scenario.network)
            guard scenario.network.observationTimesSeconds == times,
                  scenario.network.states.map(\.identifier) == stateIDs,
                  scenario.network.states.map({ $0.pathways.map(\.identifier) }) == pathIDs,
                  scenario.network.exchangeEdges.map({ $0.fromStateIdentifier+"->"+$0.toStateIdentifier }) == edgeIDs else {
                throw VivoChemistryError.invalid("joint model alternatives require the same mapped state/path/edge observables")
            }
            let c = scenario.network.states[0].pathways[0].request.replicas[0].context
            guard c.compound == baselineContext.compound, c.target == baselineContext.target,
                  c.targetVariant == baselineContext.targetVariant, c.site == baselineContext.site,
                  c.hostContext == baselineContext.hostContext, c.temperatureK == baselineContext.temperatureK,
                  c.pH == baselineContext.pH, c.ionicStrengthM == baselineContext.ionicStrengthM else {
                throw VivoChemistryError.invalid("joint scenarios differ in physical reaction conditions")
            }
            if let thermo = scenario.thermodynamics {
                _ = try VivoQMMMChemicalStateThermodynamics.calculate(thermo)
                guard thermo.states.map(\.identifier) == stateIDs,
                      thermo.temperatureK == c.temperatureK, thermo.targetPH == c.pH else {
                    throw VivoChemistryError.invalid("joint state thermodynamics differs from rate context or state order")
                }
            }
            for parameter in parameters {
                switch parameter {
                case .stateFreeEnergyShiftKJPerMol(let s):
                    guard stateIDs.indices.contains(s), scenario.thermodynamics != nil else {
                        throw VivoChemistryError.invalid("state-free-energy uncertainty requires matching thermodynamics")
                    }
                case .initialLogPopulationShift(let s):
                    guard stateIDs.indices.contains(s), scenario.thermodynamics == nil else {
                        throw VivoChemistryError.invalid("initial population perturbation cannot duplicate thermodynamic populations")
                    }
                case .exchangeLogRateShift(let e):
                    guard edgeIDs.indices.contains(e) else { throw VivoChemistryError.invalid("joint exchange edge index") }
                case .pathwayLogPMFFactor(let s,let p), .pathwayBarrierModelShiftKJPerMol(let s,let p),
                     .pathwayLogTransmissionFactor(let s,let p), .pathwayLogNuclearFactor(let s,let p):
                    guard pathIDs.indices.contains(s), pathIDs[s].indices.contains(p) else { throw VivoChemistryError.invalid("joint pathway index") }
                }
            }
            byID[scenario.modelIdentifier] = scenario
        }
        var ensemble = request.ensemble
        var retainedAssumptions = Set(ensemble.assumptions)
        for scenario in scenarios {
            for state in scenario.network.states {
                if state.populationOrigin == .assumed { retainedAssumptions.insert("Model \(scenario.modelIdentifier), state \(state.identifier): baseline initial population is assumed.") }
                for path in state.pathways where path.request.replicas.contains(where: { $0.transmissionOrigin == .assumed }) {
                    retainedAssumptions.insert("Model \(scenario.modelIdentifier), state \(state.identifier), pathway \(path.identifier): baseline classical transmission is assumed.")
                }
            }
            for edge in scenario.network.exchangeEdges where edge.rate.origin == .assumed {
                retainedAssumptions.insert("Model \(scenario.modelIdentifier): exchange rate \(edge.fromStateIdentifier) to \(edge.toStateIdentifier) is assumed.")
            }
            for state in scenario.thermodynamics?.states ?? [] where state.origin == .assumed {
                retainedAssumptions.insert("Model \(scenario.modelIdentifier), state \(state.identifier): reference semigrand free energy is assumed.")
            }
        }
        if parameters.contains(where: { if case .pathwayLogNuclearFactor = $0 { return true }; return false }) {
            retainedAssumptions.insert("Nuclear factors are supplied separable model inputs; this kinetic calculation does not evaluate quantum dynamics.")
        }
        ensemble.assumptions = retainedAssumptions.sorted()
        let observables = times.flatMap { t in
            ["t=\(t)s.survival","t=\(t)s.reacted","t=\(t)s.hazard-per-s","t=\(t)s.log-hazard-per-s"]
                + stateIDs.map { "t=\(t)s.state[\($0)].probability" }
        }
        return try VivoGlobalJointUncertainty.propagate(ensemble,observableIdentifiers: observables) { values, model, remaining in
            let scenario = byID[model]!, network = scenario.network, n = stateIDs.count
            let rt = VivoAtomicUnits.gasConstantJPerMolK*baselineContext.temperatureK/1000
            var thermo = scenario.thermodynamics
            var logits = network.states.map { $0.initialPopulation == 0 ? -Double.infinity : log($0.initialPopulation) }
            var rateLogs = network.states.map { $0.pathways.map { log($0.result.geometricMeanRatePerSecond) } }
            var kappaLogs = network.states.map { $0.pathways.map { log($0.request.replicas[0].transmissionProbability) } }
            var exchangeLogs = network.exchangeEdges.map { log($0.rate.value) }
            for (parameter,value) in zip(parameters,values) {
                switch parameter {
                case .stateFreeEnergyShiftKJPerMol(let i): thermo!.states[i].relativeSemigrandFreeEnergyKJPerMol += value
                case .initialLogPopulationShift(let i): logits[i] += value
                case .pathwayLogPMFFactor(let s,let p), .pathwayLogNuclearFactor(let s,let p): rateLogs[s][p] += value
                case .pathwayBarrierModelShiftKJPerMol(let s,let p): rateLogs[s][p] -= value/rt
                case .pathwayLogTransmissionFactor(let s,let p): rateLogs[s][p] += value; kappaLogs[s][p] += value
                case .exchangeLogRateShift(let e): exchangeLogs[e] += value
                }
            }
            guard kappaLogs.flatMap({ $0 }).allSatisfy({ $0.isFinite && $0 <= 1e-12 }) else {
                throw VivoChemistryError.invalid("joint draw makes classical kappa exceed one; do not clip the distribution")
            }
            let initial: [Double]
            if let thermo { initial = try VivoQMMMChemicalStateThermodynamics.calculate(thermo).populations.map(\.population) }
            else {
                guard let m = logits.max(), m.isFinite else { throw VivoChemistryError.invalid("joint initial population has no support") }
                let w = logits.map { exp($0-m) }, sum = w.reduce(0,+)
                guard sum.isFinite, sum > 0 else { throw VivoChemistryError.invalid("joint initial population overflow") }
                initial = w.map { $0/sum }
            }
            func rate(_ x: Double) throws -> Double {
                let v = exp(x)
                guard x.isFinite, v.isFinite, v > 0 else { throw VivoChemistryError.convergence("joint kinetic draw has unrepresentable rate") }
                return v
            }
            var q = VivoQMMatrix(n,n)
            let stateIndex = Dictionary(uniqueKeysWithValues: stateIDs.enumerated().map { ($0.element,$0.offset) })
            for (edge,x) in zip(network.exchangeEdges,exchangeLogs) {
                let i = stateIndex[edge.fromStateIdentifier]!, j = stateIndex[edge.toStateIdentifier]!
                q[i,j] = try rate(x)
            }
            for i in 0..<n {
                let chemical = try rateLogs[i].reduce(0.0) { try $0+rate($1) }
                q[i,i] = -chemical-(0..<n).reduce(0) { $0+q[i,$1] }
            }
            let trajectory = try VivoLinearKineticSensitivity.calculate(generator: q,initialProbability: initial,
                observationTimesSeconds: times,configuration: .init(maximumPrimitiveWork: min(remaining,network.maximumUniformizationWork)))
            let output: [Double?] = trajectory.observations.flatMap { point in
                [point.survivalProbability,point.reactedProbability,point.hazardPerSecond,
                 point.hazardPerSecond.flatMap { $0 > 0 ? log($0) : nil }] + point.probabilityByState.map { Optional($0) }
            }
            return .init(values: output,chargedPrimitiveWork: trajectory.chargedPrimitiveWork)
        }
    }
    public static func validate(_ result: VivoGlobalJointUncertaintyResult, request: VivoGlobalKineticUncertaintyRequest) throws {
        guard result == (try calculate(request)) else { throw VivoChemistryError.invalid("global kinetic uncertainty does not reconstruct") }
    }
}
