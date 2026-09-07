import Foundation

public enum VivoQMMMKineticSensitivityParameter: Codable, Sendable, Equatable {
    case pathway(stateIdentifier: String, pathwayIdentifier: String)
    case exchange(fromStateIdentifier: String, toStateIdentifier: String)
}
public struct VivoQMMMChemicalExchangeSensitivityRequest: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/qmmm-chemical-exchange-sensitivity/v1"
    public var schema: String
    public var network: VivoQMMMChemicalExchangeNetworkRequest
    /// Empty means every declared pathway and directed exchange edge.
    public var parameters: [VivoQMMMKineticSensitivityParameter]
    public var configuration: VivoKineticSensitivityConfiguration
    public init(network: VivoQMMMChemicalExchangeNetworkRequest, parameters: [VivoQMMMKineticSensitivityParameter] = [],
                configuration: VivoKineticSensitivityConfiguration = .init()) {
        schema = Self.schema; self.network = network; self.parameters = parameters; self.configuration = configuration
    }
}
public struct VivoQMMMKineticSensitivityBinding: Codable, Sendable, Equatable {
    public let identifier: String
    public let parameter: VivoQMMMKineticSensitivityParameter
    public let ratePerSecond: Double
    /// An individual dispersion is retained for inspection, but the adapter does
    /// not fabricate a joint independent covariance from these marginal values.
    public let marginalLogRateStandardDeviation: Double?
}
public struct VivoQMMMChemicalExchangeSensitivityResult: Codable, Sendable, Equatable {
    public let schema: String
    public let requestFingerprint: VivoFingerprint
    public let stateIdentifiers: [String]
    public let parameters: [VivoQMMMKineticSensitivityBinding]
    public let kinetics: VivoLinearKineticSensitivityResult
    public let evidenceFingerprint: VivoFingerprint
    public let interpretation: String
}
public enum VivoQMMMChemicalExchangeSensitivity {
    public static let interpretation = "Local derivatives of evidence-bound transient state populations, survival and hazard with respect to each declared conditional pathway or exchange log rate. Initial populations are fixed. Marginal rate errors do not establish a joint covariance, omitted-pathway completeness, or experimental agreement."
    public static func calculate(_ request: VivoQMMMChemicalExchangeSensitivityRequest) throws -> VivoQMMMChemicalExchangeSensitivityResult {
        guard request.schema == VivoQMMMChemicalExchangeSensitivityRequest.schema else {
            throw VivoChemistryError.invalid("kinetic sensitivity request schema")
        }
        let network = try VivoQMMMChemicalExchangeNetwork.resolve(request.network)
        var parameters = request.parameters
        if parameters.isEmpty {
            parameters = request.network.states.flatMap { state in
                state.pathways.map { .pathway(stateIdentifier: state.identifier,pathwayIdentifier: $0.identifier) }
            }
            parameters += request.network.exchangeEdges.map { .exchange(fromStateIdentifier: $0.fromStateIdentifier,toStateIdentifier: $0.toStateIdentifier) }
        }
        guard (1...256).contains(parameters.count) else { throw VivoChemistryError.resourceLimit("kinetic sensitivity parameter count") }
        for i in parameters.indices where parameters[..<i].contains(parameters[i]) {
            throw VivoChemistryError.invalid("duplicate kinetic sensitivity parameter")
        }
        let n = network.identifiers.count
        var derivatives: [VivoKineticGeneratorDerivative] = [], bindings: [VivoQMMMKineticSensitivityBinding] = []
        for parameter in parameters {
            var matrix = VivoQMMatrix(n,n)
            let value: Double, uncertainty: Double?
            switch parameter {
            case .pathway(let stateID, let pathID):
                guard let state = request.network.states.firstIndex(where: { $0.identifier == stateID }),
                      let path = request.network.states[state].pathways.first(where: { $0.identifier == pathID }) else {
                    throw VivoChemistryError.invalid("kinetic derivative references an unknown pathway")
                }
                value = path.result.geometricMeanRatePerSecond
                uncertainty = path.result.combinedConditionalLogRateStandardDeviation
                matrix[state,state] = -value
            case .exchange(let from, let to):
                guard let i = network.identifiers.firstIndex(of: from), let j = network.identifiers.firstIndex(of: to),
                      let edge = request.network.exchangeEdges.first(where: { $0.fromStateIdentifier == from && $0.toStateIdentifier == to }) else {
                    throw VivoChemistryError.invalid("kinetic derivative references an unknown exchange edge")
                }
                value = edge.rate.value
                if case .logNormal(let sigma) = edge.rate.uncertainty { uncertainty = sigma } else { uncertainty = nil }
                matrix[i,i] = -value; matrix[i,j] = value
            }
            let identifier = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(parameter)).hex
            derivatives.append(.init(identifier: identifier,generatorDerivativePerSecond: matrix))
            bindings.append(.init(identifier: identifier,parameter: parameter,ratePerSecond: value,marginalLogRateStandardDeviation: uncertainty))
        }
        let kinetics = try VivoLinearKineticSensitivity.calculate(generator: network.generator,
            initialProbability: request.network.states.map(\.initialPopulation),observationTimesSeconds: request.network.observationTimesSeconds,
            derivatives: derivatives,configuration: request.configuration)
        let requestID = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(request))
        struct Evidence: Codable { let request: VivoFingerprint; let parameters: [VivoQMMMKineticSensitivityBinding]; let kinetics: VivoLinearKineticSensitivityResult }
        let evidenceID = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(Evidence(request: requestID,parameters: bindings,kinetics: kinetics)))
        return .init(schema: "numivivo.org/qmmm-chemical-exchange-sensitivity-result/v1",requestFingerprint: requestID,
                     stateIdentifiers: network.identifiers,parameters: bindings,kinetics: kinetics,evidenceFingerprint: evidenceID,interpretation: interpretation)
    }
    public static func validate(_ result: VivoQMMMChemicalExchangeSensitivityResult,request: VivoQMMMChemicalExchangeSensitivityRequest) throws {
        guard result == (try calculate(request)) else { throw VivoChemistryError.invalid("kinetic sensitivity does not reconstruct from the original rate evidence") }
    }
}
