import Foundation

public struct VivoReactionEndpointSeedComponent: Codable, Sendable, Equatable {
    public var atomIndices:[Int]
    public var qualification:VivoNuclearQualificationRequest
    public init(atomIndices:[Int],qualification:VivoNuclearQualificationRequest) {
        self.atomIndices=atomIndices;self.qualification=qualification
    }
}
public struct VivoReactionEndpointSeed: Codable, Sendable, Equatable {
    public var identifier:String
    public var components:[VivoReactionEndpointSeedComponent]
    public init(identifier:String,components:[VivoReactionEndpointSeedComponent]) {
        self.identifier=identifier;self.components=components
    }
}
public struct VivoConnectedReactionRequest: Codable, Sendable, Equatable {
    public static let schema="numivivo.org/connected-reaction/v1"
    public var schema:String
    public var atomIdentifiers:[String]
    public var saddle:VivoNuclearQualificationRequest
    public var endpoints:[VivoReactionEndpointSeed]
    public var connectivity:VivoReactionConnectivityConfiguration
    public var reactantEndpointIdentifier:String
    public var transmission:VivoTransmissionCoefficientEvidence
    public init(atomIdentifiers:[String],saddle:VivoNuclearQualificationRequest,endpoints:[VivoReactionEndpointSeed],
                reactantEndpointIdentifier:String,connectivity:VivoReactionConnectivityConfiguration = .init(),
                transmission:VivoTransmissionCoefficientEvidence = .init()) {
        schema=Self.schema;self.atomIdentifiers=atomIdentifiers;self.saddle=saddle;self.endpoints=endpoints
        self.reactantEndpointIdentifier=reactantEndpointIdentifier;self.connectivity=connectivity;self.transmission=transmission
    }
    public func validate() throws {
        try saddle.validate();try connectivity.validate();try transmission.validate()
        let n=saddle.model.system.nuclei.count
        guard schema==Self.schema,saddle.kind == .firstOrderSaddle,atomIdentifiers.count==n,
              Set(atomIdentifiers).count==n,atomIdentifiers.allSatisfy({!$0.isEmpty && $0.utf8.count<=1024}),
              endpoints.count==2,Set(endpoints.map(\.identifier)).count==2,
              endpoints.contains(where:{$0.identifier==reactantEndpointIdentifier}) else {
            throw VivoChemistryError.invalid("connected reaction saddle, atom identities or two endpoint definitions")
        }
        let work=saddle.differences.maximumEnergyEvaluations.multipliedReportingOverflow(by:4)
        guard !work.overflow,work.partialValue<=connectivity.maximumDescentElectronicEvaluations else {
            throw VivoChemistryError.resourceLimit("connected reaction descent budget")
        }
        for endpoint in endpoints {
            guard !endpoint.identifier.isEmpty,endpoint.identifier.utf8.count<=1024,
                  !endpoint.components.isEmpty,endpoint.components.count<=min(n,16) else {
                throw VivoChemistryError.invalid("connected reaction endpoint component count")
            }
            let mapping=endpoint.components.flatMap(\.atomIndices)
            guard mapping.count==n,Set(mapping)==Set(0..<n) else {
                throw VivoChemistryError.invalid("endpoint seeds must map every saddle atom exactly once")
            }
            for component in endpoint.components {
                try component.qualification.validate()
                let q=component.qualification,m=q.model
                guard q.kind == .minimum,m.budget==saddle.model.budget,
                      component.atomIndices.count==m.system.nuclei.count,
                      m.solver==saddle.model.solver,m.scf==saddle.model.scf,m.solvent==saddle.model.solvent,
                      m.correlatedSolventConfiguration==saddle.model.correlatedSolventConfiguration,
                      m.eccFrame==nil,saddle.model.eccFrame==nil,
                      q.thermochemistry.temperatureK==saddle.thermochemistry.temperatureK,
                      q.thermochemistry.standardState==saddle.thermochemistry.standardState else {
                    throw VivoChemistryError.invalid("connected reaction endpoint model, temperature, standard state or budget differs")
                }
                for (local,global) in component.atomIndices.enumerated() {
                    guard m.system.nuclei[local].atomicNumber==saddle.model.system.nuclei[global].atomicNumber,
                          q.massesDa[local]==saddle.massesDa[global] else {
                        throw VivoChemistryError.invalid("connected reaction seed changes a mapped element or isotope")
                    }
                }
            }
        }
    }
}

/// Runs the existing stationary-point, connectivity and TST implementations in
/// their required order. Seeds remain guesses until their individual nuclear
/// qualification succeeds. Identical endpoint requests share one calculation.
public enum VivoConnectedReaction {
    public static func run(_ request:VivoConnectedReactionRequest) throws -> VivoTransitionStateTheoryResult {
        try request.validate()
        var cache:[VivoFingerprint:VivoNuclearQualifiedPoint]=[:]
        func qualify(_ seed:VivoNuclearQualificationRequest) throws -> VivoNuclearQualifiedPoint {
            let id=try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(seed))
            if let point=cache[id] { return point }
            let point=try VivoNuclearQualification.run(seed);cache[id]=point;return point
        }
        let saddle=try qualify(request.saddle)
        var endpoints:[VivoMappedReactionEndpoint]=[]
        for endpoint in request.endpoints {
            var components:[VivoMappedEndpointComponent]=[]
            for component in endpoint.components {
                components.append(.init(atomIndices:component.atomIndices,point:try qualify(component.qualification)))
            }
            endpoints.append(.init(identifier:endpoint.identifier,components:components))
        }
        let connected=try VivoReactionConnectivity.run(.init(atomIdentifiers:request.atomIdentifiers,saddle:saddle,
            endpoints:endpoints,configuration:request.connectivity))
        guard connected.converged else { throw VivoChemistryError.convergence("connected reaction did not converge at the declared path refinements") }
        return try VivoTransitionStateTheory.estimate(.init(connectivity:connected,
            reactantEndpointIdentifier:request.reactantEndpointIdentifier,transmission:request.transmission))
    }
    public static func validate(_ result:VivoTransitionStateTheoryResult,request:VivoConnectedReactionRequest) throws {
        try request.validate()
        let actual=result.request.connectivity.request
        guard actual.atomIdentifiers==request.atomIdentifiers,actual.saddle.request==request.saddle,
              actual.configuration==request.connectivity,actual.endpoints.count==request.endpoints.count,
              result.request.reactantEndpointIdentifier==request.reactantEndpointIdentifier,
              result.request.transmission==request.transmission else {
            throw VivoChemistryError.invalid("connected reaction result changed its seed or path request")
        }
        for (endpoint,seed) in zip(actual.endpoints,request.endpoints) {
            guard endpoint.identifier==seed.identifier,endpoint.components.count==seed.components.count else {
                throw VivoChemistryError.invalid("connected reaction endpoint identity changed")
            }
            for (component,input) in zip(endpoint.components,seed.components) {
                guard component.atomIndices==input.atomIndices,component.point.request==input.qualification else {
                    throw VivoChemistryError.invalid("connected reaction endpoint mapping or qualification request changed")
                }
            }
        }
        try VivoTransitionStateTheory.validate(result,request:result.request)
    }
}
