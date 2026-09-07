import Foundation

public struct VivoQMMMPathwayRate:Codable,Sendable,Equatable {
    public var identifier:String
    public var request:VivoQMMMReplicatedFreeEnergyRateRequest
    public var result:VivoQMMMReplicatedFreeEnergyRateResult
    public init(identifier:String,request:VivoQMMMReplicatedFreeEnergyRateRequest,
                result:VivoQMMMReplicatedFreeEnergyRateResult) {
        self.identifier=identifier;self.request=request;self.result=result
    }
    public func validate()throws {
        guard !identifier.trimmingCharacters(in:.whitespacesAndNewlines).isEmpty,identifier.utf8.count<=512 else {
            throw VivoKineticsError.invalid("QM/MM pathway identity")
        }
        try VivoQMMMReplicatedFreeEnergyRate.validate(result,request:request)
        guard result.converged else { throw VivoKineticsError.invalid("QM/MM pathway rate is not independently replicated/converged") }
    }
}

public struct VivoQMMMChemicalStateRate:Codable,Sendable,Equatable {
    public var identifier:String
    public var equilibriumPopulation:Double
    public var populationOrigin:VivoKineticOrigin
    public var populationEvidence:VivoKineticEvidence
    public var pathways:[VivoQMMMPathwayRate]
    public init(identifier:String,equilibriumPopulation:Double,populationOrigin:VivoKineticOrigin,
                populationEvidence:VivoKineticEvidence,pathways:[VivoQMMMPathwayRate]) {
        self.identifier=identifier;self.equilibriumPopulation=equilibriumPopulation;self.populationOrigin=populationOrigin
        self.populationEvidence=populationEvidence;self.pathways=pathways
    }
    public func validate()throws {
        guard !identifier.trimmingCharacters(in:.whitespacesAndNewlines).isEmpty,identifier.utf8.count<=512,
              equilibriumPopulation.isFinite,equilibriumPopulation>=0,equilibriumPopulation<=1,
              !pathways.isEmpty,Set(pathways.map(\.identifier)).count==pathways.count else {
            throw VivoKineticsError.invalid("QM/MM chemical-state population or pathway set")
        }
        try populationEvidence.validate(origin:populationOrigin)
        for pathway in pathways { try pathway.validate() }
    }
}

public struct VivoQMMMChemicalStateNetworkRequest:Codable,Sendable,Equatable {
    public static let schema="numivivo.org/qmmm-chemical-state-network/v1"
    public var schema:String
    public var identifier:String
    public var states:[VivoQMMMChemicalStateRate]
    /// This aggregation is physically valid only when state exchange is fast
    /// relative to chemical conversion. Otherwise an explicit kinetic exchange
    /// network is required and this request is rejected.
    public var rapidPreEquilibrium:Bool
    public init(identifier:String,states:[VivoQMMMChemicalStateRate],rapidPreEquilibrium:Bool) {
        schema=Self.schema;self.identifier=identifier;self.states=states;self.rapidPreEquilibrium=rapidPreEquilibrium
    }
}

public struct VivoQMMMPathwayContribution:Codable,Sendable,Equatable {
    public let stateIdentifier:String
    public let pathwayIdentifier:String
    public let conditionalRatePerSecond:Double
    public let statePopulation:Double
    public let contributionPerSecond:Double
}

public struct VivoQMMMChemicalStateNetworkResult:Codable,Sendable,Equatable {
    public static let schema="numivivo.org/qmmm-chemical-state-network-result/v1"
    public let schema:String
    public let requestFingerprint:VivoFingerprint
    public let effectiveRatePerSecond:Double
    public let contributions:[VivoQMMMPathwayContribution]
    public let parameter:VivoKineticParameter
    public let interpretation:String
    public let evidenceFingerprint:VivoFingerprint
}

/// Under rapid pre-equilibrium, parallel elementary pathways within one chemical
/// state add as rates, then state-specific rates are population weighted:
/// k_eff = sum_s p_s sum_j k_sj. Barriers are never averaged. This operation
/// refuses slow/intermediate state exchange because that requires explicit
/// interconversion rates in a kinetic network.
public enum VivoQMMMChemicalStateNetwork {
    public static let interpretation="Population-weighted effective first-order chemical conversion rate under an explicitly asserted rapid pre-equilibrium among bound chemical states; parallel elementary pathways add as rates. State interconversion kinetics, binding-state populations and proton chemical potentials are not inferred by this operation."

    private static func sameEnvironment(_ a:VivoKineticContext,_ b:VivoKineticContext)->Bool {
        a.compound==b.compound && a.target==b.target && a.targetVariant==b.targetVariant && a.site==b.site &&
        a.hostContext==b.hostContext && a.temperatureK==b.temperatureK && a.pH==b.pH && a.ionicStrengthM==b.ionicStrengthM
    }

    public static func calculate(_ request:VivoQMMMChemicalStateNetworkRequest)throws->VivoQMMMChemicalStateNetworkResult {
        guard request.schema==VivoQMMMChemicalStateNetworkRequest.schema,!request.identifier.isEmpty,
              request.identifier.utf8.count<=512,!request.states.isEmpty,request.states.count<=128,
              Set(request.states.map(\.identifier)).count==request.states.count,request.rapidPreEquilibrium else {
            throw VivoKineticsError.unsupported("chemical-state aggregation requires a bounded explicit state set and rapid pre-equilibrium")
        }
        for state in request.states { try state.validate() }
        let population=request.states.reduce(0.0){$0+$1.equilibriumPopulation}
        guard population.isFinite,abs(population-1)<=1e-8 else {
            throw VivoKineticsError.invalid("chemical-state equilibrium populations must sum to one")
        }
        guard let firstPath=request.states[0].pathways.first,let firstRate=firstPath.request.replicas.first else {
            throw VivoKineticsError.invalid("chemical-state network has no reference pathway")
        }
        let reference=firstRate.context
        var contributions:[VivoQMMMPathwayContribution]=[],effective=0.0,allCalculated=true
        for state in request.states {
            for pathway in state.pathways {
                guard let rateRequest=pathway.request.replicas.first,
                      sameEnvironment(reference,rateRequest.context),
                      rateRequest.context.chemicalState==state.identifier else {
                    throw VivoKineticsError.invalid("chemical-state pathway context does not match state/environment identity")
                }
                let rate=pathway.result.geometricMeanRatePerSecond
                let contribution=state.equilibriumPopulation*rate
                guard rate.isFinite,rate>0,contribution.isFinite,contribution>=0 else {
                    throw VivoKineticsError.numerical("chemical-state pathway contribution")
                }
                effective+=contribution
                contributions.append(.init(stateIdentifier:state.identifier,pathwayIdentifier:pathway.identifier,
                    conditionalRatePerSecond:rate,statePopulation:state.equilibriumPopulation,contributionPerSecond:contribution))
                if pathway.result.parameter.origin != .calculated { allCalculated=false }
            }
            if state.populationOrigin != .calculated && state.populationOrigin != .measured && state.populationOrigin != .fitted {
                allCalculated=false
            }
        }
        guard effective.isFinite,effective>0 else { throw VivoKineticsError.numerical("chemical-state effective rate") }
        let requestID=try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(request))
        struct Evidence:Codable { let schema:String;let request:VivoQMMMChemicalStateNetworkRequest;let contributions:[VivoQMMMPathwayContribution];let effective:Double }
        let evidenceData=try VivoCanonicalJSON.encode(Evidence(schema:"numivivo.org/qmmm-chemical-state-network-evidence/v1",
            request:request,contributions:contributions,effective:effective))
        let evidenceID=try VivoCanonicalJSON.fingerprint(evidenceData)
        let origin:VivoKineticOrigin = allCalculated ? .calculated:.assumed
        let evidence=VivoKineticEvidence(source:"NumiVivo explicit rapid-pre-equilibrium chemical-state/pathway network",
            locator:"sum_s population_s * sum_pathway conditional_rate",sourceFingerprint:evidenceID.hex)
        let parameter=VivoKineticParameter(value:effective,unit:.perSecond,origin:origin,uncertainty:.unknown,evidence:evidence)
        try parameter.validate(unit:.perSecond,label:"chemical-state effective rate",positive:true)
        return .init(schema:VivoQMMMChemicalStateNetworkResult.schema,requestFingerprint:requestID,
            effectiveRatePerSecond:effective,contributions:contributions,parameter:parameter,
            interpretation:interpretation,evidenceFingerprint:evidenceID)
    }

    public static func validate(_ result:VivoQMMMChemicalStateNetworkResult,
                                request:VivoQMMMChemicalStateNetworkRequest)throws {
        guard result==(try calculate(request)) else { throw VivoKineticsError.invalid("chemical-state network does not reconstruct") }
    }
}
