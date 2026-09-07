import Foundation

public enum VivoQMMMQualificationDimension:String,Codable,Sendable,CaseIterable {
    case samplingProtocol
    case reactionCoordinate
    case electronicModel
    case qmRegion
    case periodicFiniteSize
}

public struct VivoQMMMSensitivityCriterion:Codable,Sendable,Equatable {
    public var dimension:VivoQMMMQualificationDimension
    public var requiredVariants:Int
    public var maximumAbsoluteLogRateShift:Double
    public init(dimension:VivoQMMMQualificationDimension,requiredVariants:Int=1,maximumAbsoluteLogRateShift:Double=0.2) {
        self.dimension=dimension;self.requiredVariants=requiredVariants;self.maximumAbsoluteLogRateShift=maximumAbsoluteLogRateShift
    }
    public func validate()throws {
        guard (1...128).contains(requiredVariants),maximumAbsoluteLogRateShift.isFinite,maximumAbsoluteLogRateShift>0 else {
            throw VivoKineticsError.invalid("QM/MM sensitivity criterion")
        }
    }
}

public struct VivoQMMMQualificationVariant:Codable,Sendable,Equatable {
    public var identifier:String
    public var dimension:VivoQMMMQualificationDimension
    public var request:VivoQMMMReplicatedFreeEnergyRateRequest
    public var result:VivoQMMMReplicatedFreeEnergyRateResult
    public var description:String
    public init(identifier:String,dimension:VivoQMMMQualificationDimension,
                request:VivoQMMMReplicatedFreeEnergyRateRequest,result:VivoQMMMReplicatedFreeEnergyRateResult,
                description:String) {
        self.identifier=identifier;self.dimension=dimension;self.request=request;self.result=result;self.description=description
    }
    public func validate()throws {
        guard !identifier.trimmingCharacters(in:.whitespacesAndNewlines).isEmpty,identifier.utf8.count<=512,
              !description.trimmingCharacters(in:.whitespacesAndNewlines).isEmpty,description.utf8.count<=8192 else {
            throw VivoKineticsError.invalid("QM/MM qualification variant identity")
        }
        try VivoQMMMReplicatedFreeEnergyRate.validate(result,request:request)
        guard result.converged else { throw VivoKineticsError.invalid("QM/MM qualification variant is not independently replicated/converged") }
    }
}

public enum VivoQMMMExperimentalObservable:String,Codable,Sendable {
    /// A first-order chemical inactivation rate for the explicitly bound pre-reactive complex.
    case chemicalInactivationRatePerSecond
    /// Retained as external context only; it is not directly comparable to the conditional chemical rate.
    case inactivationEfficiencyPerMolarSecond
    case inhibitoryConcentrationMolar
    case targetOccupancyFraction
}

public struct VivoQMMMExternalValidationTarget:Codable,Sendable,Equatable {
    public var identifier:String
    public var observable:VivoQMMMExperimentalObservable
    public var value:Double
    public var logStandardDeviation:Double?
    public var context:VivoKineticContext
    public var source:VivoKineticEvidence
    public init(identifier:String,observable:VivoQMMMExperimentalObservable,value:Double,
                logStandardDeviation:Double?=nil,context:VivoKineticContext,source:VivoKineticEvidence) {
        self.identifier=identifier;self.observable=observable;self.value=value
        self.logStandardDeviation=logStandardDeviation;self.context=context;self.source=source
    }
    public func validate()throws {
        try context.validate();try source.validate(origin:.experimental)
        guard !identifier.isEmpty,value.isFinite,value>0 else { throw VivoKineticsError.invalid("QM/MM external validation target") }
        if let logStandardDeviation {
            guard logStandardDeviation.isFinite,logStandardDeviation>=0 else { throw VivoKineticsError.invalid("experimental log uncertainty") }
        }
    }
}

public struct VivoQMMMChemicalQualificationRequest:Codable,Sendable,Equatable {
    public static let schema="numivivo.org/qmmm-chemical-qualification/v1"
    public var schema:String
    public var identifier:String
    public var baselineRequest:VivoQMMMReplicatedFreeEnergyRateRequest
    public var baselineResult:VivoQMMMReplicatedFreeEnergyRateResult
    public var variants:[VivoQMMMQualificationVariant]
    public var criteria:[VivoQMMMSensitivityCriterion]
    public var externalTargets:[VivoQMMMExternalValidationTarget]
    public var maximumDirectValidationAbsoluteLogError:Double
    public init(identifier:String,baselineRequest:VivoQMMMReplicatedFreeEnergyRateRequest,
                baselineResult:VivoQMMMReplicatedFreeEnergyRateResult,variants:[VivoQMMMQualificationVariant],
                criteria:[VivoQMMMSensitivityCriterion],externalTargets:[VivoQMMMExternalValidationTarget]=[],
                maximumDirectValidationAbsoluteLogError:Double=0.7) {
        schema=Self.schema;self.identifier=identifier;self.baselineRequest=baselineRequest;self.baselineResult=baselineResult
        self.variants=variants;self.criteria=criteria;self.externalTargets=externalTargets
        self.maximumDirectValidationAbsoluteLogError=maximumDirectValidationAbsoluteLogError
    }
}

public struct VivoQMMMSensitivityResult:Codable,Sendable,Equatable {
    public let identifier:String
    public let dimension:VivoQMMMQualificationDimension
    public let logRateShiftFromBaseline:Double
    public let absoluteLogRateShift:Double
    public let passed:Bool
}

public struct VivoQMMMExternalValidationResult:Codable,Sendable,Equatable {
    public let identifier:String
    public let observable:VivoQMMMExperimentalObservable
    public let comparable:Bool
    public let absoluteLogError:Double?
    public let passed:Bool?
    public let reason:String
}

public struct VivoQMMMChemicalQualificationResult:Codable,Sendable,Equatable {
    public static let schema="numivivo.org/qmmm-chemical-qualification-result/v1"
    public let schema:String
    public let requestFingerprint:VivoFingerprint
    public let baselineLogRatePerSecond:Double
    public let sensitivity:[VivoQMMMSensitivityResult]
    public let externalValidation:[VivoQMMMExternalValidationResult]
    public let converged:Bool
    public let issues:[String]
    public let interpretation:String
    public let evidenceFingerprint:VivoFingerprint
}

/// This gate qualifies a declared computational protocol. It deliberately keeps
/// independent-replica dispersion, method sensitivity and external validation as
/// separate evidence categories. Alternate protonation states/pathways belong in
/// an explicit state/pathway network and are therefore not accepted as ordinary
/// sensitivity variants here.
public enum VivoQMMMChemicalQualification {
    public static let interpretation="Qualification of one explicitly declared bound-complex chemical-rate protocol against independently replicated sampling, predeclared numerical/electronic/QM-region sensitivity variants and condition-matched external rate data. Passing this gate does not establish transferability to other proteins, chemical states or mechanisms."

    private static func baselineInvariant(_ baseline:VivoQMMMFreeEnergyRateRequest,
                                          _ variant:VivoQMMMFreeEnergyRateRequest,
                                          dimension:VivoQMMMQualificationDimension)throws {
        let a=baseline.freeEnergy.provenance,b=variant.freeEnergy.provenance
        guard baseline.context==variant.context,baseline.environment==variant.environment,
              a.chemicalState==b.chemicalState,a.environment==b.environment,
              a.environmentIdentifier==b.environmentIdentifier,
              a.structureFingerprint==b.structureFingerprint else {
            throw VivoKineticsError.invalid("QM/MM sensitivity variant changes kinetic context, chemical state or host structure")
        }
        switch dimension {
        case .samplingProtocol:
            guard a.systemFingerprint==b.systemFingerprint,a.baseProviderFingerprint==b.baseProviderFingerprint,
                  baseline.freeEnergy.analysis.coordinate==variant.freeEnergy.analysis.coordinate,
                  a.reactionConnectivity==b.reactionConnectivity else {
                throw VivoKineticsError.invalid("sampling-protocol sensitivity changed Hamiltonian, coordinate or reaction")
            }
        case .reactionCoordinate:
            guard a.systemFingerprint==b.systemFingerprint,a.baseProviderFingerprint==b.baseProviderFingerprint,
                  a.reactionConnectivity.reactantEndpointIdentifier==b.reactionConnectivity.reactantEndpointIdentifier,
                  a.reactionConnectivity.productEndpointIdentifier==b.reactionConnectivity.productEndpointIdentifier else {
                throw VivoKineticsError.invalid("coordinate sensitivity changed Hamiltonian or endpoint identity")
            }
        case .electronicModel,.qmRegion:
            guard a.reactionConnectivity.reactantEndpointIdentifier==b.reactionConnectivity.reactantEndpointIdentifier,
                  a.reactionConnectivity.productEndpointIdentifier==b.reactionConnectivity.productEndpointIdentifier else {
                throw VivoKineticsError.invalid("electronic/QM-region sensitivity changed reaction endpoint identity")
            }
        case .periodicFiniteSize:
            guard baseline.freeEnergy.analysis.coordinate==variant.freeEnergy.analysis.coordinate,
                  a.reactionConnectivity.reactantEndpointIdentifier==b.reactionConnectivity.reactantEndpointIdentifier,
                  a.reactionConnectivity.productEndpointIdentifier==b.reactionConnectivity.productEndpointIdentifier else {
                throw VivoKineticsError.invalid("finite-size sensitivity changed coordinate or reaction endpoint identity")
            }
        }
    }

    public static func calculate(_ request:VivoQMMMChemicalQualificationRequest)throws->VivoQMMMChemicalQualificationResult {
        try VivoQMMMReplicatedFreeEnergyRate.validate(request.baselineResult,request:request.baselineRequest)
        guard request.schema==VivoQMMMChemicalQualificationRequest.schema,!request.identifier.isEmpty,
              request.identifier.utf8.count<=512,request.baselineResult.converged,
              request.maximumDirectValidationAbsoluteLogError.isFinite,request.maximumDirectValidationAbsoluteLogError>0,
              Set(request.variants.map(\.identifier)).count==request.variants.count,
              Set(request.criteria.map(\.dimension)).count==request.criteria.count else {
            throw VivoKineticsError.invalid("QM/MM chemical qualification identity, baseline or criteria")
        }
        for criterion in request.criteria { try criterion.validate() }
        for target in request.externalTargets { try target.validate() }
        guard let baselineRate=request.baselineRequest.replicas.first else { throw VivoKineticsError.invalid("empty replicated baseline") }
        let criteria=Dictionary(uniqueKeysWithValues:request.criteria.map{($0.dimension,$0)})
        var sensitivity:[VivoQMMMSensitivityResult]=[],issues:[String]=[]
        for variant in request.variants {
            try variant.validate()
            guard let candidate=variant.request.replicas.first,let criterion=criteria[variant.dimension] else {
                throw VivoKineticsError.invalid("qualification variant has no predeclared sensitivity criterion")
            }
            try baselineInvariant(baselineRate,candidate,dimension:variant.dimension)
            let shift=variant.result.meanLogRatePerSecond-request.baselineResult.meanLogRatePerSecond
            let passed=abs(shift)<=criterion.maximumAbsoluteLogRateShift
            sensitivity.append(.init(identifier:variant.identifier,dimension:variant.dimension,
                logRateShiftFromBaseline:shift,absoluteLogRateShift:abs(shift),passed:passed))
            if !passed { issues.append("\(variant.dimension.rawValue) variant \(variant.identifier) exceeds the predeclared log-rate sensitivity tolerance") }
        }
        for criterion in request.criteria {
            let count=sensitivity.filter{$0.dimension==criterion.dimension}.count
            if count<criterion.requiredVariants { issues.append("\(criterion.dimension.rawValue) has \(count) variants; \(criterion.requiredVariants) required") }
        }
        var external:[VivoQMMMExternalValidationResult]=[]
        for target in request.externalTargets {
            let sameContext=target.context==baselineRate.context
            switch target.observable {
            case .chemicalInactivationRatePerSecond where sameContext:
                let error=abs(request.baselineResult.meanLogRatePerSecond-log(target.value))
                let pass=error<=request.maximumDirectValidationAbsoluteLogError
                external.append(.init(identifier:target.identifier,observable:target.observable,comparable:true,
                    absoluteLogError:error,passed:pass,reason:"condition-matched first-order chemical inactivation rate"))
                if !pass { issues.append("external chemical-rate validation \(target.identifier) exceeds the predeclared absolute log-error tolerance") }
            case .chemicalInactivationRatePerSecond:
                external.append(.init(identifier:target.identifier,observable:target.observable,comparable:false,
                    absoluteLogError:nil,passed:nil,reason:"experimental kinetic context differs from the simulated bound complex"))
            default:
                external.append(.init(identifier:target.identifier,observable:target.observable,comparable:false,
                    absoluteLogError:nil,passed:nil,reason:"observable is not the same conditional first-order chemical rate"))
            }
        }
        let requestID=try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(request))
        struct Evidence:Codable { let schema:String;let request:VivoQMMMChemicalQualificationRequest;let sensitivity:[VivoQMMMSensitivityResult];let external:[VivoQMMMExternalValidationResult];let issues:[String] }
        let evidenceID=try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(Evidence(
            schema:"numivivo.org/qmmm-chemical-qualification-evidence/v1",request:request,sensitivity:sensitivity,external:external,issues:issues)))
        return .init(schema:VivoQMMMChemicalQualificationResult.schema,requestFingerprint:requestID,
            baselineLogRatePerSecond:request.baselineResult.meanLogRatePerSecond,sensitivity:sensitivity,
            externalValidation:external,converged:issues.isEmpty,issues:issues,
            interpretation:interpretation,evidenceFingerprint:evidenceID)
    }

    public static func validate(_ result:VivoQMMMChemicalQualificationResult,
                                request:VivoQMMMChemicalQualificationRequest)throws {
        guard result==(try calculate(request)) else { throw VivoKineticsError.invalid("QM/MM chemical qualification does not reconstruct") }
    }
}
