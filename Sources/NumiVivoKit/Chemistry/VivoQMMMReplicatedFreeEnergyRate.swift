import Foundation

public struct VivoQMMMReplicaAgreementConfiguration: Codable, Sendable, Equatable {
    public var minimumReplicates: Int
    /// Maximum allowed range across independent ln(k) estimates.
    public var maximumLogRateRange: Double
    /// Secondary profile diagnostic; the flux criterion remains authoritative.
    public var maximumProfileBarrierRangeKJPerMol: Double
    public init(minimumReplicates:Int=2,maximumLogRateRange:Double=0.5,
                maximumProfileBarrierRangeKJPerMol:Double=2.0) {
        self.minimumReplicates=minimumReplicates;self.maximumLogRateRange=maximumLogRateRange
        self.maximumProfileBarrierRangeKJPerMol=maximumProfileBarrierRangeKJPerMol
    }
    public func validate() throws {
        guard (2...32).contains(minimumReplicates),maximumLogRateRange.isFinite,maximumLogRateRange>0,
              maximumProfileBarrierRangeKJPerMol.isFinite,maximumProfileBarrierRangeKJPerMol>0 else {
            throw VivoKineticsError.invalid("QM/MM replica-agreement configuration")
        }
    }
}

public struct VivoQMMMReplicatedFreeEnergyRateRequest: Codable, Sendable, Equatable {
    public static let schema="numivivo.org/qmmm-replicated-free-energy-rate/v1"
    public var schema:String
    public var replicas:[VivoQMMMFreeEnergyRateRequest]
    public var agreement:VivoQMMMReplicaAgreementConfiguration
    public init(replicas:[VivoQMMMFreeEnergyRateRequest],agreement:VivoQMMMReplicaAgreementConfiguration = .init()) {
        schema=Self.schema;self.replicas=replicas;self.agreement=agreement
    }
}

public struct VivoQMMMReplicatedFreeEnergyRateResult: Codable, Sendable, Equatable {
    public static let schema="numivivo.org/qmmm-replicated-free-energy-rate-result/v1"
    public let schema:String
    public let requestFingerprint:VivoFingerprint
    public let replicaResults:[VivoQMMMFreeEnergyRateResult]
    public let geometricMeanRatePerSecond:Double
    public let meanLogRatePerSecond:Double
    /// Sample standard deviation of independent replica ln(k) values.
    public let betweenReplicaLogRateStandardDeviation:Double
    /// Quadrature combination of mean within-replica conditional variance and
    /// observed between-replica variance. It remains conditional on one model,
    /// reaction coordinate, chemical state and transmission model.
    public let combinedConditionalLogRateStandardDeviation:Double
    public let logRateRange:Double
    public let profileBarrierRangeKJPerMol:Double
    public let converged:Bool
    public let issues:[String]
    public let parameter:VivoKineticParameter
    public let evidenceFingerprint:VivoFingerprint
    public let evidenceData:Data
}

/// Independent replicas are a qualification layer over the existing single-PMF
/// result. They must use the same Hamiltonian/context/connectivity contract and
/// disjoint stochastic seeds. Agreement is scored in ln(k), because the actual
/// kinetic observable is exponentially sensitive to free energy.
public enum VivoQMMMReplicatedFreeEnergyRate {
    private struct Evidence: Codable {
        let schema:String
        let request:VivoQMMMReplicatedFreeEnergyRateRequest
        let replicaResults:[VivoQMMMFreeEnergyRateResult]
        let meanLogRatePerSecond:Double
        let betweenReplicaLogRateStandardDeviation:Double
        let combinedConditionalLogRateStandardDeviation:Double
        let issues:[String]
    }
    private enum SelfEvidence { static let schema="numivivo.org/qmmm-replicated-free-energy-rate-evidence/v1" }

    private static func environment(_ request:VivoQMMMFreeEnergyRateRequest)->VivoQMMMFreeEnergyEnvironment {
        switch request.environment { case .explicitSolution:return .explicitSolution;case .proteinEnvironment:return .proteinEnvironment }
    }

    public static func calculate(_ request:VivoQMMMReplicatedFreeEnergyRateRequest)throws->VivoQMMMReplicatedFreeEnergyRateResult {
        try request.agreement.validate()
        guard request.schema==VivoQMMMReplicatedFreeEnergyRateRequest.schema,
              request.replicas.count>=request.agreement.minimumReplicates,request.replicas.count<=32 else {
            throw VivoKineticsError.invalid("QM/MM replicated rate requires the declared number of independent replicas")
        }
        var results:[VivoQMMMFreeEnergyRateResult]=[];results.reserveCapacity(request.replicas.count)
        for replica in request.replicas { results.append(try VivoQMMMFreeEnergyRate.calculate(replica)) }
        let first=request.replicas[0],firstFE=first.freeEnergy,firstP=firstFE.provenance
        var allSeeds:[UInt64]=[]
        for replica in request.replicas {
            let fe=replica.freeEnergy,p=fe.provenance
            guard replica.context==first.context,replica.environment==first.environment,
                  replica.transmissionProbability==first.transmissionProbability,
                  replica.transmissionOrigin==first.transmissionOrigin,
                  replica.transmissionEvidence==first.transmissionEvidence,
                  fe.analysis.coordinate==firstFE.analysis.coordinate,
                  fe.analysis.temperatureK==firstFE.analysis.temperatureK,
                  fe.analysis.configuration==firstFE.analysis.configuration,
                  p.structureFingerprint==firstP.structureFingerprint,p.systemFingerprint==firstP.systemFingerprint,
                  p.baseProviderFingerprint==firstP.baseProviderFingerprint,p.dynamicsFingerprint==firstP.dynamicsFingerprint,
                  p.chemicalState==firstP.chemicalState,p.environment==firstP.environment,
                  p.environmentIdentifier==firstP.environmentIdentifier,p.reactionConnectivity==firstP.reactionConnectivity,
                  p.environment==environment(replica) else {
                throw VivoKineticsError.invalid("QM/MM replicas differ in Hamiltonian, context, reaction mapping or transmission model")
            }
            allSeeds.append(contentsOf:fe.analysis.traces.map(\.randomSeed))
        }
        guard Set(allSeeds).count==allSeeds.count else {
            throw VivoKineticsError.invalid("QM/MM independent replicas reuse stochastic window seeds")
        }
        let logs=results.map{ $0.estimate.naturalLogRatePerSecond },mean=logs.reduce(0,+)/Double(logs.count)
        let range=(logs.max() ?? mean)-(logs.min() ?? mean)
        let variance=logs.count>1 ? logs.reduce(0){$0+pow($1-mean,2)}/Double(logs.count-1):0
        let between=sqrt(max(0,variance))
        let withinVariance=results.reduce(0.0) { partial,result in
            let sd=result.estimate.conditionalLogRateStandardDeviation ?? 0
            return partial+sd*sd
        }/Double(results.count)
        let combined=sqrt(max(0,variance+withinVariance))
        let barriers=request.replicas.map{ $0.freeEnergy.analysis.activationFreeEnergyKJPerMol }
        let barrierRange=(barriers.max() ?? 0)-(barriers.min() ?? 0)
        var issues:[String]=[]
        if range>request.agreement.maximumLogRateRange { issues.append("independent replica log-rate range exceeds tolerance") }
        if barrierRange>request.agreement.maximumProfileBarrierRangeKJPerMol {
            issues.append("independent replica PMF profile-barrier range exceeds diagnostic tolerance")
        }
        let rate=exp(mean)
        guard rate.isFinite,rate>0,between.isFinite,combined.isFinite,barrierRange.isFinite else {
            throw VivoKineticsError.numerical("QM/MM replicated log-rate aggregation")
        }
        let evidenceData=try VivoCanonicalJSON.encode(Evidence(schema:SelfEvidence.schema,request:request,replicaResults:results,
            meanLogRatePerSecond:mean,betweenReplicaLogRateStandardDeviation:between,
            combinedConditionalLogRateStandardDeviation:combined,issues:issues))
        let evidenceID=try VivoCanonicalJSON.fingerprint(evidenceData)
        let evidence=VivoKineticEvidence(source:"NumiVivo independently replicated QM/MM PMF rates",
            locator:"geometric mean across disjoint-seed PMF replicas; between-replica agreement retained",
            sourceFingerprint:evidenceID.hex)
        let origin:VivoKineticOrigin = first.transmissionOrigin == .assumed ? .assumed:.calculated
        let uncertainty:VivoKineticUncertainty = combined>0 ? .logNormal(logStandardDeviation:combined):.unknown
        let parameter=VivoKineticParameter(value:rate,unit:.perSecond,origin:origin,uncertainty:uncertainty,evidence:evidence)
        try parameter.validate(unit:.perSecond,label:"replicated QM/MM inactivation",positive:true)
        return .init(schema:VivoQMMMReplicatedFreeEnergyRateResult.schema,
            requestFingerprint:try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(request)),replicaResults:results,
            geometricMeanRatePerSecond:rate,meanLogRatePerSecond:mean,betweenReplicaLogRateStandardDeviation:between,
            combinedConditionalLogRateStandardDeviation:combined,logRateRange:range,profileBarrierRangeKJPerMol:barrierRange,
            converged:issues.isEmpty,issues:issues,parameter:parameter,evidenceFingerprint:evidenceID,evidenceData:evidenceData)
    }

    public static func validate(_ result:VivoQMMMReplicatedFreeEnergyRateResult,
                                request:VivoQMMMReplicatedFreeEnergyRateRequest)throws {
        guard result == (try calculate(request)) else { throw VivoKineticsError.invalid("replicated QM/MM rate does not reconstruct") }
    }

    public static func applying(_ result:VivoQMMMReplicatedFreeEnergyRateResult,
                                request:VivoQMMMReplicatedFreeEnergyRateRequest,
                                to model:VivoCovalentKineticPack)throws->VivoCovalentKineticPack {
        try validate(result,request:request);try model.validate()
        guard result.converged,model.context==request.replicas[0].context else {
            throw VivoKineticsError.invalid("nonconverged replicated QM/MM rate or kinetic context mismatch")
        }
        let output=VivoCovalentKineticPack(identifier:model.identifier,context:model.context,
            association:model.association,dissociation:model.dissociation,inactivation:result.parameter,
            targetTurnover:model.targetTurnover,baselineTarget:model.baselineTarget,
            competitor:model.competitor,maximumUnboundDrugM:model.maximumUnboundDrugM)
        try output.validate();return output
    }
}
