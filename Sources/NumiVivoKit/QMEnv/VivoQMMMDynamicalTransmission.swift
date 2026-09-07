import Foundation
@preconcurrency import Metal

public enum VivoQMMMTransmissionOutcome:String,Codable,Sendable {
    case productCommitted,reactantRecrossed,unresolved
}

public struct VivoQMMMTransmissionAcceptance:Codable,Sendable,Equatable {
    public var minimumSurfaceCheckpoints:Int
    public var minimumEffectiveFluxSamples:Double
    public var maximumUnresolvedFluxFraction:Double
    public var maximumCoefficientStandardError:Double
    public init(minimumSurfaceCheckpoints:Int=32,minimumEffectiveFluxSamples:Double=20,
                maximumUnresolvedFluxFraction:Double=0.05,maximumCoefficientStandardError:Double=0.10) {
        self.minimumSurfaceCheckpoints=minimumSurfaceCheckpoints;self.minimumEffectiveFluxSamples=minimumEffectiveFluxSamples
        self.maximumUnresolvedFluxFraction=maximumUnresolvedFluxFraction;self.maximumCoefficientStandardError=maximumCoefficientStandardError
    }
    public func validate() throws {
        guard (4...100000).contains(minimumSurfaceCheckpoints),minimumEffectiveFluxSamples.isFinite,minimumEffectiveFluxSamples>=2,
              maximumUnresolvedFluxFraction.isFinite,(0..<1).contains(maximumUnresolvedFluxFraction),
              maximumCoefficientStandardError.isFinite,maximumCoefficientStandardError>0,maximumCoefficientStandardError<0.5 else {
            throw VivoChemistryError.invalid("QM/MM dynamical-transmission acceptance configuration")
        }
    }
}

public struct VivoQMMMDynamicalTransmissionRequest:Codable,Sendable,Equatable {
    public static let schema="numivivo.org/qmmm-dynamical-transmission/v1"
    public var schema:String
    public var freeEnergy:VivoQMMMQualifiedActivationFreeEnergy
    public var sourceDynamics:VivoMDConfiguration
    public var shootingDynamics:VivoMDConfiguration
    public var surfaceCheckpoints:[VivoMDCheckpoint]
    public var reactantCommitmentRangeNM:ClosedRange<Double>
    public var productCommitmentRangeNM:ClosedRange<Double>
    public var surfaceToleranceNM:Double
    public var maximumSteps:UInt64
    public var observeEverySteps:UInt64
    public var commitmentObservations:Int
    public var acceptance:VivoQMMMTransmissionAcceptance
    public init(freeEnergy:VivoQMMMQualifiedActivationFreeEnergy,sourceDynamics:VivoMDConfiguration,
                shootingDynamics:VivoMDConfiguration,surfaceCheckpoints:[VivoMDCheckpoint],
                reactantCommitmentRangeNM:ClosedRange<Double>,productCommitmentRangeNM:ClosedRange<Double>,
                surfaceToleranceNM:Double=0.01,maximumSteps:UInt64=5000,observeEverySteps:UInt64=5,
                commitmentObservations:Int=3,acceptance:VivoQMMMTransmissionAcceptance = .init()) {
        schema=Self.schema;self.freeEnergy=freeEnergy;self.sourceDynamics=sourceDynamics;self.shootingDynamics=shootingDynamics
        self.surfaceCheckpoints=surfaceCheckpoints;self.reactantCommitmentRangeNM=reactantCommitmentRangeNM
        self.productCommitmentRangeNM=productCommitmentRangeNM;self.surfaceToleranceNM=surfaceToleranceNM
        self.maximumSteps=maximumSteps;self.observeEverySteps=observeEverySteps
        self.commitmentObservations=commitmentObservations;self.acceptance=acceptance
    }
}

public struct VivoQMMMTransmissionTrajectory:Codable,Sendable,Equatable {
    public let sourceCheckpointFingerprint:VivoFingerprint
    public let initialCoordinateNM:Double
    public let productDirectedCoordinateVelocityNMPerPS:Double
    public let velocitiesTimeReversed:Bool
    public let outcome:VivoQMMMTransmissionOutcome
    public let committedSteps:UInt64
    public let finalCoordinateNM:Double
}

public struct VivoQMMMDynamicalTransmissionResult:Codable,Sendable,Equatable {
    public static let schema="numivivo.org/qmmm-dynamical-transmission-result/v1"
    public let schema:String
    public let requestFingerprint:VivoFingerprint
    public let retainedSystemFingerprint:VivoFingerprint
    public let providerFingerprint:VivoFingerprint
    public let trajectories:[VivoQMMMTransmissionTrajectory]
    public let transmissionCoefficient:Double
    public let coefficientStandardError:Double
    public let effectiveFluxSamples:Double
    public let unresolvedFluxFraction:Double
    public let productCommittedCount:Int
    public let reactantRecrossedCount:Int
    public let unresolvedCount:Int
    public let converged:Bool
    public let issues:[String]
    public let interpretation:String
    public let evidenceFingerprint:VivoFingerprint
}

public enum VivoQMMMDynamicalTransmission {
    public static let interpretation="Flux-weighted classical dynamical transmission coefficient from product-directed, time-reversal-symmetrized near-dividing-surface checkpoints propagated under the same unbiased Born-Oppenheimer Hamiltonian in deterministic NVE. Finite surface-band, finite-time and commitment-basin choices remain explicit. The result measures classical recrossing only; it is not a tunnelling correction or evidence for missing reaction pathways."

    private static func forceEquivalent(source:VivoMDConfiguration,shooting:VivoMDConfiguration)->Bool {
        var normalized=shooting
        normalized.ensemble=source.ensemble;normalized.thermostat=source.thermostat
        normalized.targetTemperatureK=source.targetTemperatureK;normalized.frictionPerPS=source.frictionPerPS
        normalized.randomSeed=source.randomSeed
        return normalized==source
    }

    private static func productDirection(_ request:VivoQMMMDynamicalTransmissionRequest)throws->Double {
        let surface=request.freeEnergy.analysis.configuration.dividingSurfaceNM
        let r=request.reactantCommitmentRangeNM,p=request.productCommitmentRangeNM
        guard r.lowerBound.isFinite,r.upperBound.isFinite,r.lowerBound<r.upperBound,
              p.lowerBound.isFinite,p.upperBound.isFinite,p.lowerBound<p.upperBound,
              !r.overlaps(p),!r.contains(surface),!p.contains(surface) else {
            throw VivoChemistryError.invalid("dynamical-transmission commitment ranges")
        }
        if r.upperBound<surface,p.lowerBound>surface { return 1 }
        if p.upperBound<surface,r.lowerBound>surface { return -1 }
        throw VivoChemistryError.invalid("reactant and product commitment basins must lie on opposite sides of the dividing surface")
    }

    private static func validateStatic(_ request:VivoQMMMDynamicalTransmissionRequest)throws->Double {
        try request.freeEnergy.validate();try request.sourceDynamics.validate();try request.shootingDynamics.validate();try request.acceptance.validate()
        let sourceID=try request.sourceDynamics.fingerprint()
        guard request.schema==VivoQMMMDynamicalTransmissionRequest.schema,
              request.freeEnergy.provenance.dynamicsFingerprint==sourceID,
              request.sourceDynamics.ensemble == .nvt,
              request.sourceDynamics.targetTemperatureK==request.freeEnergy.analysis.temperatureK,
              request.shootingDynamics.ensemble == .nve,request.shootingDynamics.thermostat == .none,
              request.shootingDynamics.targetTemperatureK == nil,request.shootingDynamics.frictionPerPS == nil,
              forceEquivalent(source:request.sourceDynamics,shooting:request.shootingDynamics),
              request.surfaceCheckpoints.count>=request.acceptance.minimumSurfaceCheckpoints,
              request.maximumSteps>0,request.observeEverySteps>0,request.maximumSteps>=request.observeEverySteps,
              request.commitmentObservations>0,request.commitmentObservations<=100000,
              request.surfaceToleranceNM.isFinite,request.surfaceToleranceNM>0 else {
            throw VivoChemistryError.invalid("QM/MM dynamical-transmission dynamics or execution contract")
        }
        return try productDirection(request)
    }

    private static func validateRun(_ request:VivoQMMMDynamicalTransmissionRequest,system:VivoClassicalSystem,
                                    provider:VivoMDCandidateForceProvider)throws->(VivoQMMMResolvedCoordinate,Double) {
        let direction=try validateStatic(request),systemID=try system.fingerprint(),sourceID=try request.sourceDynamics.fingerprint()
        guard request.freeEnergy.provenance.systemFingerprint==systemID,
              request.freeEnergy.provenance.structureFingerprint==system.structureFingerprint,
              request.freeEnergy.provenance.baseProviderFingerprint==provider.fingerprint,
              provider.retainedSystemFingerprint==systemID else {
            throw VivoChemistryError.invalid("QM/MM dynamical transmission differs from the qualified PMF Hamiltonian")
        }
        let coordinate=try VivoQMMMResolvedCoordinate(source:request.freeEnergy.analysis.coordinate,system:system)
        var ids=Set<VivoFingerprint>()
        for checkpoint in request.surfaceCheckpoints {
            try checkpoint.validate(particleCount:system.particles.count)
            guard checkpoint.systemFingerprint==systemID,checkpoint.configurationFingerprint==sourceID else {
                throw VivoChemistryError.invalid("surface checkpoint system or source dynamics mismatch")
            }
            guard ids.insert(try checkpoint.fingerprint()).inserted else { throw VivoChemistryError.invalid("duplicate dividing-surface checkpoint") }
            let geometry=try VivoMDCandidateGeometry(particlePositionsNM:checkpoint.positionsNM,periodicCell:checkpoint.periodicCell)
            let value=try coordinate.evaluate(geometry).valueNM
            guard abs(value-request.freeEnergy.analysis.configuration.dividingSurfaceNM)<=request.surfaceToleranceNM else {
                throw VivoChemistryError.invalid("surface checkpoint lies outside the declared dividing-surface tolerance")
            }
            try provider.validate(system:system,configuration:request.shootingDynamics,cell:checkpoint.periodicCell)
        }
        return (coordinate,direction)
    }

    private static func coordinateVelocity(_ coordinate:VivoQMMMResolvedCoordinate,checkpoint:VivoMDCheckpoint,
                                           direction:Double)throws->Double {
        let geometry=try VivoMDCandidateGeometry(particlePositionsNM:checkpoint.positionsNM,periodicCell:checkpoint.periodicCell)
        let evaluated=try coordinate.evaluate(geometry)
        var velocity=0.0
        for (particle,gradient) in evaluated.gradients {
            let v=checkpoint.velocitiesNMPerPS[Int(particle)]
            velocity += gradient.x*v.x+gradient.y*v.y+gradient.z*v.z
        }
        let value=direction*velocity
        guard value.isFinite else { throw VivoChemistryError.invalid("nonfinite dividing-surface coordinate velocity") }
        return value
    }

    private static func reverseVelocities(_ checkpoint:VivoMDCheckpoint)->VivoMDCheckpoint {
        var result=checkpoint
        result.velocitiesNMPerPS=checkpoint.velocitiesNMPerPS.map{.init(-$0.x,-$0.y,-$0.z)}
        return result
    }

    private static func transition(_ source:VivoMDCheckpoint,request:VivoQMMMDynamicalTransmissionRequest,
                                   particleCount:Int,reverse:Bool)throws->VivoMDCheckpoint {
        let staged=try VivoMDStageTransfer.prepare(checkpoint:source,source:request.sourceDynamics,
            destination:request.shootingDynamics,particleCount:particleCount,velocityInitialization:.preserve)
        let checkpoint=reverse ? reverseVelocities(staged.destinationCheckpoint):staged.destinationCheckpoint
        try checkpoint.validate(particleCount:particleCount);return checkpoint
    }

    public static func run(_ request:VivoQMMMDynamicalTransmissionRequest,system:VivoClassicalSystem,
                           baseProvider:VivoMDCandidateForceProvider,device:MTLDevice?=nil) async throws -> VivoQMMMDynamicalTransmissionResult {
        let (coordinate,direction)=try validateRun(request,system:system,provider:baseProvider)
        var trajectories:[VivoQMMMTransmissionTrajectory]=[];trajectories.reserveCapacity(request.surfaceCheckpoints.count)
        for source in request.surfaceCheckpoints {
            try Task.checkCancellation()
            let sourceGeometry=try VivoMDCandidateGeometry(particlePositionsNM:source.positionsNM,periodicCell:source.periodicCell)
            let initialCoordinate=try coordinate.evaluate(sourceGeometry).valueNM
            let rawVelocity=try coordinateVelocity(coordinate,checkpoint:source,direction:direction)
            let reverse=rawVelocity<0,weight=abs(rawVelocity),sourceID=try source.fingerprint()
            guard weight.isFinite else { throw VivoChemistryError.invalid("surface flux weight") }
            if weight==0 {
                trajectories.append(.init(sourceCheckpointFingerprint:sourceID,initialCoordinateNM:initialCoordinate,
                    productDirectedCoordinateVelocityNMPerPS:0,velocitiesTimeReversed:false,outcome:.unresolved,
                    committedSteps:0,finalCoordinateNM:initialCoordinate));continue
            }
            let start=try transition(source,request:request,particleCount:system.particles.count,reverse:reverse)
            let runtime=try await VivoMDMetalRuntime.restore(system:system,configuration:request.shootingDynamics,
                checkpoint:start,device:device,forceProvider:baseProvider)
            var productStreak=0,reactantStreak=0,outcome:VivoQMMMTransmissionOutcome = .unresolved
            var finalCoordinate=initialCoordinate,committed:UInt64=0
            while committed<request.maximumSteps,outcome == .unresolved {
                try Task.checkCancellation()
                let block=min(request.observeEverySteps,request.maximumSteps-committed)
                for _ in 0..<block {
                    guard try await runtime.step().committed else { throw VivoChemistryError.convergence("dynamical-transmission NVE shooting candidate rejected") }
                    committed+=1
                }
                let sample=try await runtime.sample(includeObservables:false)
                let geometry=try VivoMDCandidateGeometry(particlePositionsNM:sample.state.positionsNM,periodicCell:sample.state.periodicCell)
                finalCoordinate=try coordinate.evaluate(geometry).valueNM
                if request.productCommitmentRangeNM.contains(finalCoordinate) { productStreak+=1;reactantStreak=0 }
                else if request.reactantCommitmentRangeNM.contains(finalCoordinate) { reactantStreak+=1;productStreak=0 }
                else { productStreak=0;reactantStreak=0 }
                if productStreak>=request.commitmentObservations { outcome = .productCommitted }
                else if reactantStreak>=request.commitmentObservations { outcome = .reactantRecrossed }
            }
            trajectories.append(.init(sourceCheckpointFingerprint:sourceID,initialCoordinateNM:initialCoordinate,
                productDirectedCoordinateVelocityNMPerPS:weight,velocitiesTimeReversed:reverse,outcome:outcome,
                committedSteps:committed,finalCoordinateNM:finalCoordinate))
        }
        return try analyze(request:request,systemFingerprint:try system.fingerprint(),providerFingerprint:baseProvider.fingerprint,
                           trajectories:trajectories)
    }

    public static func analyze(request:VivoQMMMDynamicalTransmissionRequest,systemFingerprint:VivoFingerprint,
                               providerFingerprint:VivoFingerprint,trajectories:[VivoQMMMTransmissionTrajectory]) throws -> VivoQMMMDynamicalTransmissionResult {
        _=try validateStatic(request)
        guard systemFingerprint==request.freeEnergy.provenance.systemFingerprint,
              providerFingerprint==request.freeEnergy.provenance.baseProviderFingerprint else {
            throw VivoChemistryError.invalid("dynamical-transmission evidence is not bound to the qualified PMF Hamiltonian")
        }
        let sourceID=try request.sourceDynamics.fingerprint(),surface=request.freeEnergy.analysis.configuration.dividingSurfaceNM
        var expected=Set<VivoFingerprint>()
        for checkpoint in request.surfaceCheckpoints {
            let id=try checkpoint.fingerprint()
            guard checkpoint.systemFingerprint==systemFingerprint,checkpoint.configurationFingerprint==sourceID,
                  expected.insert(id).inserted else { throw VivoChemistryError.invalid("dynamical-transmission source checkpoint identity") }
        }
        let actual=Set(trajectories.map(\.sourceCheckpointFingerprint))
        guard trajectories.count==request.surfaceCheckpoints.count,actual==expected,
              trajectories.allSatisfy({ trajectory in
                  trajectory.initialCoordinateNM.isFinite && abs(trajectory.initialCoordinateNM-surface)<=request.surfaceToleranceNM &&
                  trajectory.productDirectedCoordinateVelocityNMPerPS.isFinite && trajectory.productDirectedCoordinateVelocityNMPerPS>=0 &&
                  trajectory.committedSteps<=request.maximumSteps && trajectory.finalCoordinateNM.isFinite &&
                  (trajectory.outcome != .productCommitted || request.productCommitmentRangeNM.contains(trajectory.finalCoordinateNM)) &&
                  (trajectory.outcome != .reactantRecrossed || request.reactantCommitmentRangeNM.contains(trajectory.finalCoordinateNM))
              }) else { throw VivoChemistryError.invalid("dynamical-transmission trajectory evidence shape or commitment identity") }
        let weights=trajectories.map(\.productDirectedCoordinateVelocityNMPerPS),total=weights.reduce(0,+),square=weights.reduce(0){$0+$1*$1}
        guard total.isFinite,total>0,square.isFinite,square>0 else { throw VivoChemistryError.convergence("dynamical-transmission surface flux is zero") }
        var productWeight=0.0,unresolvedWeight=0.0
        for (trajectory,weight) in zip(trajectories,weights) {
            if trajectory.outcome == .productCommitted { productWeight+=weight }
            else if trajectory.outcome == .unresolved { unresolvedWeight+=weight }
        }
        let coefficient=productWeight/total,effective=total*total/square,unresolved=unresolvedWeight/total
        let standardError=sqrt(max(0,(coefficient*(1-coefficient)+1/(4*effective))/effective))
        guard coefficient.isFinite,coefficient>=0,coefficient<=1,effective.isFinite,standardError.isFinite,unresolved.isFinite else {
            throw VivoChemistryError.convergence("dynamical-transmission estimator is nonfinite")
        }
        var issues:[String]=[]
        if effective<request.acceptance.minimumEffectiveFluxSamples { issues.append("effective flux-weighted surface sample count is below threshold") }
        if unresolved>request.acceptance.maximumUnresolvedFluxFraction { issues.append("unresolved finite-time flux fraction exceeds threshold") }
        if standardError>request.acceptance.maximumCoefficientStandardError { issues.append("transmission-coefficient standard error exceeds threshold") }
        if coefficient<=0 { issues.append("no product-committing positive flux was observed") }
        let requestID=try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(request))
        struct Evidence:Codable {
            let schema:String;let requestFingerprint:VivoFingerprint;let system:VivoFingerprint;let provider:VivoFingerprint
            let trajectories:[VivoQMMMTransmissionTrajectory];let coefficient:Double;let standardError:Double;let effective:Double;let unresolved:Double
        }
        let evidenceID=try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(Evidence(
            schema:"numivivo.org/qmmm-dynamical-transmission-evidence/v1",requestFingerprint:requestID,
            system:systemFingerprint,provider:providerFingerprint,trajectories:trajectories,coefficient:coefficient,
            standardError:standardError,effective:effective,unresolved:unresolved)))
        return .init(schema:VivoQMMMDynamicalTransmissionResult.schema,requestFingerprint:requestID,
            retainedSystemFingerprint:systemFingerprint,providerFingerprint:providerFingerprint,trajectories:trajectories,
            transmissionCoefficient:coefficient,coefficientStandardError:standardError,effectiveFluxSamples:effective,
            unresolvedFluxFraction:unresolved,productCommittedCount:trajectories.filter{$0.outcome == .productCommitted}.count,
            reactantRecrossedCount:trajectories.filter{$0.outcome == .reactantRecrossed}.count,
            unresolvedCount:trajectories.filter{$0.outcome == .unresolved}.count,converged:issues.isEmpty,issues:issues,
            interpretation:interpretation,evidenceFingerprint:evidenceID)
    }

    public static func validate(_ result:VivoQMMMDynamicalTransmissionResult,
                                request:VivoQMMMDynamicalTransmissionRequest)throws {
        let rebuilt=try analyze(request:request,systemFingerprint:result.retainedSystemFingerprint,
                                providerFingerprint:result.providerFingerprint,trajectories:result.trajectories)
        guard rebuilt==result else { throw VivoChemistryError.invalid("dynamical-transmission result does not reconstruct") }
    }

    public static func applying(_ result:VivoQMMMDynamicalTransmissionResult,
                                transmissionRequest:VivoQMMMDynamicalTransmissionRequest,
                                to rateRequest:VivoQMMMFreeEnergyRateRequest)throws->VivoQMMMFreeEnergyRateRequest {
        try validate(result,request:transmissionRequest)
        guard result.converged,result.transmissionCoefficient>0,
              rateRequest.freeEnergy.evidenceFingerprint==transmissionRequest.freeEnergy.evidenceFingerprint else {
            throw VivoKineticsError.invalid("nonconverged dynamical transmission or PMF evidence mismatch")
        }
        var output=rateRequest
        output.transmissionProbability=result.transmissionCoefficient;output.transmissionOrigin=.calculated
        output.transmissionEvidence=VivoKineticEvidence(source:"NumiVivo classical dividing-surface shooting",
            locator:"flux-weighted NVE recrossing coefficient; SE=\(result.coefficientStandardError); effectiveFluxSamples=\(result.effectiveFluxSamples)",
            sourceFingerprint:result.evidenceFingerprint.hex)
        return output
    }
}
