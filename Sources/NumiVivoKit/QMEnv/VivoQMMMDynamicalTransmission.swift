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
    public var minimumPlateauObservations:Int
    public var maximumPlateauRange:Double
    public init(minimumSurfaceCheckpoints:Int=32,minimumEffectiveFluxSamples:Double=20,
                maximumUnresolvedFluxFraction:Double=0.05,maximumCoefficientStandardError:Double=0.10,
                minimumPlateauObservations:Int=4,maximumPlateauRange:Double=0.10) {
        self.minimumSurfaceCheckpoints=minimumSurfaceCheckpoints;self.minimumEffectiveFluxSamples=minimumEffectiveFluxSamples
        self.maximumUnresolvedFluxFraction=maximumUnresolvedFluxFraction;self.maximumCoefficientStandardError=maximumCoefficientStandardError
        self.minimumPlateauObservations=minimumPlateauObservations;self.maximumPlateauRange=maximumPlateauRange
    }
    public func validate() throws {
        guard (4...100000).contains(minimumSurfaceCheckpoints),minimumEffectiveFluxSamples.isFinite,minimumEffectiveFluxSamples>=2,
              maximumUnresolvedFluxFraction.isFinite,(0..<1).contains(maximumUnresolvedFluxFraction),
              maximumCoefficientStandardError.isFinite,maximumCoefficientStandardError>0,maximumCoefficientStandardError<0.5,
              (2...10000).contains(minimumPlateauObservations),maximumPlateauRange.isFinite,
              maximumPlateauRange>0,maximumPlateauRange<=1 else {
            throw VivoChemistryError.invalid("QM/MM dynamical-transmission acceptance configuration")
        }
    }
}

public struct VivoQMMMDynamicalTransmissionRequest:Codable,Sendable,Equatable {
    public static let schema="numivivo.org/qmmm-dynamical-transmission/v3"
    public var schema:String
    public var surfaceRequest:VivoQMMMSurfaceEnsembleRequest
    public var surfaceEnsemble:VivoQMMMSurfaceEnsembleResult
    public var shootingDynamics:VivoMDConfiguration
    public var reactantCommitmentRangeNM:ClosedRange<Double>
    public var productCommitmentRangeNM:ClosedRange<Double>
    public var maximumSteps:UInt64
    public var observeEverySteps:UInt64
    public var commitmentObservations:Int
    public var acceptance:VivoQMMMTransmissionAcceptance
    public var freeEnergy:VivoQMMMQualifiedActivationFreeEnergy { surfaceRequest.freeEnergy }
    public init(surfaceRequest:VivoQMMMSurfaceEnsembleRequest,surfaceEnsemble:VivoQMMMSurfaceEnsembleResult,
                shootingDynamics:VivoMDConfiguration,reactantCommitmentRangeNM:ClosedRange<Double>,
                productCommitmentRangeNM:ClosedRange<Double>,maximumSteps:UInt64=5000,observeEverySteps:UInt64=5,
                commitmentObservations:Int=3,acceptance:VivoQMMMTransmissionAcceptance = .init()) {
        schema=Self.schema;self.surfaceRequest=surfaceRequest;self.surfaceEnsemble=surfaceEnsemble
        self.shootingDynamics=shootingDynamics;self.reactantCommitmentRangeNM=reactantCommitmentRangeNM
        self.productCommitmentRangeNM=productCommitmentRangeNM;self.maximumSteps=maximumSteps
        self.observeEverySteps=observeEverySteps;self.commitmentObservations=commitmentObservations;self.acceptance=acceptance
    }
}

public struct VivoQMMMTransmissionObservation:Codable,Sendable,Equatable {
    public let step:UInt64
    public let coordinateNM:Double
    public let inProduct:Bool
    public let inReactant:Bool
}

public struct VivoQMMMTransmissionBranch:Codable,Sendable,Equatable {
    public let velocitiesTimeReversed:Bool
    public let observations:[VivoQMMMTransmissionObservation]
    public let outcome:VivoQMMMTransmissionOutcome
    public let finalCoordinateNM:Double
}

public struct VivoQMMMTransmissionPair:Codable,Sendable,Equatable {
    public let sourceStateFingerprint:VivoFingerprint
    public let initialCoordinateNM:Double
    public let positiveCoordinateVelocityNMPerPS:Double
    public let positiveBranch:VivoQMMMTransmissionBranch
    public let negativeBranch:VivoQMMMTransmissionBranch
}

public struct VivoQMMMTransmissionPlateauPoint:Codable,Sendable,Equatable {
    public let step:UInt64
    public let coefficient:Double
}

public struct VivoQMMMDynamicalTransmissionResult:Codable,Sendable,Equatable {
    public static let schema="numivivo.org/qmmm-dynamical-transmission-result/v3"
    public let schema:String
    public let requestFingerprint:VivoFingerprint
    public let retainedSystemFingerprint:VivoFingerprint
    public let providerFingerprint:VivoFingerprint
    public let pairs:[VivoQMMMTransmissionPair]
    public let coefficientHistory:[VivoQMMMTransmissionPlateauPoint]
    public let transmissionCoefficient:Double
    public let plateauRange:Double
    public let coefficientStandardError:Double
    public let effectiveFluxSamples:Double
    public let unresolvedFluxFraction:Double
    public let positiveProductCommittedCount:Int
    public let negativeReactantCommittedCount:Int
    public let unresolvedPairCount:Int
    public let converged:Bool
    public let issues:[String]
    public let interpretation:String
    public let evidenceFingerprint:VivoFingerprint
}

/// Bennett-Chandler-style paired reactive-flux estimator. Every canonical
/// near-surface state produces a product-directed velocity and its exact time
/// reverse. The signed score h_P(+v,t)-h_P(-v,t) is weighted by the finite-band
/// surface importance weight and positive coordinate speed. The reported kappa
/// is the mean of a required late-time plateau, not an eventual commitment ratio.
public enum VivoQMMMDynamicalTransmission {
    public static let interpretation="Paired signed classical reactive-flux transmission coefficient from a statistically qualified, finite-band-reweighted dividing-surface ensemble. Exact +v/-v NVE trajectories are propagated under the same unbiased Born-Oppenheimer Hamiltonian and kappa is accepted only from a stable late-time plateau. This measures classical recrossing only; tunnelling, alternate mechanisms and unresolved chemical-state populations remain separate."

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

    private static func expectedObservationSteps(_ request:VivoQMMMDynamicalTransmissionRequest)->[UInt64] {
        var values:[UInt64]=[],step:UInt64=0
        while step<request.maximumSteps {
            step+=min(request.observeEverySteps,request.maximumSteps-step);values.append(step)
        }
        return values
    }

    private static func validateStatic(_ request:VivoQMMMDynamicalTransmissionRequest)throws->Double {
        try VivoQMMMSurfaceEnsemble.validate(request.surfaceEnsemble,request:request.surfaceRequest)
        try request.shootingDynamics.validate();try request.acceptance.validate()
        let source=request.surfaceRequest.dynamics,observationCount=expectedObservationSteps(request).count
        guard request.schema==VivoQMMMDynamicalTransmissionRequest.schema,request.surfaceEnsemble.converged,
              request.surfaceEnsemble.states.count>=request.acceptance.minimumSurfaceCheckpoints,
              request.shootingDynamics.ensemble == .nve,request.shootingDynamics.thermostat == .none,
              request.shootingDynamics.targetTemperatureK == nil,request.shootingDynamics.frictionPerPS == nil,
              forceEquivalent(source:source,shooting:request.shootingDynamics),
              request.maximumSteps>0,request.observeEverySteps>0,request.maximumSteps>=request.observeEverySteps,
              request.commitmentObservations>0,request.commitmentObservations<=observationCount,
              request.acceptance.minimumPlateauObservations<=observationCount else {
            throw VivoChemistryError.invalid("QM/MM dynamical-transmission surface evidence, dynamics or execution contract")
        }
        return try productDirection(request)
    }

    private static func coordinateVelocity(_ coordinate:VivoQMMMResolvedCoordinate,state:VivoQMMMSurfaceState,
                                           direction:Double)throws->Double {
        let geometry=try VivoMDCandidateGeometry(particlePositionsNM:state.positionsNM,periodicCell:state.periodicCell)
        let evaluated=try coordinate.evaluate(geometry)
        var velocity=0.0
        for (particle,gradient) in evaluated.gradients {
            let v=state.velocitiesNMPerPS[Int(particle)]
            velocity += gradient.x*v.x+gradient.y*v.y+gradient.z*v.z
        }
        let value=direction*velocity
        guard value.isFinite else { throw VivoChemistryError.invalid("nonfinite dividing-surface coordinate velocity") }
        return value
    }

    private static func validateRun(_ request:VivoQMMMDynamicalTransmissionRequest,system:VivoClassicalSystem,
                                    provider:VivoMDCandidateForceProvider)throws->(VivoQMMMResolvedCoordinate,Double) {
        let direction=try validateStatic(request),systemID=try system.fingerprint(),freeEnergy=request.freeEnergy
        guard freeEnergy.provenance.systemFingerprint==systemID,
              freeEnergy.provenance.structureFingerprint==system.structureFingerprint,
              freeEnergy.provenance.baseProviderFingerprint==provider.fingerprint,
              provider.retainedSystemFingerprint==systemID,
              request.surfaceEnsemble.retainedSystemFingerprint==systemID,
              request.surfaceEnsemble.baseProviderFingerprint==provider.fingerprint else {
            throw VivoChemistryError.invalid("QM/MM dynamical transmission differs from the qualified PMF Hamiltonian")
        }
        let coordinate=try VivoQMMMResolvedCoordinate(source:freeEnergy.analysis.coordinate,system:system)
        var sourceDynamics=request.surfaceRequest.dynamics;sourceDynamics.randomSeed=request.surfaceRequest.randomSeed
        let surfaceProvider=try VivoQMMMUmbrellaBias.provider(base:provider,system:system,coordinate:freeEnergy.analysis.coordinate,
                                                              window:request.surfaceRequest.surfaceWindow)
        let surfaceExecution=try VivoMDCandidateForceProvider.executionFingerprint(configuration:sourceDynamics,provider:surfaceProvider)
        guard surfaceProvider.fingerprint==request.surfaceEnsemble.surfaceProviderFingerprint,
              surfaceExecution==request.surfaceEnsemble.surfaceExecutionFingerprint else {
            throw VivoChemistryError.invalid("surface ensemble provider/execution identity differs from the requested umbrella Hamiltonian")
        }
        for state in request.surfaceEnsemble.states {
            guard state.positionsNM.count==system.particles.count,state.velocitiesNMPerPS.count==system.particles.count else {
                throw VivoChemistryError.invalid("surface state particle shape")
            }
            let geometry=try VivoMDCandidateGeometry(particlePositionsNM:state.positionsNM,periodicCell:state.periodicCell)
            let value=try coordinate.evaluate(geometry).valueNM
            let velocity=try coordinateVelocity(coordinate,state:state,direction:1)
            guard abs(value-state.coordinateNM)<=1e-10,
                  abs(velocity-state.coordinateVelocityNMPerPS)<=max(1e-10,abs(velocity)*1e-9) else {
                throw VivoChemistryError.invalid("surface state coordinate/Jacobian velocity differs from mapped geometry")
            }
            try provider.validate(system:system,configuration:request.shootingDynamics,cell:state.periodicCell)
        }
        return (coordinate,direction)
    }

    private static func outcome(_ observations:[VivoQMMMTransmissionObservation],commitment:Int)->VivoQMMMTransmissionOutcome {
        guard observations.count>=commitment else { return .unresolved }
        let tail=observations.suffix(commitment)
        if tail.allSatisfy(\.inProduct) { return .productCommitted }
        if tail.allSatisfy(\.inReactant) { return .reactantRecrossed }
        return .unresolved
    }

    private static func propagate(state:VivoQMMMSurfaceState,velocities:[VivoVector3D],reversed:Bool,
                                  request:VivoQMMMDynamicalTransmissionRequest,coordinate:VivoQMMMResolvedCoordinate,
                                  system:VivoClassicalSystem,provider:VivoMDCandidateForceProvider,systemID:VivoFingerprint,
                                  device:MTLDevice?) async throws->VivoQMMMTransmissionBranch {
        let initial=VivoClassicalInitialState(systemFingerprint:systemID,positionsNM:state.positionsNM,
                                               periodicCell:state.periodicCell,sourceTimePS:state.timePS)
        let runtime=try await VivoMDMetalRuntime.make(system:system,initialState:initial,configuration:request.shootingDynamics,
            initialVelocitiesNMPerPS:velocities,device:device,forceProvider:provider)
        var observations:[VivoQMMMTransmissionObservation]=[],committed:UInt64=0,final=state.coordinateNM
        while committed<request.maximumSteps {
            try Task.checkCancellation()
            let block=min(request.observeEverySteps,request.maximumSteps-committed)
            for _ in 0..<block {
                guard try await runtime.step().committed else { throw VivoChemistryError.convergence("dynamical-transmission NVE shooting candidate rejected") }
                committed+=1
            }
            let sample=try await runtime.sample(includeObservables:false)
            let geometry=try VivoMDCandidateGeometry(particlePositionsNM:sample.state.positionsNM,periodicCell:sample.state.periodicCell)
            final=try coordinate.evaluate(geometry).valueNM
            observations.append(.init(step:committed,coordinateNM:final,
                inProduct:request.productCommitmentRangeNM.contains(final),inReactant:request.reactantCommitmentRangeNM.contains(final)))
        }
        return .init(velocitiesTimeReversed:reversed,observations:observations,
                     outcome:outcome(observations,commitment:request.commitmentObservations),finalCoordinateNM:final)
    }

    public static func run(_ request:VivoQMMMDynamicalTransmissionRequest,system:VivoClassicalSystem,
                           baseProvider:VivoMDCandidateForceProvider,device:MTLDevice?=nil) async throws -> VivoQMMMDynamicalTransmissionResult {
        let (coordinate,direction)=try validateRun(request,system:system,provider:baseProvider),systemID=try system.fingerprint()
        var pairs:[VivoQMMMTransmissionPair]=[];pairs.reserveCapacity(request.surfaceEnsemble.states.count)
        for state in request.surfaceEnsemble.states {
            try Task.checkCancellation()
            let directed=try coordinateVelocity(coordinate,state:state,direction:direction)
            let speed=abs(directed)
            guard speed.isFinite else { throw VivoChemistryError.invalid("surface positive flux speed") }
            let original=state.velocitiesNMPerPS,reversedVelocities=original.map{VivoVector3D(-$0.x,-$0.y,-$0.z)}
            let positiveVelocities=directed>=0 ? original:reversedVelocities
            let negativeVelocities=directed>=0 ? reversedVelocities:original
            let positive=try await propagate(state:state,velocities:positiveVelocities,reversed:directed<0,request:request,
                coordinate:coordinate,system:system,provider:baseProvider,systemID:systemID,device:device)
            let negative=try await propagate(state:state,velocities:negativeVelocities,reversed:directed>=0,request:request,
                coordinate:coordinate,system:system,provider:baseProvider,systemID:systemID,device:device)
            pairs.append(.init(sourceStateFingerprint:state.stateFingerprint,initialCoordinateNM:state.coordinateNM,
                positiveCoordinateVelocityNMPerPS:speed,positiveBranch:positive,negativeBranch:negative))
        }
        return try analyze(request:request,systemFingerprint:systemID,providerFingerprint:baseProvider.fingerprint,pairs:pairs)
    }

    public static func analyze(request:VivoQMMMDynamicalTransmissionRequest,systemFingerprint:VivoFingerprint,
                               providerFingerprint:VivoFingerprint,pairs:[VivoQMMMTransmissionPair]) throws -> VivoQMMMDynamicalTransmissionResult {
        _=try validateStatic(request)
        guard systemFingerprint==request.freeEnergy.provenance.systemFingerprint,
              providerFingerprint==request.freeEnergy.provenance.baseProviderFingerprint else {
            throw VivoChemistryError.invalid("dynamical-transmission evidence is not bound to the qualified PMF Hamiltonian")
        }
        let states=Dictionary(uniqueKeysWithValues:request.surfaceEnsemble.states.map{($0.stateFingerprint,$0)})
        let expected=Set(states.keys),actual=Set(pairs.map(\.sourceStateFingerprint)),steps=expectedObservationSteps(request)
        guard pairs.count==request.surfaceEnsemble.states.count,actual==expected else {
            throw VivoChemistryError.invalid("dynamical-transmission pair/source identity")
        }
        for pair in pairs {
            guard let state=states[pair.sourceStateFingerprint],pair.initialCoordinateNM.isFinite,
                  abs(pair.initialCoordinateNM-state.coordinateNM)<=1e-10,
                  pair.positiveCoordinateVelocityNMPerPS.isFinite,pair.positiveCoordinateVelocityNMPerPS>=0 else {
                throw VivoChemistryError.invalid("dynamical-transmission initial pair evidence")
            }
            for branch in [pair.positiveBranch,pair.negativeBranch] {
                guard branch.observations.map(\.step)==steps,branch.observations.allSatisfy({ observation in
                    observation.coordinateNM.isFinite && observation.inProduct==request.productCommitmentRangeNM.contains(observation.coordinateNM) &&
                    observation.inReactant==request.reactantCommitmentRangeNM.contains(observation.coordinateNM) &&
                    !(observation.inProduct && observation.inReactant)
                }),branch.finalCoordinateNM==branch.observations.last?.coordinateNM,
                  branch.outcome==outcome(branch.observations,commitment:request.commitmentObservations) else {
                    throw VivoChemistryError.invalid("dynamical-transmission branch history or commitment evidence")
                }
            }
            guard pair.positiveBranch.velocitiesTimeReversed != pair.negativeBranch.velocitiesTimeReversed else {
                throw VivoChemistryError.invalid("dynamical-transmission velocity pair is not an exact time-reversal pair")
            }
        }
        let weights=pairs.map { pair -> Double in
            let state=states[pair.sourceStateFingerprint]!
            return state.statisticalWeight*pair.positiveCoordinateVelocityNMPerPS
        }
        let total=weights.reduce(0,+),square=weights.reduce(0){$0+$1*$1}
        guard total.isFinite,total>0,square.isFinite,square>0 else { throw VivoChemistryError.convergence("dynamical-transmission surface flux is zero") }
        var history:[VivoQMMMTransmissionPlateauPoint]=[];history.reserveCapacity(steps.count)
        for index in steps.indices {
            var numerator=0.0
            for (pair,weight) in zip(pairs,weights) {
                let plus=pair.positiveBranch.observations[index].inProduct ? 1.0:0.0
                let minus=pair.negativeBranch.observations[index].inProduct ? 1.0:0.0
                numerator+=weight*(plus-minus)
            }
            history.append(.init(step:steps[index],coefficient:numerator/total))
        }
        let plateau=Array(history.suffix(request.acceptance.minimumPlateauObservations)),values=plateau.map(\.coefficient)
        let coefficient=values.reduce(0,+)/Double(values.count),plateauRange=(values.max() ?? 0)-(values.min() ?? 0)
        let effective=total*total/square
        var weightedVariance=0.0,unresolvedWeight=0.0,unresolvedPairs=0
        for (pair,weight) in zip(pairs,weights) {
            let tail=0..<request.acceptance.minimumPlateauObservations
            var score=0.0
            let pObs=Array(pair.positiveBranch.observations.suffix(request.acceptance.minimumPlateauObservations))
            let nObs=Array(pair.negativeBranch.observations.suffix(request.acceptance.minimumPlateauObservations))
            for i in tail { score+=(pObs[i].inProduct ? 1.0:0.0)-(nObs[i].inProduct ? 1.0:0.0) }
            score/=Double(request.acceptance.minimumPlateauObservations)
            weightedVariance+=(weight/total)*pow(score-coefficient,2)
            if pair.positiveBranch.outcome != .productCommitted || pair.negativeBranch.outcome != .reactantRecrossed {
                unresolvedWeight+=weight;unresolvedPairs+=1
            }
        }
        let standardError=effective>1 ? sqrt(max(0,weightedVariance)/effective):Double.greatestFiniteMagnitude
        let unresolved=unresolvedWeight/total
        guard coefficient.isFinite,plateauRange.isFinite,effective.isFinite,standardError.isFinite,unresolved.isFinite else {
            throw VivoChemistryError.convergence("dynamical-transmission reactive-flux estimator is nonfinite")
        }
        var issues:[String]=[]
        if coefficient<=0 || coefficient>1 { issues.append("signed reactive-flux plateau does not yield a physical positive transmission coefficient") }
        if plateauRange>request.acceptance.maximumPlateauRange { issues.append("late-time reactive-flux coefficient has no accepted plateau") }
        if effective<request.acceptance.minimumEffectiveFluxSamples { issues.append("effective flux-weighted surface sample count is below threshold") }
        if unresolved>request.acceptance.maximumUnresolvedFluxFraction { issues.append("unresolved finite-time paired flux fraction exceeds threshold") }
        if standardError>request.acceptance.maximumCoefficientStandardError { issues.append("transmission-coefficient standard error exceeds threshold") }
        let requestID=try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(request))
        struct Evidence:Codable {
            let schema:String;let requestFingerprint:VivoFingerprint;let system:VivoFingerprint;let provider:VivoFingerprint
            let surfaceEvidence:VivoFingerprint;let pairs:[VivoQMMMTransmissionPair];let history:[VivoQMMMTransmissionPlateauPoint]
            let coefficient:Double;let plateauRange:Double;let standardError:Double;let effective:Double;let unresolved:Double
        }
        let evidenceID=try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(Evidence(
            schema:"numivivo.org/qmmm-dynamical-transmission-evidence/v3",requestFingerprint:requestID,
            system:systemFingerprint,provider:providerFingerprint,surfaceEvidence:request.surfaceEnsemble.evidenceFingerprint,
            pairs:pairs,history:history,coefficient:coefficient,plateauRange:plateauRange,
            standardError:standardError,effective:effective,unresolved:unresolved)))
        return .init(schema:VivoQMMMDynamicalTransmissionResult.schema,requestFingerprint:requestID,
            retainedSystemFingerprint:systemFingerprint,providerFingerprint:providerFingerprint,pairs:pairs,
            coefficientHistory:history,transmissionCoefficient:coefficient,plateauRange:plateauRange,
            coefficientStandardError:standardError,effectiveFluxSamples:effective,unresolvedFluxFraction:unresolved,
            positiveProductCommittedCount:pairs.filter{$0.positiveBranch.outcome == .productCommitted}.count,
            negativeReactantCommittedCount:pairs.filter{$0.negativeBranch.outcome == .reactantRecrossed}.count,
            unresolvedPairCount:unresolvedPairs,converged:issues.isEmpty,issues:issues,
            interpretation:interpretation,evidenceFingerprint:evidenceID)
    }

    public static func validate(_ result:VivoQMMMDynamicalTransmissionResult,
                                request:VivoQMMMDynamicalTransmissionRequest)throws {
        let rebuilt=try analyze(request:request,systemFingerprint:result.retainedSystemFingerprint,
                                providerFingerprint:result.providerFingerprint,pairs:result.pairs)
        guard rebuilt==result else { throw VivoChemistryError.invalid("dynamical-transmission result does not reconstruct") }
    }

    public static func applying(_ result:VivoQMMMDynamicalTransmissionResult,
                                transmissionRequest:VivoQMMMDynamicalTransmissionRequest,
                                to rateRequest:VivoQMMMFreeEnergyRateRequest)throws->VivoQMMMFreeEnergyRateRequest {
        try validate(result,request:transmissionRequest)
        guard result.converged,result.transmissionCoefficient>0,result.transmissionCoefficient<=1,
              rateRequest.freeEnergy.evidenceFingerprint==transmissionRequest.freeEnergy.evidenceFingerprint else {
            throw VivoKineticsError.invalid("nonconverged dynamical transmission or PMF evidence mismatch")
        }
        var output=rateRequest
        output.sampledFluxNormalization=try VivoQMMMSurfaceFluxNormalization.make(
            surface:transmissionRequest.surfaceEnsemble,request:transmissionRequest.surfaceRequest)
        output.transmissionProbability=result.transmissionCoefficient; output.transmissionOrigin = .calculated
        output.transmissionEvidence=VivoKineticEvidence(source:"NumiVivo paired signed classical reactive-flux shooting",
            locator:"late-time Bennett-Chandler-style plateau; range=\(result.plateauRange); SE=\(result.coefficientStandardError); effectiveFluxSamples=\(result.effectiveFluxSamples)",
            sourceFingerprint:result.evidenceFingerprint.hex)
        return output
    }
}
