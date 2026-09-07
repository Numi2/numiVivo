import Foundation
@preconcurrency import Metal

public struct VivoQMMMSurfaceEnsembleRequest: Codable, Sendable, Equatable {
    public static let schema="numivivo.org/qmmm-surface-ensemble/v2"
    public var schema:String
    public var freeEnergy:VivoQMMMQualifiedActivationFreeEnergy
    /// Base NVT numerical configuration from the PMF protocol. randomSeed is
    /// rebound to `randomSeed` below for this independent surface trajectory.
    public var dynamics:VivoMDConfiguration
    public var initialState:VivoClassicalInitialState
    public var randomSeed:UInt64
    public var surfaceWindow:VivoQMMMUmbrellaWindow
    public var equilibrationSteps:UInt64
    public var productionSteps:UInt64
    public var candidateEverySteps:UInt64
    public var minimumAcceptedStepSeparation:UInt64
    public var surfaceToleranceNM:Double
    /// Gaussian delta-kernel width used to reweight the finite surface band.
    public var surfaceKernelBandwidthNM:Double
    public var desiredStates:Int
    public var minimumEffectiveSurfaceSamples:Double
    public var maximumAutocorrelationLag:Int
    public init(freeEnergy:VivoQMMMQualifiedActivationFreeEnergy,dynamics:VivoMDConfiguration,
                initialState:VivoClassicalInitialState,randomSeed:UInt64,surfaceWindow:VivoQMMMUmbrellaWindow,
                equilibrationSteps:UInt64,productionSteps:UInt64,candidateEverySteps:UInt64=10,
                minimumAcceptedStepSeparation:UInt64=100,surfaceToleranceNM:Double=0.01,
                surfaceKernelBandwidthNM:Double=0.003,desiredStates:Int=64,
                minimumEffectiveSurfaceSamples:Double=20,maximumAutocorrelationLag:Int=4096) {
        schema=Self.schema;self.freeEnergy=freeEnergy;self.dynamics=dynamics;self.initialState=initialState
        self.randomSeed=randomSeed;self.surfaceWindow=surfaceWindow;self.equilibrationSteps=equilibrationSteps
        self.productionSteps=productionSteps;self.candidateEverySteps=candidateEverySteps
        self.minimumAcceptedStepSeparation=minimumAcceptedStepSeparation;self.surfaceToleranceNM=surfaceToleranceNM
        self.surfaceKernelBandwidthNM=surfaceKernelBandwidthNM;self.desiredStates=desiredStates
        self.minimumEffectiveSurfaceSamples=minimumEffectiveSurfaceSamples;self.maximumAutocorrelationLag=maximumAutocorrelationLag
    }
}

public struct VivoQMMMSurfaceState: Codable, Sendable, Equatable {
    /// Fingerprint of the stored state payload. This is reconstructable without
    /// pretending the biased checkpoint is an unbiased restart artifact.
    public let stateFingerprint:VivoFingerprint
    /// Exact biased MD checkpoint fingerprint retained only as execution provenance.
    public let biasedCheckpointFingerprint:VivoFingerprint
    public let acceptedStep:UInt64
    public let timePS:Double
    public let coordinateNM:Double
    public let coordinateVelocityNMPerPS:Double
    /// Normalized finite-band importance weight; all retained state weights sum to one.
    public let statisticalWeight:Double
    public let positionsNM:[VivoVector3D]
    public let velocitiesNMPerPS:[VivoVector3D]
    public let periodicCell:VivoPeriodicCell?
}

public struct VivoQMMMSurfaceEnsembleResult: Codable, Sendable, Equatable {
    public static let schema="numivivo.org/qmmm-surface-ensemble-result/v2"
    public let schema:String
    public let requestFingerprint:VivoFingerprint
    public let freeEnergyEvidenceFingerprint:VivoFingerprint
    public let retainedSystemFingerprint:VivoFingerprint
    public let baseProviderFingerprint:VivoFingerprint
    public let surfaceProviderFingerprint:VivoFingerprint
    public let surfaceExecutionFingerprint:VivoFingerprint
    public let candidateObservations:Int
    public let nearSurfaceObservations:Int
    public let coordinateAutocorrelationStride:Int
    public let effectiveAcceptedStepSeparation:UInt64
    public let states:[VivoQMMMSurfaceState]
    public let effectiveSurfaceSamples:Double
    /// <dot(xi) theta(dot(xi))> at the dividing surface. Paired velocity symmetry
    /// is used as 0.5*<|dot(xi)|>, which is valid for the canonical velocity law.
    public let positiveCoordinateVelocityNMPerPS:Double
    public let positiveCoordinateVelocityStandardErrorNMPerPS:Double
    public let converged:Bool
    public let issues:[String]
    public let evidenceFingerprint:VivoFingerprint
}

/// Generates a canonical finite-band dividing-surface ensemble under an
/// independent harmonic window. Configurations are reweighted back toward the
/// unbiased dividing surface using exp(+beta U_bias) times a Gaussian delta
/// kernel. Correlation is estimated from the complete production coordinate
/// trace before step-separated states are selected across the trajectory span.
/// Subsequent shooting starts fresh unbiased NVE trajectories from the retained
/// positions/velocities; the biased checkpoint fingerprint remains provenance.
public enum VivoQMMMSurfaceEnsemble {
    private static let gasConstantKJ=0.00831446261815324

    private struct StatePayload:Codable {
        let acceptedStep:UInt64;let timePS:Double;let coordinateNM:Double;let coordinateVelocityNMPerPS:Double
        let positionsNM:[VivoVector3D];let velocitiesNMPerPS:[VivoVector3D];let periodicCell:VivoPeriodicCell?
    }
    private struct Evidence:Codable {
        let schema:String;let requestFingerprint:VivoFingerprint;let freeEnergyEvidenceFingerprint:VivoFingerprint
        let retainedSystemFingerprint:VivoFingerprint;let baseProviderFingerprint:VivoFingerprint
        let surfaceProviderFingerprint:VivoFingerprint;let surfaceExecutionFingerprint:VivoFingerprint
        let candidateObservations:Int;let nearSurfaceObservations:Int;let coordinateAutocorrelationStride:Int
        let effectiveAcceptedStepSeparation:UInt64;let states:[VivoQMMMSurfaceState];let effectiveSurfaceSamples:Double
        let positiveCoordinateVelocityNMPerPS:Double;let positiveCoordinateVelocityStandardErrorNMPerPS:Double
        let converged:Bool;let issues:[String]
    }
    private struct Candidate {
        let checkpoint:VivoMDCheckpoint;let coordinate:Double;let coordinateVelocity:Double
    }

    private static func validateRequest(_ request:VivoQMMMSurfaceEnsembleRequest,
                                        system:VivoClassicalSystem,
                                        baseProvider:VivoMDCandidateForceProvider)throws->(VivoQMMMResolvedCoordinate,VivoMDConfiguration,VivoMDCandidateForceProvider) {
        try request.freeEnergy.validate();try request.dynamics.validate();try request.surfaceWindow.validate()
        let systemID=try system.fingerprint(),dynamicsID=try request.dynamics.fingerprint()
        let surface=request.freeEnergy.analysis.configuration.dividingSurfaceNM
        guard request.schema==VivoQMMMSurfaceEnsembleRequest.schema,
              request.freeEnergy.provenance.systemFingerprint==systemID,
              request.freeEnergy.provenance.structureFingerprint==system.structureFingerprint,
              request.freeEnergy.provenance.baseProviderFingerprint==baseProvider.fingerprint,
              request.freeEnergy.provenance.dynamicsFingerprint==dynamicsID,
              baseProvider.retainedSystemFingerprint==systemID,
              request.initialState.systemFingerprint==systemID,
              request.dynamics.ensemble == .nvt,request.dynamics.thermostat == .langevinMiddle,
              request.dynamics.targetTemperatureK==request.freeEnergy.analysis.temperatureK,
              request.dynamics.frictionPerPS.map({$0>0}) == true,
              request.surfaceWindow.centerNM==surface,
              request.equilibrationSteps>0,request.productionSteps>0,request.candidateEverySteps>0,
              request.productionSteps>=request.candidateEverySteps,
              request.minimumAcceptedStepSeparation>=request.candidateEverySteps,
              request.surfaceToleranceNM.isFinite,request.surfaceToleranceNM>0,
              request.surfaceKernelBandwidthNM.isFinite,request.surfaceKernelBandwidthNM>0,
              request.surfaceKernelBandwidthNM<=request.surfaceToleranceNM,
              (4...100000).contains(request.desiredStates),
              request.minimumEffectiveSurfaceSamples.isFinite,request.minimumEffectiveSurfaceSamples>=2,
              request.maximumAutocorrelationLag>=3 else {
            throw VivoChemistryError.invalid("QM/MM surface-ensemble Hamiltonian, schedule or surface identity")
        }
        try request.initialState.validate(particleCount:system.particles.count)
        let coordinate=try VivoQMMMResolvedCoordinate(source:request.freeEnergy.analysis.coordinate,system:system)
        var dynamics=request.dynamics;dynamics.randomSeed=request.randomSeed
        let provider=try VivoQMMMUmbrellaBias.provider(base:baseProvider,system:system,
            coordinate:request.freeEnergy.analysis.coordinate,window:request.surfaceWindow)
        try provider.validate(system:system,configuration:dynamics,cell:request.initialState.periodicCell)
        return (coordinate,dynamics,provider)
    }

    private static func coordinateVelocity(_ coordinate:VivoQMMMResolvedCoordinate,
                                           checkpoint:VivoMDCheckpoint)throws->Double {
        let geometry=try VivoMDCandidateGeometry(particlePositionsNM:checkpoint.positionsNM,periodicCell:checkpoint.periodicCell)
        let evaluated=try coordinate.evaluate(geometry)
        var value=0.0
        for (particle,gradient) in evaluated.gradients {
            let velocity=checkpoint.velocitiesNMPerPS[Int(particle)]
            value += gradient.x*velocity.x+gradient.y*velocity.y+gradient.z*velocity.z
        }
        guard value.isFinite else { throw VivoChemistryError.convergence("surface coordinate velocity is nonfinite") }
        return value
    }

    private static func autocorrelationStride(_ x:[Double],maximumLag:Int)->Int {
        guard x.count>3 else { return 1 }
        let mean=x.reduce(0,+)/Double(x.count),variance=x.reduce(0){$0+pow($1-mean,2)}/Double(x.count)
        guard variance>0,variance.isFinite else { return x.count }
        var tau=1.0
        for lag in 1...min(maximumLag,x.count-2) {
            var covariance=0.0
            for i in 0..<(x.count-lag) { covariance+=(x[i]-mean)*(x[i+lag]-mean) }
            let rho=covariance/Double(x.count-lag)/variance
            if !rho.isFinite || rho<=0 { break }
            tau+=2*rho
        }
        return max(1,min(x.count,Int(ceil(tau))))
    }

    private static func selectedCandidates(_ eligible:[Candidate],desired:Int,separation:UInt64)->[Candidate] {
        var thinned:[Candidate]=[]
        for candidate in eligible {
            if let last=thinned.last, candidate.checkpoint.acceptedStep-last.checkpoint.acceptedStep<separation { continue }
            thinned.append(candidate)
        }
        guard thinned.count>desired else { return thinned }
        if desired==1 { return [thinned[thinned.count/2]] }
        var output:[Candidate]=[];output.reserveCapacity(desired)
        var used=Set<Int>()
        for i in 0..<desired {
            let raw=Double(i)*Double(thinned.count-1)/Double(desired-1)
            var index=Int(raw.rounded())
            while used.contains(index),index+1<thinned.count { index+=1 }
            while used.contains(index),index>0 { index-=1 }
            used.insert(index);output.append(thinned[index])
        }
        return output.sorted{$0.checkpoint.acceptedStep<$1.checkpoint.acceptedStep}
    }

    private static func normalizedWeights(_ candidates:[Candidate],request:VivoQMMMSurfaceEnsembleRequest)throws->[Double] {
        guard let temperature=request.dynamics.targetTemperatureK,temperature.isFinite,temperature>0 else {
            throw VivoChemistryError.invalid("surface ensemble has no positive temperature")
        }
        let beta=1/(gasConstantKJ*temperature),center=request.surfaceWindow.centerNM,band=request.surfaceKernelBandwidthNM
        let logWeights=candidates.map { candidate -> Double in
            let z=(candidate.coordinate-center)/band
            return beta*request.surfaceWindow.biasKJPerMol(candidate.coordinate)-0.5*z*z
        }
        guard let maximum=logWeights.max(),maximum.isFinite else { throw VivoChemistryError.convergence("surface importance weights are empty or nonfinite") }
        let raw=logWeights.map{exp($0-maximum)},sum=raw.reduce(0,+)
        guard sum.isFinite,sum>0 else { throw VivoChemistryError.convergence("surface importance normalization") }
        let weights=raw.map{$0/sum}
        guard weights.allSatisfy({$0.isFinite && $0>=0}) else { throw VivoChemistryError.convergence("surface normalized weights") }
        return weights
    }

    private static func stateFingerprint(_ candidate:Candidate)throws->VivoFingerprint {
        try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(StatePayload(acceptedStep:candidate.checkpoint.acceptedStep,
            timePS:candidate.checkpoint.timePS,coordinateNM:candidate.coordinate,coordinateVelocityNMPerPS:candidate.coordinateVelocity,
            positionsNM:candidate.checkpoint.positionsNM,velocitiesNMPerPS:candidate.checkpoint.velocitiesNMPerPS,
            periodicCell:candidate.checkpoint.periodicCell)))
    }

    public static func run(_ request:VivoQMMMSurfaceEnsembleRequest,system:VivoClassicalSystem,
                           baseProvider:VivoMDCandidateForceProvider,device:MTLDevice?=nil) async throws -> VivoQMMMSurfaceEnsembleResult {
        let (coordinate,dynamics,surfaceProvider)=try validateRequest(request,system:system,baseProvider:baseProvider)
        let runtime=try await VivoMDMetalRuntime.make(system:system,initialState:request.initialState,
            configuration:dynamics,device:device,forceProvider:surfaceProvider)
        _ = try await runtime.thermalize(temperatureK:dynamics.targetTemperatureK!,
            seed:request.randomSeed^0x5355524641434558)
        for _ in 0..<request.equilibrationSteps {
            try Task.checkCancellation()
            guard try await runtime.step().committed else { throw VivoChemistryError.convergence("surface-ensemble equilibration candidate rejected") }
        }
        var eligible:[Candidate]=[],allCoordinates:[Double]=[],candidateObservations=0,production:UInt64=0
        while production<request.productionSteps {
            try Task.checkCancellation()
            let block=min(request.candidateEverySteps,request.productionSteps-production)
            for _ in 0..<block {
                guard try await runtime.step().committed else { throw VivoChemistryError.convergence("surface-ensemble production candidate rejected") }
                production+=1
            }
            candidateObservations+=1
            let checkpoint=try await runtime.checkpoint()
            let geometry=try VivoMDCandidateGeometry(particlePositionsNM:checkpoint.positionsNM,periodicCell:checkpoint.periodicCell)
            let value=try coordinate.evaluate(geometry).valueNM
            guard value.isFinite else { throw VivoChemistryError.convergence("surface-ensemble reaction coordinate is nonfinite") }
            allCoordinates.append(value)
            if abs(value-request.surfaceWindow.centerNM)<=request.surfaceToleranceNM {
                eligible.append(.init(checkpoint:checkpoint,coordinate:value,
                    coordinateVelocity:try coordinateVelocity(coordinate,checkpoint:checkpoint)))
            }
        }
        let stride=autocorrelationStride(allCoordinates,maximumLag:request.maximumAutocorrelationLag)
        let statisticalSeparation=UInt64(stride)*request.candidateEverySteps
        let separation=max(request.minimumAcceptedStepSeparation,statisticalSeparation)
        let selected=selectedCandidates(eligible,desired:request.desiredStates,separation:separation)
        let weights=try normalizedWeights(selected,request:request)
        var states:[VivoQMMMSurfaceState]=[];states.reserveCapacity(selected.count)
        for (candidate,weight) in zip(selected,weights) {
            states.append(.init(stateFingerprint:try stateFingerprint(candidate),
                biasedCheckpointFingerprint:try candidate.checkpoint.fingerprint(),acceptedStep:candidate.checkpoint.acceptedStep,
                timePS:candidate.checkpoint.timePS,coordinateNM:candidate.coordinate,
                coordinateVelocityNMPerPS:candidate.coordinateVelocity,statisticalWeight:weight,
                positionsNM:candidate.checkpoint.positionsNM,velocitiesNMPerPS:candidate.checkpoint.velocitiesNMPerPS,
                periodicCell:candidate.checkpoint.periodicCell))
        }
        let sumSquare=weights.reduce(0){$0+$1*$1},effective=sumSquare>0 ? 1/sumSquare:0
        let positiveVelocity=0.5*zip(weights,selected).reduce(0.0){$0+$1.0*abs($1.1.coordinateVelocity)}
        let weightedVariance=zip(weights,selected).reduce(0.0) { partial,pair in
            let positive=0.5*abs(pair.1.coordinateVelocity)
            return partial+pair.0*pow(positive-positiveVelocity,2)
        }
        let velocitySE=effective>1 ? sqrt(max(0,weightedVariance)/effective):Double.greatestFiniteMagnitude
        let requestID=try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(request)),systemID=try system.fingerprint()
        let executionID=try VivoMDCandidateForceProvider.executionFingerprint(configuration:dynamics,provider:surfaceProvider)
        var issues:[String]=[]
        if selected.count<request.desiredStates { issues.append("surface trajectory did not yield the requested number of correlation-separated near-surface states") }
        if effective<request.minimumEffectiveSurfaceSamples { issues.append("finite-band reweighted surface effective sample count is below threshold") }
        if positiveVelocity<=0 || !positiveVelocity.isFinite { issues.append("surface ensemble has no finite positive reaction-coordinate flux") }
        if velocitySE == Double.greatestFiniteMagnitude || !velocitySE.isFinite { issues.append("surface positive-velocity uncertainty is unresolved") }
        let converged=issues.isEmpty
        let evidence=Evidence(schema:"numivivo.org/qmmm-surface-ensemble-evidence/v2",requestFingerprint:requestID,
            freeEnergyEvidenceFingerprint:request.freeEnergy.evidenceFingerprint,retainedSystemFingerprint:systemID,
            baseProviderFingerprint:baseProvider.fingerprint,surfaceProviderFingerprint:surfaceProvider.fingerprint,
            surfaceExecutionFingerprint:executionID,candidateObservations:candidateObservations,nearSurfaceObservations:eligible.count,
            coordinateAutocorrelationStride:stride,effectiveAcceptedStepSeparation:separation,states:states,
            effectiveSurfaceSamples:effective,positiveCoordinateVelocityNMPerPS:positiveVelocity,
            positiveCoordinateVelocityStandardErrorNMPerPS:velocitySE,converged:converged,issues:issues)
        let evidenceID=try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(evidence))
        return .init(schema:VivoQMMMSurfaceEnsembleResult.schema,requestFingerprint:requestID,
            freeEnergyEvidenceFingerprint:request.freeEnergy.evidenceFingerprint,retainedSystemFingerprint:systemID,
            baseProviderFingerprint:baseProvider.fingerprint,surfaceProviderFingerprint:surfaceProvider.fingerprint,
            surfaceExecutionFingerprint:executionID,candidateObservations:candidateObservations,nearSurfaceObservations:eligible.count,
            coordinateAutocorrelationStride:stride,effectiveAcceptedStepSeparation:separation,states:states,
            effectiveSurfaceSamples:effective,positiveCoordinateVelocityNMPerPS:positiveVelocity,
            positiveCoordinateVelocityStandardErrorNMPerPS:velocitySE,converged:converged,issues:issues,evidenceFingerprint:evidenceID)
    }

    public static func validate(_ result:VivoQMMMSurfaceEnsembleResult,
                                request:VivoQMMMSurfaceEnsembleRequest)throws {
        try request.freeEnergy.validate();try request.dynamics.validate();try request.surfaceWindow.validate()
        let requestID=try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(request))
        guard result.schema==VivoQMMMSurfaceEnsembleResult.schema,result.requestFingerprint==requestID,
              result.freeEnergyEvidenceFingerprint==request.freeEnergy.evidenceFingerprint,
              result.retainedSystemFingerprint==request.freeEnergy.provenance.systemFingerprint,
              result.baseProviderFingerprint==request.freeEnergy.provenance.baseProviderFingerprint,
              result.candidateObservations>=result.nearSurfaceObservations,
              result.nearSurfaceObservations>=result.states.count,result.states.count<=request.desiredStates,
              result.coordinateAutocorrelationStride>=1,
              result.effectiveAcceptedStepSeparation>=request.minimumAcceptedStepSeparation,
              result.states.map(\.stateFingerprint).count==Set(result.states.map(\.stateFingerprint)).count,
              result.states.map(\.biasedCheckpointFingerprint).count==Set(result.states.map(\.biasedCheckpointFingerprint)).count,
              result.effectiveSurfaceSamples.isFinite,result.effectiveSurfaceSamples>=0,
              result.positiveCoordinateVelocityNMPerPS.isFinite,result.positiveCoordinateVelocityNMPerPS>=0,
              result.positiveCoordinateVelocityStandardErrorNMPerPS.isFinite,
              result.states.allSatisfy({ state in
                  state.timePS.isFinite && state.timePS>=0 && state.coordinateNM.isFinite &&
                  abs(state.coordinateNM-request.surfaceWindow.centerNM)<=request.surfaceToleranceNM &&
                  state.coordinateVelocityNMPerPS.isFinite && state.statisticalWeight.isFinite && state.statisticalWeight>=0 &&
                  state.positionsNM.count==state.velocitiesNMPerPS.count && !state.positionsNM.isEmpty &&
                  state.positionsNM.allSatisfy(\.isFinite) && state.velocitiesNMPerPS.allSatisfy(\.isFinite)
              }) else {
            throw VivoChemistryError.invalid("QM/MM surface-ensemble result does not match its immutable request")
        }
        if !result.states.isEmpty {
            let sum=result.states.reduce(0.0){$0+$1.statisticalWeight}
            guard abs(sum-1)<=1e-10 else { throw VivoChemistryError.invalid("surface state statistical weights do not normalize") }
        }
        for i in result.states.indices {
            let state=result.states[i]
            let payload=StatePayload(acceptedStep:state.acceptedStep,timePS:state.timePS,coordinateNM:state.coordinateNM,
                coordinateVelocityNMPerPS:state.coordinateVelocityNMPerPS,positionsNM:state.positionsNM,
                velocitiesNMPerPS:state.velocitiesNMPerPS,periodicCell:state.periodicCell)
            guard state.stateFingerprint==(try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(payload))) else {
                throw VivoChemistryError.invalid("surface state payload fingerprint mismatch")
            }
            if i>0 {
                guard state.acceptedStep>result.states[i-1].acceptedStep,
                      state.acceptedStep-result.states[i-1].acceptedStep>=result.effectiveAcceptedStepSeparation else {
                    throw VivoChemistryError.invalid("QM/MM surface states violate correlation-aware accepted-step separation")
                }
            }
        }
        let sumSquare=result.states.reduce(0.0){$0+$1.statisticalWeight*$1.statisticalWeight}
        let effective=sumSquare>0 ? 1/sumSquare:0
        let positive=0.5*result.states.reduce(0.0){$0+$1.statisticalWeight*abs($1.coordinateVelocityNMPerPS)}
        let variance=result.states.reduce(0.0){$0+$1.statisticalWeight*pow(0.5*abs($1.coordinateVelocityNMPerPS)-positive,2)}
        let se=effective>1 ? sqrt(max(0,variance)/effective):Double.greatestFiniteMagnitude
        guard abs(effective-result.effectiveSurfaceSamples)<=max(1e-10,abs(effective)*1e-10),
              abs(positive-result.positiveCoordinateVelocityNMPerPS)<=max(1e-12,abs(positive)*1e-10),
              se==result.positiveCoordinateVelocityStandardErrorNMPerPS || abs(se-result.positiveCoordinateVelocityStandardErrorNMPerPS)<=max(1e-12,abs(se)*1e-10) else {
            throw VivoChemistryError.invalid("surface ensemble effective count or velocity normalization does not reconstruct")
        }
        let expectedConverged=result.states.count==request.desiredStates && effective>=request.minimumEffectiveSurfaceSamples && positive>0 && se.isFinite
        guard result.converged==(result.issues.isEmpty && expectedConverged) else {
            throw VivoChemistryError.invalid("surface ensemble convergence flag differs from acceptance contract")
        }
        let evidence=Evidence(schema:"numivivo.org/qmmm-surface-ensemble-evidence/v2",requestFingerprint:requestID,
            freeEnergyEvidenceFingerprint:result.freeEnergyEvidenceFingerprint,retainedSystemFingerprint:result.retainedSystemFingerprint,
            baseProviderFingerprint:result.baseProviderFingerprint,surfaceProviderFingerprint:result.surfaceProviderFingerprint,
            surfaceExecutionFingerprint:result.surfaceExecutionFingerprint,candidateObservations:result.candidateObservations,
            nearSurfaceObservations:result.nearSurfaceObservations,coordinateAutocorrelationStride:result.coordinateAutocorrelationStride,
            effectiveAcceptedStepSeparation:result.effectiveAcceptedStepSeparation,states:result.states,
            effectiveSurfaceSamples:result.effectiveSurfaceSamples,positiveCoordinateVelocityNMPerPS:result.positiveCoordinateVelocityNMPerPS,
            positiveCoordinateVelocityStandardErrorNMPerPS:result.positiveCoordinateVelocityStandardErrorNMPerPS,
            converged:result.converged,issues:result.issues)
        guard result.evidenceFingerprint==(try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(evidence))) else {
            throw VivoChemistryError.invalid("QM/MM surface-ensemble evidence fingerprint mismatch")
        }
    }
}
