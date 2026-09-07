import Foundation
@preconcurrency import Metal

public struct VivoQMMMSurfaceEnsembleRequest: Codable, Sendable, Equatable {
    public static let schema="numivivo.org/qmmm-surface-ensemble/v1"
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
    public var desiredStates:Int
    public init(freeEnergy:VivoQMMMQualifiedActivationFreeEnergy,dynamics:VivoMDConfiguration,
                initialState:VivoClassicalInitialState,randomSeed:UInt64,surfaceWindow:VivoQMMMUmbrellaWindow,
                equilibrationSteps:UInt64,productionSteps:UInt64,candidateEverySteps:UInt64=10,
                minimumAcceptedStepSeparation:UInt64=100,surfaceToleranceNM:Double=0.01,desiredStates:Int=64) {
        schema=Self.schema;self.freeEnergy=freeEnergy;self.dynamics=dynamics;self.initialState=initialState
        self.randomSeed=randomSeed;self.surfaceWindow=surfaceWindow;self.equilibrationSteps=equilibrationSteps
        self.productionSteps=productionSteps;self.candidateEverySteps=candidateEverySteps
        self.minimumAcceptedStepSeparation=minimumAcceptedStepSeparation;self.surfaceToleranceNM=surfaceToleranceNM
        self.desiredStates=desiredStates
    }
}

public struct VivoQMMMSurfaceState: Codable, Sendable, Equatable {
    public let checkpointFingerprint:VivoFingerprint
    public let acceptedStep:UInt64
    public let timePS:Double
    public let coordinateNM:Double
    public let positionsNM:[VivoVector3D]
    public let velocitiesNMPerPS:[VivoVector3D]
    public let periodicCell:VivoPeriodicCell?
}

public struct VivoQMMMSurfaceEnsembleResult: Codable, Sendable, Equatable {
    public static let schema="numivivo.org/qmmm-surface-ensemble-result/v1"
    public let schema:String
    public let requestFingerprint:VivoFingerprint
    public let freeEnergyEvidenceFingerprint:VivoFingerprint
    public let retainedSystemFingerprint:VivoFingerprint
    public let baseProviderFingerprint:VivoFingerprint
    public let surfaceProviderFingerprint:VivoFingerprint
    public let surfaceExecutionFingerprint:VivoFingerprint
    public let candidateObservations:Int
    public let states:[VivoQMMMSurfaceState]
    public let converged:Bool
    public let issues:[String]
    public let evidenceFingerprint:VivoFingerprint
}

/// Generates a bounded canonical near-surface ensemble under an independent
/// harmonic window centered exactly at the qualified PMF dividing surface. The
/// output retains positions and canonical velocities; subsequent shooting starts
/// a new unbiased NVE runtime instead of weakening checkpoint restart identity.
public enum VivoQMMMSurfaceEnsemble {
    private struct Evidence:Codable {
        let schema:String
        let requestFingerprint:VivoFingerprint
        let freeEnergyEvidenceFingerprint:VivoFingerprint
        let retainedSystemFingerprint:VivoFingerprint
        let baseProviderFingerprint:VivoFingerprint
        let surfaceProviderFingerprint:VivoFingerprint
        let surfaceExecutionFingerprint:VivoFingerprint
        let candidateObservations:Int
        let states:[VivoQMMMSurfaceState]
        let converged:Bool
        let issues:[String]
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
              (4...100000).contains(request.desiredStates) else {
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
        var states:[VivoQMMMSurfaceState]=[],candidateObservations=0,lastAccepted:UInt64?
        states.reserveCapacity(request.desiredStates)
        var production:UInt64=0
        while production<request.productionSteps && states.count<request.desiredStates {
            try Task.checkCancellation()
            let block=min(request.candidateEverySteps,request.productionSteps-production)
            for _ in 0..<block {
                guard try await runtime.step().committed else { throw VivoChemistryError.convergence("surface-ensemble production candidate rejected") }
                production+=1
            }
            candidateObservations+=1
            let checkpoint=try await runtime.checkpoint()
            if let lastAccepted,checkpoint.acceptedStep-lastAccepted<request.minimumAcceptedStepSeparation { continue }
            let geometry=try VivoMDCandidateGeometry(particlePositionsNM:checkpoint.positionsNM,periodicCell:checkpoint.periodicCell)
            let value=try coordinate.evaluate(geometry).valueNM
            guard value.isFinite else { throw VivoChemistryError.convergence("surface-ensemble reaction coordinate is nonfinite") }
            if abs(value-request.surfaceWindow.centerNM)<=request.surfaceToleranceNM {
                let id=try checkpoint.fingerprint()
                states.append(.init(checkpointFingerprint:id,acceptedStep:checkpoint.acceptedStep,timePS:checkpoint.timePS,
                    coordinateNM:value,positionsNM:checkpoint.positionsNM,velocitiesNMPerPS:checkpoint.velocitiesNMPerPS,
                    periodicCell:checkpoint.periodicCell))
                lastAccepted=checkpoint.acceptedStep
            }
        }
        let requestID=try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(request))
        let systemID=try system.fingerprint()
        let executionID=try VivoMDCandidateForceProvider.executionFingerprint(configuration:dynamics,provider:surfaceProvider)
        var issues:[String]=[]
        if states.count<request.desiredStates { issues.append("surface trajectory did not yield the requested number of step-separated near-surface states") }
        let evidence=Evidence(schema:"numivivo.org/qmmm-surface-ensemble-evidence/v1",requestFingerprint:requestID,
            freeEnergyEvidenceFingerprint:request.freeEnergy.evidenceFingerprint,retainedSystemFingerprint:systemID,
            baseProviderFingerprint:baseProvider.fingerprint,surfaceProviderFingerprint:surfaceProvider.fingerprint,
            surfaceExecutionFingerprint:executionID,candidateObservations:candidateObservations,states:states,
            converged:issues.isEmpty,issues:issues)
        let evidenceID=try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(evidence))
        return .init(schema:VivoQMMMSurfaceEnsembleResult.schema,requestFingerprint:requestID,
            freeEnergyEvidenceFingerprint:request.freeEnergy.evidenceFingerprint,retainedSystemFingerprint:systemID,
            baseProviderFingerprint:baseProvider.fingerprint,surfaceProviderFingerprint:surfaceProvider.fingerprint,
            surfaceExecutionFingerprint:executionID,candidateObservations:candidateObservations,states:states,
            converged:issues.isEmpty,issues:issues,evidenceFingerprint:evidenceID)
    }

    public static func validate(_ result:VivoQMMMSurfaceEnsembleResult,
                                request:VivoQMMMSurfaceEnsembleRequest)throws {
        try request.freeEnergy.validate();try request.dynamics.validate();try request.surfaceWindow.validate()
        let requestID=try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(request))
        guard result.schema==VivoQMMMSurfaceEnsembleResult.schema,result.requestFingerprint==requestID,
              result.freeEnergyEvidenceFingerprint==request.freeEnergy.evidenceFingerprint,
              result.retainedSystemFingerprint==request.freeEnergy.provenance.systemFingerprint,
              result.baseProviderFingerprint==request.freeEnergy.provenance.baseProviderFingerprint,
              result.candidateObservations>=result.states.count,
              result.states.count<=request.desiredStates,
              result.states.map(\.checkpointFingerprint).count==Set(result.states.map(\.checkpointFingerprint)).count,
              result.states.allSatisfy({ state in
                  state.timePS.isFinite && state.timePS>=0 && state.coordinateNM.isFinite &&
                  abs(state.coordinateNM-request.surfaceWindow.centerNM)<=request.surfaceToleranceNM &&
                  state.positionsNM.count==state.velocitiesNMPerPS.count && !state.positionsNM.isEmpty &&
                  state.positionsNM.allSatisfy(\.isFinite) && state.velocitiesNMPerPS.allSatisfy(\.isFinite)
              }),result.converged==(result.issues.isEmpty && result.states.count==request.desiredStates) else {
            throw VivoChemistryError.invalid("QM/MM surface-ensemble result does not match its immutable request")
        }
        for i in 1..<result.states.count {
            guard result.states[i].acceptedStep>result.states[i-1].acceptedStep,
                  result.states[i].acceptedStep-result.states[i-1].acceptedStep>=request.minimumAcceptedStepSeparation else {
                throw VivoChemistryError.invalid("QM/MM surface states violate declared accepted-step separation")
            }
        }
        let evidence=Evidence(schema:"numivivo.org/qmmm-surface-ensemble-evidence/v1",requestFingerprint:requestID,
            freeEnergyEvidenceFingerprint:result.freeEnergyEvidenceFingerprint,retainedSystemFingerprint:result.retainedSystemFingerprint,
            baseProviderFingerprint:result.baseProviderFingerprint,surfaceProviderFingerprint:result.surfaceProviderFingerprint,
            surfaceExecutionFingerprint:result.surfaceExecutionFingerprint,candidateObservations:result.candidateObservations,
            states:result.states,converged:result.converged,issues:result.issues)
        guard result.evidenceFingerprint==(try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(evidence))) else {
            throw VivoChemistryError.invalid("QM/MM surface-ensemble evidence fingerprint mismatch")
        }
    }
}
