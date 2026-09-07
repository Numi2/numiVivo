import Foundation
@preconcurrency import Metal
import Testing
@testable import NumiVivoKit

@Suite(.serialized) struct QMMMTransmissionTests {
    private func fingerprint(_ text:String) throws -> VivoFingerprint {
        try VivoCanonicalJSON.fingerprint(Data(text.utf8))
    }

    private func qualifiedPMF(source:VivoMDConfiguration,systemID:VivoFingerprint,providerID:VivoFingerprint) throws -> VivoQMMMQualifiedActivationFreeEnergy {
        let coordinate=VivoQMMMReactionCoordinate(identifier:"synthetic-transfer",kind:.distanceDifference,atomIndices:[0,1,2,3])
        let windows=[VivoQMMMUmbrellaWindow(identifier:"left",centerNM:-0.15,forceConstantKJPerMolNM2:10),
                     VivoQMMMUmbrellaWindow(identifier:"right",centerNM:0.15,forceConstantKJPerMolNM2:10)]
        var values:[Double]=[]
        for _ in 0..<20 { for i in 0...60 {
            let x=Double(i)*0.01;values.append(-0.30+x);values.append(0.30-x)
        } }
        let traces=windows.enumerated().map { index,window in
            VivoQMMMUmbrellaTrace(window:window,randomSeed:UInt64(index+1),coordinateNM:values,
                                  potentialEnergyKJPerMol:[Double](repeating:0,count:values.count))
        }
        let config=VivoQMMMFreeEnergyAnalysisConfiguration(bins:81,kernelBandwidthNM:0.025,
            reactantRangeNM:-0.28 ... -0.12,dividingSurfaceNM:0.29,
            minimumDecorrelatedSamplesPerWindow:20,minimumAdjacentOverlap:0)
        let analysis=try VivoQMMMFreeEnergy.analyze(coordinate:coordinate,temperatureK:300,traces:traces,configuration:config)
        #expect(analysis.converged)
        let protocolID=try fingerprint("transmission-test-protocol")
        let provenance=VivoQMMMFreeEnergyProvenance(
            structureFingerprint:try fingerprint("structure"),systemFingerprint:systemID,
            baseProviderFingerprint:providerID,dynamicsFingerprint:try source.fingerprint(),
            samplingExecution:.init(requestFingerprint:try fingerprint("execution"),replicaProtocolFingerprint:protocolID),
            chemicalState:"state",environment:.proteinEnvironment,environmentIdentifier:"host",
            methodDescription:"synthetic transmission estimator fixture",
            reactionConnectivity:.init(connectivityFingerprint:try fingerprint("connectivity"),
                reactantEndpointIdentifier:"R",productEndpointIdentifier:"P",mappedReactionAtomIndices:coordinate.atomIndices))
        return try VivoQMMMQualifiedActivationFreeEnergy(analysis:analysis,provenance:provenance,
            fluxNormalization:.init(surfaceToReactantDensityPerNM:try VivoQMMMFreeEnergyQualification.surfaceToReactantDensityPerNM(analysis),
                                    inverseMassMetricPerDa:1))
    }

    private func syntheticSurface(freeEnergy:VivoQMMMQualifiedActivationFreeEnergy,source:VivoMDConfiguration,
                                  systemID:VivoFingerprint,providerID:VivoFingerprint)
    throws -> (VivoQMMMSurfaceEnsembleRequest,VivoQMMMSurfaceEnsembleResult,[VivoFingerprint]) {
        let initial=VivoClassicalInitialState(systemFingerprint:systemID,positionsNM:[.zero],periodicCell:nil)
        let window=VivoQMMMUmbrellaWindow(identifier:"surface",centerNM:0.29,forceConstantKJPerMolNM2:20)
        let request=VivoQMMMSurfaceEnsembleRequest(freeEnergy:freeEnergy,dynamics:source,initialState:initial,
            randomSeed:991,surfaceWindow:window,equilibrationSteps:1,productionSteps:4,candidateEverySteps:1,
            minimumAcceptedStepSeparation:1,surfaceToleranceNM:0.01,desiredStates:4)
        let ids=try (0..<4).map { try fingerprint("surface-state-\($0)") }
        let states=(0..<4).map { i in
            VivoQMMMSurfaceState(checkpointFingerprint:ids[i],acceptedStep:UInt64(i+1),timePS:Double(i+1)*0.001,
                coordinateNM:0.29,positionsNM:[.zero],velocitiesNMPerPS:[.zero],periodicCell:nil)
        }
        let requestID=try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(request))
        let surfaceProviderID=try fingerprint("surface-provider"),executionID=try fingerprint("surface-execution")
        struct Evidence:Codable {
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
        let evidence=Evidence(schema:"numivivo.org/qmmm-surface-ensemble-evidence/v1",requestFingerprint:requestID,
            freeEnergyEvidenceFingerprint:freeEnergy.evidenceFingerprint,retainedSystemFingerprint:systemID,
            baseProviderFingerprint:providerID,surfaceProviderFingerprint:surfaceProviderID,
            surfaceExecutionFingerprint:executionID,candidateObservations:4,states:states,converged:true,issues:[])
        let evidenceID=try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(evidence))
        let result=VivoQMMMSurfaceEnsembleResult(schema:VivoQMMMSurfaceEnsembleResult.schema,requestFingerprint:requestID,
            freeEnergyEvidenceFingerprint:freeEnergy.evidenceFingerprint,retainedSystemFingerprint:systemID,
            baseProviderFingerprint:providerID,surfaceProviderFingerprint:surfaceProviderID,
            surfaceExecutionFingerprint:executionID,candidateObservations:4,states:states,converged:true,issues:[],
            evidenceFingerprint:evidenceID)
        try VivoQMMMSurfaceEnsemble.validate(result,request:request)
        return (request,result,ids)
    }

    @Test func fluxWeightedEstimatorReconstructsAndProducesCalculatedTransmissionEvidence() throws {
        let source=VivoMDConfiguration(timeStepPS:0.001,electrostatics:.cutoff,ensemble:.nvt,thermostat:.langevinMiddle,
            targetTemperatureK:300,frictionPerPS:1,neighborListEnabled:false)
        let shooting=VivoMDConfiguration(timeStepPS:0.001,electrostatics:.cutoff,ensemble:.nve,thermostat:.none,
            targetTemperatureK:nil,frictionPerPS:nil,neighborListEnabled:false)
        let systemID=try fingerprint("dummy-system"),providerID=try fingerprint("dummy-provider")
        let freeEnergy=try qualifiedPMF(source:source,systemID:systemID,providerID:providerID)
        let (surfaceRequest,surfaceEnsemble,ids)=try syntheticSurface(freeEnergy:freeEnergy,source:source,
            systemID:systemID,providerID:providerID)
        let request=VivoQMMMDynamicalTransmissionRequest(surfaceRequest:surfaceRequest,surfaceEnsemble:surfaceEnsemble,
            shootingDynamics:shooting,reactantCommitmentRangeNM:-0.5 ... -0.35,productCommitmentRangeNM:0.35 ... 0.5,
            maximumSteps:100,observeEverySteps:1,commitmentObservations:1,
            acceptance:.init(minimumSurfaceCheckpoints:4,minimumEffectiveFluxSamples:3.5,
                             maximumUnresolvedFluxFraction:0.1,maximumCoefficientStandardError:0.3))
        let trajectories=[
            VivoQMMMTransmissionTrajectory(sourceCheckpointFingerprint:ids[0],initialCoordinateNM:0.29,
                productDirectedCoordinateVelocityNMPerPS:1,velocitiesTimeReversed:false,outcome:.productCommitted,committedSteps:10,finalCoordinateNM:0.4),
            VivoQMMMTransmissionTrajectory(sourceCheckpointFingerprint:ids[1],initialCoordinateNM:0.29,
                productDirectedCoordinateVelocityNMPerPS:1,velocitiesTimeReversed:true,outcome:.productCommitted,committedSteps:12,finalCoordinateNM:0.4),
            VivoQMMMTransmissionTrajectory(sourceCheckpointFingerprint:ids[2],initialCoordinateNM:0.29,
                productDirectedCoordinateVelocityNMPerPS:1,velocitiesTimeReversed:false,outcome:.productCommitted,committedSteps:9,finalCoordinateNM:0.4),
            VivoQMMMTransmissionTrajectory(sourceCheckpointFingerprint:ids[3],initialCoordinateNM:0.29,
                productDirectedCoordinateVelocityNMPerPS:1,velocitiesTimeReversed:true,outcome:.reactantRecrossed,committedSteps:7,finalCoordinateNM:-0.4)]
        let result=try VivoQMMMDynamicalTransmission.analyze(request:request,systemFingerprint:systemID,
            providerFingerprint:providerID,trajectories:trajectories)
        #expect(result.converged)
        #expect(abs(result.transmissionCoefficient-0.75)<1e-12)
        #expect(result.effectiveFluxSamples==4)
        #expect(result.unresolvedFluxFraction==0)
        #expect(result.coefficientStandardError<=0.3)
        try VivoQMMMDynamicalTransmission.validate(result,request:request)

        let context=VivoKineticContext(compound:"synthetic",target:"protein",targetVariant:"reference",site:"site",
            chemicalState:"state",hostContext:"host",temperatureK:300,pH:7,ionicStrengthM:0.15)
        let rate=VivoQMMMFreeEnergyRateRequest(context:context,environment:.proteinEnvironment,freeEnergy:freeEnergy,
            transmissionProbability:1,transmissionOrigin:.assumed,
            transmissionEvidence:.init(source:"assumption",locator:"test"),samplingDescription:"fixture")
        let computed=try VivoQMMMDynamicalTransmission.applying(result,transmissionRequest:request,to:rate)
        #expect(computed.transmissionProbability==0.75)
        #expect(computed.transmissionOrigin == .calculated)
        #expect(computed.transmissionEvidence.sourceFingerprint==result.evidenceFingerprint.hex)

        var unresolved=trajectories
        unresolved[0] = .init(sourceCheckpointFingerprint:ids[0],initialCoordinateNM:0.29,
            productDirectedCoordinateVelocityNMPerPS:1,velocitiesTimeReversed:false,outcome:.unresolved,committedSteps:100,finalCoordinateNM:0.1)
        let failed=try VivoQMMMDynamicalTransmission.analyze(request:request,systemFingerprint:systemID,
            providerFingerprint:providerID,trajectories:unresolved)
        #expect(!failed.converged)
        #expect(failed.unresolvedFluxFraction==0.25)
        #expect(throws:(any Error).self) { _ = try VivoQMMMDynamicalTransmission.applying(failed,transmissionRequest:request,to:rate) }
        #expect(throws:(any Error).self) {
            _ = try VivoQMMMDynamicalTransmission.analyze(request:request,
                systemFingerprint:try fingerprint("wrong-system"),providerFingerprint:providerID,trajectories:trajectories)
        }
    }

    @Test func metalSurfaceCollectionHandsExactStatesToFreshUnbiasedNVEShooting() async throws {
        let structureID=try fingerprint("structure")
        let particles=(0..<4).map { i in VivoClassicalParticle(index:UInt32(i),atomIndex:UInt32(i),
            typeIdentifier:"H",massDa:1,chargeE:0,sigmaNM:0,epsilonKJPerMol:0) }
        let system=VivoClassicalSystem(identifier:"transmission-metal-fixture",structureFingerprint:structureID,particles:particles)
        let systemID=try system.fingerprint(),providerID=try fingerprint("transmission-zero-provider")
        let provider=try VivoMDCandidateForceProvider(fingerprint:providerID,retainedSystemFingerprint:systemID,
            boundary:.finiteCluster,supportsCellMoves:false,maximumAcceptedResidual:1e-8,molecularConnectivitySystem:system) { geometry in
            try VivoMDCandidateForceEvaluation(providerFingerprint:providerID,geometry:geometry,additionalEnergyKJPerMol:0,
                physicalParticleForcesKJPerMolNM:Array(repeating:.zero,count:4),derivativeMethod:"zero synthetic BO correction",
                convergenceResidual:0,requiredResidual:1e-8)
        }
        let source=VivoMDConfiguration(timeStepPS:0.0001,electrostatics:.cutoff,ensemble:.nvt,thermostat:.langevinMiddle,
            targetTemperatureK:300,frictionPerPS:1,neighborListEnabled:false)
        let freeEnergy=try qualifiedPMF(source:source,systemID:systemID,providerID:providerID)
        let positions:[VivoVector3D]=[.init(0,0,0),.init(0.49,0,0),.init(0,1,0),.init(0.20,1,0)]
        let initial=VivoClassicalInitialState(systemFingerprint:systemID,positionsNM:positions,periodicCell:nil)
        let surfaceRequest=VivoQMMMSurfaceEnsembleRequest(freeEnergy:freeEnergy,dynamics:source,initialState:initial,
            randomSeed:223,surfaceWindow:.init(identifier:"surface",centerNM:0.29,forceConstantKJPerMolNM2:10_000),
            equilibrationSteps:2,productionSteps:16,candidateEverySteps:1,minimumAcceptedStepSeparation:2,
            surfaceToleranceNM:0.25,desiredStates:4)
        let device=try VivoMetalDeviceSelector.productionDevice()
        #expect(device.hasUnifiedMemory)
        let surface=try await VivoQMMMSurfaceEnsemble.run(surfaceRequest,system:system,baseProvider:provider,device:device)
        #expect(surface.converged)
        #expect(surface.states.count==4)
        try VivoQMMMSurfaceEnsemble.validate(surface,request:surfaceRequest)

        let shooting=VivoMDConfiguration(timeStepPS:0.0001,electrostatics:.cutoff,ensemble:.nve,thermostat:.none,
            targetTemperatureK:nil,frictionPerPS:nil,neighborListEnabled:false)
        let transmission=VivoQMMMDynamicalTransmissionRequest(surfaceRequest:surfaceRequest,surfaceEnsemble:surface,
            shootingDynamics:shooting,reactantCommitmentRangeNM:-2 ... 0,productCommitmentRangeNM:0.6 ... 2,
            maximumSteps:2,observeEverySteps:1,commitmentObservations:1,
            acceptance:.init(minimumSurfaceCheckpoints:4,minimumEffectiveFluxSamples:2,
                             maximumUnresolvedFluxFraction:0.99,maximumCoefficientStandardError:0.49))
        let result=try await VivoQMMMDynamicalTransmission.run(transmission,system:system,baseProvider:provider,device:device)
        #expect(result.trajectories.count==surface.states.count)
        #expect(Set(result.trajectories.map(\.sourceCheckpointFingerprint))==Set(surface.states.map(\.checkpointFingerprint)))
        #expect(result.trajectories.allSatisfy{$0.committedSteps<=2})
        #expect(result.trajectories.contains{$0.committedSteps>0})
        try VivoQMMMDynamicalTransmission.validate(result,request:transmission)
    }
}
