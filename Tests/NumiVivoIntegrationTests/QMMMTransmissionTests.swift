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

    private struct SurfaceStatePayload:Codable {
        let acceptedStep:UInt64
        let timePS:Double
        let coordinateNM:Double
        let coordinateVelocityNMPerPS:Double
        let positionsNM:[VivoVector3D]
        let velocitiesNMPerPS:[VivoVector3D]
        let periodicCell:VivoPeriodicCell?
    }

    private func syntheticSurface(freeEnergy:VivoQMMMQualifiedActivationFreeEnergy,source:VivoMDConfiguration,
                                  systemID:VivoFingerprint,providerID:VivoFingerprint)
    throws -> (VivoQMMMSurfaceEnsembleRequest,VivoQMMMSurfaceEnsembleResult,[VivoFingerprint]) {
        let initial=VivoClassicalInitialState(systemFingerprint:systemID,positionsNM:[.zero],periodicCell:nil)
        let window=VivoQMMMUmbrellaWindow(identifier:"surface",centerNM:0.29,forceConstantKJPerMolNM2:20)
        let request=VivoQMMMSurfaceEnsembleRequest(freeEnergy:freeEnergy,dynamics:source,initialState:initial,
            randomSeed:991,surfaceWindow:window,equilibrationSteps:1,productionSteps:4,candidateEverySteps:1,
            minimumAcceptedStepSeparation:1,surfaceToleranceNM:0.01,surfaceKernelBandwidthNM:0.003,
            desiredStates:4,minimumEffectiveSurfaceSamples:2,maximumAutocorrelationLag:3)
        var states:[VivoQMMMSurfaceState]=[],ids:[VivoFingerprint]=[]
        for i in 0..<4 {
            let step=UInt64(i+1),time=Double(i+1)*0.001,positions=[VivoVector3D.zero],velocities=[VivoVector3D(1,0,0)]
            let payload=SurfaceStatePayload(acceptedStep:step,timePS:time,coordinateNM:0.29,coordinateVelocityNMPerPS:1,
                positionsNM:positions,velocitiesNMPerPS:velocities,periodicCell:nil)
            let stateID=try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(payload));ids.append(stateID)
            states.append(.init(stateFingerprint:stateID,biasedCheckpointFingerprint:try fingerprint("biased-checkpoint-\(i)"),
                acceptedStep:step,timePS:time,coordinateNM:0.29,coordinateVelocityNMPerPS:1,statisticalWeight:0.25,
                positionsNM:positions,velocitiesNMPerPS:velocities,periodicCell:nil))
        }
        let requestID=try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(request))
        let surfaceProviderID=try fingerprint("surface-provider"),executionID=try fingerprint("surface-execution")
        struct Evidence:Codable {
            let schema:String;let requestFingerprint:VivoFingerprint;let freeEnergyEvidenceFingerprint:VivoFingerprint
            let retainedSystemFingerprint:VivoFingerprint;let baseProviderFingerprint:VivoFingerprint
            let surfaceProviderFingerprint:VivoFingerprint;let surfaceExecutionFingerprint:VivoFingerprint
            let candidateObservations:Int;let nearSurfaceObservations:Int;let coordinateAutocorrelationStride:Int
            let effectiveAcceptedStepSeparation:UInt64;let states:[VivoQMMMSurfaceState];let effectiveSurfaceSamples:Double
            let positiveCoordinateVelocityNMPerPS:Double;let positiveCoordinateVelocityStandardErrorNMPerPS:Double
            let converged:Bool;let issues:[String]
        }
        let evidence=Evidence(schema:"numivivo.org/qmmm-surface-ensemble-evidence/v2",requestFingerprint:requestID,
            freeEnergyEvidenceFingerprint:freeEnergy.evidenceFingerprint,retainedSystemFingerprint:systemID,
            baseProviderFingerprint:providerID,surfaceProviderFingerprint:surfaceProviderID,surfaceExecutionFingerprint:executionID,
            candidateObservations:4,nearSurfaceObservations:4,coordinateAutocorrelationStride:1,
            effectiveAcceptedStepSeparation:1,states:states,effectiveSurfaceSamples:4,
            positiveCoordinateVelocityNMPerPS:0.5,positiveCoordinateVelocityStandardErrorNMPerPS:0,
            converged:true,issues:[])
        let evidenceID=try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(evidence))
        let result=VivoQMMMSurfaceEnsembleResult(schema:VivoQMMMSurfaceEnsembleResult.schema,requestFingerprint:requestID,
            freeEnergyEvidenceFingerprint:freeEnergy.evidenceFingerprint,retainedSystemFingerprint:systemID,
            baseProviderFingerprint:providerID,surfaceProviderFingerprint:surfaceProviderID,surfaceExecutionFingerprint:executionID,
            candidateObservations:4,nearSurfaceObservations:4,coordinateAutocorrelationStride:1,effectiveAcceptedStepSeparation:1,
            states:states,effectiveSurfaceSamples:4,positiveCoordinateVelocityNMPerPS:0.5,
            positiveCoordinateVelocityStandardErrorNMPerPS:0,converged:true,issues:[],evidenceFingerprint:evidenceID)
        try VivoQMMMSurfaceEnsemble.validate(result,request:request)
        return (request,result,ids)
    }

    private func observations(product:Bool,reactant:Bool,coordinate:Double,steps:Int=4)->[VivoQMMMTransmissionObservation] {
        (1...steps).map { .init(step:UInt64($0),coordinateNM:coordinate,inProduct:product,inReactant:reactant) }
    }

    @Test func pairedReactiveFluxPlateauReconstructsAndProducesCalculatedTransmissionEvidence() throws {
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
            maximumSteps:4,observeEverySteps:1,commitmentObservations:1,
            acceptance:.init(minimumSurfaceCheckpoints:4,minimumEffectiveFluxSamples:3.5,
                             maximumUnresolvedFluxFraction:0.3,maximumCoefficientStandardError:0.3,
                             minimumPlateauObservations:4,maximumPlateauRange:0.01))
        let plusProduct=VivoQMMMTransmissionBranch(velocitiesTimeReversed:false,
            observations:observations(product:true,reactant:false,coordinate:0.4),outcome:.productCommitted,finalCoordinateNM:0.4)
        let plusReactant=VivoQMMMTransmissionBranch(velocitiesTimeReversed:false,
            observations:observations(product:false,reactant:true,coordinate:-0.4),outcome:.reactantRecrossed,finalCoordinateNM:-0.4)
        let minusReactant=VivoQMMMTransmissionBranch(velocitiesTimeReversed:true,
            observations:observations(product:false,reactant:true,coordinate:-0.4),outcome:.reactantRecrossed,finalCoordinateNM:-0.4)
        var pairs=ids.enumerated().map { index,id in
            VivoQMMMTransmissionPair(sourceStateFingerprint:id,initialCoordinateNM:0.29,positiveCoordinateVelocityNMPerPS:1,
                positiveBranch:index<3 ? plusProduct:plusReactant,negativeBranch:minusReactant)
        }
        let result=try VivoQMMMDynamicalTransmission.analyze(request:request,systemFingerprint:systemID,
            providerFingerprint:providerID,pairs:pairs)
        #expect(result.converged)
        #expect(abs(result.transmissionCoefficient-0.75)<1e-12)
        #expect(result.plateauRange==0)
        #expect(result.coefficientHistory.count==4)
        #expect(result.effectiveFluxSamples==4)
        #expect(result.unresolvedFluxFraction==0.25)
        #expect(result.positiveProductCommittedCount==3)
        #expect(result.negativeReactantCommittedCount==4)
        #expect(result.unresolvedPairCount==1)
        try VivoQMMMDynamicalTransmission.validate(result,request:request)

        let analysisDocument=VivoQMMMDynamicalTransmissionAnalysisRequest(transmissionRequest:request,
            retainedSystemFingerprint:systemID,providerFingerprint:providerID,pairs:pairs)
        #expect(try analysisDocument.calculate()==result)

        let context=VivoKineticContext(compound:"synthetic",target:"protein",targetVariant:"reference",site:"site",
            chemicalState:"state",hostContext:"host",temperatureK:300,pH:7,ionicStrengthM:0.15)
        let rate=VivoQMMMFreeEnergyRateRequest(context:context,environment:.proteinEnvironment,freeEnergy:freeEnergy,
            transmissionProbability:1,transmissionOrigin:.assumed,
            transmissionEvidence:.init(source:"assumption",locator:"test"),samplingDescription:"fixture")
        let application=VivoQMMMComputedTransmissionApplicationRequest(analysisRequest:analysisDocument,result:result,rateRequest:rate)
        let computed=try application.calculate()
        #expect(computed.transmissionProbability==0.75)
        #expect(computed.transmissionOrigin == .calculated)
        #expect(computed.transmissionEvidence.sourceFingerprint==result.evidenceFingerprint.hex)
        #expect(computed.sampledFluxNormalization != nil)

        let unstable=VivoQMMMTransmissionBranch(velocitiesTimeReversed:false,
            observations:[.init(step:1,coordinateNM:0.4,inProduct:true,inReactant:false),
                          .init(step:2,coordinateNM:0.1,inProduct:false,inReactant:false),
                          .init(step:3,coordinateNM:0.4,inProduct:true,inReactant:false),
                          .init(step:4,coordinateNM:0.1,inProduct:false,inReactant:false)],
            outcome:.unresolved,finalCoordinateNM:0.1)
        pairs[0]=.init(sourceStateFingerprint:ids[0],initialCoordinateNM:0.29,positiveCoordinateVelocityNMPerPS:1,
            positiveBranch:unstable,negativeBranch:minusReactant)
        let failed=try VivoQMMMDynamicalTransmission.analyze(request:request,systemFingerprint:systemID,
            providerFingerprint:providerID,pairs:pairs)
        #expect(!failed.converged)
        #expect(failed.plateauRange>request.acceptance.maximumPlateauRange)
        #expect(throws:(any Error).self) { _ = try VivoQMMMDynamicalTransmission.applying(failed,transmissionRequest:request,to:rate) }
        #expect(throws:(any Error).self) {
            _ = try VivoQMMMDynamicalTransmission.analyze(request:request,
                systemFingerprint:try fingerprint("wrong-system"),providerFingerprint:providerID,pairs:pairs)
        }
    }

    @Test func metalSurfaceCollectionHandsExactStatesToPairedFreshUnbiasedNVEShooting() async throws {
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
            equilibrationSteps:2,productionSteps:128,candidateEverySteps:1,minimumAcceptedStepSeparation:1,
            surfaceToleranceNM:0.25,surfaceKernelBandwidthNM:0.10,desiredStates:4,
            minimumEffectiveSurfaceSamples:2,maximumAutocorrelationLag:3)
        let device=try VivoMetalDeviceSelector.productionDevice()
        #expect(device.hasUnifiedMemory)
        let surface=try await VivoQMMMSurfaceEnsemble.run(surfaceRequest,system:system,baseProvider:provider,device:device)
        #expect(surface.converged)
        #expect(surface.states.count==4)
        #expect(surface.effectiveSurfaceSamples>=2)
        #expect(surface.states.allSatisfy{$0.coordinateVelocityNMPerPS.isFinite && $0.statisticalWeight>0})
        try VivoQMMMSurfaceEnsemble.validate(surface,request:surfaceRequest)

        let shooting=VivoMDConfiguration(timeStepPS:0.0001,electrostatics:.cutoff,ensemble:.nve,thermostat:.none,
            targetTemperatureK:nil,frictionPerPS:nil,neighborListEnabled:false)
        let transmission=VivoQMMMDynamicalTransmissionRequest(surfaceRequest:surfaceRequest,surfaceEnsemble:surface,
            shootingDynamics:shooting,reactantCommitmentRangeNM:-2 ... 0,productCommitmentRangeNM:0.6 ... 2,
            maximumSteps:2,observeEverySteps:1,commitmentObservations:1,
            acceptance:.init(minimumSurfaceCheckpoints:4,minimumEffectiveFluxSamples:2,
                             maximumUnresolvedFluxFraction:0.99,maximumCoefficientStandardError:0.49,
                             minimumPlateauObservations:2,maximumPlateauRange:1))
        let result=try await VivoQMMMDynamicalTransmission.run(transmission,system:system,baseProvider:provider,device:device)
        #expect(result.pairs.count==surface.states.count)
        #expect(Set(result.pairs.map(\.sourceStateFingerprint))==Set(surface.states.map(\.stateFingerprint)))
        #expect(result.pairs.allSatisfy{$0.positiveBranch.observations.count==2 && $0.negativeBranch.observations.count==2})
        #expect(result.pairs.allSatisfy{$0.positiveBranch.velocitiesTimeReversed != $0.negativeBranch.velocitiesTimeReversed})
        try VivoQMMMDynamicalTransmission.validate(result,request:transmission)
    }
}
