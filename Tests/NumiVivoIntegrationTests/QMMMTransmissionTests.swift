import Foundation
import Testing
@testable import NumiVivoKit

@Suite(.serialized) struct QMMMTransmissionTests {
    private func qualifiedPMF() throws -> VivoQMMMQualifiedActivationFreeEnergy {
        let coordinate=VivoQMMMReactionCoordinate(identifier:"synthetic-transfer",kind:.distanceDifference,atomIndices:[0,1,2,3])
        let windows=[VivoQMMMUmbrellaWindow(identifier:"left",centerNM:-0.15,forceConstantKJPerMolNM2:10),
                     VivoQMMMUmbrellaWindow(identifier:"right",centerNM:0.15,forceConstantKJPerMolNM2:10)]
        var values:[Double]=[]
        for _ in 0..<20 { for i in 0...60 {
            let x=Double(i)*0.01
            values.append(-0.30+x);values.append(0.30-x)
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
        let protocolID=try VivoCanonicalJSON.fingerprint(Data("transmission-test-protocol".utf8))
        let provenance=VivoQMMMFreeEnergyProvenance(
            structureFingerprint:try VivoCanonicalJSON.fingerprint(Data("structure".utf8)),
            systemFingerprint:try VivoCanonicalJSON.fingerprint(Data("system".utf8)),
            baseProviderFingerprint:try VivoCanonicalJSON.fingerprint(Data("provider".utf8)),
            dynamicsFingerprint:try VivoCanonicalJSON.fingerprint(Data("source-dynamics".utf8)),
            samplingExecution:.init(requestFingerprint:try VivoCanonicalJSON.fingerprint(Data("execution".utf8)),
                                    replicaProtocolFingerprint:protocolID),
            chemicalState:"state",environment:.proteinEnvironment,environmentIdentifier:"host",
            methodDescription:"synthetic transmission estimator fixture",
            reactionConnectivity:.init(connectivityFingerprint:try VivoCanonicalJSON.fingerprint(Data("connectivity".utf8)),
                reactantEndpointIdentifier:"R",productEndpointIdentifier:"P",mappedReactionAtomIndices:coordinate.atomIndices))
        return try VivoQMMMQualifiedActivationFreeEnergy(analysis:analysis,provenance:provenance,
            fluxNormalization:.init(surfaceToReactantDensityPerNM:try VivoQMMMFreeEnergyQualification.surfaceToReactantDensityPerNM(analysis),
                                    inverseMassMetricPerDa:1))
    }

    @Test func fluxWeightedEstimatorReconstructsAndProducesCalculatedTransmissionEvidence() throws {
        let freeEnergy=try qualifiedPMF()
        let source=VivoMDConfiguration(timeStepPS:0.001,electrostatics:.cutoff,ensemble:.nvt,thermostat:.langevinMiddle,
            targetTemperatureK:300,frictionPerPS:1,neighborListEnabled:false)
        let shooting=VivoMDConfiguration(timeStepPS:0.001,electrostatics:.cutoff,ensemble:.nve,thermostat:.none,
            targetTemperatureK:nil,frictionPerPS:nil,neighborListEnabled:false)
        let sourceID=try source.fingerprint(),systemID=try VivoCanonicalJSON.fingerprint(Data("dummy-system".utf8))
        let checkpoints=(0..<4).map { i in VivoMDCheckpoint(systemFingerprint:systemID,configurationFingerprint:sourceID,
            acceptedStep:UInt64(i),timePS:Double(i)*0.001,positionsNM:[.zero],velocitiesNMPerPS:[.zero],periodicCell:nil) }
        let request=VivoQMMMDynamicalTransmissionRequest(freeEnergy:freeEnergy,sourceDynamics:source,shootingDynamics:shooting,
            surfaceCheckpoints:checkpoints,reactantCommitmentRangeNM:-0.5 ... -0.35,productCommitmentRangeNM:0.35 ... 0.5,
            maximumSteps:100,observeEverySteps:1,commitmentObservations:1,
            acceptance:.init(minimumSurfaceCheckpoints:4,minimumEffectiveFluxSamples:3.5,
                             maximumUnresolvedFluxFraction:0.1,maximumCoefficientStandardError:0.3))
        let ids=try checkpoints.map{$0.fingerprint()}
        let trajectories=[
            VivoQMMMTransmissionTrajectory(sourceCheckpointFingerprint:ids[0],initialCoordinateNM:0.29,
                productDirectedCoordinateVelocityNMPerPS:1,velocitiesTimeReversed:false,outcome:.productCommitted,committedSteps:10,finalCoordinateNM:0.4),
            VivoQMMMTransmissionTrajectory(sourceCheckpointFingerprint:ids[1],initialCoordinateNM:0.29,
                productDirectedCoordinateVelocityNMPerPS:1,velocitiesTimeReversed:true,outcome:.productCommitted,committedSteps:12,finalCoordinateNM:0.4),
            VivoQMMMTransmissionTrajectory(sourceCheckpointFingerprint:ids[2],initialCoordinateNM:0.29,
                productDirectedCoordinateVelocityNMPerPS:1,velocitiesTimeReversed:false,outcome:.productCommitted,committedSteps:9,finalCoordinateNM:0.4),
            VivoQMMMTransmissionTrajectory(sourceCheckpointFingerprint:ids[3],initialCoordinateNM:0.29,
                productDirectedCoordinateVelocityNMPerPS:1,velocitiesTimeReversed:true,outcome:.reactantRecrossed,committedSteps:7,finalCoordinateNM:-0.4)]
        let providerID=try VivoCanonicalJSON.fingerprint(Data("dummy-provider".utf8))
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
        unresolved[0]=.init(sourceCheckpointFingerprint:ids[0],initialCoordinateNM:0.29,
            productDirectedCoordinateVelocityNMPerPS:1,velocitiesTimeReversed:false,outcome:.unresolved,committedSteps:100,finalCoordinateNM:0.1)
        let failed=try VivoQMMMDynamicalTransmission.analyze(request:request,systemFingerprint:systemID,
            providerFingerprint:providerID,trajectories:unresolved)
        #expect(!failed.converged)
        #expect(failed.unresolvedFluxFraction==0.25)
        #expect(throws:(any Error).self) { _ = try VivoQMMMDynamicalTransmission.applying(failed,transmissionRequest:request,to:rate) }
    }
}
