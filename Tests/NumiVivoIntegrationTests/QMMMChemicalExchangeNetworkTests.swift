import Foundation
import Testing
@testable import NumiVivoKit

@Suite(.serialized) struct QMMMChemicalExchangeNetworkTests {
    private func fingerprint(_ value:String) throws -> VivoFingerprint {
        try VivoCanonicalJSON.fingerprint(Data(value.utf8))
    }

    func replicatedPathway(state:String,executionPrefix:String,seedBase:UInt64) throws -> VivoQMMMPathwayRate {
        let coordinate=VivoQMMMReactionCoordinate(identifier:"exchange-transfer",kind:.distance,atomIndices:[0,1])
        let windows=[VivoQMMMUmbrellaWindow(identifier:"left",centerNM:-0.10,forceConstantKJPerMolNM2:20),
                     VivoQMMMUmbrellaWindow(identifier:"right",centerNM:0.10,forceConstantKJPerMolNM2:20)]
        var values:[Double]=[]
        for _ in 0..<20 { for i in 0...40 { values.append(-0.20+Double(i)*0.01);values.append(0.20-Double(i)*0.01) } }
        let traces=windows.enumerated().map { index,window in
            VivoQMMMUmbrellaTrace(window:window,randomSeed:seedBase+UInt64(index),coordinateNM:values,
                potentialEnergyKJPerMol:[Double](repeating:0,count:values.count))
        }
        let config=VivoQMMMFreeEnergyAnalysisConfiguration(bins:61,kernelBandwidthNM:0.02,
            reactantRangeNM:-0.18 ... -0.08,dividingSurfaceNM:0.18,
            minimumDecorrelatedSamplesPerWindow:20,minimumAdjacentOverlap:0)
        let analysis=try VivoQMMMFreeEnergy.analyze(coordinate:coordinate,temperatureK:300,traces:traces,configuration:config)
        #expect(analysis.converged)
        let protocolID=try fingerprint("exchange-shared-protocol")
        func request(replica:Int,seeds:[UInt64]) throws -> VivoQMMMFreeEnergyRateRequest {
            let replicaTraces=zip(traces,seeds).map { trace,seed in
                VivoQMMMUmbrellaTrace(window:trace.window,randomSeed:seed,coordinateNM:trace.coordinateNM,
                    potentialEnergyKJPerMol:trace.potentialEnergyKJPerMol)
            }
            let replicaAnalysis=try VivoQMMMFreeEnergy.analyze(coordinate:coordinate,temperatureK:300,
                traces:replicaTraces,configuration:config)
            let provenance=VivoQMMMFreeEnergyProvenance(
                structureFingerprint:try fingerprint("exchange-structure"),
                systemFingerprint:try fingerprint("exchange-system"),
                baseProviderFingerprint:try fingerprint("exchange-provider"),
                dynamicsFingerprint:try fingerprint("exchange-dynamics"),
                samplingExecution:.init(requestFingerprint:try fingerprint("\(executionPrefix)-\(replica)"),
                    replicaProtocolFingerprint:protocolID),
                chemicalState:state,environment:.proteinEnvironment,environmentIdentifier:"exchange-host",
                methodDescription:"synthetic exchange-network qualification fixture",
                reactionConnectivity:.init(connectivityFingerprint:try fingerprint("exchange-connectivity"),
                    reactantEndpointIdentifier:"R",productEndpointIdentifier:"P",mappedReactionAtomIndices:coordinate.atomIndices))
            let qualified=try VivoQMMMQualifiedActivationFreeEnergy(analysis:replicaAnalysis,provenance:provenance,
                fluxNormalization:.init(surfaceToReactantDensityPerNM:try VivoQMMMFreeEnergyQualification.surfaceToReactantDensityPerNM(replicaAnalysis),
                    inverseMassMetricPerDa:1e-12))
            let context=VivoKineticContext(compound:"synthetic",target:"protein",targetVariant:"reference",site:"site",
                chemicalState:state,hostContext:"exchange-host",temperatureK:300,pH:7,ionicStrengthM:0.15)
            return VivoQMMMFreeEnergyRateRequest(context:context,environment:.proteinEnvironment,freeEnergy:qualified,
                transmissionProbability:1,transmissionOrigin:.assumed,
                transmissionEvidence:.init(source:"assumption",locator:"exchange-network test"),
                samplingDescription:"independent synthetic exchange-state replica")
        }
        let first=try request(replica:0,seeds:[seedBase+100,seedBase+101])
        let second=try request(replica:1,seeds:[seedBase+200,seedBase+201])
        let replicatedRequest=VivoQMMMReplicatedFreeEnergyRateRequest(replicas:[first,second],
            agreement:.init(maximumLogRateRange:1e-8,maximumProfileBarrierRangeKJPerMol:1e-8))
        let result=try VivoQMMMReplicatedFreeEnergyRate.calculate(replicatedRequest)
        #expect(result.converged)
        return .init(identifier:"path",request:replicatedRequest,result:result)
    }

    private func assumedRate(_ value:Double,_ label:String)->VivoKineticParameter {
        .init(value:value,unit:.perSecond,origin:.assumed,uncertainty:.unknown,
              evidence:.init(source:"assumption",locator:label))
    }

    private func state(_ identifier:String,_ population:Double,_ pathway:VivoQMMMPathwayRate)->VivoQMMMChemicalExchangeState {
        .init(identifier:identifier,initialPopulation:population,populationOrigin:.assumed,
              populationEvidence:.init(source:"assumption",locator:"initial population \(identifier)"),pathways:[pathway])
    }

    @Test func equalChemicalRatesRemainSingleExponentialUnderArbitraryExchange() throws {
        let a=try replicatedPathway(state:"A",executionPrefix:"A-execution",seedBase:10)
        let b=try replicatedPathway(state:"B",executionPrefix:"B-execution",seedBase:1000)
        let k=a.result.geometricMeanRatePerSecond
        #expect(abs(k-b.result.geometricMeanRatePerSecond)/k<1e-12)
        let request=VivoQMMMChemicalExchangeNetworkRequest(identifier:"equal-rate-exchange",states:[
            state("A",0.25,a),state("B",0.75,b)],exchangeEdges:[
                .init(fromStateIdentifier:"A",toStateIdentifier:"B",rate:assumedRate(5*k,"A to B")),
                .init(fromStateIdentifier:"B",toStateIdentifier:"A",rate:assumedRate(7*k,"B to A"))],
            observationTimesSeconds:[0,0.1/k,1/k])
        let result=try VivoQMMMChemicalExchangeNetwork.calculate(request)
        try VivoQMMMChemicalExchangeNetwork.validate(result,request:request)
        #expect(result.stateChemicalRatesPerSecond.count==2)
        for observation in result.observations {
            let expected=exp(-k*observation.timeSeconds)
            #expect(abs(observation.survivalProbability-expected)<=max(1e-11,expected*2e-10))
            #expect(abs(observation.instantaneousHazardPerSecond-k)<=max(1e-10,k*2e-10))
            if observation.timeSeconds>0 {
                #expect(abs((observation.apparentFirstOrderRatePerSecond ?? -1)-k)<=max(1e-10,k*2e-10))
            }
            #expect(abs(observation.unreactedProbabilityByState.reduce(0,+)-observation.survivalProbability)<1e-12)
        }

        let time=request.observationTimesSeconds[2],predicted=result.observations[2].survivalProbability
        let measured=VivoQMMMChemicalExchangeValidationTarget(identifier:"time-course-1",timeSeconds:time,
            survivalProbability:predicted+0.01,standardDeviationProbability:0.02,
            context:a.request.replicas[0].context,
            source:.init(source:"synthetic measured time course",locator:"row 1",
                         sourceFingerprint:(try fingerprint("exchange-measured-time-course")).hex))
        let validationRequest=VivoQMMMChemicalExchangeValidationRequest(identifier:"exchange-validation",
            networkRequest:request,networkResult:result,targets:[measured],
            maximumAbsoluteProbabilityError:0.02,maximumStandardizedResidual:1)
        let validation=try VivoQMMMChemicalExchangeValidation.calculate(validationRequest)
        #expect(validation.converged)
        #expect(validation.points.count==1 && validation.points[0].comparable && validation.points[0].passed==true)
        #expect(abs((validation.points[0].standardizedResidual ?? -1)-0.5)<1e-10)
        try VivoQMMMChemicalExchangeValidation.validate(validation,request:validationRequest)

        let offGrid=VivoQMMMChemicalExchangeValidationTarget(identifier:"off-grid",timeSeconds:0.5/k,
            survivalProbability:0.5,context:a.request.replicas[0].context,
            source:.init(source:"synthetic measured time course",locator:"row 2",
                         sourceFingerprint:(try fingerprint("exchange-measured-off-grid")).hex))
        let offGridResult=try VivoQMMMChemicalExchangeValidation.calculate(.init(identifier:"off-grid-validation",
            networkRequest:request,networkResult:result,targets:[offGrid]))
        #expect(!offGridResult.converged)
        #expect(offGridResult.points[0].comparable==false)
    }

    @Test func unequalStateRatesProduceTimeDependentHazardAndInvalidEdgesReject() throws {
        let a=try replicatedPathway(state:"A",executionPrefix:"A2-execution",seedBase:3000)
        var b=try replicatedPathway(state:"B",executionPrefix:"B2-execution",seedBase:5000)
        var modified=b.request
        modified.replicas=modified.replicas.map { request in
            var copy=request
            copy.transmissionProbability=0.5
            return copy
        }
        b = .init(identifier:"path",request:modified,result:try VivoQMMMReplicatedFreeEnergyRate.calculate(modified))
        let k=a.result.geometricMeanRatePerSecond
        let request=VivoQMMMChemicalExchangeNetworkRequest(identifier:"unequal-rate-exchange",states:[
            state("A",0.5,a),state("B",0.5,b)],exchangeEdges:[
                .init(fromStateIdentifier:"A",toStateIdentifier:"B",rate:assumedRate(0.2*k,"A to B")),
                .init(fromStateIdentifier:"B",toStateIdentifier:"A",rate:assumedRate(0.1*k,"B to A"))],
            observationTimesSeconds:[0,0.1/k,1/k,3/k])
        let result=try VivoQMMMChemicalExchangeNetwork.calculate(request)
        let hazards=result.observations.map(\.instantaneousHazardPerSecond)
        #expect((hazards.max() ?? 0)-(hazards.min() ?? 0)>k*1e-4)
        #expect(result.observations.last!.reactedProbability>result.observations[1].reactedProbability)

        let duplicate=VivoQMMMChemicalExchangeNetworkRequest(identifier:"duplicate-edge",states:request.states,
            exchangeEdges:[request.exchangeEdges[0],request.exchangeEdges[0]],observationTimesSeconds:[0])
        #expect(throws:(any Error).self) { _ = try VivoQMMMChemicalExchangeNetwork.calculate(duplicate) }
        let unknown=VivoQMMMChemicalExchangeNetworkRequest(identifier:"unknown-edge",states:request.states,
            exchangeEdges:[.init(fromStateIdentifier:"A",toStateIdentifier:"missing",rate:assumedRate(k,"bad"))],
            observationTimesSeconds:[0])
        #expect(throws:(any Error).self) { _ = try VivoQMMMChemicalExchangeNetwork.calculate(unknown) }
    }
}
