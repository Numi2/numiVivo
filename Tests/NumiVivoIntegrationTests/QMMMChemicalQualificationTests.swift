import Foundation
import Testing
@testable import NumiVivoKit

@Suite(.serialized) struct QMMMChemicalQualificationTests {
    private func fingerprint(_ text:String)throws->VivoFingerprint {
        try VivoCanonicalJSON.fingerprint(Data(text.utf8))
    }

    private func replicated(transmission:Double,transmissionID:String,
                            origin:VivoKineticOrigin = .calculated) throws -> VivoQMMMReplicatedFreeEnergyRateRequest {
        let coordinate=VivoQMMMReactionCoordinate(identifier:"qualification-transfer",kind:.distance,atomIndices:[0,1])
        let windows=[VivoQMMMUmbrellaWindow(identifier:"left",centerNM:-0.10,forceConstantKJPerMolNM2:20),
                     VivoQMMMUmbrellaWindow(identifier:"right",centerNM:0.10,forceConstantKJPerMolNM2:20)]
        var samples:[Double]=[]
        for _ in 0..<20 { for i in 0...40 { samples.append(-0.20+Double(i)*0.01);samples.append(0.20-Double(i)*0.01) } }
        let config=VivoQMMMFreeEnergyAnalysisConfiguration(bins:61,kernelBandwidthNM:0.02,
            reactantRangeNM:-0.18 ... -0.08,dividingSurfaceNM:0.18,
            minimumDecorrelatedSamplesPerWindow:20,minimumAdjacentOverlap:0)
        let protocolID=try fingerprint("qualification-shared-protocol")
        var requests:[VivoQMMMFreeEnergyRateRequest]=[]
        for replica in 0..<2 {
            let seeds=[UInt64(100+replica*10),UInt64(101+replica*10)]
            let traces=zip(windows,seeds).map { window,seed in
                VivoQMMMUmbrellaTrace(window:window,randomSeed:seed,coordinateNM:samples,
                    potentialEnergyKJPerMol:[Double](repeating:0,count:samples.count))
            }
            let analysis=try VivoQMMMFreeEnergy.analyze(coordinate:coordinate,temperatureK:300,traces:traces,configuration:config)
            #expect(analysis.converged)
            let provenance=VivoQMMMFreeEnergyProvenance(
                structureFingerprint:try fingerprint("qualification-structure"),
                systemFingerprint:try fingerprint("qualification-system"),
                baseProviderFingerprint:try fingerprint("qualification-provider"),
                dynamicsFingerprint:try fingerprint("qualification-dynamics"),
                samplingExecution:.init(requestFingerprint:try fingerprint("qualification-execution-\(replica)"),
                    replicaProtocolFingerprint:protocolID),chemicalState:"state",
                environment:.proteinEnvironment,environmentIdentifier:"host",
                methodDescription:"synthetic transmission-sensitivity qualification fixture",
                reactionConnectivity:.init(connectivityFingerprint:try fingerprint("qualification-connectivity"),
                    reactantEndpointIdentifier:"R",productEndpointIdentifier:"P",mappedReactionAtomIndices:coordinate.atomIndices))
            let qualified=try VivoQMMMQualifiedActivationFreeEnergy(analysis:analysis,provenance:provenance,
                fluxNormalization:.init(surfaceToReactantDensityPerNM:try VivoQMMMFreeEnergyQualification.surfaceToReactantDensityPerNM(analysis),
                    inverseMassMetricPerDa:1e-12))
            let context=VivoKineticContext(compound:"synthetic",target:"protein",targetVariant:"reference",site:"site",
                chemicalState:"state",hostContext:"host",temperatureK:300,pH:7,ionicStrengthM:0.15)
            let evidence:VivoKineticEvidence
            if origin == .assumed {
                evidence=.init(source:"assumption",locator:transmissionID)
            } else {
                evidence=.init(source:"synthetic calculated transmission",locator:transmissionID,
                    sourceFingerprint:(try fingerprint(transmissionID)).hex)
            }
            requests.append(.init(context:context,environment:.proteinEnvironment,freeEnergy:qualified,
                transmissionProbability:transmission,transmissionOrigin:origin,transmissionEvidence:evidence,
                samplingDescription:"synthetic independently replicated PMF"))
        }
        return .init(replicas:requests,
            agreement:.init(maximumLogRateRange:1e-8,maximumProfileBarrierRangeKJPerMol:1e-8))
    }

    @Test func transmissionProtocolSensitivityHoldsPMFFixedAndMeasuresLogRateShift() throws {
        let baselineRequest=try replicated(transmission:1,transmissionID:"transmission-baseline")
        let baselineResult=try VivoQMMMReplicatedFreeEnergyRate.calculate(baselineRequest)
        #expect(baselineResult.converged)
        var variantRequest=baselineRequest
        variantRequest.replicas=try baselineRequest.replicas.map { request in
            var copy=request
            copy.transmissionProbability=0.9
            copy.transmissionOrigin = .calculated
            copy.transmissionEvidence=.init(source:"synthetic calculated transmission",locator:"transmission-variant",
                sourceFingerprint:(try fingerprint("transmission-variant")).hex)
            return copy
        }
        let variantResult=try VivoQMMMReplicatedFreeEnergyRate.calculate(variantRequest)
        #expect(variantResult.converged)
        let variant=VivoQMMMQualificationVariant(identifier:"shorter-shooting-horizon",
            dimension:.transmissionProtocol,request:variantRequest,result:variantResult,
            description:"synthetic perturbation representing an alternate qualified surface/shooting protocol")
        let qualificationRequest=VivoQMMMChemicalQualificationRequest(identifier:"transmission-protocol-test",
            baselineRequest:baselineRequest,baselineResult:baselineResult,variants:[variant],
            criteria:[.init(dimension:.transmissionProtocol,requiredVariants:1,maximumAbsoluteLogRateShift:0.2)])
        let result=try VivoQMMMChemicalQualification.calculate(qualificationRequest)
        #expect(result.converged)
        #expect(result.sensitivity.count==1)
        #expect(result.sensitivity[0].dimension == .transmissionProtocol)
        #expect(abs(result.sensitivity[0].logRateShiftFromBaseline-log(0.9))<1e-12)
        try VivoQMMMChemicalQualification.validate(result,request:qualificationRequest)

        let strict=VivoQMMMChemicalQualificationRequest(identifier:"transmission-protocol-strict",
            baselineRequest:baselineRequest,baselineResult:baselineResult,variants:[variant],
            criteria:[.init(dimension:.transmissionProtocol,requiredVariants:1,maximumAbsoluteLogRateShift:0.05)])
        let failed=try VivoQMMMChemicalQualification.calculate(strict)
        #expect(!failed.converged)
        #expect(failed.sensitivity[0].passed==false)
    }

    @Test func transmissionProtocolDimensionRejectsAssumedTransmission() throws {
        let baselineRequest=try replicated(transmission:1,transmissionID:"calculated-baseline")
        let baselineResult=try VivoQMMMReplicatedFreeEnergyRate.calculate(baselineRequest)
        var assumedRequest=baselineRequest
        assumedRequest.replicas=baselineRequest.replicas.map { request in
            var copy=request
            copy.transmissionProbability=0.9
            copy.transmissionOrigin = .assumed
            copy.transmissionEvidence=.init(source:"assumption",locator:"not dynamical evidence")
            return copy
        }
        let assumedResult=try VivoQMMMReplicatedFreeEnergyRate.calculate(assumedRequest)
        let variant=VivoQMMMQualificationVariant(identifier:"invalid-assumed-transmission",
            dimension:.transmissionProtocol,request:assumedRequest,result:assumedResult,
            description:"must not qualify as a dynamical-transmission protocol variant")
        let qualification=VivoQMMMChemicalQualificationRequest(identifier:"reject-assumed-transmission",
            baselineRequest:baselineRequest,baselineResult:baselineResult,variants:[variant],
            criteria:[.init(dimension:.transmissionProtocol)])
        #expect(throws:(any Error).self) { _ = try VivoQMMMChemicalQualification.calculate(qualification) }
    }
}
