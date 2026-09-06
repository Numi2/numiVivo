import Foundation
import Testing
@testable import NumiVivoKit

@Suite(.serialized) struct QMMMFreeEnergyTests {
    private func system() throws -> VivoClassicalSystem {
        let structureID=try VivoCanonicalJSON.fingerprint(Data("qmmm-free-energy-fixture".utf8))
        return VivoClassicalSystem(identifier:"free-energy-fixture",structureFingerprint:structureID,particles:[
            .init(index:0,atomIndex:0,typeIdentifier:"A",massDa:12,chargeE:0,sigmaNM:0,epsilonKJPerMol:0),
            .init(index:1,atomIndex:1,typeIdentifier:"B",massDa:12,chargeE:0,sigmaNM:0,epsilonKJPerMol:0)])
    }

    @Test func umbrellaProviderIsAnalyticConservativeAndAdditive() async throws {
        let system=try system(),systemID=try system.fingerprint()
        let baseID=try VivoCanonicalJSON.fingerprint(Data("base-provider".utf8))
        let base=try VivoMDCandidateForceProvider(fingerprint:baseID,retainedSystemFingerprint:systemID,boundary:.finiteCluster,
            supportsCellMoves:false,maximumAcceptedResidual:1e-6,molecularConnectivitySystem:system) { geometry in
            try VivoMDCandidateForceEvaluation(providerFingerprint:baseID,geometry:geometry,additionalEnergyKJPerMol:7,
                physicalParticleForcesKJPerMolNM:[.zero,.zero],derivativeMethod:"synthetic",convergenceResidual:0,requiredResidual:1e-6)
        }
        let coordinate=VivoQMMMReactionCoordinate(identifier:"forming",kind:.distance,atomIndices:[0,1])
        let window=VivoQMMMUmbrellaWindow(identifier:"w",centerNM:0.2,forceConstantKJPerMolNM2:100)
        let provider=try VivoQMMMUmbrellaBias.provider(base:base,system:system,coordinate:coordinate,window:window)
        func evaluate(_ x:Double) async throws -> VivoMDCandidateForceEvaluation {
            let geometry=try VivoMDCandidateGeometry(particlePositionsNM:[.init(x,0,0),.zero],periodicCell:nil)
            return try await provider.evaluate(geometry)
        }
        let value=try await evaluate(0.3)
        #expect(abs(value.additionalEnergyKJPerMol-7.5)<1e-12)
        #expect(abs(value.physicalParticleForcesKJPerMolNM[0].x+10)<1e-12)
        #expect(abs(value.physicalParticleForcesKJPerMolNM[1].x-10)<1e-12)
        let h=1e-6,plus=try await evaluate(0.3+h),minus=try await evaluate(0.3-h)
        let finiteDifference=-(plus.additionalEnergyKJPerMol-minus.additionalEnergyKJPerMol)/(2*h)
        #expect(abs(finiteDifference-value.physicalParticleForcesKJPerMolNM[0].x)<1e-7)
        #expect((value.physicalParticleForcesKJPerMolNM[0]+value.physicalParticleForcesKJPerMolNM[1]).norm<1e-12)
    }

    @Test func unbinnedMBARReconstructsBarrierAndProteinRate() throws {
        let temperature=300.0,beta=1/(0.00831446261815324*temperature)
        let windows=[-0.24,-0.12,0,0.12,0.24].enumerated().map {
            VivoQMMMUmbrellaWindow(identifier:"w\($0.offset)",centerNM:$0.element,forceConstantKJPerMolNM2:180)
        }
        func unbiased(_ x:Double)->Double {
            // Symmetric double well: U(±0.2)=0, U(0)=20 kJ/mol.
            20*pow(1-pow(x/0.2,2),2)
        }
        struct RNG { var x:UInt64;mutating func uniform()->Double{x = x&*6364136223846793005&+1442695040888963407;return Double(x>>11)/Double(UInt64(1)<<53)} }
        var traces:[VivoQMMMUmbrellaTrace]=[]
        for (index,window) in windows.enumerated() {
            var rng=RNG(x:UInt64(1234+index)),x=window.centerNM,values:[Double]=[],energies:[Double]=[]
            for step in 0..<140000 {
                let proposal=x+(rng.uniform()-0.5)*0.035
                let old=unbiased(x)+window.biasKJPerMol(x),new=unbiased(proposal)+window.biasKJPerMol(proposal)
                if log(max(rng.uniform(),Double.leastNonzeroMagnitude)) < -beta*(new-old) { x=proposal }
                if step>=20000 && step%20==0 { values.append(x);energies.append(unbiased(x)+window.biasKJPerMol(x)) }
            }
            traces.append(.init(window:window,randomSeed:UInt64(index+1),coordinateNM:values,potentialEnergyKJPerMol:energies))
        }
        let coordinate=VivoQMMMReactionCoordinate(identifier:"symmetric-transfer",kind:.distanceDifference,atomIndices:[0,1,2,3])
        let cfg=VivoQMMMFreeEnergyAnalysisConfiguration(bins:121,kernelBandwidthNM:0.012,reactantRangeNM:-0.27 ... -0.12,
            dividingSurfaceNM:0,minimumDecorrelatedSamplesPerWindow:60,minimumAdjacentOverlap:0.005)
        let result=try VivoQMMMFreeEnergy.analyze(coordinate:coordinate,temperatureK:temperature,traces:traces,configuration:cfg)
        #expect(result.converged)
        #expect(result.activationFreeEnergyKJPerMol>15 && result.activationFreeEnergyKJPerMol<25)
        #expect(result.conditionalStandardDeviationKJPerMol.isFinite)
        try VivoQMMMFreeEnergy.validate(result)
        let context=VivoKineticContext(compound:"synthetic",target:"protein",targetVariant:"reference",site:"reactive-site",
            chemicalState:"prepared-reactive-state",hostContext:"protein-pocket",temperatureK:temperature,pH:7.4,ionicStrengthM:0.15)
        let transmission=VivoKineticEvidence(source:"assumption",locator:"unit test")
        let request=VivoQMMMFreeEnergyRateRequest(context:context,environment:.proteinEnvironment,freeEnergy:result,
            transmissionProbability:1,transmissionOrigin:.assumed,transmissionEvidence:transmission,
            samplingDescription:"synthetic converged protein-environment PMF fixture")
        let rate=try VivoQMMMFreeEnergyRate.calculate(request)
        #expect(rate.estimate.ratePerSecond.isFinite && rate.estimate.ratePerSecond>0)
        #expect(rate.barrier.origin == .calculated)
        #expect(rate.parameter.origin == .assumed)
    }
}
