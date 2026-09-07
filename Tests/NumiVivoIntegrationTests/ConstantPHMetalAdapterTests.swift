import Foundation
import Testing
@testable import NumiVivoKit

@Suite(.serialized) struct ConstantPHMetalAdapterTests {
    private func evidence(_ value:String)->VivoKineticEvidence { .init(source:"synthetic Metal constant-pH fixture",locator:value) }
    private func system(_ charge:Double)throws->VivoClassicalSystem {
        let structureID=try VivoCanonicalJSON.fingerprint(Data("same-latent-one-atom-structure".utf8))
        return .init(identifier:"constant-pH-metal-\(charge)",structureFingerprint:structureID,
            particles:[.init(index:0,atomIndex:0,typeIdentifier:"X",massDa:12,chargeE:charge,sigmaNM:0,epsilonKJPerMol:0)])
    }
    private func configuration()->VivoMDConfiguration {
        .init(timeStepPS:0.001,cutoffNM:1,neighborSkinNM:0,
              electrostatics:.cutoff,ensemble:.nve,thermostat:.none,targetTemperatureK:nil,
              frictionPerPS:nil,neighborListEnabled:false,randomSeed:991)
    }

    @Test func nativeFactoryPropagatesAndEvaluatesAcceptedMetalState() async throws {
        let config=configuration(),specifications=[
            VivoConstantPHMetalStateSpecification(identifier:"neutral",system:try system(0),configuration:config),
            VivoConstantPHMetalStateSpecification(identifier:"charged",system:try system(1),configuration:config)]
        let executable=try VivoConstantPHMetalStateFactory.make(specifications:specifications,samplingTemperatureK:300)
        #expect(executable.count==2)
        #expect(Set(executable.map(\.physicalManifoldFingerprint).map(\.hex)).count==1)
        #expect(Set(executable.map(\.hamiltonianFingerprint).map(\.hex)).count==2)
        let physical=try VivoConstantPHPhysicalState(stepIndex:0,timePS:0,
            positionsNM:[.zero],velocitiesNMPerPS:[.init(0.001,0,0)],periodicCell:nil)
        let e0=try await executable[0].potentialEnergyKJPerMol(physical)
        let e1=try await executable[1].potentialEnergyKJPerMol(physical)
        #expect(abs(e0)<1e-12 && abs(e1)<1e-12)
        let states=[
            VivoConstantPHStateDefinition(identifier:"neutral",boundProtonOffset:0,referenceSemigrandBiasKJPerMol:0,
                origin:.assumed,evidence:evidence("neutral"),neighbors:["charged"]),
            VivoConstantPHStateDefinition(identifier:"charged",boundProtonOffset:1,referenceSemigrandBiasKJPerMol:0,
                origin:.assumed,evidence:evidence("charged"),neighbors:["neutral"])]
        let sampler=try VivoConstantPHDiscreteMD(configuration:.init(identifier:"real-metal-cph",
            temperatureK:300,referencePH:7,targetPH:7,states:states,initialStateIdentifier:"neutral",
            mdStepsPerAttempt:2,attemptCount:2,seed:123),initialPhysicalState:physical,executableStates:executable)
        let result=try await sampler.run()
        #expect(result.finalCheckpoint.physicalState.stepIndex==4)
        #expect(abs(result.finalCheckpoint.physicalState.timePS-0.004)<1e-12)
        #expect(result.finalCheckpoint.attempts.count==2)
    }

    @Test func factoryRejectsNPTAndDifferentMassManifold() throws {
        var npt=configuration();npt.ensemble = .npt;npt.thermostat = .langevinMiddle;npt.targetTemperatureK=300
        npt.frictionPerPS=1;npt.barostat = .monteCarloIsotropic;npt.targetPressureBar=1
        #expect(throws:(any Error).self) {
            _=try VivoConstantPHMetalStateFactory.make(specifications:[
                .init(identifier:"a",system:try system(0),configuration:npt),
                .init(identifier:"b",system:try system(1),configuration:npt)],samplingTemperatureK:300)
        }
        var changed=try system(1);changed.particles[0].massDa=13
        #expect(throws:(any Error).self) {
            _=try VivoConstantPHMetalStateFactory.make(specifications:[
                .init(identifier:"a",system:try system(0),configuration:configuration()),
                .init(identifier:"b",system:changed,configuration:configuration())],samplingTemperatureK:300)
        }
    }
}