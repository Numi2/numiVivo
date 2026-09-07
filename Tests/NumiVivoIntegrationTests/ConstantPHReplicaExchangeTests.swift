import Foundation
import Testing
@testable import NumiVivoKit

@Suite(.serialized) struct ConstantPHReplicaExchangeTests {
    private func fingerprint(_ value:String)throws->VivoFingerprint { try VivoCanonicalJSON.fingerprint(Data(value.utf8)) }
    private func evidence(_ value:String)->VivoKineticEvidence { .init(source:"synthetic replica-exchange regression",locator:value) }
    private func physical(_ index:Int)throws->VivoConstantPHPhysicalState {
        try .init(stepIndex:0,timePS:0,positionsNM:[.init(Double(index),0,0)],velocitiesNMPerPS:[.zero],periodicCell:nil)
    }
    private func executable(_ id:String,_ manifold:VivoFingerprint)throws->VivoConstantPHExecutableState {
        .init(identifier:id,physicalManifoldFingerprint:manifold,hamiltonianFingerprint:try fingerprint("hx-\(id)"),
              potentialEnergyKJPerMol:{_ in 0},propagate:{state,steps in
            try .init(stepIndex:state.stepIndex+steps,timePS:state.timePS+Double(steps)*0.001,
                      positionsNM:state.positionsNM,velocitiesNMPerPS:state.velocitiesNMPerPS,periodicCell:state.periodicCell)
        })
    }
    private func configuration(cycles:Int=20_000)->VivoConstantPHReplicaExchangeConfiguration {
        let states=[
            VivoConstantPHStateDefinition(identifier:"protonated",boundProtonOffset:1,referenceSemigrandBiasKJPerMol:0,
                origin:.assumed,evidence:evidence("protonated"),neighbors:["deprotonated"]),
            VivoConstantPHStateDefinition(identifier:"deprotonated",boundProtonOffset:0,referenceSemigrandBiasKJPerMol:0,
                origin:.assumed,evidence:evidence("deprotonated"),neighbors:["protonated"])]
        return .init(identifier:"two-lane-pH-rex",temperatureK:300,referencePH:7,pHLadder:[7,8],states:states,
                     initialStateIdentifier:"protonated",mdStepsPerCycle:1,cycleCount:cycles,seed:314159)
    }

    @Test func pHLanesRecoverTheirSemigrandPopulations() async throws {
        let cfg=configuration(),manifold=try fingerprint("rex-manifold")
        let sampler=try VivoConstantPHReplicaExchange(configuration:cfg,
            initialPhysicalStates:[physical(0),physical(1)],
            executableStates:[executable("protonated",manifold),executable("deprotonated",manifold)])
        let result=try await sampler.run()
        func population(_ pH:Double,_ state:String)throws->Double {
            try #require(result.populations.first{$0.targetPH==pH && $0.stateIdentifier==state}).population
        }
        #expect(abs(try population(7,"protonated")-0.5)<0.025)
        #expect(abs(try population(8,"deprotonated")-10.0/11.0)<0.025)
        #expect(result.swapAcceptanceFraction>0 && result.swapAcceptanceFraction<1)
        #expect(result.finalCheckpoint.chemicalMoveAttempts==cfg.cycleCount*cfg.pHLadder.count)
        #expect(result.finalCheckpoint.lanes.allSatisfy{$0.physicalState.stepIndex==UInt64(cfg.cycleCount)})
        let unfavorable=try #require(result.finalCheckpoint.swapAttempts.first{
            $0.lowerStateIdentifier=="protonated" && $0.upperStateIdentifier=="deprotonated"})
        #expect(abs(unfavorable.logAcceptanceProbability + log(10)) < 1e-12)
        let favorable=try #require(result.finalCheckpoint.swapAttempts.first{
            $0.lowerStateIdentifier=="deprotonated" && $0.upperStateIdentifier=="protonated"})
        #expect(favorable.logAcceptanceProbability == 0)
    }

    @Test func replayIsBitwiseDeterministicAndWalkersExchangeLanes() async throws {
        let cfg=configuration(cycles:2000),manifold=try fingerprint("replay-rex-manifold")
        let states=try [executable("protonated",manifold),executable("deprotonated",manifold)]
        let initial=try [physical(0),physical(1)]
        let first=try await VivoConstantPHReplicaExchange(configuration:cfg,initialPhysicalStates:initial,executableStates:states).run()
        let second=try await VivoConstantPHReplicaExchange(configuration:cfg,initialPhysicalStates:initial,executableStates:states).run()
        #expect(first==second)
        #expect(first.finalCheckpoint.swapAttempts.contains(\.accepted))
        #expect(Set(first.finalCheckpoint.lanes.map(\.walkerIdentifier))==Set([0,1]))
    }

    @Test func ladderMustBeStrictAndExecutablesShareOneManifold() throws {
        var bad=configuration(cycles:2);bad.pHLadder=[7,7]
        #expect(throws:(any Error).self){try bad.validate()}
        let cfg=configuration(cycles:2)
        let a=try executable("protonated",fingerprint("m1")),b=try executable("deprotonated",fingerprint("m2"))
        #expect(throws:(any Error).self){
            _=try VivoConstantPHReplicaExchange(configuration:cfg,initialPhysicalStates:[physical(0),physical(1)],executableStates:[a,b])
        }
    }
}