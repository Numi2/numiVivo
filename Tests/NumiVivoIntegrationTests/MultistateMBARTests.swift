import Foundation
import Testing
@testable import NumiVivoKit

@Suite(.serialized) struct MultistateMBARTests {
    @Test func constantEnergyOffsetsRecoverExactFreeEnergiesAndBootstrapZero() throws {
        var samples:[VivoMBARReducedPotentialSample]=[]
        let offsets=[0.0,1.5,-0.7]
        for origin in 0..<3 {
            for n in 0..<30 {
                let base=0.03*Double(n)+0.11*Double(origin)
                samples.append(.init(identifier:"\(origin)-\(n)",originStateIndex:origin,
                    reducedPotentials:offsets.map{base+$0}))
            }
        }
        let result=try VivoMultistateMBAR.solve(samples:samples,stateCount:3,
            configuration:.init(minimumTargetEffectiveSamples:20,bootstrapReplicates:32,bootstrapSeed:7))
        #expect(result.converged && result.overlapGraphConnected)
        #expect(abs(result.estimates[0].relativeFreeEnergy)<1e-12)
        #expect(abs(result.estimates[1].relativeFreeEnergy-1.5)<1e-10)
        #expect(abs(result.estimates[2].relativeFreeEnergy+0.7)<1e-10)
        #expect((result.estimates[1].bootstrapStandardDeviation ?? 1)<1e-10)
        #expect(result.pairwiseOverlap.allSatisfy{$0.symmetricAcceptanceOverlap>0.999999})
        #expect(result.pairwiseOverlap.allSatisfy{$0.reweightingESSFractionAToB>0.999999 && $0.reweightingESSFractionBToA>0.999999})
    }

    @Test func mutuallyInaccessibleWellsAreDisconnectedDespitePerfectWithinDirectionESS() throws {
        var samples:[VivoMBARReducedPotentialSample]=[]
        for n in 0..<40 {
            samples.append(.init(identifier:"a-\(n)",originStateIndex:0,reducedPotentials:[0,100]))
            samples.append(.init(identifier:"b-\(n)",originStateIndex:1,reducedPotentials:[100,0]))
        }
        let result=try VivoMultistateMBAR.solve(samples:samples,stateCount:2,
            configuration:.init(minimumBidirectionalAcceptanceOverlap:0.01,
                                minimumTargetEffectiveSamples:2,bootstrapReplicates:0))
        #expect(!result.converged && !result.overlapGraphConnected)
        let edge=try #require(result.pairwiseOverlap.first)
        #expect(edge.reweightingESSFractionAToB>0.999999 && edge.reweightingESSFractionBToA>0.999999)
        #expect(edge.symmetricAcceptanceOverlap<1e-20)
        #expect(result.issues.contains{$0.contains("disconnected")})
    }

    @Test func intermediateStateConnectsOtherwiseSeparatedEndpoints() throws {
        var samples:[VivoMBARReducedPotentialSample]=[]
        for n in 0..<50 {
            let jitter=0.01*Double(n%7)
            samples.append(.init(identifier:"a-\(n)",originStateIndex:0,reducedPotentials:[jitter,1+jitter,100+jitter]))
            samples.append(.init(identifier:"m-\(n)",originStateIndex:1,reducedPotentials:[1+jitter,jitter,1+jitter]))
            samples.append(.init(identifier:"b-\(n)",originStateIndex:2,reducedPotentials:[100+jitter,1+jitter,jitter]))
        }
        let result=try VivoMultistateMBAR.solve(samples:samples,stateCount:3,
            configuration:.init(minimumBidirectionalAcceptanceOverlap:0.1,
                                minimumTargetEffectiveSamples:20,bootstrapReplicates:0))
        #expect(result.converged && result.overlapGraphConnected)
        let direct=try #require(result.pairwiseOverlap.first{$0.stateA==0 && $0.stateB==2})
        #expect(!direct.connected)
        #expect(result.pairwiseOverlap.first{$0.stateA==0 && $0.stateB==1}?.connected==true)
        #expect(result.pairwiseOverlap.first{$0.stateA==1 && $0.stateB==2}?.connected==true)
        #expect(abs(result.estimates[0].relativeFreeEnergy)<1e-10)
        #expect(abs(result.estimates[2].relativeFreeEnergy)<1e-10)
    }

    @Test func deterministicBootstrapAndCapacityFailures() throws {
        var samples:[VivoMBARReducedPotentialSample]=[]
        for origin in 0..<2 { for n in 0..<25 {
            let x=Double(n-12)/10
            let u0=0.5*x*x
            let u1=0.5*(x-0.3)*(x-0.3)+0.2
            samples.append(.init(identifier:"\(origin)-\(n)",originStateIndex:origin,reducedPotentials:[u0,u1]))
        } }
        let cfg=VivoMBARConfiguration(minimumBidirectionalAcceptanceOverlap:0.01,
            minimumTargetEffectiveSamples:10,bootstrapReplicates:24,bootstrapSeed:1234,maximumWorkElements:10_000)
        let a=try VivoMultistateMBAR.solve(samples:samples,stateCount:2,configuration:cfg)
        let b=try VivoMultistateMBAR.solve(samples:samples,stateCount:2,configuration:cfg)
        #expect(a==b)
        #expect((a.estimates[1].bootstrapStandardDeviation ?? 0)>0)
        #expect(throws:(any Error).self) {
            _=try VivoMultistateMBAR.solve(samples:samples,stateCount:2,
                configuration:.init(minimumTargetEffectiveSamples:2,bootstrapReplicates:24,maximumWorkElements:100))
        }
        let missing=samples.filter{$0.originStateIndex==0}
        #expect(throws:(any Error).self) { _=try VivoMultistateMBAR.solve(samples:missing,stateCount:2) }
    }
}