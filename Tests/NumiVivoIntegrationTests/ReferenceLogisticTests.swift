import Foundation
import Testing
@testable import NumiVivoKit

@Suite struct ReferenceLogisticTests {
    @Test func weightedMultinomialGradientMatchesFiniteDifferencesIncludingIntercepts() throws {
        let x = [[-2.0,1],[1,3],[0.5,-1],[2,0]], labels = [0,1,2,2], weights = [4.0/3,4.0/3,2.0/3]
        let p = [0.2,-0.3,0.4,-0.1,0.5,-0.2,0.3,-0.1,0.1]
        let result = try VivoReferenceLogistic.objective(p,x: x,labels: labels,classes: 3,weights: weights,penalty: 1)
        for j in p.indices {
            var plus=p,minus=p;plus[j]+=1e-5;minus[j]-=1e-5
            let difference = try (VivoReferenceLogistic.objective(plus,x: x,labels: labels,classes: 3,weights: weights,penalty: 1).value - VivoReferenceLogistic.objective(minus,x: x,labels: labels,classes: 3,weights: weights,penalty: 1).value)/2e-5
            #expect(abs(difference-result.gradient[j])<1e-9)
        }
    }
    @Test func balancedConstantDataHaveEqualProbabilitiesAndUnpenalizedIntercepts() throws {
        let model=try VivoReferenceLogistic.fit(scores: Array(repeating: [2.0,5],count: 10),labels: [0,1,1,2,2,2,2,2,2,2],classes: 3,options: .init())
        #expect(model.classCounts == [1,2,7]);#expect(model.scales == [1,1])
        let p=try VivoReferenceLogistic.probabilities([2,5],model: model)
        for v in p { #expect(abs(v-1.0/3)<1e-14) }
        #expect(model.gradientMaximum<=1e-7)
    }
    @Test func solverConvergesAndScalingIsTrainingOnly() throws {
        let x=[[-3.0,-2],[-2,-3],[3,2],[2,3],[0,5],[1,6]],labels=[0,0,1,1,2,2]
        let model=try VivoReferenceLogistic.fit(scores: x,labels: labels,classes: 3,options: .init())
        #expect(model.gradientMaximum<=1e-7 && model.iterations>0)
        #expect(zip(model.objectives.dropFirst(),model.objectives).allSatisfy { $0 <= $1 })
        for i in x.indices {
            let p=try VivoReferenceLogistic.probabilities(x[i],model: model)
            #expect(abs(p.reduce(0,+)-1)<1e-14);#expect(p.firstIndex(of:p.max()!)==labels[i])
        }
        let shifted=try VivoReferenceLogistic.fit(scores: x.map { [$0[0]*7+11,$0[1]*3-9] },labels: labels,classes: 3,options: .init())
        let a=try VivoReferenceLogistic.probabilities([0.1,1.2],model: model)
        let b=try VivoReferenceLogistic.probabilities([11.7,-5.4],model: shifted)
        for j in a.indices { #expect(abs(a[j]-b[j])<1e-10) }
    }
    @Test func invalidInputsAndUnconvergedOrOverBudgetFitsReject() throws {
        let x=[[-2.0,1],[1,3],[0.5,-1],[2,0]],labels=[0,1,2,2]
        var limited=VivoReferenceLogisticOptions();limited.maximumIterations=1
        #expect(throws:(any Error).self) { try VivoReferenceLogistic.fit(scores:x,labels:labels,classes:3,options:limited) }
        limited.maximumWork=1
        #expect(throws:(any Error).self) { try VivoReferenceLogistic.fit(scores:x,labels:labels,classes:3,options:limited) }
        #expect(throws:(any Error).self) { try VivoReferenceLogistic.fit(scores:x,labels:[0,0,1,1],classes:3,options:.init()) }
        #expect(throws:(any Error).self) { try VivoReferenceLogistic.fit(scores:[[.nan],[1]],labels:[0,1],classes:2,options:.init()) }
        #expect(throws:(any Error).self) { try VivoCanonicalJSON.decode(VivoReferenceLogisticOptions.self,from: Data("{\"unknown\":true}".utf8)) }
    }
}
