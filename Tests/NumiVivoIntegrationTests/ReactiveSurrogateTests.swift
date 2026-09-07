import Foundation
import Testing
@testable import NumiVivoKit

@Suite(.serialized) struct ReactiveSurrogateTests {
    @Test func groupedEnergyForceFitAndNoHeldOutLeakage() async throws {
        let (model,a,b,labels) = try await PrecisionSamplingFixtures.trainedPair()
        #expect(model.payload.qualification.passed)
        #expect(model.payload.qualification.heldOutGroups == ["held"])
        #expect(model.payload.qualification.trainingGroupCount == 21)
        let again = try VivoReactiveDeltaSurrogate.train(authority: a.definition,baseline: b.definition,labels: labels,
            heldOutGroups: ["held"],configuration: model.payload.configuration)
        #expect(again == model)
        #expect(throws: (any Error).self) {
            try VivoReactiveDeltaSurrogate.train(authority: a.definition,baseline: b.definition,labels: labels+[labels[0]],
                heldOutGroups: ["held"],configuration: model.payload.configuration)
        }
        #expect(throws: (any Error).self) {
            try VivoReactiveDeltaSurrogate.train(authority: a.definition,baseline: b.definition,labels: labels,
                heldOutGroups: labels.map(\.sourceGroup),configuration: model.payload.configuration)
        }
    }
    @Test func learnedForcesAreEnergyDerivativesAndConserveForceTorque() async throws {
        let (model,_,_,_) = try await PrecisionSamplingFixtures.trainedPair()
        let q: [VivoVector3D] = [.init(-0.4,0.1,0.1),.init(0.5,0.3,-0.2)]
        let result = try VivoReactiveDeltaSurrogate.predict(model,positionsNM: q), h = 1e-6
        for atom in 0..<2 { for axis in 0..<3 {
            func moved(_ delta: Double) -> [VivoVector3D] {
                var out = q
                if axis == 0 { out[atom].x += delta }
                else if axis == 1 { out[atom].y += delta }
                else { out[atom].z += delta }
                return out
            }
            let plus = try VivoReactiveDeltaSurrogate.predict(model,positionsNM: moved(h)).deltaEnergyKJPerMol
            let minus = try VivoReactiveDeltaSurrogate.predict(model,positionsNM: moved(-h)).deltaEnergyKJPerMol
            let f = result.deltaForcesKJPerMolNM[atom]
            let analytic = axis == 0 ? f.x : (axis == 1 ? f.y : f.z)
            #expect(abs(analytic+(plus-minus)/(2*h)) < 1e-5)
        } }
        #expect(result.deltaForcesKJPerMolNM.reduce(.zero,+).norm < 1e-12)
        #expect(zip(q,result.deltaForcesKJPerMolNM).reduce(VivoVector3D.zero) { $0+$1.0.cross($1.1) }.norm < 1e-12)
        let translated = try VivoReactiveDeltaSurrogate.predict(model,positionsNM: q.map { $0 + .init(1,2,3) })
        #expect(abs(translated.deltaEnergyKJPerMol-result.deltaEnergyKJPerMol) < 1e-10)
    }
    @Test func outOfSupportIsNotAnAuthoritativePrediction() async throws {
        let (model,_,_,_) = try await PrecisionSamplingFixtures.trainedPair()
        let p = try VivoReactiveDeltaSurrogate.predict(model,positionsNM: PrecisionSamplingFixtures.pairPositions(4))
        #expect(!p.eligible && !p.reasons.isEmpty)
        var object = try #require(JSONSerialization.jsonObject(with: VivoCanonicalJSON.encode(model)) as? [String:Any])
        var payload = try #require(object["payload"] as? [String:Any]); payload["centersNM"] = [[4.0]]; object["payload"] = payload
        let altered = try VivoCanonicalJSON.decode(VivoReactiveSurrogateModel.self,from: JSONSerialization.data(withJSONObject: object))
        #expect(throws: (any Error).self) { try altered.validate() }
    }
    @Test func correctedSamplingUsesAuthorityEndpointsAndRestartsExactly() async throws {
        let (model,a,b,_) = try await PrecisionSamplingFixtures.trainedPair()
        let cfg = VivoReactiveSamplingConfiguration(temperatureK: 300,timeStepPS: 0.001,integrationSteps: 8)
        let cp = try VivoReactiveSamplingCheckpoint(model: model,configuration: cfg,chainIdentifier: "new-independent-chain",seed: 0,
            positionsNM: PrecisionSamplingFixtures.pairPositions(1))
        let whole = try await VivoReactiveSurrogateSampling.run(model: model,authority: a,baseline: b,checkpoint: cp,sweeps: 12)
        let first = try await VivoReactiveSurrogateSampling.run(model: model,authority: a,baseline: b,checkpoint: cp,sweeps: 5)
        let second = try await VivoReactiveSurrogateSampling.run(model: model,authority: a,baseline: b,checkpoint: first.end,sweeps: 7)
        try whole.validate(model: model); try second.validate(model: model)
        #expect(whole.end == second.end && whole.observations == first.observations+second.observations)
        #expect(whole.authorityEvaluations == 13 && whole.baselineEvaluations > whole.authorityEvaluations)
        #expect(whole.authorityLabels.count == 12 && whole.observations.count == 12)
        for observation in whole.observations {
            #expect(try await a.checked(observation.positionsNM).energyKJPerMol == observation.authorityEnergyKJPerMol)
        }
        #expect(Set(whole.authorityLabels.map(\.sourceGroup)) == ["new-independent-chain"])
    }
    @Test func domainRejectionRetainsStateAndEmitsAcquisition() async throws {
        let (model,a,b,_) = try await PrecisionSamplingFixtures.trainedPair()
        let q = PrecisionSamplingFixtures.pairPositions(4)
        let cp = try VivoReactiveSamplingCheckpoint(model: model,configuration: .init(temperatureK: 300,authorityRefreshProbability: 1e-9),chainIdentifier: "domain",seed: 2,positionsNM: q)
        let run = try await VivoReactiveSurrogateSampling.run(model: model,authority: a,baseline: b,checkpoint: cp,sweeps: 3)
        #expect(run.observations.allSatisfy { !$0.accepted && $0.positionsNM == q })
        #expect(run.acquisitions.count == 3 && run.authorityEvaluations == 1)
        #expect(run.authorityLabels.isEmpty)
    }
    @Test func authorityRefreshCanMoveBeyondSurrogateSupport() async throws {
        let (model,a,b,_) = try await PrecisionSamplingFixtures.trainedPair()
        let q = PrecisionSamplingFixtures.pairPositions(4)
        let cfg = VivoReactiveSamplingConfiguration(temperatureK: 300,authorityRefreshProbability: 1,authorityRefreshLengthNM: 0.01)
        let cp = try VivoReactiveSamplingCheckpoint(model: model,configuration: cfg,chainIdentifier: "refresh",seed: 91,positionsNM: q)
        let run = try await VivoReactiveSurrogateSampling.run(model: model,authority: a,baseline: b,checkpoint: cp,sweeps: 8)
        try run.validate(model: model)
        #expect(run.observations.allSatisfy { $0.proposalKind == .authorityRandomWalk })
        #expect(run.authorityEvaluations == 9 && run.authorityLabels.count == 8)
        #expect(run.observations.contains { $0.accepted && $0.positionsNM != q })
    }

    @Test func independentGroupedCoverageCannotHideUnsupportedRegions() async throws {
        let (model,authority,baseline,training) = try await PrecisionSamplingFixtures.trainedPair()
        var labels: [VivoReactiveTrainingLabel] = []
        let firstDistances=(0..<8).map { 0.81+0.02*Double($0) }
        let secondDistances=(0..<8).map { 1.03+0.02*Double($0) }
        let coverageGroups:[(String,[Double])] = [
            ("coverage-a",firstDistances),("coverage-b",secondDistances)]
        for (group,distances) in coverageGroups {
            for (index,distance) in distances.enumerated() {
                let positions=PrecisionSamplingFixtures.pairPositions(distance)
                labels.append(.init(identifier:"\(group)-\(index)",sourceGroup:group,
                    authority:try await authority.checked(positions),baseline:try await baseline.checked(positions)))
            }
        }
        let source = try PrecisionSamplingFixtures.id("independent-coverage-source")
        let request = VivoReactiveSurrogateCoverageRequest(campaignIdentifier:"pair-coverage",
            model:model,labels:labels,evidenceSourceFingerprint:source,
            independenceDeclaration:"Coverage trajectories were fixed before the model and were not used for fit or threshold selection.",
            domainDescription:"Two independent pair-distance source groups spanning 0.81 through 1.17 nm.")
        let result = try VivoReactiveSurrogateCoverage.assess(request)
        try VivoReactiveSurrogateCoverage.validate(result,request:request)
        #expect(result.passed)
        #expect(result.groups.count == 2 && result.overallEligibleFraction == 1)
        #expect(result.groups.allSatisfy { $0.passed && $0.eligibleCount == 8 })

        var unsupported:[VivoReactiveTrainingLabel]=[]
        for group in ["far-a","far-b"] { for index in 0..<8 {
            let positions=PrecisionSamplingFixtures.pairPositions(2+0.02*Double(index)+(group=="far-b" ? 0.3:0))
            unsupported.append(.init(identifier:"\(group)-\(index)",sourceGroup:group,
                authority:try await authority.checked(positions),baseline:try await baseline.checked(positions)))
        } }
        let failed = try VivoReactiveSurrogateCoverage.assess(campaignIdentifier:"unsupported-control",
            model:model,labels:labels+unsupported,evidenceSourceFingerprint:source,
            independenceDeclaration:"Deliberate external-domain negative control.",
            domainDescription:"Distances outside the frozen descriptor support.",
            configuration:.init(minimumOverallEligibleFraction:0.5,minimumPerGroupEligibleFraction:0.5))
        #expect(!failed.passed)
        #expect(failed.overallEligibleFraction >= 0.5)
        #expect(failed.groups.contains { $0.passed })
        #expect(failed.groups.contains { !$0.passed && $0.outsideDescriptorSupportCount == 8 })
        #expect(throws:(any Error).self) {
            try VivoReactiveSurrogateCoverage.assess(campaignIdentifier:"source-reuse",model:model,
                labels:labels,evidenceSourceFingerprint:model.payload.trainingDataFingerprint,
                independenceDeclaration:"A reused development-source identity must be rejected.",
                domainDescription:"Invalid source provenance.")
        }
        #expect(throws:(any Error).self) {
            try VivoReactiveSurrogateCoverage.assess(campaignIdentifier:"training-reuse",model:model,
                labels:training,evidenceSourceFingerprint:source,
                independenceDeclaration:"This false declaration must not bypass exact training-data identity.",
                domainDescription:"Invalid reuse of fitting and held-out labels.")
        }
    }

}
