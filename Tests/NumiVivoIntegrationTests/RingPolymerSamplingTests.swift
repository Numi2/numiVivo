import Foundation
import Testing
@testable import NumiVivoKit

@Suite(.serialized) struct RingPolymerSamplingTests {
    @Test func correlatedDiagnosticsDetectDisplacedChains() throws {
        var stable: [[Double]] = []
        for seed in 1...4 {
            var rng = VivoSplitMix64(state: UInt64(seed))
            stable.append((0..<256).map { _ in rng.normal() })
        }
        let mixed = try VivoCorrelatedSamplingAnalysis.calculate(chains: stable)
        #expect(mixed.splitRHat < 1.1)
        #expect(mixed.autocorrelationEffectiveSampleSize > 100)
        let displaced = try VivoCorrelatedSamplingAnalysis.calculate(chains: [
            [Double](repeating: -1,count: 64),[Double](repeating: 1,count: 64)])
        #expect(displaced.splitRHat == Double.greatestFiniteMagnitude)
        #expect(displaced.autocorrelationEffectiveSampleSize == 4)
        let trending = [Array(0..<64).map { Double($0) },Array(1...64).map { Double($0) }]
        let bounded = try VivoCorrelatedSamplingAnalysis.calculate(
            chains:trending,maximumAutocorrelationProducts:1)
        #expect(bounded.autocorrelationSequenceTruncated)
        #expect(bounded.autocorrelationEffectiveSampleSize == 4)
    }

    @Test func oneBeadLimitAndIsotopeNormalization() throws {
        let q: [[VivoVector3D]] = [[.init(0.2,0,0)]], rt = VivoAtomicUnits.gasConstantJPerMolK*0.3
        #expect(try VivoRingPolymerSampling.springEnergy(positions: q,massesDa: [1],temperatureK: 300) == 0)
        let a = try VivoRingPolymerSampling.reducedConfigurationalPotential(positions: q,massesDa: [1],temperatureK: 300,beadPotentialEnergiesKJPerMol: [3])
        #expect(abs(a-3/rt) < 1e-12)
        let beads = Array(repeating: q[0],count: 8), energies = Array(repeating: 3.0,count: 8)
        let light = try VivoRingPolymerSampling.reducedConfigurationalPotential(positions: beads,massesDa: [1],temperatureK: 300,beadPotentialEnergiesKJPerMol: energies)
        let heavy = try VivoRingPolymerSampling.reducedConfigurationalPotential(positions: beads,massesDa: [2],temperatureK: 300,beadPotentialEnergiesKJPerMol: energies)
        #expect(abs(heavy-light+12*log(2)) < 1e-12)
    }
    @Test func restartPreservesEverySweepAndSeedZero() async throws {
        let potential = try PrecisionSamplingFixtures.harmonic()
        let cfg = VivoRingPolymerConfiguration(temperatureK: 300,beadCount: 4,timeStepPS: 0.0005,integrationSteps: 3)
        let checkpoint = try VivoRingPolymerCheckpoint(definition: potential.definition,configuration: cfg,seed: 0,
            beadPositionsNM: Array(repeating: [.zero],count: 4))
        let whole = try await VivoRingPolymerSampling.run(potential: potential,checkpoint: checkpoint,sweeps: 12)
        let first = try await VivoRingPolymerSampling.run(potential: potential,checkpoint: checkpoint,sweeps: 5)
        let decoded = try VivoCanonicalJSON.decode(VivoRingPolymerCheckpoint.self,from: VivoCanonicalJSON.encode(first.end))
        let second = try await VivoRingPolymerSampling.run(potential: potential,checkpoint: decoded,sweeps: 7)
        try whole.validate(); try second.validate()
        #expect(whole.end == second.end)
        #expect(whole.observations == first.observations+second.observations)
        #expect(whole.observations.count == 12 && whole.forceEvaluations == 4*(1+12*3))
        #expect(whole.observations.allSatisfy { $0.logAcceptanceProbability <= 0 })
    }
    @Test func harmonicQuantumDistributionMatchesFiniteBeadSolution() async throws {
        let k = 200_000.0, temperature = 300.0, p = 8
        let potential = try PrecisionSamplingFixtures.harmonic(k: k)
        let cfg = VivoRingPolymerConfiguration(temperatureK: temperature,beadCount: p,timeStepPS: 0.0007,integrationSteps: 5)
        let initial = try VivoRingPolymerCheckpoint(definition: potential.definition,configuration: cfg,seed: 103,
            beadPositionsNM: Array(repeating: [.zero],count: p))
        let run = try await VivoRingPolymerSampling.run(potential: potential,checkpoint: initial,sweeps: 1400)
        let samples = Array(run.observations.dropFirst(200))
        let mean = samples.reduce(0.0) { $0+$1.beadPositionsNM.reduce(0) { $0+$1[0].squaredNorm }/Double(3*p) }/Double(samples.count)
        let rt = VivoAtomicUnits.gasConstantJPerMolK*temperature/1000
        let hbar = VivoAtomicUnits.planckJS*6.02214076e23*1e9/(2*Double.pi), omegaP = Double(p)*rt/hbar
        var expected = 0.0
        for j in 0..<p {
            let sine = sin(Double.pi*Double(j)/Double(p))
            expected += rt/(k+4*omegaP*omegaP*sine*sine)
        }
        #expect(abs(mean/expected-1) < 0.20)
        #expect(expected > 2*rt/k)
        #expect(samples.filter(\.accepted).count > samples.count/2)
    }
    @Test func failedAuthorityAndBudgetCannotMutateCheckpoint() async throws {
        let potential = try PrecisionSamplingFixtures.harmonic()
        let cfg = VivoRingPolymerConfiguration(temperatureK: 300,beadCount: 2,maximumForceEvaluations: 1)
        let checkpoint = try VivoRingPolymerCheckpoint(definition: potential.definition,configuration: cfg,seed: 5,beadPositionsNM: [[.zero],[.zero]])
        await #expect(throws: (any Error).self) { try await VivoRingPolymerSampling.run(potential: potential,checkpoint: checkpoint,sweeps: 1) }
        try checkpoint.validate()
        let failed = try VivoNuclearPotential(definition: potential.definition) { _ in throw VivoChemistryError.convergence("injected electronic failure") }
        let valid = try VivoRingPolymerCheckpoint(definition: potential.definition,configuration: .init(temperatureK: 300,beadCount: 1),seed: 5,beadPositionsNM: [[.zero]])
        await #expect(throws: (any Error).self) { try await VivoRingPolymerSampling.run(potential: failed,checkpoint: valid,sweeps: 1) }
        #expect(valid.sweep == 0)
    }
    @Test func periodicChartTranslatesWholeMoleculesAndAllBeadsTogether() throws {
        let cell = VivoPeriodicCell(a: .init(2,0,0),b: .init(0,2,0),c: .init(0,0,2))
        let d = try VivoNuclearPotentialDefinition(hamiltonianFingerprint: PrecisionSamplingFixtures.id("periodic-ring"),
            atomIndices: [0,1],particleIndices: [1,0],massesDa: [12,1],periodicCell: cell,periodicMoleculeGroups: [[0,1]])
        let q: [[VivoVector3D]] = [[.init(1.9,0.1,0.1),.init(2.05,0.1,0.1)], [.init(2.1,0.1,0.1),.init(2.25,0.1,0.1)]]
        let canonical = try d.canonicalBeads(q)
        #expect(abs((canonical[0][0].x+canonical[1][0].x)/2) < 1e-12)
        for b in q.indices { #expect(((q[b][1]-q[b][0])-(canonical[b][1]-canonical[b][0])).norm < 1e-12) }
        let first = try VivoRingPolymerSampling.springEnergy(positions: q,massesDa: d.massesDa,temperatureK: 300)
        let second = try VivoRingPolymerSampling.springEnergy(positions: canonical,massesDa: d.massesDa,temperatureK: 300)
        #expect(abs(first-second) < 1e-10*first)
        #expect(throws: (any Error).self) {
            try VivoNuclearPotentialDefinition(hamiltonianFingerprint: PrecisionSamplingFixtures.id("bad-chart"),atomIndices: [0,1],
                particleIndices: [0,1],massesDa: [12,1],periodicCell: cell,periodicMoleculeGroups: [[0]])
        }
    }

    @Test func retainedAcceptanceAndThinnedCheckpointReplayRejectTampering() async throws {
        let potential = try PrecisionSamplingFixtures.harmonic()
        let cfg = VivoRingPolymerConfiguration(temperatureK: 300,beadCount: 2,integrationSteps: 2)
        let cp = try VivoRingPolymerCheckpoint(definition: potential.definition,configuration: cfg,seed: 24,
            beadPositionsNM: [[.zero],[.zero]])
        let run = try await VivoRingPolymerSampling.run(potential: potential,checkpoint: cp,sweeps: 6,recordEvery: 2)
        try run.validate(recordEvery: 2)
        var json = try #require(JSONSerialization.jsonObject(with: VivoCanonicalJSON.encode(run)) as? [String:Any])
        var observations = try #require(json["observations"] as? [[String:Any]])
        observations[0]["accepted"] = !run.observations[0].accepted; json["observations"] = observations
        let altered = try VivoCanonicalJSON.decode(VivoRingPolymerRun.self,from: JSONSerialization.data(withJSONObject: json))
        #expect(throws: (any Error).self) { try altered.validate(recordEvery: 2) }
    }

    @Test func independentChainsProduceExplicitBeadConvergenceEvidence() async throws {
        let potential = try PrecisionSamplingFixtures.harmonic(k: 25_000)
        var evidence: [VivoRingPolymerConvergenceChain] = []
        for beadCount in [1,2,4] { for replica in 0..<2 {
            let cfg = VivoRingPolymerConfiguration(temperatureK:300,beadCount:beadCount,
                timeStepPS:0.0007,integrationSteps:4)
            let checkpoint = try VivoRingPolymerCheckpoint(definition:potential.definition,configuration:cfg,
                seed:UInt64(100*beadCount+replica),beadPositionsNM:Array(repeating:[.zero],count:beadCount))
            let run = try await VivoRingPolymerSampling.run(potential:potential,checkpoint:checkpoint,sweeps:160)
            evidence.append(.init(identifier:"P\(beadCount)-R\(replica)",run:run,discardedObservations:32))
        } }
        let configuration = VivoRingPolymerConvergenceConfiguration(
            observable:.primitiveTotalEnergyKJPerMol,absoluteTolerance:100,relativeTolerance:0,
            minimumEffectiveSamples:4,maximumSplitRHat:2)
        let request = VivoRingPolymerConvergenceRequest(campaignIdentifier:"harmonic-control",
            selectionProtocol:"Three bead levels and two seeds were fixed before inspecting the retained energy series.",
            chains:evidence,configuration:configuration)
        let result = try VivoRingPolymerConvergence.assess(request)
        try VivoRingPolymerConvergence.validate(result,request:request)
        #expect(result.passed)
        #expect(result.beadDiagnostics.map(\.beadCount) == [1,2,4])
        #expect(result.comparisons.count == 2)
        #expect(result.beadDiagnostics.allSatisfy { $0.statistics.retainedSamplesPerChain == 128 })
        let strict = try VivoRingPolymerConvergence.assess(campaignIdentifier:"harmonic-control-strict",
            selectionProtocol:"The deliberately impossible tolerance is a negative qualification fixture.",
            chains:evidence,configuration:.init(observable:.primitiveTotalEnergyKJPerMol,
                absoluteTolerance:1e-12,relativeTolerance:0,minimumEffectiveSamples:4,maximumSplitRHat:2))
        #expect(!strict.passed)
        var reused=evidence
        reused[1] = .init(identifier:"forged",run:evidence[0].run,discardedObservations:32)
        #expect(throws:(any Error).self) {
            try VivoRingPolymerConvergence.assess(campaignIdentifier:"reused-chain",
                selectionProtocol:"Duplicate start evidence must not be called independent.",
                chains:reused,configuration:configuration)
        }
    }

}
