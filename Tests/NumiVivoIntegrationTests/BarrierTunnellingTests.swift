import Foundation
import Testing
@testable import NumiVivoKit

@Suite(.serialized) struct BarrierTunnellingTests {
    private func request(_ model: VivoBarrierTunnellingModel = .asymmetricEckart) throws -> VivoBarrierTunnellingRequest {
        .init(model: model, temperatureK: 300, imaginaryWavenumberPerCM: 1000,
              forwardBarrierKJPerMol: 40, reverseBarrierKJPerMol: 60,
              hamiltonianFingerprint: try PrecisionSamplingFixtures.id("barrier-H"),
              stationaryPointEvidence: try PrecisionSamplingFixtures.id("barrier-saddle"), energyConvention: "electronic + specified transverse ZPE")
    }
    @Test func wignerSmallCorrectionAndRegimeRejection() throws {
        var r = try request(.wigner); r.imaginaryWavenumberPerCM = 100
        let result = try VivoBarrierTunnelling.calculate(r)
        let x = VivoAtomicUnits.planckJS*299_792_458*100*100/(VivoAtomicUnits.boltzmannJPerK*r.temperatureK)
        let factor = try #require(result.factor)
        #expect(abs(factor - (1+x*x/24)) < 1e-12)
        try VivoBarrierTunnelling.validate(result)
        r.imaginaryWavenumberPerCM = 1000
        #expect(throws: (any Error).self) { try VivoBarrierTunnelling.calculate(r) }
    }
    @Test func transmissionRespectsBothAsymptotesAndHighEnergyLimit() throws {
        var r = try request(); r.reverseBarrierKJPerMol = 20
        for e in [-1.0,0,10,20] { #expect(try VivoBarrierTunnelling.transmissionProbability(energyKJPerMol: e,request: r) == 0) }
        let high = try VivoBarrierTunnelling.transmissionProbability(energyKJPerMol: 1000,request: r)
        #expect(high > 0.999999 && high <= 1)
        let middle = try VivoBarrierTunnelling.transmissionProbability(energyKJPerMol: 35,request: r)
        #expect(middle > 0 && middle < 1)
    }
    @Test func forwardReverseDetailedBalanceAndQuadratureConvergence() throws {
        let r = try request(), a = try VivoBarrierTunnelling.calculate(r)
        var reverse = r; reverse.forwardBarrierKJPerMol = r.reverseBarrierKJPerMol; reverse.reverseBarrierKJPerMol = r.forwardBarrierKJPerMol
        let b = try VivoBarrierTunnelling.calculate(reverse)
        #expect(abs(a.logFactor-b.logFactor) < 1e-12)
        var fine = r; fine.relativeTolerance = 1e-10
        let c = try VivoBarrierTunnelling.calculate(fine)
        #expect(abs(c.logFactor-a.logFactor) < 1e-7)
        #expect(a.evaluations > 0 && a.quadratureRelativeChange <= r.relativeTolerance)
        #expect(try #require(a.logAbsoluteTailBound) - a.logFactor < log(r.relativeTolerance))
        // Independently integrated in Tools/PrecisionSampling/reference_checks.py.
        #expect(abs(a.logFactor - 1.111947262990766) < 1e-7)
    }
    @Test func shallowBarrierAnalyticContinuationAndBudgetErrors() throws {
        var r = try request(); r.forwardBarrierKJPerMol = 0.1; r.reverseBarrierKJPerMol = 0.2
        for e in [0.001,0.05,1,10] {
            let p = try VivoBarrierTunnelling.transmissionProbability(energyKJPerMol: e,request: r)
            #expect(p.isFinite && p > 0 && p <= 1)
        }
        r = try request(); r.maximumEvaluations = 256
        #expect(throws: (any Error).self) { try VivoBarrierTunnelling.calculate(r) }
        r = try request(); r.forwardBarrierKJPerMol = -1
        #expect(throws: (any Error).self) { try VivoBarrierTunnelling.calculate(r) }
    }
    @Test func rateCompositionIsSeparateBoundAndCannotDoubleCount() throws {
        let r = try request(), result = try VivoBarrierTunnelling.calculate(r), baseline = try PrecisionSamplingFixtures.id("classical-rate")
        let composed = try VivoNuclearRateComposition(baselineFingerprint: baseline, hamiltonianFingerprint: r.hamiltonianFingerprint,
            temperatureK: r.temperatureK, baselineLogRatePerSecond: -12, baselineAlreadyContainsTunnelling: false,
            correction: result, separabilityAssumption: "One-dimensional tunnelling factor separable from classical recrossing")
        #expect(composed.correctedLogRatePerSecond == -12+result.logFactor)
        #expect(throws: (any Error).self) {
            try VivoNuclearRateComposition(baselineFingerprint: baseline, hamiltonianFingerprint: r.hamiltonianFingerprint,
                temperatureK: r.temperatureK, baselineLogRatePerSecond: -12, baselineAlreadyContainsTunnelling: true,
                correction: result, separabilityAssumption: "declared")
        }
        #expect(throws: (any Error).self) {
            try VivoNuclearRateComposition(baselineFingerprint: baseline, hamiltonianFingerprint: baseline,
                temperatureK: r.temperatureK, baselineLogRatePerSecond: -12, baselineAlreadyContainsTunnelling: false,
                correction: result, separabilityAssumption: "declared")
        }
    }
}
