import Foundation
import Testing
@testable import NumiVivoKit

@Suite struct SingleCellMetalNegativeBinomialTests {
    @Test(.enabled(if: ProcessInfo.processInfo.environment["NUMIVIVO_TEST_METAL"] == "1"))
    func batchedVariableObjectiveMatchesFP64CPUOracle() throws {
        let counts: [UInt64] = [0, 3, 12, 7, 1, 19, 2, 0, 31, 5]
        let means = [0.25, 1.2, 4.5, 8.0, 2.1, 17.0, 0.75, 3.0, 29.0, 6.5]
        let dispersion = 0.2
        let metal = try VivoMetalNegativeBinomialLikelihood(counts: counts, dispersion: dispersion)
        let value = try metal.evaluate(means: means)
        let expected = zip(counts, means).reduce(0.0) { total, pair in
            let y = Double(pair.0), mu = pair.1
            return total + y * log(mu) - (y + 1 / dispersion) * log1p(dispersion * mu)
        }
        #expect(abs(value - expected) <= max(5e-5, abs(expected) * 5e-6))
    }

    @Test(.enabled(if: ProcessInfo.processInfo.environment["NUMIVIVO_TEST_METAL"] == "1"))
    func optInFitUsesMetalObjectiveAndRetainsExactReport() throws {
        let counts: [UInt64] = [11, 20, 12, 55, 43, 61]
        let design = [[1.0, 0], [1, 0], [1, 0], [1, 1], [1, 1], [1, 1]]
        let offsets = [0.0, 0.2, 0.1, 0.3, 0.1, 0.2]
        let cpu = try VivoOmicsNegativeBinomial.fit(counts: counts, design: design, offsets: offsets,
            contrast: [0, 1], dispersion: 0.1)
        let metal = try VivoOmicsNegativeBinomial.fit(counts: counts, design: design, offsets: offsets,
            contrast: [0, 1], dispersion: 0.1, backend: .metalFP32)
        #expect(metal.backend == .metalFP32)
        #expect(metal.converged)
        #expect(abs((metal.effect ?? .nan) - (cpu.effect ?? .nan)) < 1e-4)
        #expect(abs(metal.logLikelihood - cpu.logLikelihood) < 1e-6)
        let encoded = try VivoCanonicalJSON.encode(metal)
        #expect(String(decoding: encoded, as: UTF8.self).contains(#""backend":"metalFP32""#))
    }

    @Test(.enabled(if: ProcessInfo.processInfo.environment["NUMIVIVO_TEST_METAL"] == "1"))
    func objectiveRejectsCountsOutsideFP32Contract() throws {
        #expect(throws: (any Error).self) {
            _ = try VivoMetalNegativeBinomialLikelihood(counts: [16_777_217], dispersion: 0.2)
        }
    }

    @Test func backendIsAbsentFromDefaultFitEncoding() throws {
        let counts: [UInt64] = [2, 5, 10, 12, 20, 23]
        let design = counts.map { _ in [1.0] }
        let fit = try VivoOmicsNegativeBinomial.fit(counts: counts, design: design,
            offsets: Array(repeating: 0, count: counts.count), contrast: [1], dispersion: 0.2)
        #expect(fit.backend == nil)
        let encoded = try VivoCanonicalJSON.encode(fit)
        #expect(!String(decoding: encoded, as: UTF8.self).contains("backend"))
    }
}
