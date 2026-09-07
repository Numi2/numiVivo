import Foundation
@preconcurrency import Metal
import Testing
@testable import NumiVivoKit

/// Native execution evidence for the dedicated UInt32/FP32 hybrid backend.
/// Missing Apple hardware fails these checks; it is not a successful skip.
@Suite(.serialized) struct HybridRuntimeTests {
    private let seed: UInt64 = 0x4859425249445631

    private func device() throws -> MTLDevice {
        let result = try VivoMetalDeviceSelector.productionDevice()
        #expect(result.hasUnifiedMemory)
        return result
    }

    private func rejects(_ operation: () async throws -> Void) async {
        do { try await operation(); Issue.record("Expected hybrid operation rejection") }
        catch { }
    }

    private func record<T: Encodable>(_ value: T, _ name: String) throws {
        guard let directory = ProcessInfo.processInfo.environment["NUMIVIVO_TEST_ARTIFACTS"] else { return }
        let url = URL(fileURLWithPath: directory, isDirectory: true)
        try FileManager.default.createDirectory(at: url, withIntermediateDirectories: true)
        try VivoCanonicalJSON.encode(value).write(to: url.appendingPathComponent(name + ".json"), options: .atomic)
    }

    private func configuration(lanes: UInt32, dt: Float = 0.05,
                               chunk: UInt32 = 64, dispatches: UInt32 = 64) -> VivoHybridRuntimeConfiguration {
        .init(laneCount: lanes, timeStep: dt, minimumTimeStep: 1e-7, maximumTimeStep: 2,
              seed: seed, exactEventsPerDispatch: chunk, maximumExactDispatches: dispatches,
              maximumPublications: 16)
    }

    private func mixedModel() throws -> VivoExactSSAModel {
        try .init(species: ["exact-A", "exact-B", "tau-A", "tau-B", "rk-A", "rk-B", "constant"],
                  reactions: [
                    .init(id: "exact", law: .firstOrder, reactantA: 0, rate: 0.5,
                          changes: [.init(speciesIndex: 0, delta: -1), .init(speciesIndex: 1, delta: 1)]),
                    .init(id: "tau", law: .firstOrder, reactantA: 2, rate: 0.25,
                          changes: [.init(speciesIndex: 2, delta: -1), .init(speciesIndex: 3, delta: 1)]),
                    .init(id: "rk", law: .firstOrder, reactantA: 4, rate: 0.4,
                          changes: [.init(speciesIndex: 4, delta: -1), .init(speciesIndex: 5, delta: 1)])
                  ])
    }

    private func mixedPlan(_ model: VivoExactSSAModel) throws -> VivoHybridStochasticPlan {
        try VivoHybridExplicitPlanBuilder.make(model: model,
            reactionAuthorities: ["tau": .tauLeap, "rk": .deterministicRK2], maximumStep: 0.05)
    }

    private func mixedInitial(lanes: UInt32) -> [UInt32] {
        [UInt32(40), 0, 10_000, 0, 100, 0, .max].flatMap {
            Array(repeating: $0, count: Int(lanes))
        }
    }

    @Test func mixedAuthoritiesPreserveConservationAndPublicationTypes() async throws {
        // An odd lane count also exercises partial final SIMD groups.
        let lanes: UInt32 = 37, model = try mixedModel(), plan = try mixedPlan(model)
        let config = configuration(lanes: lanes)
        let runtime = try await VivoHybridReactionRuntime.make(model: model, plan: plan,
            configuration: config, initialCounts: mixedInitial(lanes: lanes), device: device())
        #expect(runtime.execution.authorities == [.exactSSA, .exactSSA, .tauLeap, .tauLeap,
                                                  .deterministicRK2, .deterministicRK2, .constant])
        let requests: [VivoHybridPublicationRequest] = [
            .init(speciesIndex: 0, laneIndex: 0), .init(speciesIndex: 2, laneIndex: lanes - 1),
            .init(speciesIndex: 4, laneIndex: 0), .init(speciesIndex: 6, laneIndex: lanes - 1)
        ]
        var publications: [VivoHybridPublication] = []
        for _ in 0..<5 {
            let accepted = try await runtime.step(publications: requests)
            try #require(accepted.disposition == .committed)
            #expect(accepted.status.complete)
            publications = accepted.publications
        }
        let state = try await runtime.snapshot()
        #expect(state.stepIndex == 5)
        #expect(state.timeSeconds == 5 * Double(config.timeStep))
        let x = Double(Float(0.4)) * Double(config.timeStep)
        let expectedRK = 100 * pow(1 - x + 0.5 * x * x, 5)
        var exactChanges: UInt64 = 0, tauChanges: UInt64 = 0
        for lane in 0..<Int(lanes) {
            #expect(UInt64(state.counts[lane]) + UInt64(state.counts[Int(lanes) + lane]) == 40)
            #expect(UInt64(state.counts[2 * Int(lanes) + lane]) + UInt64(state.counts[3 * Int(lanes) + lane]) == 10_000)
            let a = Double(state.continuousValues[4 * Int(lanes) + lane])
            let b = Double(state.continuousValues[5 * Int(lanes) + lane])
            #expect(abs(a - expectedRK) < 3e-5)
            #expect(abs(a + b - 100) < 3e-5)
            #expect(state.counts[6 * Int(lanes) + lane] == UInt32.max)
            exactChanges += UInt64(state.counts[Int(lanes) + lane])
            tauChanges += UInt64(state.counts[3 * Int(lanes) + lane])
            for species in 0..<7 {
                let index = species * Int(lanes) + lane
                if species == 4 || species == 5 { #expect(state.counts[index] == 0) }
                else { #expect(state.continuousValues[index] == 0) }
            }
        }
        #expect(exactChanges > 0 && tauChanges > 0)
        try #require(publications.count == requests.count)
        for publication in publications {
            #expect(publication.value == state.value(species: publication.speciesIndex, lane: publication.laneIndex))
            #expect((publication.exactCount == nil) == (publication.authority == .deterministicRK2))
        }
        #expect(publications[3].exactCount == UInt32.max)
        #expect(publications[3].value == Double(UInt32.max))
        try record(state, "hybrid-mixed-authorities")
    }

    @Test func exactTrajectoryIsInvariantToDispatchChunkSize() async throws {
        let model = try VivoExactSSAModel(species: ["population"], reactions: [
            .init(id: "birth", law: .zeroOrder, rate: 240, changes: [.init(speciesIndex: 0, delta: 1)]),
            .init(id: "death", law: .firstOrder, reactantA: 0, rate: 2,
                  changes: [.init(speciesIndex: 0, delta: -1)])
        ])
        let plan = try VivoHybridExplicitPlanBuilder.make(model: model), lanes: UInt32 = 53
        let initial = [UInt32](repeating: 50, count: Int(lanes)), gpu = try device()
        let small = try await VivoHybridReactionRuntime.make(model: model, plan: plan,
            configuration: configuration(lanes: lanes, dt: 0.1, chunk: 1, dispatches: 512),
            initialCounts: initial, device: gpu)
        let large = try await VivoHybridReactionRuntime.make(model: model, plan: plan,
            configuration: configuration(lanes: lanes, dt: 0.1, chunk: 128, dispatches: 4),
            initialCounts: initial, device: gpu)
        for _ in 0..<3 {
            let a = try await small.step(), b = try await large.step()
            try #require(a.disposition == .committed && b.disposition == .committed)
            #expect(a.exactDispatches > 1)
            #expect(a.status.maximumExactEvents == b.status.maximumExactEvents)
            #expect(try await small.snapshot() == large.snapshot())
        }
        #expect(try await small.checkpoint() == large.checkpoint())
        try record(try await small.checkpoint(), "hybrid-exact-chunk-invariance")
    }

    private func changed(_ checkpoint: VivoHybridCheckpoint,
                         _ transform: (inout [String: Any]) -> Void) throws -> VivoHybridCheckpoint {
        var object = try #require(JSONSerialization.jsonObject(with: JSONEncoder().encode(checkpoint)) as? [String: Any])
        transform(&object)
        return try JSONDecoder().decode(VivoHybridCheckpoint.self, from: JSONSerialization.data(withJSONObject: object))
    }

    @Test func checkpointContinuationAndFailedRestoreAreAtomic() async throws {
        let lanes: UInt32 = 19, model = try mixedModel(), plan = try mixedPlan(model), gpu = try device()
        let config = configuration(lanes: lanes)
        let original = try await VivoHybridReactionRuntime.make(model: model, plan: plan,
            configuration: config, initialCounts: mixedInitial(lanes: lanes), device: gpu)
        for _ in 0..<2 { try #require(try await original.step().disposition == .committed) }
        let accepted = try await original.checkpoint()
        let resumed = try await VivoHybridReactionRuntime.make(model: model, plan: plan,
            configuration: config, initialCounts: mixedInitial(lanes: lanes), device: gpu)
        try await resumed.restore(accepted)
        #expect(try await resumed.checkpoint() == accepted)
        let invalid: [VivoHybridCheckpoint] = try [
            changed(accepted) { $0["seed"] = NSNumber(value: seed ^ 1) },
            changed(accepted) { $0["modelFingerprint"] = String(repeating: "0", count: 64) },
            changed(accepted) { $0["planFingerprint"] = String(repeating: "0", count: 64) },
            changed(accepted) { $0["numericalABIVersion"] = 2 },
            changed(accepted) { $0["countsUInt32LE"] = Data(accepted.countsUInt32LE.dropLast()).base64EncodedString() },
            changed(accepted) {
                var inactive = accepted.continuousFP32LE
                inactive.replaceSubrange(0..<4, with: [0, 0, 128, 63]) // 1.0 in an exact-owned slot.
                $0["continuousFP32LE"] = inactive.base64EncodedString()
            },
            changed(accepted) {
                var inactive = accepted.countsUInt32LE
                let offset = 4 * Int(lanes) * 4 // First RK2 species has no count representation.
                inactive.replaceSubrange(offset..<(offset + 4), with: [1, 0, 0, 0])
                $0["countsUInt32LE"] = inactive.base64EncodedString()
            }
        ]
        for checkpoint in invalid {
            await rejects { try await resumed.restore(checkpoint) }
            #expect(try await resumed.checkpoint() == accepted)
        }
        let pending = try await resumed.prepareStep()
        try #require(pending.canCommit)
        await rejects { try await resumed.restore(accepted) }
        await rejects { _ = try await resumed.snapshot() }
        try await resumed.discardPreparedStep(transactionID: pending.transactionID)
        #expect(try await resumed.checkpoint() == accepted)
        for _ in 0..<4 {
            try #require(try await original.step().disposition == .committed)
            try #require(try await resumed.step().disposition == .committed)
            #expect(try await original.checkpoint() == resumed.checkpoint())
        }
        try record(try await resumed.checkpoint(), "hybrid-resumed-continuation")
    }

    private func rejectingModel() throws -> VivoExactSSAModel {
        try .init(species: ["exact", "continuous"], reactions: [
            .init(id: "birth", law: .zeroOrder, rate: 64, changes: [.init(speciesIndex: 0, delta: 1)]),
            .init(id: "decay", law: .firstOrder, reactantA: 1, rate: 0.1,
                  changes: [.init(speciesIndex: 1, delta: -1)])
        ])
    }

    @Test func rejectedExactCandidateDoesNotCommitOtherAuthorities() async throws {
        let model = try rejectingModel()
        let plan = try VivoHybridExplicitPlanBuilder.make(model: model, reactionAuthorities: ["decay": .deterministicRK2])
        let runtime = try await VivoHybridReactionRuntime.make(model: model, plan: plan,
            configuration: configuration(lanes: 1, dt: 1), initialCounts: [.max, 100], device: device())
        let before = try await runtime.checkpoint()
        // The FP32 open-interval uniform bounds imply the first wait is < 1 s
        // at rate 64. Overflow is therefore guaranteed, independent of the seed.
        let rejected = try await runtime.prepareStep(publications: [.init(speciesIndex: 1, laneIndex: 0)])
        #expect(rejected.disposition == .rejected)
        #expect(rejected.status.flags & 4 != 0)
        #expect(rejected.publications.isEmpty)
        #expect(await runtime.hasPendingTransaction() == false)
        await rejects { _ = try await runtime.commitPreparedStep(transactionID: rejected.transactionID) }
        #expect(try await runtime.checkpoint() == before)
    }

    @Test func exactWorkExhaustionDoesNotCommitOrShortenTheHorizon() async throws {
        let model = try rejectingModel()
        let plan = try VivoHybridExplicitPlanBuilder.make(model: model, reactionAuthorities: ["decay": .deterministicRK2])
        let runtime = try await VivoHybridReactionRuntime.make(model: model, plan: plan,
            configuration: configuration(lanes: 1, dt: 1, chunk: 1, dispatches: 1),
            initialCounts: [0, 100], device: device())
        let before = try await runtime.checkpoint()
        let exhausted = try await runtime.step(publications: [.init(speciesIndex: 1, laneIndex: 0)])
        #expect(exhausted.disposition == .exactWorkBudgetExceeded)
        #expect(exhausted.status.valid && exhausted.status.unfinishedExactLanes == 1)
        #expect(exhausted.exactDispatches == 1 && exhausted.status.maximumExactEvents == 1)
        #expect(exhausted.timeBefore == 0 && exhausted.timeAfter == 0 && exhausted.requestedTimeStep == 1)
        #expect(exhausted.publications.isEmpty)
        #expect(await runtime.hasPendingTransaction() == false)
        #expect(try await runtime.checkpoint() == before)
    }

    private struct Moments: Encodable {
        let authority: VivoHybridExecutionMode
        let distribution: String
        let sampleCount: Int
        let expectedMean: Double
        let observedMean: Double
        let expectedVariance: Double
        let observedVariance: Double
        let meanTolerance: Double
        let varianceTolerance: Double
    }

    private func moments(_ counts: [UInt32], authority: VivoHybridExecutionMode, distribution: String,
                         mean: Double, variance: Double, fourthCentralMoment: Double) -> Moments {
        let n = Double(counts.count), values = counts.map(Double.init)
        let observedMean = values.reduce(0, +) / n
        let observedVariance = values.reduce(0) { $0 + ($1 - observedMean) * ($1 - observedMean) } / (n - 1)
        // Bounds are fixed from analytical sampling errors before observing data.
        // Six standard errors are deliberately conservative for fixed-seed CI.
        let meanTolerance = 6 * sqrt(variance / n)
        let varianceTolerance = 6 * sqrt((fourthCentralMoment - ((n - 3) / (n - 1)) * variance * variance) / n)
        #expect(abs(observedMean - mean) <= meanTolerance)
        #expect(abs(observedVariance - variance) <= varianceTolerance)
        return .init(authority: authority, distribution: distribution, sampleCount: counts.count,
                     expectedMean: mean, observedMean: observedMean, expectedVariance: variance,
                     observedVariance: observedVariance, meanTolerance: meanTolerance, varianceTolerance: varianceTolerance)
    }

    @Test func zeroOrderBirthsMatchPoissonMoments() async throws {
        let lanes: UInt32 = 8_192, gpu = try device()
        var reports: [Moments] = []
        // Lambda=4 covers SSA and inversion; lambda=64 covers tau PTRS.
        for (authority, rate) in [(VivoHybridExecutionMode.exactSSA, Float(8)), (.tauLeap, Float(8)), (.tauLeap, Float(128))] {
            let model = try VivoExactSSAModel(species: ["births"], reactions: [
                .init(id: "birth", law: .zeroOrder, rate: rate, changes: [.init(speciesIndex: 0, delta: 1)])
            ])
            let plan = try VivoHybridExplicitPlanBuilder.make(model: model, reactionAuthorities: ["birth": authority])
            let runtime = try await VivoHybridReactionRuntime.make(model: model, plan: plan,
                configuration: configuration(lanes: lanes, dt: 0.5),
                initialCounts: .init(repeating: 0, count: Int(lanes)), device: gpu)
            try #require(try await runtime.step().disposition == .committed)
            let state = try await runtime.snapshot(), lambda = Double(rate) * 0.5
            reports.append(moments(state.counts, authority: authority, distribution: "Poisson",
                mean: lambda, variance: lambda, fourthCentralMoment: lambda + 3 * lambda * lambda))
            if lambda == 4 {
                let probability = exp(-lambda), expectedZeros = Double(lanes) * probability
                let zeros = Double(state.counts.filter { $0 == 0 }.count)
                #expect(abs(zeros - expectedZeros) <= 6 * sqrt(Double(lanes) * probability * (1 - probability)))
            }
        }
        try record(reports, "hybrid-poisson-moments")
    }

    @Test func exactFirstOrderDecayMatchesBinomialMoments() async throws {
        let lanes: UInt32 = 8_192, initial: UInt32 = 40
        let model = try VivoExactSSAModel(species: ["survivors"], reactions: [
            .init(id: "death", law: .firstOrder, reactantA: 0, rate: 0.5,
                  changes: [.init(speciesIndex: 0, delta: -1)])
        ])
        let plan = try VivoHybridExplicitPlanBuilder.make(model: model)
        let runtime = try await VivoHybridReactionRuntime.make(model: model, plan: plan,
            configuration: configuration(lanes: lanes, dt: 1),
            initialCounts: .init(repeating: initial, count: Int(lanes)), device: device())
        try #require(try await runtime.step().disposition == .committed)
        let state = try await runtime.snapshot(), survival = exp(-0.5)
        #expect(state.counts.allSatisfy { $0 <= initial })
        let variance = Double(initial) * survival * (1 - survival)
        let fourth = 3 * variance * variance + variance * (1 - 6 * survival * (1 - survival))
        let report = moments(state.counts, authority: .exactSSA, distribution: "Binomial",
            mean: Double(initial) * survival, variance: variance, fourthCentralMoment: fourth)
        try record(report, "hybrid-exact-death-moments")
    }
}
