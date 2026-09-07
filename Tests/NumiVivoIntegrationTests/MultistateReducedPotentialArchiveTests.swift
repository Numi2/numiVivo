import Foundation
import Testing
@testable import NumiVivoKit

private actor MultistateArchiveFailureGate {
    private var calls = 0
    let failOnCall: Int
    init(failOnCall: Int) { self.failOnCall = failOnCall }
    func check() throws {
        calls += 1
        if calls == failOnCall { throw VivoChemistryError.convergence("synthetic interrupted origin propagation") }
    }
}

@Suite(.serialized) struct MultistateReducedPotentialArchiveTests {
    private func fingerprint(_ text:String) throws -> VivoFingerprint {
        try VivoCanonicalJSON.fingerprint(Data(text.utf8))
    }

    private func directory() throws -> URL {
        let root = URL(fileURLWithPath:NSTemporaryDirectory(), isDirectory:true)
            .appendingPathComponent("numivivo-multistate-\(UUID().uuidString)", isDirectory:true)
        try FileManager.default.createDirectory(at:root, withIntermediateDirectories:true)
        return root
    }

    private func physical() throws -> VivoConstantPHPhysicalState {
        try .init(stepIndex:0, timePS:0, positionsNM:[.zero],
                  velocitiesNMPerPS:[.init(0.001,0,0)], periodicCell:nil)
    }

    private func executables(manifold:VivoFingerprint,
                             hamiltonians:[VivoFingerprint],
                             temperatureK:Double,
                             failureGate:MultistateArchiveFailureGate? = nil) -> [VivoConstantPHExecutableState] {
        let rt = 0.00831446261815324 * temperatureK
        return hamiltonians.indices.map { index in
            VivoConstantPHExecutableState(identifier:"state-\(index)",
                physicalManifoldFingerprint:manifold,
                hamiltonianFingerprint:hamiltonians[index],
                potentialEnergyKJPerMol:{ _ in Double(index) * rt },
                propagate:{ state, steps in
                    if index == 0, let failureGate { try await failureGate.check() }
                    return try VivoConstantPHPhysicalState(stepIndex:state.stepIndex + steps,
                        timePS:state.timePS + Double(steps) * 0.001,
                        positionsNM:state.positionsNM,
                        velocitiesNMPerPS:state.velocitiesNMPerPS,
                        periodicCell:state.periodicCell)
                })
        }
    }

    private func request(manifold:VivoFingerprint,
                         hamiltonians:[VivoFingerprint],
                         initials:[VivoConstantPHPhysicalState],
                         temperatureK:Double = 300) throws -> VivoMultistateReducedPotentialArchiveRequest {
        try .init(identifier:"two-state-archive", temperatureK:temperatureK,
                  stateIdentifiers:["state-0","state-1"],
                  physicalManifoldFingerprint:manifold,
                  hamiltonianFingerprints:hamiltonians,
                  initialPhysicalStateFingerprints:initials.map { try $0.fingerprint() },
                  equilibrationSteps:0, productionSteps:8, sampleEvery:2,
                  checkpointBlockSteps:4, maximumMatrixElements:10_000,
                  mbar:.init(residualTolerance:1e-12, maximumIterations:10_000,
                             minimumBidirectionalAcceptanceOverlap:0.5,
                             minimumTargetEffectiveSamples:2,
                             bootstrapReplicates:4, bootstrapSeed:73,
                             maximumWorkElements:10_000))
    }

    @Test func archiveEvaluatesEveryHamiltonianAndSolvesMBAR() async throws {
        let root = try directory(); defer { try? FileManager.default.removeItem(at:root) }
        let store = try VivoArtifactStore(rootURL:root)
        let manifold = try fingerprint("common-manifold")
        let hamiltonians = try [fingerprint("h0"), fingerprint("h1")]
        let initials = try [physical(), physical()]
        let req = try request(manifold:manifold, hamiltonians:hamiltonians, initials:initials)
        let result = try await VivoMultistateReducedPotentialArchiveRunner.run(req,
            initialPhysicalStates:initials,
            executableStates:executables(manifold:manifold, hamiltonians:hamiltonians, temperatureK:req.temperatureK),
            store:store)
        #expect(result.status == .converged)
        #expect(result.completedOrigins == 2)
        #expect(result.chunkFingerprints.count == 4)
        let mbar = try #require(result.mbar)
        #expect(mbar.sampleCount == 8)
        #expect(mbar.converged)
        #expect(abs(mbar.estimates[0].relativeFreeEnergy) < 1e-12)
        #expect(abs(mbar.estimates[1].relativeFreeEnergy - 1.0) < 1e-9)
    }

    @Test func rejectedPartialBlockResumesFromLastDurableBlockWithoutDuplicatingSamples() async throws {
        let root = try directory(); defer { try? FileManager.default.removeItem(at:root) }
        let store = try VivoArtifactStore(rootURL:root)
        let manifold = try fingerprint("resume-manifold")
        let hamiltonians = try [fingerprint("resume-h0"), fingerprint("resume-h1")]
        let initials = try [physical(), physical()]
        let req = try request(manifold:manifold, hamiltonians:hamiltonians, initials:initials)
        let gate = MultistateArchiveFailureGate(failOnCall:3)
        let interrupted = try await VivoMultistateReducedPotentialArchiveRunner.run(req,
            initialPhysicalStates:initials,
            executableStates:executables(manifold:manifold, hamiltonians:hamiltonians,
                                         temperatureK:req.temperatureK, failureGate:gate),
            store:store)
        #expect(interrupted.status == .rejected)
        #expect(interrupted.completedOrigins == 0)
        #expect(interrupted.currentOriginProductionSteps == 4)
        #expect(interrupted.chunkFingerprints.count == 1)

        let resumed = try await VivoMultistateReducedPotentialArchiveRunner.run(req,
            initialPhysicalStates:initials,
            executableStates:executables(manifold:manifold, hamiltonians:hamiltonians,
                                         temperatureK:req.temperatureK),
            store:store, resumeFrom:interrupted.checkpointFingerprint)
        #expect(resumed.status == .converged)
        #expect(resumed.completedOrigins == 2)
        #expect(resumed.chunkFingerprints.count == 4)
        #expect(resumed.mbar?.sampleCount == 8)
        #expect(abs((resumed.mbar?.estimates[1].relativeFreeEnergy ?? -99) - 1.0) < 1e-9)
    }

    @Test func requestRejectsMatrixWorkBeyondBothArchiveAndMBARBudgets() throws {
        let manifold = try fingerprint("budget-manifold")
        let hamiltonians = try [fingerprint("budget-h0"), fingerprint("budget-h1")]
        let initials = try [physical(), physical()]
        var req = try request(manifold:manifold, hamiltonians:hamiltonians, initials:initials)
        req.maximumMatrixElements = 4
        #expect(throws:(any Error).self) { try req.validate() }
        req.maximumMatrixElements = 10_000
        req.mbar.maximumWorkElements = 4
        #expect(throws:(any Error).self) { try req.validate() }
    }
}
