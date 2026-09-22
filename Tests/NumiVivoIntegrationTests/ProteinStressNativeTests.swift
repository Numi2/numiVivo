import Foundation
import Testing
@testable import NumiVivoKit

@Suite("Protein stress native integration")
struct ProteinStressNativeTests {
    @Test func immutableProviderMatchesAnalyticSpring() async throws {
        let request = try VivoProteinStressExample.harmonicFixture()
        let compiled = try VivoProteinStressCompilation(request)
        let provider = try #require(try compiled.forceProvider(stage: 1))
        let geometry = try VivoMDCandidateGeometry(particlePositionsNM: request.sourceCheckpoint.positionsNM, periodicCell: nil)
        let evaluation = try await provider.evaluate(geometry)
        #expect(abs(evaluation.additionalEnergyKJPerMol - 0.0625) < 1e-12)
        #expect(abs(evaluation.physicalParticleForcesKJPerMolNM[1].x - 2.5) < 1e-12)
        #expect(evaluation.physicalParticleForcesKJPerMolNM.reduce(.zero, +).norm < 1e-12)
        #expect(try compiled.configurationFingerprint(stage: 0) != request.sourceCheckpoint.configurationFingerprint)
    }

    @Test func transitionPreservesStateAndAccountsForWork() throws {
        let request = try VivoProteinStressExample.harmonicFixture(), compiled = try VivoProteinStressCompilation(request)
        let first = try compiled.transfer(request.sourceCheckpoint, from: nil, to: 0)
        #expect(first.workKJPerMol == 0)
        let second = try compiled.transfer(first.checkpoint, from: 0, to: 1)
        #expect(abs(second.workKJPerMol - 0.0625) < 1e-12)
        #expect(second.checkpoint.positionsNM == request.sourceCheckpoint.positionsNM)
        #expect(second.checkpoint.velocitiesNMPerPS == request.sourceCheckpoint.velocitiesNMPerPS)
        #expect(second.checkpoint.acceptedStep == request.sourceCheckpoint.acceptedStep)
        #expect(second.checkpoint.timePS == request.sourceCheckpoint.timePS)
        #expect(throws: (any Error).self) { try compiled.transfer(request.sourceCheckpoint, from: 0, to: 1) }
    }

    @Test func malformedSchemaAndDuplicateReplicaSeedsAreRejected() throws {
        var request = try VivoProteinStressExample.harmonicFixture()
        request.schema = "invalid"
        #expect(throws: (any Error).self) { try VivoProteinStressCompilation(request) }
        request = try VivoProteinStressExample.harmonicFixture()
        let campaign = VivoProteinStressCampaignRequest(candidates: [.init(identifier: "duplicate", replicas: [request, request])],
                                                        comparisonDescription: "synthetic matching attachments")
        #expect(throws: (any Error).self) { try VivoProteinStressCampaign.validate(campaign) }
    }

    @Test func journalVerifiesMetricsAndRejectsForgedSwitchWork() async throws {
        let directory = FileManager.default.temporaryDirectory.appendingPathComponent(UUID().uuidString)
        defer { try? FileManager.default.removeItem(at: directory) }
        let store = try VivoArtifactStore(rootURL: directory)
        let request = try VivoProteinStressExample.harmonicFixture(), compiled = try VivoProteinStressCompilation(request)
        let requestArtifact = try await VivoProteinStressRunner.put(request, kind: "protein-stress-request", store: store)
        let checkpoint = try await VivoProteinStressRunner.put(request.sourceCheckpoint, kind: "md-checkpoint", store: store)
        let source = VivoProteinStressJournalEntry(request: requestArtifact.fingerprint, previous: nil, action: .source,
            stageIndex: -1, stepsInStage: 0, checkpoint: checkpoint.fingerprint, cumulativeProtocolWorkKJPerMol: 0,
            frame: try compiled.frame(checkpoint: request.sourceCheckpoint, stage: nil), rejection: nil)
        let root = try await VivoProteinStressRunner.put(source, kind: "protein-stress-journal", store: store)
        let transfer = try compiled.transfer(request.sourceCheckpoint, from: nil, to: 0)
        let destination = try await VivoProteinStressRunner.put(transfer.checkpoint, kind: "md-checkpoint", store: store)
        let valid = VivoProteinStressJournalEntry(request: requestArtifact.fingerprint, previous: root.fingerprint, action: .stageStart,
            stageIndex: 0, stepsInStage: 0, checkpoint: destination.fingerprint, cumulativeProtocolWorkKJPerMol: 0,
            frame: try compiled.frame(checkpoint: transfer.checkpoint, stage: 0), rejection: nil)
        let validArtifact = try await VivoProteinStressRunner.put(valid, kind: "protein-stress-journal", store: store)
        #expect(try await VivoProteinStressRunner.verify(request, store: store, journalTail: validArtifact.fingerprint).count == 2)
        let forged = VivoProteinStressJournalEntry(request: requestArtifact.fingerprint, previous: root.fingerprint, action: .stageStart,
            stageIndex: 0, stepsInStage: 0, checkpoint: destination.fingerprint, cumulativeProtocolWorkKJPerMol: 123,
            frame: valid.frame, rejection: nil)
        let forgedArtifact = try await VivoProteinStressRunner.put(forged, kind: "protein-stress-journal", store: store)
        await #expect(throws: (any Error).self) {
            try await VivoProteinStressRunner.verify(request, store: store, journalTail: forgedArtifact.fingerprint)
        }
    }

    /// Requires an actual supported Metal device. No skip or CPU fallback is
    /// accepted as a successful Apple-runtime qualification.
    @Test func metalExecutionAndExactAcceptedPrefixResume() async throws {
        let directory = FileManager.default.temporaryDirectory.appendingPathComponent(UUID().uuidString)
        defer { try? FileManager.default.removeItem(at: directory) }
        let store = try VivoArtifactStore(rootURL: directory), request = try VivoProteinStressExample.harmonicFixture()
        let full = try await VivoProteinStressRunner.run(request, store: store)
        #expect(full.disposition == .completed, "\(full.diagnostic ?? "no diagnostic")")
        let tail = try #require(full.journalTail)
        let verified = try await VivoProteinStressRunner.verify(request, store: store, journalTail: tail)
        let last = try #require(verified.last)
        #expect(last.acceptedStep == 8)
        #expect(last.entry.stageIndex == 1 && last.entry.stepsInStage == 4)
        let prefix = try #require(verified.first { $0.entry.action == .sample && $0.entry.stageIndex == 0 })
        let resumed = try await VivoProteinStressRunner.run(request, store: store, resumeFrom: prefix.fingerprint)
        #expect(resumed.disposition == .completed, "\(resumed.diagnostic ?? "no diagnostic")")
        #expect(resumed.journalTail == full.journalTail)
    }
}
