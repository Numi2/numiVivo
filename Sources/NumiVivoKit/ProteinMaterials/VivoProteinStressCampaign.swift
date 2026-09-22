import Foundation

public struct VivoProteinStressCandidate: Codable, Sendable, Equatable {
    public let identifier: String
    public let replicas: [VivoProteinStressRequest]
    public init(identifier: String, replicas: [VivoProteinStressRequest]) { self.identifier = identifier; self.replicas = replicas }
}

public struct VivoProteinStressCampaignRequest: Codable, Sendable, Equatable {
    public let candidates: [VivoProteinStressCandidate]
    /// Declared correspondence of attachment sites and selected structures across
    /// candidates; this is a caller assertion, not an automatic anatomical match.
    public let comparisonDescription: String
    public init(candidates: [VivoProteinStressCandidate], comparisonDescription: String) {
        self.candidates = candidates; self.comparisonDescription = comparisonDescription
    }
}

public struct VivoProteinStressCandidateResult: Codable, Sendable, Equatable {
    public let identifier: String
    public let runs: [VivoProteinStressReceipt]
    public let allReplicasCompleted: Bool
    public let sampledPeakTensileForcePN: VivoProteinReplicaStatistics?
    public let finalRetainedContactFraction: VivoProteinReplicaStatistics?
    public let exclusion: String?
}

public struct VivoProteinStressCampaignReport: Codable, Sendable, Equatable {
    public var schema: String = "numivivo.org/protein-stress-campaign/v1"
    public let request: VivoFingerprint
    public let candidates: [VivoProteinStressCandidateResult]
    /// Nondominated complete candidates on sampled peak force and final contact
    /// retention. This is not a ranking of true strength or thermal stability.
    public let conditionalParetoCandidates: [String]
    public let interpretation: String
}

public enum VivoProteinStressCampaign {
    public static func validate(_ request: VivoProteinStressCampaignRequest) throws {
        guard !request.candidates.isEmpty, request.candidates.count <= 128,
              !request.comparisonDescription.isEmpty, request.comparisonDescription.utf8.count <= 4096,
              Set(request.candidates.map(\.identifier)).count == request.candidates.count,
              let baseline = request.candidates.first?.replicas.first else {
            throw VivoProteinStressError.invalid("empty or oversized candidate campaign")
        }
        // Native topology/force-field provenance stays with each request. Matching
        // numerical protocols alone does not establish comparable model accuracy.
        var baselineConfiguration = baseline.sourceConfiguration; baselineConfiguration.randomSeed = 0
        var totalSteps: UInt64 = 0
        for candidate in request.candidates {
            guard !candidate.identifier.isEmpty, candidate.identifier.utf8.count <= 128,
                  (2...32).contains(candidate.replicas.count),
                  Set(candidate.replicas.map(\.replicaID)).count == candidate.replicas.count,
                  Set(candidate.replicas.map(\.randomSeed)).count == candidate.replicas.count else {
                throw VivoProteinStressError.invalid("candidate requires 2–32 distinct replica identities and seeds")
            }
            let candidateSource = candidate.replicas[0]
            for replica in candidate.replicas {
                _ = try VivoProteinStressCompilation(replica)
                guard replica.system == candidateSource.system, replica.pull == candidateSource.pull,
                      replica.selection == candidateSource.selection, replica.hydrogenBonds == candidateSource.hydrogenBonds,
                      replica.nativeContacts == candidateSource.nativeContacts else {
                    throw VivoProteinStressError.invalid("replicas of one candidate changed the system or analysis definition")
                }
                var configuration = replica.sourceConfiguration; configuration.randomSeed = 0
                guard configuration == baselineConfiguration, replica.stages == baseline.stages,
                      replica.sampleEvery == baseline.sampleEvery,
                      replica.pull?.stiffnessKJPerMolNM2 == baseline.pull?.stiffnessKJPerMolNM2,
                      replica.pull?.projectionAxis == baseline.pull?.projectionAxis,
                      (replica.pull == nil) == (baseline.pull == nil),
                      replica.hydrogenBondCriteria == baseline.hydrogenBondCriteria,
                      replica.maximumContactDistanceRatio == baseline.maximumContactDistanceRatio else {
                    throw VivoProteinStressError.invalid("candidate protocols differ; do not pool incomparable measurements")
                }
                totalSteps += replica.stages.reduce(0) { $0 + $1.steps }
            }
        }
        guard totalSteps <= 10_000_000 else { throw VivoProteinStressError.invalid("campaign step budget exceeds 10 million") }
    }

    public static func run(_ request: VivoProteinStressCampaignRequest, store: VivoArtifactStore) async throws -> VivoProteinStressCampaignReport {
        try validate(request)
        let source = try await VivoProteinStressRunner.put(request, kind: "protein-stress-campaign-request", store: store)
        var results: [VivoProteinStressCandidateResult] = []
        for candidate in request.candidates {
            var runs: [VivoProteinStressReceipt] = [], peaks: [Double] = [], contacts: [Double] = []
            for replica in candidate.replicas {
                try Task.checkCancellation()
                let receipt = try await VivoProteinStressRunner.run(replica, store: store)
                _ = try await VivoProteinStressRunner.put(receipt, kind: "protein-stress-receipt", store: store)
                runs.append(receipt)
                guard receipt.disposition == .completed, let tail = receipt.journalTail else { continue }
                let verified = try await VivoProteinStressRunner.verify(replica, store: store, journalTail: tail)
                // Include entry frames: the protocol may begin with a preloaded
                // restraint. Never silently discard a high initial force.
                let values = verified.compactMap { $0.entry.frame.tensileForcePN }
                if let peak = values.max() { peaks.append(peak) }
                if let value = verified.last?.entry.frame.structure.retainedContactFraction { contacts.append(value) }
            }
            let complete = runs.count == candidate.replicas.count && runs.allSatisfy { $0.disposition == .completed }
            results.append(.init(identifier: candidate.identifier, runs: runs, allReplicasCompleted: complete,
                sampledPeakTensileForcePN: complete && peaks.count == runs.count ? try VivoProteinReplicaStatistics.calculate(peaks) : nil,
                finalRetainedContactFraction: complete && contacts.count == runs.count ? try VivoProteinReplicaStatistics.calculate(contacts) : nil,
                exclusion: complete ? nil : "one or more declared replicas failed, were rejected, or were cancelled; no successful-subset score"))
        }
        let eligible = results.filter { $0.allReplicasCompleted && $0.sampledPeakTensileForcePN != nil && $0.finalRetainedContactFraction != nil }
        let frontier = eligible.filter { candidate in
            !eligible.contains { other in
                guard other.identifier != candidate.identifier else { return false }
                let f = candidate.sampledPeakTensileForcePN!.mean, q = candidate.finalRetainedContactFraction!.mean
                let of = other.sampledPeakTensileForcePN!.mean, oq = other.finalRetainedContactFraction!.mean
                return of >= f && oq >= q && (of > f || oq > q)
            }
        }.map(\.identifier).sorted()
        let report = VivoProteinStressCampaignReport(request: source.fingerprint, candidates: results, conditionalParetoCandidates: frontier,
            interpretation: "Descriptive complete-replica screening under matching declared protocols. Attachment/selection correspondence is caller-declared. Force-field accuracy, independent conformational exploration, protein unfolding, experimental thermostability and material-scale transfer are not established. No sequence-generation model is invoked.")
        _ = try await VivoProteinStressRunner.put(report, kind: "protein-stress-campaign-report", store: store)
        return report
    }
}
