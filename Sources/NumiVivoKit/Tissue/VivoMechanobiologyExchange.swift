import Foundation

/// Interface for future Matter ↔ biology exchange. No constitutive law or
/// committed state is invented by this contract. Runtime participants must use
/// native prepare/release semantics (as in VivoMolecularPhysiologyCoordinator).
public struct VivoMechanobiologyChannel: Codable, Sendable {
    public let identifier: String
    public let direction: String // mechanicsToBiology | biologyToMechanics
    public let quantity: String // geometry, stress, strain, transport, permeability, force, growth, adhesion, materialState
    public let unit: String
    public let entityIDs: [String]
    public let modelFingerprint: String
    public let evidence: VivoTissueEvidence
}
public struct VivoMechanobiologyProposal: Codable, Sendable {
    public let transactionID: UUID
    public let specimenFingerprint: String
    public let mechanicsCheckpoint: String
    public let biologyCheckpoint: String
    public let timeBeforeSeconds: Double
    public let proposedTimeAfterSeconds: Double
    public let channels: [VivoMechanobiologyChannel]
    public let convergenceTolerance: Double
}
public struct VivoMechanobiologyPreparedState: Sendable {
    public let transactionID: UUID
    public let owner: String
    public let previousCheckpoint: String
    public let candidateCheckpoint: String
    public let acceptedTimeSeconds: Double
    public let converged: Bool
    public let residual: Double
}
/// Prepared states are private until both owners accept identical time and
/// lineage. A failed release must quarantine the experiment; rollback cannot
/// be claimed unless both native owners verify restoration of exact checkpoints.
public protocol VivoMechanobiologyParticipant: Actor {
    func prepare(_ proposal: VivoMechanobiologyProposal) async throws -> VivoMechanobiologyPreparedState
    func release(transactionID: UUID, peer: VivoMechanobiologyPreparedState) async throws
    func abort(transactionID: UUID) async throws -> String
}
