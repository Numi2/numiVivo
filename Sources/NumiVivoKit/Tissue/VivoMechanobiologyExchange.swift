import Foundation

/// Exchange contract, not a second physics/biology solver. Runtime participants
/// retain native prepare/release authority and checkpoint ownership.
public struct VivoMechanobiologyChannel: Codable, Sendable {
    public let identifier: String
    public let direction: String // mechanicsToBiology | biologyToMechanics
    public let quantity: String
    public let unit: String
    public let entityIDs: [String]
    /// Row-major components per entity; units and component names are explicit.
    public let components: [String]
    public let values: [Double]
    public let modelFingerprint: String
    public let evidence: VivoTissueEvidence
    public init(identifier: String, direction: String, quantity: String, unit: String,
                entityIDs: [String], components: [String], values: [Double],
                modelFingerprint: String, evidence: VivoTissueEvidence) {
        self.identifier=identifier; self.direction=direction; self.quantity=quantity
        self.unit=unit; self.entityIDs=entityIDs; self.components=components
        self.values=values; self.modelFingerprint=modelFingerprint; self.evidence=evidence
    }
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
    public init(transactionID: UUID, specimenFingerprint: String, mechanicsCheckpoint: String,
                biologyCheckpoint: String, timeBeforeSeconds: Double, proposedTimeAfterSeconds: Double,
                channels: [VivoMechanobiologyChannel], convergenceTolerance: Double) {
        self.transactionID=transactionID; self.specimenFingerprint=specimenFingerprint
        self.mechanicsCheckpoint=mechanicsCheckpoint; self.biologyCheckpoint=biologyCheckpoint
        self.timeBeforeSeconds=timeBeforeSeconds; self.proposedTimeAfterSeconds=proposedTimeAfterSeconds
        self.channels=channels; self.convergenceTolerance=convergenceTolerance
    }
    private func need(_ condition: Bool, _ reason: String) throws {
        if !condition { throw NSError(domain:"NumiVivo.Mechanobiology",code:1,userInfo:[NSLocalizedDescriptionKey:reason]) }
    }
    public func validate() throws {
        try need(!specimenFingerprint.isEmpty && !mechanicsCheckpoint.isEmpty && !biologyCheckpoint.isEmpty,"missing lineage")
        try need(timeBeforeSeconds.isFinite && proposedTimeAfterSeconds.isFinite && timeBeforeSeconds>=0 && proposedTimeAfterSeconds>timeBeforeSeconds,"invalid clock interval")
        try need(convergenceTolerance.isFinite && convergenceTolerance>=0,"invalid convergence tolerance")
        try need(!channels.isEmpty && Set(channels.map(\.identifier)).count==channels.count,"channel identity")
        let directions=["mechanicsToBiology":Set(["geometry","stress","strain","transport"]),
                        "biologyToMechanics":Set(["permeability","force","growth","adhesion","materialState"])]
        for c in channels {
            try need(!c.identifier.isEmpty && directions[c.direction]?.contains(c.quantity)==true,"unsupported channel direction/quantity")
            try need(!c.unit.isEmpty && !c.modelFingerprint.isEmpty && !c.entityIDs.isEmpty && Set(c.entityIDs).count==c.entityIDs.count && c.entityIDs.allSatisfy{!$0.isEmpty},"channel units/model/entities")
            try need(!c.components.isEmpty && Set(c.components).count==c.components.count && c.components.allSatisfy{!$0.isEmpty},"channel components")
            let (size,overflow)=c.entityIDs.count.multipliedReportingOverflow(by:c.components.count)
            try need(!overflow && c.values.count==size && c.values.allSatisfy{$0.isFinite},"channel shape/value")
            try need(c.evidence.state != .unavailable && !c.evidence.sourceID.isEmpty && !c.evidence.uncertainty.isEmpty,"unavailable channel cannot supply values")
            if [.inference,.validatedDomain,.outOfDistribution].contains(c.evidence.state) {
                try need(c.evidence.modelID?.isEmpty==false,"model evidence identity")
            }
            if c.evidence.state == .validatedDomain { try need(c.evidence.validationDomain?.isEmpty==false,"validation domain") }
        }
    }
    /// Structural admission before release. This cannot establish native rollback,
    /// convergence correctness, physical validity or atomic distributed commit.
    public func validatePrepared(mechanics: VivoMechanobiologyPreparedState,
                                 biology: VivoMechanobiologyPreparedState) throws {
        try validate()
        for (state,owner,checkpoint) in [(mechanics,"mechanics",mechanicsCheckpoint),(biology,"biology",biologyCheckpoint)] {
            try need(state.transactionID==transactionID && state.owner==owner && state.previousCheckpoint==checkpoint,"prepared lineage mismatch")
            try need(state.acceptedTimeSeconds==proposedTimeAfterSeconds,"prepared clock mismatch")
            try need(state.converged && state.residual.isFinite && state.residual>=0 && state.residual<=convergenceTolerance,"participant not converged")
            try need(!state.candidateCheckpoint.isEmpty,"missing candidate checkpoint")
        }
    }
}
public struct VivoMechanobiologyPreparedState: Sendable {
    public let transactionID: UUID
    public let owner: String
    public let previousCheckpoint: String
    public let candidateCheckpoint: String
    public let acceptedTimeSeconds: Double
    public let converged: Bool
    public let residual: Double
    public init(transactionID: UUID, owner: String, previousCheckpoint: String,
                candidateCheckpoint: String, acceptedTimeSeconds: Double, converged: Bool, residual: Double) {
        self.transactionID=transactionID; self.owner=owner; self.previousCheckpoint=previousCheckpoint
        self.candidateCheckpoint=candidateCheckpoint; self.acceptedTimeSeconds=acceptedTimeSeconds
        self.converged=converged; self.residual=residual
    }
}
/// Prepared states remain private until both owners accept the same time and
/// lineage. A failed release must quarantine the experiment. Rollback requires
/// both native owners to verify restoration of exact checkpoints.
public protocol VivoMechanobiologyParticipant: Actor {
    func prepare(_ proposal: VivoMechanobiologyProposal) async throws -> VivoMechanobiologyPreparedState
    func release(transactionID: UUID, peer: VivoMechanobiologyPreparedState) async throws
    func abort(transactionID: UUID) async throws -> String
}
