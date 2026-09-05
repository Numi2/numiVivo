import Foundation

public enum VivoReactionAtlasStateRole: String, Codable, Sendable {
    case reactant, precomplex, transitionState, intermediate, product
}
public struct VivoReactionAtlasState: Codable, Sendable, Equatable {
    public let identifier: String
    public let structureFingerprint: VivoFingerprint
    public let molecularCharge: Int
    public let spinMultiplicity: Int
    public let role: VivoReactionAtlasStateRole
    public let labels: [String:String]
    public init(identifier: String, structureFingerprint: VivoFingerprint, molecularCharge: Int,
                spinMultiplicity: Int, role: VivoReactionAtlasStateRole, labels: [String:String] = [:]) {
        self.identifier=identifier;self.structureFingerprint=structureFingerprint;self.molecularCharge=molecularCharge
        self.spinMultiplicity=spinMultiplicity;self.role=role;self.labels=labels
    }
}
public struct VivoReactionAtlasRate: Codable, Sendable, Equatable {
    public let request: VivoTransitionStateRateRequest
    public let estimate: VivoTransitionStateRateEstimate
}
public struct VivoReactionAtlasEdge: Codable, Sendable, Equatable {
    public let identifier: String
    public let reactantStateIdentifiers: [String]
    public let productStateIdentifiers: [String]
    public let mechanism: String
    public let barrier: VivoReactionBarrier?
    public let rate: VivoReactionAtlasRate?
    public let evidenceFingerprints: [VivoFingerprint]
    public let confidence: Double?
    public let reversible: Bool
    public init(identifier: String, reactantStateIdentifiers: [String], productStateIdentifiers: [String],
                mechanism: String, barrier: VivoReactionBarrier? = nil, rate: VivoReactionAtlasRate? = nil,
                evidenceFingerprints: [VivoFingerprint] = [], confidence: Double? = nil, reversible: Bool = false) {
        self.identifier=identifier;self.reactantStateIdentifiers=reactantStateIdentifiers.sorted()
        self.productStateIdentifiers=productStateIdentifiers.sorted();self.mechanism=mechanism
        self.barrier=barrier;self.rate=rate;self.evidenceFingerprints=evidenceFingerprints.sorted{$0.hex<$1.hex}
        self.confidence=confidence;self.reversible=reversible
    }
}
public struct VivoReactionAtlas: Codable, Sendable, Equatable {
    public static let schema="numivivo.org/reaction-atlas/v1"
    public let schema:String
    public let identifier:String
    public let states:[VivoReactionAtlasState]
    public let edges:[VivoReactionAtlasEdge]
    public let provenance:[String:String]
    public init(identifier:String,states:[VivoReactionAtlasState],edges:[VivoReactionAtlasEdge],provenance:[String:String]=[:]) {
        schema=Self.schema;self.identifier=identifier;self.states=states.sorted{$0.identifier<$1.identifier}
        self.edges=edges.sorted{$0.identifier<$1.identifier};self.provenance=provenance
    }
    public func validate() throws {
        guard schema==Self.schema,!identifier.isEmpty,!states.isEmpty,states.count<=1_000_000,edges.count<=4_000_000,
              Set(states.map(\.identifier)).count==states.count,Set(edges.map(\.identifier)).count==edges.count else {
            throw VivoChemistryError.invalid("reaction atlas schema/identity/capacity")
        }
        let ids=Set(states.map(\.identifier))
        for state in states {
            guard !state.identifier.isEmpty,state.identifier.utf8.count<=1024,state.spinMultiplicity>0,
                  state.labels.count<=128,state.labels.allSatisfy({$0.key.utf8.count<=256 && $0.value.utf8.count<=4096}) else {
                throw VivoChemistryError.invalid("reaction atlas state")
            }
        }
        for edge in edges {
            guard !edge.identifier.isEmpty,!edge.mechanism.isEmpty,!edge.reactantStateIdentifiers.isEmpty,
                  !edge.productStateIdentifiers.isEmpty,Set(edge.reactantStateIdentifiers).count==edge.reactantStateIdentifiers.count,
                  Set(edge.productStateIdentifiers).count==edge.productStateIdentifiers.count,
                  Set(edge.reactantStateIdentifiers).isSubset(of:ids),Set(edge.productStateIdentifiers).isSubset(of:ids),
                  edge.confidence?.isFinite != false,edge.confidence.map({(0...1).contains($0)}) != false else {
                throw VivoChemistryError.invalid("reaction atlas edge")
            }
            try edge.barrier?.validate()
            if let rate=edge.rate {
                let expected=try VivoTransitionStateRateEstimator.estimate(rate.request)
                guard expected==rate.estimate,rate.estimate.ratePerSecond.isFinite,rate.estimate.ratePerSecond>0,
                      rate.estimate.naturalLogRatePerSecond.isFinite else { throw VivoChemistryError.invalid("atlas rate estimate/request mismatch") }
                guard edge.reactantStateIdentifiers.count==1,edge.productStateIdentifiers.count==1,
                      let barrier=edge.barrier,barrier.precomplex.context.energyKind == .gibbs else {
                    throw VivoChemistryError.invalid("unimolecular TST rate may only annotate a 1->1 edge with its Gibbs barrier")
                }
            }
        }
    }
    public func fingerprint() throws -> VivoFingerprint {
        try validate();return try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(self))
    }
}

public enum VivoReactionAtlasBuilder {
    public static func calculatedUnimolecularEdge(identifier:String,reactantState:String,productState:String,
                                                  mechanism:String,barrier:VivoReactionBarrier,
                                                  kineticContext:VivoKineticContext,barrierEvidence:VivoKineticEvidence,
                                                  transmissionProbability:Double=1,
                                                  transmissionOrigin:VivoKineticOrigin,
                                                  transmissionEvidence:VivoKineticEvidence,
                                                  conditionalBarrierStandardDeviationKJPerMol:Double?=nil,
                                                  evidenceFingerprints:[VivoFingerprint]=[],confidence:Double?=nil) throws -> VivoReactionAtlasEdge {
        let activation=try VivoReactionEnergies.activationBarrier(barrier,kineticContext:kineticContext,
            evidence:barrierEvidence,conditionalStandardDeviationKJPerMol:conditionalBarrierStandardDeviationKJPerMol)
        let request=VivoTransitionStateRateRequest(barrier:activation,transmissionProbability:transmissionProbability,
            transmissionOrigin:transmissionOrigin,transmissionEvidence:transmissionEvidence)
        let estimate=try VivoTransitionStateRateEstimator.estimate(request)
        return .init(identifier:identifier,reactantStateIdentifiers:[reactantState],productStateIdentifiers:[productState],
                     mechanism:mechanism,barrier:barrier,rate:.init(request:request,estimate:estimate),
                     evidenceFingerprints:evidenceFingerprints,confidence:confidence,reversible:false)
    }
}

public enum VivoReactionAtlasKineticsAdapter {
    public static func inactivationParameter(from edge:VivoReactionAtlasEdge,
                                             evidence:VivoKineticEvidence) throws -> VivoKineticParameter {
        guard edge.reactantStateIdentifiers.count==1,edge.productStateIdentifiers.count==1,let rate=edge.rate,
              rate.estimate.ratePerSecond.isFinite,rate.estimate.ratePerSecond>0 else {
            throw VivoKineticsError.invalid("atlas edge is not a characterized unimolecular inactivation step")
        }
        try evidence.validate(origin:.calculated)
        let uncertainty:VivoKineticUncertainty
        if let sd=rate.estimate.conditionalLogRateStandardDeviation,sd>0 { uncertainty = .logNormal(logStandardDeviation:sd) }
        else { uncertainty = .unknown }
        return .init(value:rate.estimate.ratePerSecond,unit:.perSecond,origin:.calculated,
                     uncertainty:uncertainty,evidence:evidence)
    }

    public static func covalentPack(identifier:String,edge:VivoReactionAtlasEdge,
                                    association:VivoKineticParameter,dissociation:VivoKineticParameter,
                                    targetTurnover:VivoKineticParameter,baselineTarget:VivoKineticParameter,
                                    maximumUnboundDrugM:Double,inactivationEvidence:VivoKineticEvidence,
                                    competitor:VivoKineticCompetitor?=nil) throws -> VivoCovalentKineticPack {
        guard let context=edge.rate?.request.barrier.context else { throw VivoKineticsError.invalid("atlas edge has no kinetic context") }
        let inactivation=try inactivationParameter(from:edge,evidence:inactivationEvidence)
        let pack=VivoCovalentKineticPack(identifier:identifier,context:context,association:association,dissociation:dissociation,
            inactivation:inactivation,targetTurnover:targetTurnover,baselineTarget:baselineTarget,
            competitor:competitor,maximumUnboundDrugM:maximumUnboundDrugM)
        try pack.validate();return pack
    }
}
