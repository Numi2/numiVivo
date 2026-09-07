import Foundation

/// Refinement meaning is more specific than the downstream qualification
/// dimension. Fragment and bath changes are electronic-model variants at the
/// conditional-rate boundary; QM-region changes retain the dedicated dimension.
public enum VivoAdaptiveQMMMModelSpaceKind: String, Codable, Sendable {
    case activeElectronicSpace
    case embeddingFragment
    case embeddingBath
    case qmRegion
    case electronicModel

    var qualificationDimension: VivoQMMMQualificationDimension {
        switch self {
        case .qmRegion: return .qmRegion
        case .activeElectronicSpace,.embeddingFragment,.embeddingBath,.electronicModel: return .electronicModel
        }
    }
    var requiresNestedComponents: Bool {
        switch self {
        case .activeElectronicSpace,.embeddingFragment,.embeddingBath,.qmRegion: return true
        case .electronicModel: return false
        }
    }
}

public enum VivoAdaptiveModelSpaceRole: String, Codable, Sendable { case discovery, confirmation }

/// One complete independently replicated candidate Hamiltonian/rate protocol.
/// componentIdentifiers are semantic orbital/fragment/bath/atom identifiers used
/// only to prove monotone nesting; the scientific calculation remains `variant`.
public struct VivoAdaptiveQMMMModelSpaceCandidate: Codable, Sendable, Equatable {
    public var identifier: String
    public var kind: VivoAdaptiveQMMMModelSpaceKind
    public var role: VivoAdaptiveModelSpaceRole
    public var parentIdentifier: String?
    public var componentIdentifiers: [String]
    public var variant: VivoQMMMQualificationVariant
    public var costClass: String
    public var declaredWorkUnits: Int
    public var initialEstimatedSeconds: Double
    public var expectedMetricReduction: Double
    public init(identifier: String,kind: VivoAdaptiveQMMMModelSpaceKind,role: VivoAdaptiveModelSpaceRole = .discovery,
                parentIdentifier: String? = nil,componentIdentifiers: [String],variant: VivoQMMMQualificationVariant,
                costClass: String,declaredWorkUnits: Int,initialEstimatedSeconds: Double,
                expectedMetricReduction: Double) {
        self.identifier=identifier; self.kind=kind; self.role=role; self.parentIdentifier=parentIdentifier
        self.componentIdentifiers=componentIdentifiers; self.variant=variant; self.costClass=costClass
        self.declaredWorkUnits=declaredWorkUnits; self.initialEstimatedSeconds=initialEstimatedSeconds
        self.expectedMetricReduction=expectedMetricReduction
    }
}

public struct VivoAdaptiveQMMMModelSpaceFamily: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/adaptive-qmmm-model-space-family/v1"
    public var schema: String
    public var baselineRequest: VivoQMMMReplicatedFreeEnergyRateRequest
    public var baselineResult: VivoQMMMReplicatedFreeEnergyRateResult
    public var criterionIdentifier: String
    public var metricUnit: String
    public var observableFingerprint: String
    public var candidates: [VivoAdaptiveQMMMModelSpaceCandidate]
    public var standardDeviationMultiplier: Double
    public init(baselineRequest: VivoQMMMReplicatedFreeEnergyRateRequest,
                baselineResult: VivoQMMMReplicatedFreeEnergyRateResult,
                criterionIdentifier: String,metricUnit: String = "statistically-guarded-absolute-log-rate-shift",
                observableFingerprint: String,candidates: [VivoAdaptiveQMMMModelSpaceCandidate],
                standardDeviationMultiplier: Double = 2) {
        schema=Self.schema; self.baselineRequest=baselineRequest; self.baselineResult=baselineResult
        self.criterionIdentifier=criterionIdentifier; self.metricUnit=metricUnit; self.observableFingerprint=observableFingerprint
        self.candidates=candidates; self.standardDeviationMultiplier=standardDeviationMultiplier
    }

    public func validate() throws {
        guard schema==Self.schema,!criterionIdentifier.isEmpty,criterionIdentifier.utf8.count<=256,
              metricUnit=="statistically-guarded-absolute-log-rate-shift",
              observableFingerprint.count==64,observableFingerprint.utf8.allSatisfy({ (48...57).contains($0) || (97...102).contains($0) }),
              (1...128).contains(candidates.count),Set(candidates.map(\.identifier)).count==candidates.count,
              standardDeviationMultiplier.isFinite,standardDeviationMultiplier>=0,standardDeviationMultiplier<=10 else {
            throw VivoChemistryError.invalid("adaptive QM/MM model-space family identity")
        }
        try VivoQMMMReplicatedFreeEnergyRate.validate(baselineResult,request:baselineRequest)
        guard baselineResult.converged,
              observableFingerprint == (try VivoAdaptiveMetricContext.rate(baselineRequest).hex) else {
            throw VivoChemistryError.invalid("adaptive QM/MM model-space baseline or observable binding")
        }
        let byID=Dictionary(uniqueKeysWithValues:candidates.map{($0.identifier,$0)})
        for candidate in candidates {
            try candidate.variant.validate()
            guard !candidate.identifier.isEmpty,candidate.identifier.utf8.count<=256,
                  candidate.variant.dimension==candidate.kind.qualificationDimension,
                  !candidate.costClass.isEmpty,candidate.costClass.utf8.count<=256,
                  candidate.declaredWorkUnits>0,candidate.declaredWorkUnits<=1_000_000_000_000,
                  candidate.initialEstimatedSeconds.isFinite,candidate.initialEstimatedSeconds>0,
                  candidate.expectedMetricReduction.isFinite,candidate.expectedMetricReduction>=0,
                  !candidate.componentIdentifiers.isEmpty,
                  Set(candidate.componentIdentifiers).count==candidate.componentIdentifiers.count,
                  candidate.componentIdentifiers.allSatisfy({ !$0.isEmpty && $0.utf8.count<=512 }) else {
                throw VivoChemistryError.invalid("adaptive QM/MM model-space candidate")
            }
            _ = try VivoQMMMVariantSensitivity.calculate(.init(baselineRequest:baselineRequest,baselineResult:baselineResult,
                variant:candidate.variant,standardDeviationMultiplier:standardDeviationMultiplier))
            if let parentID=candidate.parentIdentifier {
                guard let parent=byID[parentID],parent.kind==candidate.kind,parent.identifier != candidate.identifier else {
                    throw VivoChemistryError.invalid("adaptive model-space parent identity/kind")
                }
                if candidate.kind.requiresNestedComponents {
                    let a=Set(parent.componentIdentifiers),b=Set(candidate.componentIdentifiers)
                    guard a.isStrictSubset(of:b) else {
                        throw VivoChemistryError.invalid("adaptive model-space expansion must be a strict component superset")
                    }
                }
            }
        }
        for candidate in candidates {
            var seen=Set<String>(),cursor:VivoAdaptiveQMMMModelSpaceCandidate?=candidate
            while let value=cursor,let parent=value.parentIdentifier {
                guard seen.insert(parent).inserted,let next=byID[parent] else {
                    throw VivoChemistryError.invalid("adaptive model-space parent cycle")
                }
                cursor=next
            }
            if candidate.role == .confirmation {
                guard !candidates.contains(where:{$0.parentIdentifier==candidate.identifier}) else {
                    throw VivoChemistryError.invalid("held-out model-space confirmation cannot be a refinement parent")
                }
            }
        }
    }

    /// Proposals plug directly into VivoAdaptiveCampaign. Each workflow recipe
    /// must independently produce the matching native metric evidence; this
    /// function does not fabricate a metric from candidate metadata.
    public func proposals() throws -> [VivoRefinementActionProposal] {
        try validate()
        return candidates.map { candidate in
            VivoRefinementActionProposal(identifier:candidate.identifier,criterionIdentifier:criterionIdentifier,
                role:candidate.role == .discovery ? .discovery:.confirmation,
                prerequisites:candidate.parentIdentifier.map{[$0]} ?? [],costClass:candidate.costClass,
                declaredWorkUnits:candidate.declaredWorkUnits,initialEstimatedSeconds:candidate.initialEstimatedSeconds,
                expectedMetricReduction:candidate.expectedMetricReduction)
        }
    }

    public func sensitivityRequest(for identifier:String) throws -> VivoQMMMVariantSensitivityRequest {
        try validate()
        guard let candidate=candidates.first(where:{$0.identifier==identifier}) else {
            throw VivoChemistryError.invalid("unknown adaptive QM/MM model-space candidate")
        }
        return .init(baselineRequest:baselineRequest,baselineResult:baselineResult,variant:candidate.variant,
                     standardDeviationMultiplier:standardDeviationMultiplier)
    }

    /// Direct bridge to the already allowlisted chemical-sensitivity workflow.
    /// The tolerance here only controls that workflow's own `converged` flag;
    /// adaptive campaign acceptance still uses its separately declared criterion.
    public func qualificationRequest(for identifier:String,
                                     maximumAbsoluteLogRateShift:Double=1e12) throws -> VivoQMMMChemicalQualificationRequest {
        try validate()
        guard maximumAbsoluteLogRateShift.isFinite,maximumAbsoluteLogRateShift>0,
              let candidate=candidates.first(where:{$0.identifier==identifier}) else {
            throw VivoChemistryError.invalid("adaptive model-space qualification request")
        }
        return .init(identifier:"adaptive-model-space:"+candidate.identifier,
            baselineRequest:baselineRequest,baselineResult:baselineResult,variants:[candidate.variant],
            criteria:[.init(dimension:candidate.variant.dimension,requiredVariants:1,
                            maximumAbsoluteLogRateShift:maximumAbsoluteLogRateShift)],
            externalTargets:[],maximumDirectValidationAbsoluteLogError:1e12)
    }
}
