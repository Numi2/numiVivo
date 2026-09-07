import Foundation

public struct VivoReactiveSurrogateCoverageConfiguration: Codable, Sendable, Equatable {
    public let minimumSourceGroups: Int
    public let minimumLabelsPerGroup: Int
    public let minimumOverallEligibleFraction: Double
    public let minimumPerGroupEligibleFraction: Double
    public let maximumPrimitiveWork: Int
    public init(minimumSourceGroups: Int = 2, minimumLabelsPerGroup: Int = 8,
                minimumOverallEligibleFraction: Double = 0.95,
                minimumPerGroupEligibleFraction: Double = 0.8,
                maximumPrimitiveWork: Int = 100_000_000) {
        self.minimumSourceGroups=minimumSourceGroups;self.minimumLabelsPerGroup=minimumLabelsPerGroup
        self.minimumOverallEligibleFraction=minimumOverallEligibleFraction
        self.minimumPerGroupEligibleFraction=minimumPerGroupEligibleFraction
        self.maximumPrimitiveWork=maximumPrimitiveWork
    }
    public func validate() throws {
        guard (2...1024).contains(minimumSourceGroups),(2...1_000_000).contains(minimumLabelsPerGroup),
              minimumOverallEligibleFraction.isFinite,(0...1).contains(minimumOverallEligibleFraction),
              minimumPerGroupEligibleFraction.isFinite,(0...1).contains(minimumPerGroupEligibleFraction),
              maximumPrimitiveWork > 0 else {
            throw VivoChemistryError.invalid("reactive surrogate coverage criteria")
        }
    }
}

public struct VivoReactiveSurrogateCoverageRequest: Codable, Sendable, Equatable {
    public let campaignIdentifier: String
    public let model: VivoReactiveSurrogateModel
    public let labels: [VivoReactiveTrainingLabel]
    public let evidenceSourceFingerprint: VivoFingerprint
    public let independenceDeclaration: String
    public let domainDescription: String
    public let configuration: VivoReactiveSurrogateCoverageConfiguration
    public init(campaignIdentifier: String, model: VivoReactiveSurrogateModel,
                labels: [VivoReactiveTrainingLabel], evidenceSourceFingerprint: VivoFingerprint,
                independenceDeclaration: String, domainDescription: String,
                configuration: VivoReactiveSurrogateCoverageConfiguration = .init()) {
        self.campaignIdentifier=campaignIdentifier;self.model=model;self.labels=labels
        self.evidenceSourceFingerprint=evidenceSourceFingerprint
        self.independenceDeclaration=independenceDeclaration;self.domainDescription=domainDescription
        self.configuration=configuration
    }
}

public struct VivoReactiveSurrogateCoverageGroup: Codable, Sendable, Equatable {
    public let sourceGroup: String
    public let labelCount: Int
    public let eligibleCount: Int
    public let eligibleFraction: Double
    public let eligibleEnergyRMSEKJPerMol: Double?
    public let eligibleMaximumEnergyErrorKJPerMol: Double?
    public let eligibleForceRMSEKJPerMolNM: Double?
    public let eligibleMaximumForceErrorKJPerMolNM: Double?
    public let outsideDescriptorSupportCount: Int
    public let committeeDisagreementCount: Int
    public let passed: Bool
}

public struct VivoReactiveSurrogateCoverageResult: Codable, Sendable, Equatable {
    public let campaignIdentifier: String
    public let modelFingerprint: VivoFingerprint
    public let coverageDataFingerprint: VivoFingerprint
    public let evidenceSourceFingerprint: VivoFingerprint
    public let configuration: VivoReactiveSurrogateCoverageConfiguration
    public let groups: [VivoReactiveSurrogateCoverageGroup]
    public let overallEligibleFraction: Double
    public let passed: Bool
    public let independenceDeclaration: String
    public let domainDescription: String
    public let evidenceFingerprint: VivoFingerprint
    public let interpretation: String
}

private struct VivoReactiveCoverageEvaluation {
    let label: VivoReactiveTrainingLabel
    let prediction: VivoReactiveSurrogatePrediction
    let energyError: Double
    let forceSquaredError: Double
    let forceVectorCount: Int
    let maximumForceError: Double
}

/// Evaluates a frozen model on grouped authoritative labels that were not used
/// for training or threshold selection. The caller's independence declaration
/// is retained; a fingerprint cannot authenticate how the data were collected.
public enum VivoReactiveSurrogateCoverage {
    public static let interpretation = "Grouped external-domain coverage of a frozen delta surrogate. Eligibility uses the frozen descriptor-distance and committee thresholds; energy and force errors use the authoritative complete Hamiltonian. Every source group must meet coverage and eligible-error limits so a large correlated trajectory cannot hide a failed region. The retained independence declaration is not authenticated, coverage is not posterior calibration, and passing does not make surrogate forces authoritative or establish sampling efficiency."

    public static func assess(campaignIdentifier: String, model: VivoReactiveSurrogateModel,
                              labels: [VivoReactiveTrainingLabel],
                              evidenceSourceFingerprint: VivoFingerprint,
                              independenceDeclaration: String, domainDescription: String,
                              configuration: VivoReactiveSurrogateCoverageConfiguration = .init()) throws -> VivoReactiveSurrogateCoverageResult {
        try model.validate();try configuration.validate()
        func validText(_ value:String)->Bool {
            !value.trimmingCharacters(in:.whitespacesAndNewlines).isEmpty && value.utf8.count <= 4096
        }
        guard validText(campaignIdentifier),validText(independenceDeclaration),validText(domainDescription),
              labels.count <= 1_000_000,Set(labels.map(\.identifier)).count==labels.count,
              labels.allSatisfy({ validText($0.identifier) && validText($0.sourceGroup) }) else {
            throw VivoChemistryError.invalid("reactive surrogate coverage campaign or labels")
        }
        let authority=model.payload.authorityDefinition,baseline=model.payload.baselineDefinition
        let work=Double(labels.count)*Double(model.payload.configuration.features.count)
            * Double(model.payload.centersNM.count)
            * (Double(model.payload.configuration.committeeSize)
                * (Double(authority.atomIndices.count)+1)+1)
        guard work <= Double(configuration.maximumPrimitiveWork) else {
            throw VivoChemistryError.resourceLimit("reactive coverage analysis budget")
        }
        var geometryFingerprints=Set<VivoFingerprint>()
        var sourceGroupCounts:[String:Int]=[:]
        for label in labels {
            let positions=label.authority.requestedPositionsNM
            try label.authority.validate(definition:authority,positionsNM:positions)
            try label.baseline.validate(definition:baseline,positionsNM:positions)
            guard label.baseline.requestedPositionsNM==positions else {
                throw VivoChemistryError.invalid("reactive coverage authority/baseline geometry mismatch")
            }
            let geometryID=try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(positions))
            guard geometryFingerprints.insert(geometryID).inserted else {
                throw VivoChemistryError.invalid("reactive coverage duplicate geometry")
            }
            sourceGroupCounts[label.sourceGroup,default:0] += 1
            guard sourceGroupCounts.count <= 1024 else {
                throw VivoChemistryError.resourceLimit("reactive coverage source-group budget")
            }
        }
        let dataID=try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(labels))
        let sourceGroups=sourceGroupCounts.keys.sorted()
        guard dataID != model.payload.trainingDataFingerprint,
              evidenceSourceFingerprint != model.payload.trainingDataFingerprint,
              sourceGroups.count >= configuration.minimumSourceGroups,sourceGroups.count <= 1024,
              Set(sourceGroups).isDisjoint(with:Set(model.payload.qualification.heldOutGroups)),
              sourceGroups.allSatisfy({ sourceGroupCounts[$0,default:0] >= configuration.minimumLabelsPerGroup }) else {
            throw VivoChemistryError.invalid("reactive coverage must be grouped and distinct from training/selection data")
        }
        let indexed=try labels.map { label -> VivoReactiveCoverageEvaluation in
            let prediction=try VivoReactiveDeltaSurrogate.predict(model,positionsNM:label.authority.requestedPositionsNM)
            let observedEnergy=label.authority.energyKJPerMol-label.baseline.energyKJPerMol
            let energyError=prediction.deltaEnergyKJPerMol-observedEnergy
            let forceErrors=zip(label.authority.forcesKJPerMolNM,
                zip(label.baseline.forcesKJPerMolNM,prediction.deltaForcesKJPerMolNM))
                .map { item in item.1.1-(item.0-item.1.0) }
            let forceSquaredError=forceErrors.reduce(0) { $0+$1.squaredNorm }
            let maximumForceError=forceErrors.map(\.norm).max() ?? 0
            guard energyError.isFinite,forceSquaredError.isFinite,maximumForceError.isFinite else {
                throw VivoChemistryError.convergence("reactive coverage error overflow")
            }
            return .init(label:label,prediction:prediction,energyError:energyError,
                forceSquaredError:forceSquaredError,forceVectorCount:forceErrors.count,
                maximumForceError:maximumForceError)
        }
        let cfg=model.payload.configuration
        let valuesByGroup=Dictionary(grouping:indexed) { $0.label.sourceGroup }
        var results:[VivoReactiveSurrogateCoverageGroup]=[]
        for group in sourceGroups {
            guard let values=valuesByGroup[group] else {
                throw VivoChemistryError.invalid("reactive coverage group aggregation")
            }
            let eligible=values.filter { $0.prediction.eligible }
            let fraction=Double(eligible.count)/Double(values.count)
            let energyRMSE=eligible.isEmpty ? nil : sqrt(eligible.reduce(0) {
                $0+$1.energyError*$1.energyError
            }/Double(eligible.count))
            let maximumEnergy=eligible.map { abs($0.energyError) }.max()
            let forceVectorCount=eligible.reduce(0) { $0+$1.forceVectorCount }
            let forceRMSE=forceVectorCount==0 ? nil : sqrt(eligible.reduce(0) {
                $0+$1.forceSquaredError
            }/Double(forceVectorCount))
            let maximumForce=eligible.map(\.maximumForceError).max()
            let outside=values.filter { $0.prediction.reasons.contains("outside training-descriptor support") }.count
            let disagreement=values.filter { $0.prediction.reasons.contains("committee disagreement requires authority labels") }.count
            let passed=fraction >= configuration.minimumPerGroupEligibleFraction
                && maximumEnergy.map { $0 <= cfg.maximumHeldOutEnergyErrorKJPerMol } == true
                && maximumForce.map { $0 <= cfg.maximumHeldOutForceErrorKJPerMolNM } == true
            results.append(.init(sourceGroup:group,labelCount:values.count,eligibleCount:eligible.count,
                eligibleFraction:fraction,eligibleEnergyRMSEKJPerMol:energyRMSE,
                eligibleMaximumEnergyErrorKJPerMol:maximumEnergy,eligibleForceRMSEKJPerMolNM:forceRMSE,
                eligibleMaximumForceErrorKJPerMolNM:maximumForce,outsideDescriptorSupportCount:outside,
                committeeDisagreementCount:disagreement,passed:passed))
        }
        let overall=Double(indexed.filter { $0.prediction.eligible }.count)/Double(indexed.count)
        let passed=overall >= configuration.minimumOverallEligibleFraction && results.allSatisfy(\.passed)
        struct Evidence:Codable {
            let campaignIdentifier:String;let model:VivoFingerprint;let data:VivoFingerprint
            let source:VivoFingerprint;let configuration:VivoReactiveSurrogateCoverageConfiguration
            let groups:[VivoReactiveSurrogateCoverageGroup];let overall:Double;let passed:Bool
            let independenceDeclaration:String;let domainDescription:String
        }
        let evidenceID=try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(Evidence(
            campaignIdentifier:campaignIdentifier,model:model.fingerprint,data:dataID,source:evidenceSourceFingerprint,
            configuration:configuration,groups:results,overall:overall,passed:passed,
            independenceDeclaration:independenceDeclaration,domainDescription:domainDescription)))
        return .init(campaignIdentifier:campaignIdentifier,modelFingerprint:model.fingerprint,
            coverageDataFingerprint:dataID,evidenceSourceFingerprint:evidenceSourceFingerprint,
            configuration:configuration,groups:results,overallEligibleFraction:overall,passed:passed,
            independenceDeclaration:independenceDeclaration,domainDescription:domainDescription,
            evidenceFingerprint:evidenceID,interpretation:interpretation)
    }

    public static func assess(_ request: VivoReactiveSurrogateCoverageRequest) throws -> VivoReactiveSurrogateCoverageResult {
        try assess(campaignIdentifier:request.campaignIdentifier,model:request.model,labels:request.labels,
            evidenceSourceFingerprint:request.evidenceSourceFingerprint,
            independenceDeclaration:request.independenceDeclaration,domainDescription:request.domainDescription,
            configuration:request.configuration)
    }
    public static func validate(_ result: VivoReactiveSurrogateCoverageResult,
                                request: VivoReactiveSurrogateCoverageRequest) throws {
        guard result == (try assess(request)) else {
            throw VivoChemistryError.invalid("reactive surrogate coverage evidence does not reconstruct")
        }
    }
}
