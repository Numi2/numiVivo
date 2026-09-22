import Foundation

/// Minimum independently sourced biological coverage required before a
/// contextual response learner may fit a residual. These are admission gates,
/// not model hyperparameters. A matched context or a technical view never
/// counts as another biological unit.
public struct VivoCellResponseCohortRequirements: Codable, Sendable, Equatable {
    public let minimumTrainingBiologicalUnitsPerTarget: Int
    public let minimumValidationBiologicalUnitsPerTarget: Int
    public let minimumTestBiologicalUnitsPerTarget: Int
    public let minimumBiologicalUnits: Int

    public init(minimumTrainingBiologicalUnitsPerTarget: Int = 8,
                minimumValidationBiologicalUnitsPerTarget: Int = 2,
                minimumTestBiologicalUnitsPerTarget: Int = 2,
                minimumBiologicalUnits: Int = 3) {
        self.minimumTrainingBiologicalUnitsPerTarget = minimumTrainingBiologicalUnitsPerTarget
        self.minimumValidationBiologicalUnitsPerTarget = minimumValidationBiologicalUnitsPerTarget
        self.minimumTestBiologicalUnitsPerTarget = minimumTestBiologicalUnitsPerTarget
        self.minimumBiologicalUnits = minimumBiologicalUnits
    }

    private enum CodingKeys: String, CodingKey {
        case minimumTrainingBiologicalUnitsPerTarget, minimumValidationBiologicalUnitsPerTarget,
             minimumTestBiologicalUnitsPerTarget, minimumBiologicalUnits
    }

    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: [
            "minimumTrainingBiologicalUnitsPerTarget", "minimumValidationBiologicalUnitsPerTarget",
            "minimumTestBiologicalUnitsPerTarget", "minimumBiologicalUnits"
        ])
        let values = try decoder.container(keyedBy: CodingKeys.self)
        minimumTrainingBiologicalUnitsPerTarget = try values.decodeIfPresent(Int.self, forKey: .minimumTrainingBiologicalUnitsPerTarget) ?? 8
        minimumValidationBiologicalUnitsPerTarget = try values.decodeIfPresent(Int.self, forKey: .minimumValidationBiologicalUnitsPerTarget) ?? 2
        minimumTestBiologicalUnitsPerTarget = try values.decodeIfPresent(Int.self, forKey: .minimumTestBiologicalUnitsPerTarget) ?? 2
        minimumBiologicalUnits = try values.decodeIfPresent(Int.self, forKey: .minimumBiologicalUnits) ?? 3
    }

    public func validate() throws {
        guard (1...1_024).contains(minimumTrainingBiologicalUnitsPerTarget),
              (1...1_024).contains(minimumValidationBiologicalUnitsPerTarget),
              (1...1_024).contains(minimumTestBiologicalUnitsPerTarget),
              (3...1_024).contains(minimumBiologicalUnits) else {
            throw VivoOmicsError.invalid("cell-response cohort requirements")
        }
    }
}

/// One receipt-verified prepared corpus and the immutable original source from
/// which it was made. The corpus binds plan, row assignments, and source-store
/// receipt; the source scopes biological identities across separately prepared
/// views of the same raw material.
public struct VivoCellResponseCohortSource: Codable, Sendable, Equatable {
    public let corpus: VivoFingerprint
    public let plan: VivoFingerprint
    public let source: VivoFingerprint

    private enum CodingKeys: String, CodingKey { case corpus, plan, source }

    public init(corpus: VivoFingerprint, plan: VivoFingerprint, source: VivoFingerprint) {
        self.corpus = corpus
        self.plan = plan
        self.source = source
    }

    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["corpus", "plan", "source"])
        let values = try decoder.container(keyedBy: CodingKeys.self)
        corpus = try values.decode(VivoFingerprint.self, forKey: .corpus)
        plan = try values.decode(VivoFingerprint.self, forKey: .plan)
        source = try values.decode(VivoFingerprint.self, forKey: .source)
    }
}

/// Exact treated samples for one target in a matched context. The admission
/// retains these rows instead of collapsing them to target labels, so a future
/// learner can replay precisely the treated/control membership it was given.
public struct VivoCellResponseCohortTargetSamples: Codable, Sendable, Equatable {
    public let targetID: String
    public let sampleIDs: [String]

    private enum CodingKeys: String, CodingKey { case targetID, sampleIDs }

    public init(targetID: String, sampleIDs: [String]) {
        self.targetID = targetID
        self.sampleIDs = sampleIDs
    }

    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["targetID", "sampleIDs"])
        let values = try decoder.container(keyedBy: CodingKeys.self)
        targetID = try values.decode(String.self, forKey: .targetID)
        sampleIDs = try values.decode([String].self, forKey: .sampleIDs)
    }
}

/// One source-bound matched treated/control context. donorID, when present,
/// is the biological split authority; otherwise the declared biological
/// replicate is scoped to the original raw source. `studyID` is retained as
/// provenance only, because relabeling a study must never create an apparent
/// independent biological unit. A corpus receipt is deliberately not used as
/// a biological identity because a source can be prepared into more than one
/// corpus artifact.
public struct VivoCellResponseCohortContext: Codable, Sendable, Equatable {
    public let corpus: VivoFingerprint
    public let studyID: String
    public let biologicalReplicateID: String
    public let donorID: String?
    public let contextID: String
    public let pairID: String
    public let modality: String
    public let partition: VivoCellResponsePartition
    public let treatedSamples: [VivoCellResponseCohortTargetSamples]
    public let controlSampleIDs: [String]

    private enum CodingKeys: String, CodingKey {
        case corpus, studyID, biologicalReplicateID, donorID, contextID, pairID, modality,
             partition, treatedSamples, controlSampleIDs
    }

    public init(corpus: VivoFingerprint, studyID: String, biologicalReplicateID: String,
                donorID: String?, contextID: String, pairID: String, modality: String,
                partition: VivoCellResponsePartition,
                treatedSamples: [VivoCellResponseCohortTargetSamples],
                controlSampleIDs: [String]) {
        self.corpus = corpus
        self.studyID = studyID
        self.biologicalReplicateID = biologicalReplicateID
        self.donorID = donorID
        self.contextID = contextID
        self.pairID = pairID
        self.modality = modality
        self.partition = partition
        self.treatedSamples = treatedSamples
        self.controlSampleIDs = controlSampleIDs
    }

    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: [
            "corpus", "studyID", "biologicalReplicateID", "donorID", "contextID", "pairID", "modality",
            "partition", "treatedSamples", "controlSampleIDs"
        ])
        let values = try decoder.container(keyedBy: CodingKeys.self)
        corpus = try values.decode(VivoFingerprint.self, forKey: .corpus)
        studyID = try values.decode(String.self, forKey: .studyID)
        biologicalReplicateID = try values.decode(String.self, forKey: .biologicalReplicateID)
        donorID = try values.decodeIfPresent(String.self, forKey: .donorID)
        contextID = try values.decode(String.self, forKey: .contextID)
        pairID = try values.decode(String.self, forKey: .pairID)
        modality = try values.decode(String.self, forKey: .modality)
        partition = try values.decode(VivoCellResponsePartition.self, forKey: .partition)
        treatedSamples = try values.decode([VivoCellResponseCohortTargetSamples].self, forKey: .treatedSamples)
        controlSampleIDs = try values.decode([String].self, forKey: .controlSampleIDs)
    }
}

/// Data-volume diagnostics for a target and split. matchedContexts is shown
/// explicitly, but only biologicalUnits can satisfy admission requirements.
public struct VivoCellResponseCohortCoverage: Codable, Sendable, Equatable {
    public let targetID: String
    public let partition: VivoCellResponsePartition
    public let biologicalUnits: Int
    public let matchedContexts: Int

    private enum CodingKeys: String, CodingKey { case targetID, partition, biologicalUnits, matchedContexts }

    public init(targetID: String, partition: VivoCellResponsePartition,
                biologicalUnits: Int, matchedContexts: Int) {
        self.targetID = targetID
        self.partition = partition
        self.biologicalUnits = biologicalUnits
        self.matchedContexts = matchedContexts
    }

    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["targetID", "partition", "biologicalUnits", "matchedContexts"])
        let values = try decoder.container(keyedBy: CodingKeys.self)
        targetID = try values.decode(String.self, forKey: .targetID)
        partition = try values.decode(VivoCellResponsePartition.self, forKey: .partition)
        biologicalUnits = try values.decode(Int.self, forKey: .biologicalUnits)
        matchedContexts = try values.decode(Int.self, forKey: .matchedContexts)
    }
}

/// A canonical, source-replayable biological split contract. Decoding checks
/// its structural invariants; verify against readers additionally replays all
/// memberships from live receipt-verified corpus readers.
public struct VivoCellResponseCohortAdmission: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let format: String
    public let requirements: VivoCellResponseCohortRequirements
    public let sources: [VivoCellResponseCohortSource]
    public let featureAxis: VivoCellResponseFeatureAxis
    public let descriptorSource: VivoFingerprint
    public let targets: [VivoCellResponseTarget]
    public let contexts: [VivoCellResponseCohortContext]
    public let coverage: [VivoCellResponseCohortCoverage]

    private enum CodingKeys: String, CodingKey {
        case schemaVersion, format, requirements, sources, featureAxis, descriptorSource, targets, contexts, coverage
    }

    public init(requirements: VivoCellResponseCohortRequirements,
                sources: [VivoCellResponseCohortSource],
                featureAxis: VivoCellResponseFeatureAxis,
                descriptorSource: VivoFingerprint,
                targets: [VivoCellResponseTarget],
                contexts: [VivoCellResponseCohortContext],
                coverage: [VivoCellResponseCohortCoverage]) {
        schemaVersion = 2
        format = "numivivo-cell-response-cohort-admission/v2"
        self.requirements = requirements
        self.sources = sources
        self.featureAxis = featureAxis
        self.descriptorSource = descriptorSource
        self.targets = targets
        self.contexts = contexts
        self.coverage = coverage
    }

    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: [
            "schemaVersion", "format", "requirements", "sources", "featureAxis", "descriptorSource", "targets", "contexts", "coverage"
        ])
        let values = try decoder.container(keyedBy: CodingKeys.self)
        schemaVersion = try values.decode(Int.self, forKey: .schemaVersion)
        format = try values.decode(String.self, forKey: .format)
        requirements = try values.decode(VivoCellResponseCohortRequirements.self, forKey: .requirements)
        sources = try values.decode([VivoCellResponseCohortSource].self, forKey: .sources)
        featureAxis = try values.decode(VivoCellResponseFeatureAxis.self, forKey: .featureAxis)
        descriptorSource = try values.decode(VivoFingerprint.self, forKey: .descriptorSource)
        targets = try values.decode([VivoCellResponseTarget].self, forKey: .targets)
        contexts = try values.decode([VivoCellResponseCohortContext].self, forKey: .contexts)
        coverage = try values.decode([VivoCellResponseCohortCoverage].self, forKey: .coverage)
    }
}

/// Admission gate for contextual residual learning. It binds every declared
/// row to a source receipt and refuses split leakage. It does not infer
/// biological independence from counts; qualification applies a fixed
/// coverage floor separately from caller-selected development thresholds.
public enum VivoCellResponseCohort {
    /// A composite reader preserves source row identity and canonical ordering
    /// across these independently receipt-verified corpus views.
    private static let maximumReaders = 64
    private static let maximumAssignments = 2_000_000
    private static let maximumContexts = 2_000_000
    private static let maximumSamplesPerContext = 1_048_576

    /// The non-negotiable biological floor for a model-evaluation claim.
    /// Callers may admit smaller source-bound development cohorts, but they
    /// cannot qualify an evaluation until every target has this independent
    /// coverage and the split contains at least twelve distinct units.
    public static let qualificationRequirements = VivoCellResponseCohortRequirements(
        minimumTrainingBiologicalUnitsPerTarget: 8,
        minimumValidationBiologicalUnitsPerTarget: 2,
        minimumTestBiologicalUnitsPerTarget: 2,
        minimumBiologicalUnits: 12
    )

    private struct ContextKey: Hashable {
        let corpus: String
        let studyID: String
        let biologicalReplicateID: String
        let donorID: String
        let contextID: String
        let pairID: String
        let modality: String
    }

    private struct BiologicalKey: Hashable {
        let source: String
        let identityKind: String
        let identity: String
    }

    /// A raw-source-scoped replicate label has exactly one donor declaration.
    /// Otherwise an omitted donor could silently make the same biological unit
    /// appear independent in another split.
    private struct SourceReplicateKey: Hashable {
        let source: String
        let biologicalReplicateID: String
    }

    /// Qualification must retain held-out evidence for each target in every
    /// raw source where that target was observed. Corpus receipts may be
    /// separate preparations of one raw source, so this key deliberately uses
    /// the source fingerprint rather than a prepared-corpus fingerprint.
    private struct SourceTargetKey: Hashable {
        let source: String
        let targetID: String
    }

    private struct CollectedContext {
        let corpus: VivoFingerprint
        let studyID: String
        let biologicalReplicateID: String
        let donorID: String?
        let contextID: String
        let pairID: String
        let modality: String
        var partition: VivoCellResponsePartition
        var treated: [String: Set<String>]
        var controls: Set<String>
    }

    private struct BoundReader {
        let reader: VivoCellResponseCorpusReader
        let source: VivoCellResponseCohortSource
    }

    private static func corpusFingerprint(_ reader: VivoCellResponseCorpusReader) throws -> VivoFingerprint {
        try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(reader.receipt))
    }

    private static func boundReaders(_ readers: [VivoCellResponseCorpusReader]) throws -> [BoundReader] {
        guard (1...maximumReaders).contains(readers.count) else {
            throw VivoOmicsError.invalid("cell-response cohort readers")
        }
        // Build the executable view up front. Admission must not accept a
        // source set whose axes/targets cannot later be replayed by the
        // learner's canonical global-row namespace.
        _ = try VivoCellResponseCompositeCorpusReader(readers: readers)
        let result = try readers.map { reader in
            try reader.plan.validate()
            return BoundReader(reader: reader, source: .init(corpus: try corpusFingerprint(reader),
                                                               plan: reader.receipt.plan,
                                                               source: reader.receipt.sourceStore.source))
        }.sorted { $0.source.corpus.hex < $1.source.corpus.hex }
        guard Set(result.map(\.source.corpus.hex)).count == result.count else {
            throw VivoOmicsError.invalid("cell-response cohort duplicate corpus receipt")
        }
        return result
    }

    private static func contextKey(_ context: VivoCellResponseCohortContext) -> ContextKey {
        .init(corpus: context.corpus.hex, studyID: context.studyID,
              biologicalReplicateID: context.biologicalReplicateID, donorID: context.donorID ?? "",
              contextID: context.contextID, pairID: context.pairID, modality: context.modality)
    }

    private static func contextSortKey(_ context: VivoCellResponseCohortContext) -> [String] {
        [context.corpus.hex, context.studyID, context.biologicalReplicateID, context.donorID ?? "",
         context.contextID, context.pairID, context.modality]
    }

    private static func biologicalKey(source: VivoFingerprint,
                                      biologicalReplicateID: String, donorID: String?) -> BiologicalKey {
        if let donorID {
            return .init(source: source.hex, identityKind: "donor", identity: donorID)
        }
        return .init(source: source.hex, identityKind: "biological-replicate", identity: biologicalReplicateID)
    }

    private static func coverage(targets: [VivoCellResponseTarget],
                                 sourcesByCorpus: [String: VivoCellResponseCohortSource],
                                 contexts: [VivoCellResponseCohortContext]) throws -> [VivoCellResponseCohortCoverage] {
        var units: [String: [String: Set<BiologicalKey>]] = [:]
        var matchedContexts: [String: [String: Set<ContextKey>]] = [:]
        for context in contexts {
            guard let source = sourcesByCorpus[context.corpus.hex] else {
                throw VivoOmicsError.invalid("cell-response cohort source membership")
            }
            let unit = biologicalKey(source: source.source, biologicalReplicateID: context.biologicalReplicateID,
                                     donorID: context.donorID)
            let matched = contextKey(context)
            for target in context.treatedSamples.map(\.targetID) {
                units[context.partition.rawValue, default: [:]][target, default: []].insert(unit)
                matchedContexts[context.partition.rawValue, default: [:]][target, default: []].insert(matched)
            }
        }
        return targets.flatMap { target in
            [VivoCellResponsePartition.training, .validation, .test].map { partition in
                .init(targetID: target.id, partition: partition,
                      biologicalUnits: units[partition.rawValue]?[target.id]?.count ?? 0,
                      matchedContexts: matchedContexts[partition.rawValue]?[target.id]?.count ?? 0)
            }
        }
    }

    private static func requirementsSatisfied(_ requirements: VivoCellResponseCohortRequirements,
                                              coverage: [VivoCellResponseCohortCoverage],
                                              sourcesByCorpus: [String: VivoCellResponseCohortSource],
                                              contexts: [VivoCellResponseCohortContext]) throws {
        var donorByReplicate: [SourceReplicateKey: String] = [:]
        var partitions: [BiologicalKey: Set<String>] = [:]
        for context in contexts {
            guard let source = sourcesByCorpus[context.corpus.hex] else {
                throw VivoOmicsError.invalid("cell-response cohort source membership")
            }
            let replicate = SourceReplicateKey(source: source.source.hex,
                                                biologicalReplicateID: context.biologicalReplicateID)
            let donor = context.donorID ?? ""
            if let existing = donorByReplicate[replicate], existing != donor {
                throw VivoOmicsError.invalid("cell-response cohort inconsistent donor identity")
            }
            donorByReplicate[replicate] = donor
            let unit = biologicalKey(source: source.source, biologicalReplicateID: context.biologicalReplicateID,
                                     donorID: context.donorID)
            partitions[unit, default: []].insert(context.partition.rawValue)
        }
        guard partitions.count >= requirements.minimumBiologicalUnits,
              partitions.values.allSatisfy({ $0.count == 1 }) else {
            throw VivoOmicsError.invalid("cell-response cohort biological unit split")
        }
        for entry in coverage {
            let minimum: Int
            switch entry.partition {
            case .training: minimum = requirements.minimumTrainingBiologicalUnitsPerTarget
            case .validation: minimum = requirements.minimumValidationBiologicalUnitsPerTarget
            case .test: minimum = requirements.minimumTestBiologicalUnitsPerTarget
            }
            guard entry.biologicalUnits >= minimum else {
                throw VivoOmicsError.invalid("cell-response cohort target biological coverage")
            }
        }
    }

    private static func qualificationSourceTargetTestCoverage(
        sourcesByCorpus: [String: VivoCellResponseCohortSource],
        contexts: [VivoCellResponseCohortContext]
    ) throws {
        var eligible: Set<SourceTargetKey> = []
        var heldOutUnits: [SourceTargetKey: Set<BiologicalKey>] = [:]
        for context in contexts {
            guard let source = sourcesByCorpus[context.corpus.hex] else {
                throw VivoOmicsError.invalid("cell-response cohort source membership")
            }
            let unit = biologicalKey(source: source.source,
                                     biologicalReplicateID: context.biologicalReplicateID,
                                     donorID: context.donorID)
            for targetID in context.treatedSamples.map(\.targetID) {
                let key = SourceTargetKey(source: source.source.hex, targetID: targetID)
                eligible.insert(key)
                if context.partition == .test {
                    heldOutUnits[key, default: []].insert(unit)
                }
            }
        }
        guard eligible.allSatisfy({
            (heldOutUnits[$0]?.count ?? 0) >= qualificationRequirements.minimumTestBiologicalUnitsPerTarget
        }) else {
            throw VivoOmicsError.invalid("cell-response qualification lacks held-out source target coverage")
        }
    }

    /// Create a canonical admission only from live receipt-verified readers.
    /// Every declared treated and control sample must resolve to at least one
    /// source row before it can enter the contract.
    public static func admit(readers: [VivoCellResponseCorpusReader],
                             requirements: VivoCellResponseCohortRequirements = .init()) throws -> VivoCellResponseCohortAdmission {
        try requirements.validate()
        let bound = try boundReaders(readers)
        guard let first = bound.first else { throw VivoOmicsError.invalid("cell-response cohort readers") }
        let expectedAxis = first.reader.plan.featureAxis
        let expectedDescriptor = first.reader.plan.descriptorSource
        let expectedTargets = first.reader.plan.targets
        guard expectedTargets.map(\.id) == expectedTargets.map(\.id).sorted() else {
            throw VivoOmicsError.invalid("cell-response cohort target order")
        }
        var collected: [ContextKey: CollectedContext] = [:]
        var assignments = 0

        for input in bound {
            let plan = input.reader.plan
            guard plan.featureAxis == expectedAxis,
                  plan.descriptorSource == expectedDescriptor,
                  plan.targets == expectedTargets else {
                throw VivoOmicsError.invalid("cell-response cohort axes or descriptors differ")
            }
            for assignment in plan.assignments {
                assignments += 1
                guard assignments <= maximumAssignments else { throw VivoOmicsError.limit("cell-response cohort assignments") }
                guard !(try input.reader.rows(forSampleID: assignment.sampleID)).isEmpty else {
                    throw VivoOmicsError.invalid("cell-response cohort assignment has no source rows")
                }
                let key = ContextKey(corpus: input.source.corpus.hex, studyID: assignment.studyID,
                                     biologicalReplicateID: assignment.sourceSample.biologicalReplicateID,
                                     donorID: assignment.sourceSample.donorID ?? "", contextID: assignment.contextID,
                                     pairID: assignment.pairID, modality: assignment.modality)
                if var existing = collected[key] {
                    guard existing.partition == assignment.partition else {
                        throw VivoOmicsError.invalid("cell-response cohort splits one matched context")
                    }
                    switch assignment.role {
                    case .control:
                        existing.controls.insert(assignment.sampleID)
                    case .perturbed:
                        guard let targetID = assignment.targetID else {
                            throw VivoOmicsError.invalid("cell-response cohort perturbed target")
                        }
                        existing.treated[targetID, default: []].insert(assignment.sampleID)
                    }
                    guard existing.controls.count + existing.treated.values.reduce(0, { $0 + $1.count }) <= maximumSamplesPerContext else {
                        throw VivoOmicsError.limit("cell-response cohort context samples")
                    }
                    collected[key] = existing
                } else {
                    var treated: [String: Set<String>] = [:]
                    var controls: Set<String> = []
                    switch assignment.role {
                    case .control:
                        controls.insert(assignment.sampleID)
                    case .perturbed:
                        guard let targetID = assignment.targetID else {
                            throw VivoOmicsError.invalid("cell-response cohort perturbed target")
                        }
                        treated[targetID] = [assignment.sampleID]
                    }
                    collected[key] = .init(corpus: input.source.corpus, studyID: assignment.studyID,
                                           biologicalReplicateID: assignment.sourceSample.biologicalReplicateID,
                                           donorID: assignment.sourceSample.donorID, contextID: assignment.contextID,
                                           pairID: assignment.pairID, modality: assignment.modality,
                                           partition: assignment.partition, treated: treated, controls: controls)
                    guard collected.count <= maximumContexts else { throw VivoOmicsError.limit("cell-response cohort contexts") }
                }
            }
        }

        let contexts = try collected.values.map { value -> VivoCellResponseCohortContext in
            guard !value.treated.isEmpty, !value.controls.isEmpty else {
                throw VivoOmicsError.invalid("cell-response cohort context lacks matched rows")
            }
            let treated = value.treated.keys.sorted().map { targetID in
                VivoCellResponseCohortTargetSamples(targetID: targetID, sampleIDs: value.treated[targetID]!.sorted())
            }
            return .init(corpus: value.corpus, studyID: value.studyID,
                         biologicalReplicateID: value.biologicalReplicateID, donorID: value.donorID,
                         contextID: value.contextID, pairID: value.pairID, modality: value.modality,
                         partition: value.partition, treatedSamples: treated,
                         controlSampleIDs: value.controls.sorted())
        }.sorted { contextSortKey($0).lexicographicallyPrecedes(contextSortKey($1)) }
        let sources = bound.map(\.source)
        let sourceMap = Dictionary(uniqueKeysWithValues: sources.map { ($0.corpus.hex, $0) })
        let diagnostics = try coverage(targets: expectedTargets, sourcesByCorpus: sourceMap, contexts: contexts)
        try requirementsSatisfied(requirements, coverage: diagnostics, sourcesByCorpus: sourceMap, contexts: contexts)
        let admission = VivoCellResponseCohortAdmission(requirements: requirements, sources: sources,
                                                         featureAxis: expectedAxis, descriptorSource: expectedDescriptor,
                                                         targets: expectedTargets, contexts: contexts, coverage: diagnostics)
        try validate(admission)
        return admission
    }

    /// Verify canonical structure without claiming that an arbitrary decoded
    /// admission names real rows. Use reader-bound verification before training
    /// or qualification to replay the source binding.
    public static func validate(_ admission: VivoCellResponseCohortAdmission) throws {
        try admission.requirements.validate()
        try admission.featureAxis.validate()
        let contextsAreCanonical = admission.contexts.indices.dropFirst().allSatisfy {
            contextSortKey(admission.contexts[$0 - 1]).lexicographicallyPrecedes(contextSortKey(admission.contexts[$0]))
        }
        guard admission.schemaVersion == 2,
              admission.format == "numivivo-cell-response-cohort-admission/v2",
              !admission.sources.isEmpty, admission.sources.count <= maximumReaders,
              admission.sources.map(\.corpus.hex) == admission.sources.map(\.corpus.hex).sorted(),
              Set(admission.sources.map(\.corpus.hex)).count == admission.sources.count,
              !admission.targets.isEmpty, admission.targets.count <= 100_000,
              admission.targets.map(\.id) == admission.targets.map(\.id).sorted(),
              Set(admission.targets.map(\.id)).count == admission.targets.count,
              !admission.contexts.isEmpty, admission.contexts.count <= maximumContexts,
              contextsAreCanonical else {
            throw VivoOmicsError.invalid("cell-response cohort admission")
        }
        let descriptorCount = admission.targets.first?.descriptors.count
        for target in admission.targets { try target.validate(expectedDescriptorCount: descriptorCount) }
        let sourceMap = Dictionary(uniqueKeysWithValues: admission.sources.map { ($0.corpus.hex, $0) })
        let targetIDs = Set(admission.targets.map(\.id))
        var seenContexts: Set<ContextKey> = []
        for context in admission.contexts {
            let treatedIDs = context.treatedSamples.map(\.targetID)
            let treatedSampleIDs = context.treatedSamples.flatMap(\.sampleIDs)
            let treatedSamplesAreCanonical = context.treatedSamples.allSatisfy { samples in
                !samples.sampleIDs.isEmpty &&
                    samples.sampleIDs == samples.sampleIDs.sorted() &&
                    samples.sampleIDs.allSatisfy(vivoOmicsID)
            }
            guard sourceMap[context.corpus.hex] != nil,
                  vivoOmicsID(context.studyID), vivoOmicsID(context.biologicalReplicateID),
                  context.donorID.map(vivoOmicsID) ?? true,
                  vivoOmicsID(context.contextID), vivoOmicsID(context.pairID), vivoOmicsID(context.modality),
                  !context.treatedSamples.isEmpty, !context.controlSampleIDs.isEmpty,
                  treatedIDs == treatedIDs.sorted(), Set(treatedIDs).count == treatedIDs.count,
                  treatedIDs.allSatisfy(vivoOmicsID), Set(treatedIDs).isSubset(of: targetIDs),
                  treatedSamplesAreCanonical,
                  context.controlSampleIDs == context.controlSampleIDs.sorted(),
                  context.controlSampleIDs.allSatisfy(vivoOmicsID),
                  Set(context.controlSampleIDs).count == context.controlSampleIDs.count,
                  Set(treatedSampleIDs).count == treatedSampleIDs.count,
                  Set(treatedSampleIDs).isDisjoint(with: Set(context.controlSampleIDs)),
                  context.controlSampleIDs.count + treatedSampleIDs.count <= maximumSamplesPerContext,
                  seenContexts.insert(contextKey(context)).inserted else {
                throw VivoOmicsError.invalid("cell-response cohort admission context")
            }
        }
        let expectedCoverage = try coverage(targets: admission.targets, sourcesByCorpus: sourceMap, contexts: admission.contexts)
        guard admission.coverage == expectedCoverage else {
            throw VivoOmicsError.invalid("cell-response cohort admission coverage")
        }
        try requirementsSatisfied(admission.requirements, coverage: admission.coverage,
                                  sourcesByCorpus: sourceMap, contexts: admission.contexts)
    }

    /// Require the fixed biological floor for reporting a model evaluation as
    /// qualified. This evaluates actual replayable coverage rather than the
    /// caller-selected development thresholds stored in the admission.
    public static func requireQualificationCoverage(_ admission: VivoCellResponseCohortAdmission) throws {
        try validate(admission)
        let sourceMap = Dictionary(uniqueKeysWithValues: admission.sources.map { ($0.corpus.hex, $0) })
        try requirementsSatisfied(qualificationRequirements, coverage: admission.coverage,
                                  sourcesByCorpus: sourceMap, contexts: admission.contexts)
        try qualificationSourceTargetTestCoverage(sourcesByCorpus: sourceMap, contexts: admission.contexts)
        // Aggregate target coverage is insufficient for a composite model:
        // one study could supply all reported held-out strata while another
        // source has none. A qualified evaluation needs at least one replayable
        // treated test stratum from every receipt-bound source.
        let sourcesWithHeldOutTreatment = Set(admission.contexts.compactMap { context -> String? in
            guard context.partition == .test, !context.treatedSamples.isEmpty else { return nil }
            return context.corpus.hex
        })
        guard sourcesWithHeldOutTreatment == Set(admission.sources.map(\.corpus.hex)) else {
            throw VivoOmicsError.invalid("cell-response qualification lacks held-out source coverage")
        }
    }

    /// Recreate the admission from verified readers and require canonical exact
    /// equality. This proves every corpus receipt, plan, target, control row,
    /// and treated row still matches the declared contract.
    public static func verify(_ admission: VivoCellResponseCohortAdmission,
                              readers: [VivoCellResponseCorpusReader]) throws {
        try validate(admission)
        let rebuilt = try admit(readers: readers, requirements: admission.requirements)
        guard rebuilt == admission,
              try VivoCanonicalJSON.encode(rebuilt) == VivoCanonicalJSON.encode(admission) else {
            throw VivoOmicsError.invalid("cell-response cohort source binding")
        }
    }
}
