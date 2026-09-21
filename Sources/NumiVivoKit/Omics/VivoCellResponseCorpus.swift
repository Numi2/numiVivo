import Foundation

/// The role of a measured sample in a perturbation-response experiment.
/// The role is always supplied by the source plan; it is never inferred from a
/// condition label or a target identifier.
public enum VivoCellResponseRole: String, Codable, Sendable, Equatable {
    case control
    case perturbed
}

/// Frozen split membership for a source sample. A split is assigned before a
/// learner is created, so held-out outcomes cannot be selected after scoring.
public enum VivoCellResponsePartition: String, Codable, Sendable, Equatable {
    case training
    case validation
    case test
}

/// A target vocabulary entry. Descriptors are supplied by the caller and may
/// represent annotations or precomputed biological features. They are retained
/// verbatim rather than inferred from expression data.
public struct VivoCellResponseTarget: Codable, Sendable, Equatable {
    public let id: String
    public let descriptors: [Double]

    private enum CodingKeys: String, CodingKey { case id, descriptors }

    public init(id: String, descriptors: [Double]) {
        self.id = id
        self.descriptors = descriptors
    }

    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["id", "descriptors"])
        let values = try decoder.container(keyedBy: CodingKeys.self)
        id = try values.decode(String.self, forKey: .id)
        descriptors = try values.decode([Double].self, forKey: .descriptors)
    }

    /// Check the supplied identity and descriptor vector before it is used in
    /// corpus construction or a prediction query.
    public func validate(expectedDescriptorCount: Int? = nil) throws {
        guard vivoOmicsID(id), !descriptors.isEmpty, descriptors.count <= 4_096,
              descriptors.allSatisfy({ $0.isFinite && Float($0).isFinite }),
              expectedDescriptorCount == nil || descriptors.count == expectedDescriptorCount else {
            throw VivoOmicsError.invalid("cell-response target identity or descriptor vector")
        }
    }

    /// A canonical identity for this exact target and its ordered descriptor
    /// vector. It prevents an inference query from combining a fitted target
    /// embedding with altered descriptors.
    public func fingerprint() throws -> VivoFingerprint {
        try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(self))
    }
}

/// An immutable, external descriptor document. Its target vectors must exactly
/// match the corpus plan, while its own fingerprints identify the annotation,
/// sequence, or other source used to produce them.
public struct VivoCellResponseDescriptorSource: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let format: String
    public let sourceDescription: String
    public let targets: [VivoCellResponseTarget]
    public let sourceArtifacts: [VivoFingerprint]

    private enum CodingKeys: String, CodingKey { case schemaVersion, format, sourceDescription, targets, sourceArtifacts }

    public init(sourceDescription: String, targets: [VivoCellResponseTarget], sourceArtifacts: [VivoFingerprint]) {
        schemaVersion = 1
        format = "numivivo-cell-response-descriptor-source/v1"
        self.sourceDescription = sourceDescription
        self.targets = targets
        self.sourceArtifacts = sourceArtifacts
    }

    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["schemaVersion", "format", "sourceDescription", "targets", "sourceArtifacts"])
        let values = try decoder.container(keyedBy: CodingKeys.self)
        schemaVersion = try values.decode(Int.self, forKey: .schemaVersion)
        format = try values.decode(String.self, forKey: .format)
        sourceDescription = try values.decode(String.self, forKey: .sourceDescription)
        targets = try values.decode([VivoCellResponseTarget].self, forKey: .targets)
        sourceArtifacts = try values.decode([VivoFingerprint].self, forKey: .sourceArtifacts)
    }

    public func validate() throws {
        guard schemaVersion == 1, format == "numivivo-cell-response-descriptor-source/v1",
              !sourceDescription.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty,
              sourceDescription.utf8.count <= 16_384, !targets.isEmpty, targets.count <= 100_000,
              Set(targets.map(\.id)).count == targets.count, !sourceArtifacts.isEmpty,
              sourceArtifacts.count <= 100_000 else {
            throw VivoOmicsError.invalid("cell-response descriptor source")
        }
        let descriptorCount = targets.first?.descriptors.count
        for target in targets { try target.validate(expectedDescriptorCount: descriptorCount) }
    }
}

/// The subset of target vocabulary entries whose embeddings received training
/// gradients. Held-out targets deliberately have no binding and must use their
/// supplied descriptors without a learned target embedding.
public struct VivoCellResponseTrainedTargetBinding: Codable, Sendable, Equatable {
    public let id: String
    public let descriptorFingerprint: VivoFingerprint

    private enum CodingKeys: String, CodingKey { case id, descriptorFingerprint }

    public init(id: String, descriptorFingerprint: VivoFingerprint) {
        self.id = id
        self.descriptorFingerprint = descriptorFingerprint
    }

    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["id", "descriptorFingerprint"])
        let values = try decoder.container(keyedBy: CodingKeys.self)
        id = try values.decode(String.self, forKey: .id)
        descriptorFingerprint = try values.decode(VivoFingerprint.self, forKey: .descriptorFingerprint)
        guard vivoOmicsID(id) else { throw VivoOmicsError.invalid("cell-response trained target identity") }
    }
}

/// Explicit experimental identity for every source sample. `pairID` binds a
/// treated sample to controls observed in the same declared context; it is not
/// a computationally inferred cell match.
public struct VivoCellResponseAssignment: Codable, Sendable, Equatable {
    public let sampleID: String
    /// Exact source metadata copied from the receipt-verified count store.
    /// It pins the direct guide and batch labels used by this plan to the
    /// source sample rather than allowing them to be inferred from cell data.
    public let sourceSample: VivoOmicsSample
    public let role: VivoCellResponseRole
    public let targetID: String?
    public let guideID: String?
    public let modality: String
    public let studyID: String
    public let contextID: String
    public let pairID: String
    public let partition: VivoCellResponsePartition

    private enum CodingKeys: String, CodingKey {
        case sampleID, sourceSample, role, targetID, guideID, modality, studyID, contextID, pairID, partition
    }

    public init(sampleID: String, sourceSample: VivoOmicsSample, role: VivoCellResponseRole, targetID: String? = nil, guideID: String? = nil,
                modality: String, studyID: String, contextID: String, pairID: String,
                partition: VivoCellResponsePartition) {
        self.sampleID = sampleID
        self.sourceSample = sourceSample
        self.role = role
        self.targetID = targetID
        self.guideID = guideID
        self.modality = modality
        self.studyID = studyID
        self.contextID = contextID
        self.pairID = pairID
        self.partition = partition
    }

    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["sampleID", "sourceSample", "role", "targetID", "guideID", "modality", "studyID", "contextID", "pairID", "partition"])
        let values = try decoder.container(keyedBy: CodingKeys.self)
        sampleID = try values.decode(String.self, forKey: .sampleID)
        sourceSample = try values.decode(VivoOmicsSample.self, forKey: .sourceSample)
        role = try values.decode(VivoCellResponseRole.self, forKey: .role)
        targetID = try values.decodeIfPresent(String.self, forKey: .targetID)
        guideID = try values.decodeIfPresent(String.self, forKey: .guideID)
        modality = try values.decode(String.self, forKey: .modality)
        studyID = try values.decode(String.self, forKey: .studyID)
        contextID = try values.decode(String.self, forKey: .contextID)
        pairID = try values.decode(String.self, forKey: .pairID)
        partition = try values.decode(VivoCellResponsePartition.self, forKey: .partition)
    }

    func validate(targets: Set<String>) throws {
        try sourceSample.validate()
        guard let guideID, sourceSample.id == sampleID,
              guideID == sourceSample.condition, contextID == sourceSample.batchID,
              [sampleID, modality, studyID, contextID, pairID, guideID].allSatisfy(vivoOmicsID) else {
            throw VivoOmicsError.invalid("cell-response assignment identity")
        }
        switch role {
        case .control:
            guard targetID == nil else { throw VivoOmicsError.invalid("control assignment has a target") }
        case .perturbed:
            guard let targetID, targetID == guideID, vivoOmicsID(targetID), targets.contains(targetID) else {
                throw VivoOmicsError.invalid("perturbed assignment target")
            }
        }
    }
}

/// Ordered prediction axis. A missing feature stays distinct from a measured
/// zero through `missingValue`, which must be negative because log-normalized
/// measured counts are nonnegative.
public struct VivoCellResponseFeatureAxis: Codable, Sendable, Equatable {
    public let featureIDs: [String]
    public let normalizationTarget: Double
    public let missingValue: Double

    private enum CodingKeys: String, CodingKey { case featureIDs, normalizationTarget, missingValue }

    public init(featureIDs: [String], normalizationTarget: Double = 10_000, missingValue: Double = -1) {
        self.featureIDs = featureIDs
        self.normalizationTarget = normalizationTarget
        self.missingValue = missingValue
    }

    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["featureIDs", "normalizationTarget", "missingValue"])
        let values = try decoder.container(keyedBy: CodingKeys.self)
        featureIDs = try values.decode([String].self, forKey: .featureIDs)
        normalizationTarget = try values.decodeIfPresent(Double.self, forKey: .normalizationTarget) ?? 10_000
        missingValue = try values.decodeIfPresent(Double.self, forKey: .missingValue) ?? -1
    }

    public func validate() throws {
        guard !featureIDs.isEmpty, featureIDs.count <= 200_000,
              Set(featureIDs).count == featureIDs.count, featureIDs.allSatisfy(vivoOmicsID),
              normalizationTarget.isFinite, normalizationTarget > 0, normalizationTarget <= 1_000_000_000,
              missingValue.isFinite, missingValue < 0, missingValue >= -100 else {
            throw VivoOmicsError.invalid("cell-response feature axis")
        }
    }
}

/// A source-bound plan for one prepared count store. Multi-source training is
/// composed from independently verified corpora so no source identity is lost.
public struct VivoCellResponseCorpusPlan: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let id: String
    public let sourceDescription: String
    public let featureAxis: VivoCellResponseFeatureAxis
    public let targets: [VivoCellResponseTarget]
    /// Canonical fingerprint of the external descriptor document copied into
    /// the prepared corpus and checked against the target vectors below.
    public let descriptorSource: VivoFingerprint
    public let assignments: [VivoCellResponseAssignment]
    public let heldOutTargetIDs: [String]
    public let heldOutContextIDs: [String]

    private enum CodingKeys: String, CodingKey {
        case schemaVersion, id, sourceDescription, featureAxis, targets, descriptorSource, assignments, heldOutTargetIDs, heldOutContextIDs
    }

    public init(id: String, sourceDescription: String, featureAxis: VivoCellResponseFeatureAxis,
                targets: [VivoCellResponseTarget], descriptorSource: VivoFingerprint,
                assignments: [VivoCellResponseAssignment],
                heldOutTargetIDs: [String] = [], heldOutContextIDs: [String] = []) {
        schemaVersion = 3
        self.id = id
        self.sourceDescription = sourceDescription
        self.featureAxis = featureAxis
        self.targets = targets
        self.descriptorSource = descriptorSource
        self.assignments = assignments
        self.heldOutTargetIDs = heldOutTargetIDs
        self.heldOutContextIDs = heldOutContextIDs
    }

    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["schemaVersion", "id", "sourceDescription", "featureAxis", "targets", "descriptorSource", "assignments", "heldOutTargetIDs", "heldOutContextIDs"])
        let values = try decoder.container(keyedBy: CodingKeys.self)
        schemaVersion = try values.decode(Int.self, forKey: .schemaVersion)
        id = try values.decode(String.self, forKey: .id)
        sourceDescription = try values.decode(String.self, forKey: .sourceDescription)
        featureAxis = try values.decode(VivoCellResponseFeatureAxis.self, forKey: .featureAxis)
        targets = try values.decode([VivoCellResponseTarget].self, forKey: .targets)
        descriptorSource = try values.decode(VivoFingerprint.self, forKey: .descriptorSource)
        assignments = try values.decode([VivoCellResponseAssignment].self, forKey: .assignments)
        heldOutTargetIDs = try values.decodeIfPresent([String].self, forKey: .heldOutTargetIDs) ?? []
        heldOutContextIDs = try values.decodeIfPresent([String].self, forKey: .heldOutContextIDs) ?? []
    }

    public func validate() throws {
        try featureAxis.validate()
        guard schemaVersion == 3, vivoOmicsID(id), !sourceDescription.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty,
              sourceDescription.utf8.count <= 16_384, !targets.isEmpty, targets.count <= 100_000,
              !assignments.isEmpty, assignments.count <= 2_000_000,
              Set(targets.map(\.id)).count == targets.count,
              Set(assignments.map(\.sampleID)).count == assignments.count,
              Set(heldOutTargetIDs).count == heldOutTargetIDs.count, Set(heldOutContextIDs).count == heldOutContextIDs.count,
              heldOutTargetIDs.allSatisfy(vivoOmicsID), heldOutContextIDs.allSatisfy(vivoOmicsID) else {
            throw VivoOmicsError.invalid("cell-response corpus plan identity or bounds")
        }
        let descriptorCount = targets.first?.descriptors.count
        for target in targets { try target.validate(expectedDescriptorCount: descriptorCount) }
        let targetIDs = Set(targets.map(\.id))
        guard Set(heldOutTargetIDs).isSubset(of: targetIDs) else { throw VivoOmicsError.invalid("unknown held-out target") }
        for assignment in assignments { try assignment.validate(targets: targetIDs) }
        for assignment in assignments where assignment.partition == .training {
            if let target = assignment.targetID, heldOutTargetIDs.contains(target) {
                throw VivoOmicsError.invalid("held-out target appears in training")
            }
            if heldOutContextIDs.contains(assignment.contextID) {
                throw VivoOmicsError.invalid("held-out context appears in training")
            }
        }
        for target in heldOutTargetIDs {
            guard assignments.filter({ $0.targetID == target }).allSatisfy({ $0.partition != .training }) else {
                throw VivoOmicsError.invalid("held-out target outcome appears in training")
            }
        }
        for context in heldOutContextIDs {
            guard assignments.filter({ $0.contextID == context }).allSatisfy({ $0.partition != .training }) else {
                throw VivoOmicsError.invalid("held-out context appears in training")
            }
        }
        for assignment in assignments where assignment.role == .perturbed {
            guard assignments.contains(where: {
                $0.role == .control && $0.pairID == assignment.pairID &&
                    $0.contextID == assignment.contextID && $0.studyID == assignment.studyID &&
                    $0.modality == assignment.modality && $0.partition == assignment.partition
            }) else { throw VivoOmicsError.invalid("perturbed assignment has no declared matched control") }
        }
    }
}

/// Immutable, replayable corpus receipt. The source count store remains the
/// authority; the corpus adds only a verified row index and sample labels.
public struct VivoCellResponseCorpusReceipt: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let format: String
    public let sourceStore: VivoH5ADCountStoreReceipt
    public let plan: VivoFingerprint
    public let descriptorSource: VivoFingerprint
    public let metadata: VivoFingerprint
    public let quality: VivoFingerprint
    public let rowOffsets: VivoFingerprint
    public let rowAssignments: VivoFingerprint
    public let implementation: VivoFingerprint
}

/// One sparse source row. Values are exact raw counts; callers choose an
/// explicit transform rather than receiving silently normalized values.
public struct VivoCellResponseSparseRow: Sendable, Equatable {
    public let sourceRow: Int
    public let featureIndices: [Int]
    public let counts: [UInt64]
}

/// A treated row and its declared control context. The row identities stay
/// visible so callers can audit every training example.
public struct VivoCellResponseTrainingExample: Sendable, Equatable {
    public let targetRow: Int
    public let contextRows: [Int]
    public let targetIndex: Int
    /// One declared perturbed sample and its matched control context. The
    /// learner samples these strata uniformly before selecting a target cell.
    public let stratumIndex: Int
    public let partition: VivoCellResponsePartition
}

public final class VivoCellResponseCorpusReader {
    public let plan: VivoCellResponseCorpusPlan
    public let receipt: VivoCellResponseCorpusReceipt
    public let metadata: VivoSingleCellCountMetadata
    public let quality: VivoCountStoreQuality
    public let featureMask: [Float]

    // Retain the verified count-store lease for the full reader lifetime. Row
    // reads must never fall back to the mutable source-store path.
    private let sourceSnapshot: VivoH5ADCountStoreSnapshot
    private let records: VivoWindowedCountRecords
    private let rowOffsets: [UInt64]
    private let rowAssignments: [Int32]
    private let featureToAxis: [Int]
    private let assignmentRows: [[Int]]
    private let assignmentIndicesBySampleID: [String: Int]

    init(plan: VivoCellResponseCorpusPlan, receipt: VivoCellResponseCorpusReceipt,
         metadata: VivoSingleCellCountMetadata, quality: VivoCountStoreQuality,
         sourceSnapshot: VivoH5ADCountStoreSnapshot, rowOffsets: [UInt64], rowAssignments: [Int32]) throws {
        self.plan = plan
        self.receipt = receipt
        self.metadata = metadata
        self.quality = quality
        self.sourceSnapshot = sourceSnapshot
        self.records = sourceSnapshot.records
        self.rowOffsets = rowOffsets
        self.rowAssignments = rowAssignments
        let sourceFeatureIndex = Dictionary(uniqueKeysWithValues: metadata.features.enumerated().map { ($0.element.id, $0.offset) })
        let axisIndex = Dictionary(uniqueKeysWithValues: plan.featureAxis.featureIDs.enumerated().map { ($0.element, $0.offset) })
        featureToAxis = metadata.features.map { axisIndex[$0.id] ?? -1 }
        featureMask = plan.featureAxis.featureIDs.map { sourceFeatureIndex[$0] == nil ? Float(0) : Float(1) }
        var collected = Array(repeating: [Int](), count: plan.assignments.count)
        for (row, assignment) in rowAssignments.enumerated() {
            guard assignment >= 0, Int(assignment) < collected.count else { throw VivoOmicsError.invalid("cell-response row assignment") }
            collected[Int(assignment)].append(row)
        }
        assignmentRows = collected
        assignmentIndicesBySampleID = Dictionary(uniqueKeysWithValues: plan.assignments.enumerated().map {
            ($0.element.sampleID, $0.offset)
        })
    }

    public var featureCount: Int { plan.featureAxis.featureIDs.count }
    public var targetCount: Int { plan.targets.count }
    public var descriptorCount: Int { plan.targets[0].descriptors.count }

    /// Return only targets that have at least one declared, observed training
    /// row. The order follows the corpus vocabulary and is therefore stable
    /// for checkpointed embedding indices.
    public func trainingTargetBindings() throws -> [VivoCellResponseTrainedTargetBinding] {
        var targetIndices: Set<Int> = []
        let lookup = Dictionary(uniqueKeysWithValues: plan.targets.enumerated().map { ($0.element.id, $0.offset) })
        for assignmentIndex in plan.assignments.indices {
            let assignment = plan.assignments[assignmentIndex]
            guard assignment.role == .perturbed, assignment.partition == .training,
                  let targetID = assignment.targetID, let targetIndex = lookup[targetID] else {
                continue
            }
            guard !assignmentRows[assignmentIndex].isEmpty else {
                throw VivoOmicsError.invalid("cell-response training assignment has no rows")
            }
            targetIndices.insert(targetIndex)
        }
        guard !targetIndices.isEmpty else { throw VivoOmicsError.invalid("cell-response corpus has no training targets") }
        return try plan.targets.enumerated().compactMap { index, target in
            guard targetIndices.contains(index) else { return nil }
            return .init(id: target.id, descriptorFingerprint: try target.fingerprint())
        }
    }

    public func rows(for assignmentIndex: Int) throws -> [Int] {
        guard assignmentIndex >= 0, assignmentIndex < assignmentRows.count else { throw VivoOmicsError.invalid("cell-response assignment index") }
        return assignmentRows[assignmentIndex]
    }

    /// Resolve a declared source sample to its cell rows without inferring
    /// condition identity from the count matrix.
    public func rows(forSampleID sampleID: String) throws -> [Int] {
        guard let assignment = assignmentIndicesBySampleID[sampleID] else {
            throw VivoOmicsError.invalid("cell-response source sample")
        }
        return try rows(for: assignment)
    }

    public func assignment(forSampleID sampleID: String) throws -> VivoCellResponseAssignment {
        guard let index = assignmentIndicesBySampleID[sampleID] else {
            throw VivoOmicsError.invalid("cell-response source sample")
        }
        return plan.assignments[index]
    }

    public func sparseRow(_ sourceRow: Int) throws -> VivoCellResponseSparseRow {
        guard sourceRow >= 0, sourceRow + 1 < rowOffsets.count else { throw VivoOmicsError.invalid("cell-response source row") }
        let start64 = rowOffsets[sourceRow], end64 = rowOffsets[sourceRow + 1]
        guard start64 <= end64, end64 <= UInt64(records.count), start64 <= UInt64(Int.max), end64 <= UInt64(Int.max) else {
            throw VivoOmicsError.invalid("cell-response row offset")
        }
        var features: [Int] = []
        var counts: [UInt64] = []
        features.reserveCapacity(Int(end64 - start64)); counts.reserveCapacity(Int(end64 - start64))
        for recordIndex in Int(start64)..<Int(end64) {
            let record = try records.record(recordIndex)
            guard record.row == sourceRow, record.feature >= 0, record.feature < metadata.features.count, record.bits > 0 else {
                throw VivoOmicsError.invalid("cell-response row record")
            }
            features.append(record.feature); counts.append(record.bits)
        }
        return .init(sourceRow: sourceRow, featureIndices: features, counts: counts)
    }

    /// Log-normalized values on the explicit prediction axis. A source feature
    /// absent from that axis is ignored, while an axis feature absent from the
    /// source is kept at the declared negative missing sentinel.
    public func normalizedRow(_ sourceRow: Int) throws -> [Float] {
        guard sourceRow >= 0, sourceRow < quality.rowTotals.count, quality.rowTotals[sourceRow] > 0 else {
            throw VivoOmicsError.invalid("cell-response source row total")
        }
        let row = try sparseRow(sourceRow)
        var output = featureMask.map { $0 > 0 ? 0 : Float(plan.featureAxis.missingValue) }
        let scale = plan.featureAxis.normalizationTarget / Double(quality.rowTotals[sourceRow])
        for (feature, count) in zip(row.featureIndices, row.counts) {
            let axis = featureToAxis[feature]
            if axis >= 0 { output[axis] = Float(log1p(Double(count) * scale)) }
        }
        return output
    }

    public func examples(in partition: VivoCellResponsePartition) throws -> [VivoCellResponseTrainingExample] {
        var controls: [String: [Int]] = [:]
        for index in plan.assignments.indices {
            let assignment = plan.assignments[index]
            guard assignment.role == .control else { continue }
            let key = assignment.studyID + "\u{1f}" + assignment.modality + "\u{1f}" +
                assignment.contextID + "\u{1f}" + assignment.pairID + "\u{1f}" + assignment.partition.rawValue
            controls[key, default: []].append(contentsOf: assignmentRows[index])
        }
        var targetIndices: [String: Int] = [:]
        for (index, target) in plan.targets.enumerated() { targetIndices[target.id] = index }
        var output: [VivoCellResponseTrainingExample] = []
        for assignmentIndex in plan.assignments.indices {
            let assignment = plan.assignments[assignmentIndex]
            guard assignment.role == .perturbed, assignment.partition == partition,
                  let targetID = assignment.targetID, let targetIndex = targetIndices[targetID] else { continue }
            let key = assignment.studyID + "\u{1f}" + assignment.modality + "\u{1f}" +
                assignment.contextID + "\u{1f}" + assignment.pairID + "\u{1f}" + assignment.partition.rawValue
            guard let context = controls[key], !context.isEmpty else { throw VivoOmicsError.invalid("cell-response matched control rows") }
            for targetRow in assignmentRows[assignmentIndex] {
                output.append(.init(targetRow: targetRow, contextRows: context, targetIndex: targetIndex,
                                    stratumIndex: assignmentIndex, partition: partition))
            }
        }
        guard !output.isEmpty else { throw VivoOmicsError.invalid("cell-response partition has no treated rows") }
        return output
    }
}

public enum VivoCellResponseCorpus {
    public static let format = "numivivo-cell-response-corpus/v2"

    private static func document(_ root: URL, _ name: String, maximum: Int) throws -> Data {
        try VivoSingleCellCampaignIO.readDocument(root.appendingPathComponent(name), maximumBytes: maximum)
    }

    private static func staging(_ parent: URL) throws -> URL {
        let directory = parent.appendingPathComponent(".numivivo-cell-response-" + UUID().uuidString)
        try FileManager.default.createDirectory(at: directory, withIntermediateDirectories: false, attributes: [.posixPermissions: 0o700])
        return directory
    }

    private static func uint64Data(_ values: [UInt64]) -> Data {
        let copy = values.map(\.littleEndian)
        return copy.withUnsafeBytes { Data($0) }
    }

    private static func int32Data(_ values: [Int32]) -> Data {
        let copy = values.map(\.littleEndian)
        return copy.withUnsafeBytes { Data($0) }
    }

    private static func uint64Values(_ data: Data, count: Int) throws -> [UInt64] {
        guard count >= 0, data.count == count * MemoryLayout<UInt64>.stride else { throw VivoOmicsError.invalid("cell-response UInt64 sidecar size") }
        return stride(from: 0, to: data.count, by: MemoryLayout<UInt64>.stride).map { offset in
            data.withUnsafeBytes { UInt64(littleEndian: $0.loadUnaligned(fromByteOffset: offset, as: UInt64.self)) }
        }
    }

    private static func int32Values(_ data: Data, count: Int) throws -> [Int32] {
        guard count >= 0, data.count == count * MemoryLayout<Int32>.stride else { throw VivoOmicsError.invalid("cell-response Int32 sidecar size") }
        return stride(from: 0, to: data.count, by: MemoryLayout<Int32>.stride).map { offset in
            data.withUnsafeBytes { Int32(littleEndian: $0.loadUnaligned(fromByteOffset: offset, as: Int32.self)) }
        }
    }

    private static func buildRowOffsets(records: VivoWindowedCountRecords, cellCount: Int, featureCount: Int) throws -> [UInt64] {
        guard cellCount >= 0, featureCount >= 0 else { throw VivoOmicsError.invalid("cell-response count shape") }
        var offsets = [UInt64](repeating: 0, count: cellCount + 1)
        var nextOffset = 0, previousRow = -1, previousFeature = -1
        for index in 0..<records.count {
            try Task.checkCancellation()
            let record = try records.record(index)
            guard record.row >= 0, record.row < cellCount, record.feature >= 0, record.feature < featureCount, record.bits > 0,
                  record.row >= previousRow,
                  record.row != previousRow || record.feature > previousFeature else {
                throw VivoOmicsError.invalid("cell-response requires row-major canonical count records")
            }
            while nextOffset <= record.row { offsets[nextOffset] = UInt64(index); nextOffset += 1 }
            previousRow = record.row; previousFeature = record.feature
        }
        while nextOffset < offsets.count { offsets[nextOffset] = UInt64(records.count); nextOffset += 1 }
        return offsets
    }

    private static func assignments(_ plan: VivoCellResponseCorpusPlan, metadata: VivoSingleCellCountMetadata) throws -> [Int32] {
        let lookup = Dictionary(uniqueKeysWithValues: plan.assignments.enumerated().map { ($0.element.sampleID, $0.offset) })
        let sourceMetadata = Dictionary(uniqueKeysWithValues: metadata.samples.map { ($0.id, $0) })
        // Count-store metadata may retain source sample labels for strata that
        // have no surviving selected cells. They have no rows to learn from,
        // so a corpus must bind exactly the observed sample population rather
        // than invent a role for a metadata-only label. Every cell still has
        // one declared assignment and every assignment must own observed rows.
        let declaredSampleIDs = Set(metadata.samples.map(\.id))
        let observedSampleIDs = Set(metadata.cells.map(\.sampleID))
        guard observedSampleIDs.isSubset(of: declaredSampleIDs),
              Set(lookup.keys) == observedSampleIDs else {
            throw VivoOmicsError.invalid("cell-response assignments do not cover exactly the observed source samples")
        }
        for assignment in plan.assignments {
            guard sourceMetadata[assignment.sampleID] == assignment.sourceSample else {
                throw VivoOmicsError.invalid("cell-response assignment source metadata differs from count store")
            }
        }
        let result = try metadata.cells.map {
            guard let assignment = lookup[$0.sampleID], let stored = Int32(exactly: assignment) else {
                throw VivoOmicsError.invalid("cell-response source cell assignment")
            }
            return stored
        }
        var rowsPerAssignment = [Int](repeating: 0, count: plan.assignments.count)
        for assignment in result {
            rowsPerAssignment[Int(assignment)] += 1
        }
        guard rowsPerAssignment.allSatisfy({ $0 > 0 }) else {
            throw VivoOmicsError.invalid("cell-response assignment has no source cells")
        }
        return result
    }

    private static func validatedDescriptorSource(_ bytes: Data,
                                                  plan: VivoCellResponseCorpusPlan) throws -> VivoCellResponseDescriptorSource {
        guard bytes.count <= 67_108_864 else { throw VivoOmicsError.limit("cell-response descriptor source bytes") }
        let source = try VivoCanonicalJSON.decode(VivoCellResponseDescriptorSource.self, from: bytes)
        try source.validate()
        guard try VivoCanonicalJSON.encode(source) == bytes,
              try VivoCanonicalJSON.fingerprint(bytes) == plan.descriptorSource,
              source.targets == plan.targets else {
            throw VivoOmicsError.invalid("cell-response descriptor source differs from corpus plan")
        }
        return source
    }

    /// Fingerprint an opaque model or prediction payload with the same bounded
    /// file semantics used by the source count-store receipts.
    public static func fingerprint(_ source: URL) throws -> VivoFingerprint {
        try VivoH5ADCountStore.fingerprint(source)
    }

    public static func prepare(sourceStore: URL, plan: VivoCellResponseCorpusPlan, descriptorSourceBytes: Data,
                               implementation: VivoFingerprint, to destination: URL) throws -> VivoCellResponseCorpusReceipt {
        try plan.validate()
        _ = try validatedDescriptorSource(descriptorSourceBytes, plan: plan)
        try VivoH5ADCountStore.requireNew(destination)
        let sourceSnapshot = try VivoH5ADCountStore.openSnapshot(sourceStore, implementation: implementation)
        let sourceReceipt = sourceSnapshot.receipt
        let metadataBytes = sourceSnapshot.metadataBytes
        let qualityBytes = sourceSnapshot.qualityBytes
        let metadata = sourceSnapshot.metadata
        let quality = sourceSnapshot.quality
        guard quality.rowTotals.allSatisfy({ $0 > 0 }) else {
            throw VivoOmicsError.invalid("cell-response source count store contains an empty cell")
        }
        let sourceFeatureIDs = Set(metadata.features.map(\.id))
        guard plan.featureAxis.featureIDs.contains(where: sourceFeatureIDs.contains) else {
            throw VivoOmicsError.invalid("cell-response feature axis has no source overlap")
        }
        let rowAssignments = try assignments(plan, metadata: metadata)
        let offsets = try buildRowOffsets(records: sourceSnapshot.records, cellCount: metadata.cells.count, featureCount: metadata.features.count)
        let temporary = try staging(destination.deletingLastPathComponent())
        defer { try? FileManager.default.removeItem(at: temporary) }
        let planBytes = try VivoCanonicalJSON.encode(plan)
        guard planBytes.count <= 67_108_864 else { throw VivoOmicsError.limit("cell-response corpus plan bytes") }
        try planBytes.write(to: temporary.appendingPathComponent("plan.json"), options: .withoutOverwriting)
        try descriptorSourceBytes.write(to: temporary.appendingPathComponent("descriptor-source.json"), options: .withoutOverwriting)
        try metadataBytes.write(to: temporary.appendingPathComponent("metadata.json"), options: .withoutOverwriting)
        try qualityBytes.write(to: temporary.appendingPathComponent("quality.json"), options: .withoutOverwriting)
        let offsetsData = uint64Data(offsets), assignmentsData = int32Data(rowAssignments)
        try offsetsData.write(to: temporary.appendingPathComponent("row-offsets.bin"), options: .withoutOverwriting)
        try assignmentsData.write(to: temporary.appendingPathComponent("row-assignments.bin"), options: .withoutOverwriting)
        let receipt = try VivoCellResponseCorpusReceipt(
            schemaVersion: 2, format: format, sourceStore: sourceReceipt,
            plan: VivoCanonicalJSON.fingerprint(planBytes), descriptorSource: try VivoCanonicalJSON.fingerprint(descriptorSourceBytes),
            metadata: VivoCanonicalJSON.fingerprint(metadataBytes),
            quality: VivoCanonicalJSON.fingerprint(qualityBytes),
            rowOffsets: try VivoCanonicalJSON.fingerprint(offsetsData),
            rowAssignments: try VivoCanonicalJSON.fingerprint(assignmentsData),
            implementation: implementation)
        try VivoCanonicalJSON.encode(receipt).write(to: temporary.appendingPathComponent("receipt.json"), options: .withoutOverwriting)
        try Task.checkCancellation()
        try FileManager.default.moveItem(at: temporary, to: destination)
        return receipt
    }

    public static func open(_ corpus: URL, sourceStore: URL, implementation: VivoFingerprint) throws -> VivoCellResponseCorpusReader {
        let receiptBytes = try document(corpus, "receipt.json", maximum: 131_072)
        let receipt = try VivoCanonicalJSON.decode(VivoCellResponseCorpusReceipt.self, from: receiptBytes)
        guard receipt.schemaVersion == 2, receipt.format == format, receipt.implementation == implementation,
              try VivoCanonicalJSON.encode(receipt) == receiptBytes else { throw VivoOmicsError.invalid("cell-response corpus receipt") }
        let planBytes = try document(corpus, "plan.json", maximum: 67_108_864)
        let descriptorSourceBytes = try document(corpus, "descriptor-source.json", maximum: 67_108_864)
        let metadataBytes = try document(corpus, "metadata.json", maximum: 536_870_912)
        let qualityBytes = try document(corpus, "quality.json", maximum: 268_435_456)
        guard try VivoCanonicalJSON.fingerprint(planBytes) == receipt.plan,
              try VivoCanonicalJSON.fingerprint(descriptorSourceBytes) == receipt.descriptorSource,
              try VivoCanonicalJSON.fingerprint(metadataBytes) == receipt.metadata,
              try VivoCanonicalJSON.fingerprint(qualityBytes) == receipt.quality else {
            throw VivoOmicsError.invalid("cell-response corpus artifact fingerprints")
        }
        let plan = try VivoCanonicalJSON.decode(VivoCellResponseCorpusPlan.self, from: planBytes)
        guard try VivoCanonicalJSON.encode(plan) == planBytes else {
            throw VivoOmicsError.invalid("cell-response corpus plan is not canonical")
        }
        try plan.validate()
        _ = try validatedDescriptorSource(descriptorSourceBytes, plan: plan)
        let metadata = try VivoCanonicalJSON.decode(VivoSingleCellCountMetadata.self, from: metadataBytes)
        let quality = try VivoCanonicalJSON.decode(VivoCountStoreQuality.self, from: qualityBytes)
        guard try VivoCanonicalJSON.encode(metadata) == metadataBytes,
              try VivoCanonicalJSON.encode(quality) == qualityBytes else {
            throw VivoOmicsError.invalid("cell-response corpus metadata is not canonical")
        }
        let sourceSnapshot = try VivoH5ADCountStore.openSnapshot(sourceStore, implementation: implementation)
        guard sourceSnapshot.receipt == receipt.sourceStore,
              sourceSnapshot.metadataBytes == metadataBytes, sourceSnapshot.qualityBytes == qualityBytes,
              sourceSnapshot.metadata == metadata, sourceSnapshot.quality == quality else {
            throw VivoOmicsError.invalid("cell-response source store differs")
        }
        let rowOffsetsData = try document(corpus, "row-offsets.bin", maximum: (metadata.cells.count + 1) * MemoryLayout<UInt64>.stride)
        let rowAssignmentsData = try document(corpus, "row-assignments.bin", maximum: metadata.cells.count * MemoryLayout<Int32>.stride)
        guard try VivoCanonicalJSON.fingerprint(rowOffsetsData) == receipt.rowOffsets,
              try VivoCanonicalJSON.fingerprint(rowAssignmentsData) == receipt.rowAssignments else {
            throw VivoOmicsError.invalid("cell-response corpus sidecar fingerprints")
        }
        let rowOffsets = try uint64Values(rowOffsetsData, count: metadata.cells.count + 1)
        let rowAssignments = try int32Values(rowAssignmentsData, count: metadata.cells.count)
        let rebuiltOffsets = try buildRowOffsets(records: sourceSnapshot.records, cellCount: metadata.cells.count, featureCount: metadata.features.count)
        let rebuiltAssignments = try assignments(plan, metadata: metadata)
        guard rowOffsets == rebuiltOffsets, rowAssignments == rebuiltAssignments else {
            throw VivoOmicsError.invalid("cell-response corpus reconstruction differs")
        }
        return try .init(plan: plan, receipt: receipt, metadata: metadata, quality: quality,
                         sourceSnapshot: sourceSnapshot, rowOffsets: rowOffsets, rowAssignments: rowAssignments)
    }

    public static func verify(_ corpus: URL, sourceStore: URL, implementation: VivoFingerprint) throws -> VivoCellResponseCorpusReceipt {
        try open(corpus, sourceStore: sourceStore, implementation: implementation).receipt
    }
}
