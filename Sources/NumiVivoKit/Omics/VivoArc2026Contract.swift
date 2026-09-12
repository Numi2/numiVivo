import Foundation

/// Adapter for the published Arc 2026 task, not a replacement metric implementation.
/// Official scores are produced ONLY by the separately pinned cell-eval2 process.
public enum VivoArc2026 {
    public static let evaluatorRepository = "ArcInstitute/cell-eval2"
    public static let evaluatorCommit = "5e64833518a6603a0301cbe28185d49c30f4a986"
    public static let evaluatorVersion = "0.16.0"
    public static let ruleVersion = 3
    public static let featureCount = 18_533
    public static let targetCount = 300
    public static let contextsPerPhase = 3
    public static let control = "non-targeting"
    public static let maximumCountsPerCell: UInt64 = 1_000_000
    public static let metrics = [
        "pds_cosine", "expr_mse_unbiased_capped_norm",
        "de_wilcoxon_direction_fidelity_yield_raw", "de_wilcoxon_direction_reach_raw",
        "de_wilcoxon_sig_jaccard", "de_wilcoxon_lfc_nmae"
    ]
    // Resource admission, NOT extra competition rules. Prediction cell counts are
    // unconstrained by Arc; a larger local workload needs a reviewed storage profile.
    public static let maximumCells = 2_000_000
    public static let maximumEntries = 2_000_000_000
    public static let maximumFileBytes = 64 * 1_024 * 1_024 * 1_024
    public static let maximumJSONBytes = 16 * 1_024 * 1_024

    static func identifier(_ value: String) -> Bool {
        !value.isEmpty && value.utf8.count <= 80 && value.utf8.allSatisfy {
            (48...57).contains($0) || (65...90).contains($0) || (97...122).contains($0) || $0 == 45 || $0 == 95
        }
    }
    static func sha256(_ value: String) -> Bool {
        value.utf8.count == 64 && value.utf8.allSatisfy { (48...57).contains($0) || (97...102).contains($0) }
    }
    static func phase(_ value: String) throws {
        guard ["validation", "final-test"].contains(value) else { throw VivoOmicsError.invalid("Arc phase must be validation or final-test") }
    }
    static func contextIDs(_ values: [String]) throws {
        guard values.count == contextsPerPhase, Set(values).count == values.count,
              values.allSatisfy(identifier) else { throw VivoOmicsError.invalid("Arc phase requires three unique safe context IDs") }
    }
    static func labels(_ values: [String], count: Int) throws {
        guard values.count == count, Set(values).count == count, values.allSatisfy(vivoOmicsID) else {
            throw VivoOmicsError.invalid("Arc axis has missing, duplicate or invalid identities")
        }
    }
}

public struct VivoArc2026PrepareContext: Codable, Sendable, Equatable {
    public let id: String
    public let controlsH5AD: String
    public let targetsJSON: String
    public let perturbationColumn: String
    private enum CodingKeys: String, CodingKey { case id, controlsH5AD, targetsJSON, perturbationColumn }
    public init(id: String, controlsH5AD: String, targetsJSON: String, perturbationColumn: String) {
        self.id = id
        self.controlsH5AD = controlsH5AD
        self.targetsJSON = targetsJSON
        self.perturbationColumn = perturbationColumn
    }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["id", "controlsH5AD", "targetsJSON", "perturbationColumn"])
        let c = try decoder.container(keyedBy: CodingKeys.self)
        id = try c.decode(String.self, forKey: .id)
        controlsH5AD = try c.decode(String.self, forKey: .controlsH5AD)
        targetsJSON = try c.decode(String.self, forKey: .targetsJSON)
        perturbationColumn = try c.decode(String.self, forKey: .perturbationColumn)
    }
}

public struct VivoArc2026PreparePlan: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let phase: String
    public let contexts: [VivoArc2026PrepareContext]
    private enum CodingKeys: String, CodingKey { case schemaVersion, phase, contexts }
    public init(schemaVersion: Int, phase: String, contexts: [VivoArc2026PrepareContext]) {
        self.schemaVersion = schemaVersion
        self.phase = phase
        self.contexts = contexts
    }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["schemaVersion", "phase", "contexts"])
        let c = try decoder.container(keyedBy: CodingKeys.self)
        schemaVersion = try c.decode(Int.self, forKey: .schemaVersion)
        phase = try c.decode(String.self, forKey: .phase)
        contexts = try c.decode([VivoArc2026PrepareContext].self, forKey: .contexts)
    }
}

public struct VivoArc2026CountInspection: Codable, Sendable, Equatable {
    public let cells: Int
    public let entries: Int
    public let totalCounts: UInt64
    public let maximumCellTotal: UInt64
    public let cellsPerLabel: [String: Int]
    private enum CodingKeys: String, CodingKey { case cells, entries, totalCounts, maximumCellTotal, cellsPerLabel }
    public init(cells: Int, entries: Int, totalCounts: UInt64, maximumCellTotal: UInt64, cellsPerLabel: [String: Int]) {
        self.cells = cells
        self.entries = entries
        self.totalCounts = totalCounts
        self.maximumCellTotal = maximumCellTotal
        self.cellsPerLabel = cellsPerLabel
    }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["cells", "entries", "totalCounts", "maximumCellTotal", "cellsPerLabel"])
        let c = try decoder.container(keyedBy: CodingKeys.self)
        cells = try c.decode(Int.self, forKey: .cells)
        entries = try c.decode(Int.self, forKey: .entries)
        totalCounts = try c.decode(UInt64.self, forKey: .totalCounts)
        maximumCellTotal = try c.decode(UInt64.self, forKey: .maximumCellTotal)
        cellsPerLabel = try c.decode([String: Int].self, forKey: .cellsPerLabel)
    }
}

public struct VivoArc2026QueryContext: Codable, Sendable, Equatable {
    public let id: String
    public let perturbationColumn: String
    public let genes: [String]
    public let targets: [String]
    public let controlsSHA256: String
    public let targetsSHA256: String
    public let inspection: VivoArc2026CountInspection
    private enum CodingKeys: String, CodingKey { case id, perturbationColumn, genes, targets, controlsSHA256, targetsSHA256, inspection }
    public init(id: String, perturbationColumn: String, genes: [String], targets: [String], controlsSHA256: String, targetsSHA256: String, inspection: VivoArc2026CountInspection) {
        self.id = id
        self.perturbationColumn = perturbationColumn
        self.genes = genes
        self.targets = targets
        self.controlsSHA256 = controlsSHA256
        self.targetsSHA256 = targetsSHA256
        self.inspection = inspection
    }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["id", "perturbationColumn", "genes", "targets", "controlsSHA256", "targetsSHA256", "inspection"])
        let c = try decoder.container(keyedBy: CodingKeys.self)
        id = try c.decode(String.self, forKey: .id)
        perturbationColumn = try c.decode(String.self, forKey: .perturbationColumn)
        genes = try c.decode([String].self, forKey: .genes)
        targets = try c.decode([String].self, forKey: .targets)
        controlsSHA256 = try c.decode(String.self, forKey: .controlsSHA256)
        targetsSHA256 = try c.decode(String.self, forKey: .targetsSHA256)
        inspection = try c.decode(VivoArc2026CountInspection.self, forKey: .inspection)
    }
}

public struct VivoArc2026Query: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let phase: String
    public let evaluatorCommit: String
    public let contexts: [VivoArc2026QueryContext]
    private enum CodingKeys: String, CodingKey { case schemaVersion, phase, evaluatorCommit, contexts }
    public init(schemaVersion: Int, phase: String, evaluatorCommit: String, contexts: [VivoArc2026QueryContext]) {
        self.schemaVersion = schemaVersion
        self.phase = phase
        self.evaluatorCommit = evaluatorCommit
        self.contexts = contexts
    }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["schemaVersion", "phase", "evaluatorCommit", "contexts"])
        let c = try decoder.container(keyedBy: CodingKeys.self)
        schemaVersion = try c.decode(Int.self, forKey: .schemaVersion)
        phase = try c.decode(String.self, forKey: .phase)
        evaluatorCommit = try c.decode(String.self, forKey: .evaluatorCommit)
        contexts = try c.decode([VivoArc2026QueryContext].self, forKey: .contexts)
    }
}

public struct VivoArc2026PredictionInput: Codable, Sendable, Equatable {
    public let id: String
    public let predictionH5AD: String
    private enum CodingKeys: String, CodingKey { case id, predictionH5AD }
    public init(id: String, predictionH5AD: String) {
        self.id = id
        self.predictionH5AD = predictionH5AD
    }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["id", "predictionH5AD"])
        let c = try decoder.container(keyedBy: CodingKeys.self)
        id = try c.decode(String.self, forKey: .id)
        predictionH5AD = try c.decode(String.self, forKey: .predictionH5AD)
    }
}

public struct VivoArc2026PackPlan: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let modelID: String
    public let modelArtifact: String
    public let trainingDataManifest: String
    public let excludedPerturbedContexts: [String]
    public let contexts: [VivoArc2026PredictionInput]
    private enum CodingKeys: String, CodingKey { case schemaVersion, modelID, modelArtifact, trainingDataManifest, excludedPerturbedContexts, contexts }
    public init(schemaVersion: Int, modelID: String, modelArtifact: String, trainingDataManifest: String, excludedPerturbedContexts: [String], contexts: [VivoArc2026PredictionInput]) {
        self.schemaVersion = schemaVersion
        self.modelID = modelID
        self.modelArtifact = modelArtifact
        self.trainingDataManifest = trainingDataManifest
        self.excludedPerturbedContexts = excludedPerturbedContexts
        self.contexts = contexts
    }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["schemaVersion", "modelID", "modelArtifact", "trainingDataManifest", "excludedPerturbedContexts", "contexts"])
        let c = try decoder.container(keyedBy: CodingKeys.self)
        schemaVersion = try c.decode(Int.self, forKey: .schemaVersion)
        modelID = try c.decode(String.self, forKey: .modelID)
        modelArtifact = try c.decode(String.self, forKey: .modelArtifact)
        trainingDataManifest = try c.decode(String.self, forKey: .trainingDataManifest)
        excludedPerturbedContexts = try c.decode([String].self, forKey: .excludedPerturbedContexts)
        contexts = try c.decode([VivoArc2026PredictionInput].self, forKey: .contexts)
    }
}

public struct VivoArc2026ModelDeclaration: Codable, Sendable, Equatable {
    public let id: String
    public let artifactSHA256: String
    public let trainingManifestSHA256: String
    public let excludedPerturbedContexts: [String]
    private enum CodingKeys: String, CodingKey { case id, artifactSHA256, trainingManifestSHA256, excludedPerturbedContexts }
    public init(id: String, artifactSHA256: String, trainingManifestSHA256: String, excludedPerturbedContexts: [String]) {
        self.id = id
        self.artifactSHA256 = artifactSHA256
        self.trainingManifestSHA256 = trainingManifestSHA256
        self.excludedPerturbedContexts = excludedPerturbedContexts
    }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["id", "artifactSHA256", "trainingManifestSHA256", "excludedPerturbedContexts"])
        let c = try decoder.container(keyedBy: CodingKeys.self)
        id = try c.decode(String.self, forKey: .id)
        artifactSHA256 = try c.decode(String.self, forKey: .artifactSHA256)
        trainingManifestSHA256 = try c.decode(String.self, forKey: .trainingManifestSHA256)
        excludedPerturbedContexts = try c.decode([String].self, forKey: .excludedPerturbedContexts)
    }
}

public struct VivoArc2026PredictionReceipt: Codable, Sendable, Equatable {
    public let id: String
    public let predictionSHA256: String
    public let inspection: VivoArc2026CountInspection
    private enum CodingKeys: String, CodingKey { case id, predictionSHA256, inspection }
    public init(id: String, predictionSHA256: String, inspection: VivoArc2026CountInspection) {
        self.id = id
        self.predictionSHA256 = predictionSHA256
        self.inspection = inspection
    }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["id", "predictionSHA256", "inspection"])
        let c = try decoder.container(keyedBy: CodingKeys.self)
        id = try c.decode(String.self, forKey: .id)
        predictionSHA256 = try c.decode(String.self, forKey: .predictionSHA256)
        inspection = try c.decode(VivoArc2026CountInspection.self, forKey: .inspection)
    }
}

public struct VivoArc2026Submission: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let evidenceClass: String
    public let phase: String
    public let evaluatorCommit: String
    public let querySHA256: String
    public let model: VivoArc2026ModelDeclaration
    public let contexts: [VivoArc2026PredictionReceipt]
    private enum CodingKeys: String, CodingKey { case schemaVersion, evidenceClass, phase, evaluatorCommit, querySHA256, model, contexts }
    public init(schemaVersion: Int, evidenceClass: String, phase: String, evaluatorCommit: String, querySHA256: String, model: VivoArc2026ModelDeclaration, contexts: [VivoArc2026PredictionReceipt]) {
        self.schemaVersion = schemaVersion
        self.evidenceClass = evidenceClass
        self.phase = phase
        self.evaluatorCommit = evaluatorCommit
        self.querySHA256 = querySHA256
        self.model = model
        self.contexts = contexts
    }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["schemaVersion", "evidenceClass", "phase", "evaluatorCommit", "querySHA256", "model", "contexts"])
        let c = try decoder.container(keyedBy: CodingKeys.self)
        schemaVersion = try c.decode(Int.self, forKey: .schemaVersion)
        evidenceClass = try c.decode(String.self, forKey: .evidenceClass)
        phase = try c.decode(String.self, forKey: .phase)
        evaluatorCommit = try c.decode(String.self, forKey: .evaluatorCommit)
        querySHA256 = try c.decode(String.self, forKey: .querySHA256)
        model = try c.decode(VivoArc2026ModelDeclaration.self, forKey: .model)
        contexts = try c.decode([VivoArc2026PredictionReceipt].self, forKey: .contexts)
    }
}

extension VivoArc2026PreparePlan {
    public func validate() throws {
        guard schemaVersion == 1 else { throw VivoOmicsError.invalid("Arc preparation schema") }
        try VivoArc2026.phase(phase); try VivoArc2026.contextIDs(contexts.map(\.id))
        for context in contexts {
            guard !context.controlsH5AD.isEmpty, !context.targetsJSON.isEmpty,
                  VivoArc2026.identifier(context.perturbationColumn) else { throw VivoOmicsError.invalid("Arc source paths or label column") }
        }
    }
}
extension VivoArc2026Query {
    public func validate() throws {
        guard schemaVersion == 1, evaluatorCommit == VivoArc2026.evaluatorCommit else { throw VivoOmicsError.invalid("Arc query schema or evaluator identity") }
        try VivoArc2026.phase(phase); try VivoArc2026.contextIDs(contexts.map(\.id))
        for c in contexts {
            try VivoArc2026.labels(c.genes, count: VivoArc2026.featureCount)
            try VivoArc2026.labels(c.targets, count: VivoArc2026.targetCount)
            guard !c.targets.contains(VivoArc2026.control), VivoArc2026.identifier(c.perturbationColumn),
                  VivoArc2026.sha256(c.controlsSHA256), VivoArc2026.sha256(c.targetsSHA256),
                  c.genes == contexts[0].genes else { throw VivoOmicsError.invalid("Arc query identities, hashes or common ordered gene axis") }
        }
    }
}

/// Reuses the native exact-count reader. Entries can arrive in CSR or CSC order.
/// This object checks raw counts; it never normalizes, rounds, imputes or scores.
struct VivoArc2026CountAccumulator {
    private let labels: [String]
    private let featureCount: Int
    private var rowTotals: [UInt64]
    private var entries = 0
    init(labels: [String], featureCount: Int, expectedLabels: Set<String>) throws {
        guard !labels.isEmpty, labels.count <= VivoArc2026.maximumCells,
              featureCount > 0, Set(labels) == expectedLabels else { throw VivoOmicsError.invalid("Arc missing/extra perturbation labels or cell resource limit") }
        self.labels = labels; self.featureCount = featureCount
        rowTotals = Array(repeating: 0, count: labels.count)
    }
    mutating func append(row: Int, feature: Int, count: UInt64) throws {
        guard row >= 0, row < rowTotals.count, feature >= 0, feature < featureCount,
              count > 0, entries < VivoArc2026.maximumEntries else { throw VivoOmicsError.invalid("Arc count coordinate or entry limit") }
        let total = try vivoOmicsSum(rowTotals[row], count)
        guard total <= VivoArc2026.maximumCountsPerCell else { throw VivoOmicsError.invalid("Arc cell exceeds 1000000 raw counts") }
        rowTotals[row] = total; entries += 1
    }
    func finish() throws -> VivoArc2026CountInspection {
        var total: UInt64 = 0, perLabel: [String: Int] = [:]
        for (i, label) in labels.enumerated() {
            total = try vivoOmicsSum(total, rowTotals[i]); perLabel[label, default: 0] += 1
        }
        return .init(cells: labels.count, entries: entries, totalCounts: total,
                     maximumCellTotal: rowTotals.max() ?? 0, cellsPerLabel: perLabel)
    }
}
