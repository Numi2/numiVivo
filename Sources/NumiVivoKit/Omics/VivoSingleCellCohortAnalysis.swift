import Foundation

public struct VivoSingleCellAnalysisPlan: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public var id: String
    public var filter: VivoSingleCellFilterPolicy
    public var normalizationTarget: Double
    public var contrasts: [VivoOmicsExpressionContrast]
    public var neighbors: VivoSingleCellNeighborOptions?
    public var reduction: VivoSingleCellReductionOptions?
    public init(id: String, filter: VivoSingleCellFilterPolicy = .init(), normalizationTarget: Double = 10_000,
                contrasts: [VivoOmicsExpressionContrast] = [],reduction: VivoSingleCellReductionOptions? = nil, neighbors: VivoSingleCellNeighborOptions? = nil) {
        schemaVersion = 1; self.id = id; self.filter = filter
        self.normalizationTarget = normalizationTarget; self.contrasts = contrasts
        self.reduction = reduction; self.neighbors = neighbors
    }
    private enum CodingKeys: String, CodingKey { case schemaVersion, id, filter, normalizationTarget, contrasts, reduction, neighbors }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["schemaVersion", "id", "filter", "normalizationTarget", "contrasts", "reduction", "neighbors"])
        let c = try decoder.container(keyedBy: CodingKeys.self)
        schemaVersion = try c.decode(Int.self, forKey: .schemaVersion); id = try c.decode(String.self, forKey: .id)
        filter = try c.decodeIfPresent(VivoSingleCellFilterPolicy.self, forKey: .filter) ?? .init()
        normalizationTarget = try c.decodeIfPresent(Double.self, forKey: .normalizationTarget) ?? 10_000
        contrasts = try c.decodeIfPresent([VivoOmicsExpressionContrast].self, forKey: .contrasts) ?? []
        neighbors = try c.decodeIfPresent(VivoSingleCellNeighborOptions.self,forKey: .neighbors)
        reduction = try c.decodeIfPresent(VivoSingleCellReductionOptions.self,forKey: .reduction)
    }
    public func validate() throws {
        guard schemaVersion == 1, vivoOmicsID(id), normalizationTarget.isFinite, normalizationTarget > 0,
              contrasts.count <= 32, Set(contrasts.map(\.id)).count == contrasts.count else {
            throw VivoOmicsError.invalid("analysis plan schema, identity, normalization or contrast count")
        }
        try filter.validate()
        try reduction?.validate()
        try neighbors?.validate()
        guard neighbors == nil || reduction != nil else { throw VivoOmicsError.invalid("neighbors require PCA reduction") }
        for contrast in contrasts { try contrast.validate() }
    }
}
public struct VivoSingleCellCohortReport: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let method: String
    public let plan: VivoSingleCellAnalysisPlan
    public let processed: VivoSingleCellProcessed
    public let contrasts: [VivoOmicsExpressionResult]
    public var reduction: VivoSingleCellReductionResult? = nil
    public var neighbors: VivoSingleCellNeighborGraph? = nil
}
public enum VivoSingleCellCohortAnalysis {
    public static func run(_ dataset: VivoSingleCellDataset, plan: VivoSingleCellAnalysisPlan) throws -> VivoSingleCellCohortReport {
        try plan.validate()
        let processed = try VivoSingleCellProcessing.run(dataset, policy: plan.filter, normalizationTarget: plan.normalizationTarget)
        var results: [VivoOmicsExpressionResult] = []
        for contrast in plan.contrasts {
            try Task.checkCancellation()
            results.append(try VivoPseudobulkDifferentialExpression.evaluate(processed.dataset, bulk: processed.pseudobulk, contrast: contrast))
        }
        try Task.checkCancellation()
        let reduction=try plan.reduction.map { try VivoSingleCellReduction.run(processed,options: $0) }
        let neighbors = try plan.neighbors.map { try VivoSingleCellNeighbors.run(scores: reduction!.scores, cells: reduction!.cells, options: $0) }
        return .init(schemaVersion: 1, method: "native-count-quality-and-donor-expression-v1", plan: plan, processed: processed, contrasts: results,reduction: reduction,neighbors: neighbors)
    }
}

public struct VivoSingleCellAnalysisReceipt: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let request: VivoFingerprint
    public let result: VivoFingerprint
    public let implementation: VivoFingerprint
    public init(request: VivoFingerprint, result: VivoFingerprint, implementation: VivoFingerprint) {
        schemaVersion = 1; self.request = request; self.result = result; self.implementation = implementation
    }
}
private struct VivoSingleCellStoredAnalysisRequest: Codable, Sendable, Equatable {
    let schemaVersion: Int
    let counts: VivoSingleCellRunReceipt
    let planBytes: Data
}
private struct VivoSingleCellStoredAnalysisResult: Codable, Sendable, Equatable {
    let schemaVersion: Int
    let request: VivoFingerprint
    let implementation: VivoFingerprint
    let report: VivoSingleCellCohortReport
}
public enum VivoSingleCellAnalysisArtifacts {
    public static let requestKind = "vivo.singlecell-analysis-request-v1"
    public static let resultKind = "vivo.singlecell-analysis-result-v1"
    public static let planKind = "vivo.singlecell-analysis-plan-v1"
    public static func decodePlan(_ bytes: Data) throws -> VivoSingleCellAnalysisPlan {
        guard bytes.count <= 2 * 1_024 * 1_024 else { throw VivoOmicsError.limit("analysis plan exceeds 2 MiB") }
        let plan = try VivoCanonicalJSON.decode(VivoSingleCellAnalysisPlan.self, from: bytes)
        try plan.validate(); return plan
    }
    public static func publish(counts: VivoSingleCellRunReceipt, planBytes: Data, implementation: VivoFingerprint,
                               store: VivoArtifactStore) async throws -> VivoSingleCellAnalysisReceipt {
        let plan = try decodePlan(planBytes)
        let source = try await VivoSingleCellArtifacts.verify(receipt: counts, implementation: implementation, store: store)
        let report = try VivoSingleCellCohortAnalysis.run(source.dataset, plan: plan)
        let request = try VivoCanonicalJSON.encode(VivoSingleCellStoredAnalysisRequest(schemaVersion: 1, counts: counts, planBytes: planBytes))
        let requestID = try VivoCanonicalJSON.fingerprint(request)
        let result = try VivoCanonicalJSON.encode(VivoSingleCellStoredAnalysisResult(schemaVersion: 1, request: requestID,
                                                                                  implementation: implementation, report: report))
        try Task.checkCancellation()
        let savedRequest = try await store.put(data: request, kind: requestKind, mediaType: "application/json")
        guard savedRequest.fingerprint == requestID else { throw VivoOmicsError.invalid("stored analysis request identity differs") }
        try Task.checkCancellation()
        let savedResult = try await store.put(data: result, kind: resultKind, mediaType: "application/json")
        try Task.checkCancellation()
        return .init(request: requestID, result: savedResult.fingerprint, implementation: implementation)
    }
    public static func verify(_ receipt: VivoSingleCellAnalysisReceipt, implementation: VivoFingerprint,
                              store: VivoArtifactStore) async throws -> VivoSingleCellCohortReport {
        guard receipt.schemaVersion == 1, receipt.implementation == implementation else {
            throw VivoOmicsError.invalid("analysis receipt schema or executing implementation differs")
        }
        for (fingerprint, kind) in [(receipt.request, requestKind), (receipt.result, resultKind)] {
            let descriptor = try await store.descriptor(for: fingerprint)
            guard descriptor.kind == kind, descriptor.mediaType == "application/json" else { throw VivoOmicsError.invalid("analysis artifact kind") }
        }
        let requestBytes = try await store.data(for: receipt.request, maximumBytes: 4 * 1_024 * 1_024)
        let resultBytes = try await store.data(for: receipt.result, maximumBytes: 512 * 1_024 * 1_024)
        let request = try VivoCanonicalJSON.decode(VivoSingleCellStoredAnalysisRequest.self, from: requestBytes)
        let record = try VivoCanonicalJSON.decode(VivoSingleCellStoredAnalysisResult.self, from: resultBytes)
        guard request.schemaVersion == 1, record.schemaVersion == 1, record.request == receipt.request,
              record.implementation == implementation else { throw VivoOmicsError.invalid("analysis source binding") }
        let plan = try decodePlan(request.planBytes)
        let source = try await VivoSingleCellArtifacts.verify(receipt: request.counts, implementation: implementation, store: store)
        guard try VivoSingleCellCohortAnalysis.run(source.dataset, plan: plan) == record.report else {
            throw VivoOmicsError.invalid("analysis report differs from complete source reconstruction")
        }
        try Task.checkCancellation(); return record.report
    }
}
