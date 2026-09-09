import Foundation

public struct VivoSingleCellReferencePlan: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let mapping: VivoH5ADImportPlan
    public let reduction: VivoH5ADReductionOptions
    public let featureNamespace: String
    public let labelProvenance: String
    public let neighbors: Int
    public init(mapping: VivoH5ADImportPlan,reduction: VivoH5ADReductionOptions = .init(),featureNamespace: String,labelProvenance: String,neighbors: Int = 15) {
        schemaVersion=1;self.mapping=mapping;self.reduction=reduction;self.featureNamespace=featureNamespace;self.labelProvenance=labelProvenance;self.neighbors=neighbors
    }
    private enum CodingKeys: String,CodingKey { case schemaVersion,mapping,reduction,featureNamespace,labelProvenance,neighbors }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder,allowed: ["schemaVersion","mapping","reduction","featureNamespace","labelProvenance","neighbors"])
        let c=try decoder.container(keyedBy: CodingKeys.self)
        schemaVersion=try c.decode(Int.self,forKey: .schemaVersion);mapping=try c.decode(VivoH5ADImportPlan.self,forKey: .mapping)
        reduction=try c.decodeIfPresent(VivoH5ADReductionOptions.self,forKey: .reduction) ?? .init()
        featureNamespace=try c.decode(String.self,forKey: .featureNamespace);labelProvenance=try c.decode(String.self,forKey: .labelProvenance)
        neighbors=try c.decodeIfPresent(Int.self,forKey: .neighbors) ?? 15
    }
    func validate() throws {
        try reduction.validate()
        guard schemaVersion==1,mapping.groupColumn != nil,vivoOmicsID(featureNamespace),
              !labelProvenance.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty,labelProvenance.utf8.count<=16_384,
              (1...128).contains(neighbors) else { throw VivoOmicsError.invalid("reference schema, label source, namespace or neighbors") }
    }
    var trainingPlan: VivoH5ADPseudobulkPlan {
        var options=reduction;options.pca.retainProjectionCenters=true
        return .init(mapping: mapping,reduction: options)
    }
}
public struct VivoSingleCellReferenceQueryPlan: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let mapping: VivoH5ADImportPlan
    public let featureNamespace: String
    public let maximumDistanceOperations: Int
    public let maximumProjectionUpdates: Int
    public init(mapping: VivoH5ADImportPlan,featureNamespace: String,maximumDistanceOperations: Int = 2_000_000_000,maximumProjectionUpdates: Int = 2_000_000_000) {
        schemaVersion=1;self.mapping=mapping;self.featureNamespace=featureNamespace
        self.maximumDistanceOperations=maximumDistanceOperations;self.maximumProjectionUpdates=maximumProjectionUpdates
    }
    private enum CodingKeys: String,CodingKey { case schemaVersion,mapping,featureNamespace,maximumDistanceOperations,maximumProjectionUpdates }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder,allowed: ["schemaVersion","mapping","featureNamespace","maximumDistanceOperations","maximumProjectionUpdates"])
        let c=try decoder.container(keyedBy: CodingKeys.self)
        schemaVersion=try c.decode(Int.self,forKey: .schemaVersion);mapping=try c.decode(VivoH5ADImportPlan.self,forKey: .mapping)
        featureNamespace=try c.decode(String.self,forKey: .featureNamespace)
        maximumDistanceOperations=try c.decodeIfPresent(Int.self,forKey: .maximumDistanceOperations) ?? 2_000_000_000
        maximumProjectionUpdates=try c.decodeIfPresent(Int.self,forKey: .maximumProjectionUpdates) ?? 2_000_000_000
    }
    func validate() throws {
        guard schemaVersion==1,mapping.groupColumn==nil,vivoOmicsID(featureNamespace),
              (1...20_000_000_000).contains(maximumDistanceOperations),(1...20_000_000_000).contains(maximumProjectionUpdates) else {
            throw VivoOmicsError.invalid("query schema, namespace or work budget; query label mapping is prohibited")
        }
    }
}
public struct VivoSingleCellReferenceModel: Codable, Sendable, Equatable {
    public let method: String
    public let plan: VivoSingleCellReferencePlan
    public let source: VivoFingerprint
    public let featureIDs: [String]
    public let organism: String
    public let reduction: VivoSingleCellReductionResult
    public let referenceLabels: [String]
    public let classes: [String]
    public let qualification: String
}
public struct VivoSingleCellReferencePrediction: Codable, Sendable, Equatable {
    public let cell: VivoOmicsCellIdentity
    public let totalCounts: UInt64
    public let scores: [Double]?
    public let neighborIndices: [Int]
    public let squaredDistances: [Double]
    public let votes: [Int]
    public let candidateLabel: String?
}
public struct VivoSingleCellReferenceReport: Codable, Sendable, Equatable {
    public let method: String
    public let referenceSource: VivoFingerprint
    public let classes: [String]
    public let cells: [VivoSingleCellReferencePrediction]
    public let overlappingDonorIDs: [String]
    public let distanceOperations: Int
    public let projectionUpdates: Int
    public let hdf5Version: String
    public let qualification: String
}

public enum VivoSingleCellReference {
    static func model(_ report: VivoH5ADPseudobulkReport,plan: VivoSingleCellReferencePlan,source: VivoFingerprint) throws -> VivoSingleCellReferenceModel {
        try plan.validate()
        let labels=report.metadata.cells.compactMap(\.group),organisms=Set(report.metadata.samples.map(\.organism))
        guard let reduction=report.reduction,let centers=reduction.projectionCenters,
              labels.count==report.metadata.cells.count,labels.allSatisfy(vivoOmicsID),
              labels.count>=plan.neighbors,labels.count<=100_000,Set(labels).count<=256,
              organisms.count==1,centers.count==reduction.selectedFeatureIndices.count,
              report.quality.allSatisfy({ $0.totalCounts>0 }) else {
            throw VivoOmicsError.invalid("reference requires complete labels, one organism, positive libraries and retained centers")
        }
        return .init(method: "frozen-sparse-log-PCA-uniform-reference-kNN-v1",plan: plan,source: source,
            featureIDs: report.metadata.features.map(\.id),organism: organisms.first!,reduction: reduction,
            referenceLabels: labels,classes: Set(labels).sorted(),
            qualification: "Training-derived reference voting; candidate labels only. Uncalibrated votes, no novelty detection, no cross-study or rare-class qualification.")
    }
    static func nearest(_ score: [Double],reference: [[Double]],neighbors: Int) throws -> [VivoSingleCellNeighbors.Neighbor] {
        var heap: [VivoSingleCellNeighbors.Neighbor]=[]
        for (index,row) in reference.enumerated() {
            if index%1024==0 { try Task.checkCancellation() }
            var distance=0.0
            for j in score.indices { let delta=score[j]-row[j];distance+=delta*delta }
            guard distance.isFinite else { throw VivoOmicsError.invalid("reference distance is nonfinite") }
            VivoSingleCellNeighbors.retain(.init(index: index,squaredDistance: distance),in: &heap,capacity: neighbors)
        }
        return heap.sorted { $0.precedes($1) }
    }
    static func evaluate(snapshot: URL,plan: VivoSingleCellReferenceQueryPlan,model: VivoSingleCellReferenceModel) throws -> VivoSingleCellReferenceReport {
        try plan.validate()
        guard plan.featureNamespace==model.plan.featureNamespace,plan.mapping.countUnit==model.plan.mapping.countUnit else {
            throw VivoOmicsError.invalid("reference/query namespace or count unit mismatch")
        }
        let reduction=model.reduction,d=reduction.options.components,centers=reduction.projectionCenters!
        var metadata: VivoSingleCellCountMetadata?,totals: [UInt64]=[],local: [Int]=[]
        var operations=0
        let version=try VivoSingleCellH5AD.scanSnapshot(snapshot,plan: plan.mapping,limits: VivoH5ADPseudobulk.sourceLimits,onMetadata: { value in
            guard Set(value.samples.map(\.organism))==[model.organism],value.features.count==model.featureIDs.count,
                  Set(value.features.map(\.id))==Set(model.featureIDs) else { throw VivoOmicsError.invalid("query must match the complete reference gene universe and organism") }
            let referenceCells=Set(reduction.cells)
            guard value.cells.allSatisfy({ !referenceCells.contains(.init(sampleID: $0.sampleID,barcode: $0.barcode)) }) else {
                throw VivoOmicsError.invalid("query overlaps reference cell identities")
            }
            let n=value.cells.count
            guard n<=5_000_000/max(d,max(model.classes.count,model.plan.neighbors)),
                  n<=plan.maximumDistanceOperations/d/reduction.cells.count else { throw VivoOmicsError.limit("reference query score/output or distance-operation budget") }
            operations=n*d*reduction.cells.count
            let selectedIDs=Dictionary(uniqueKeysWithValues: reduction.selectedFeatureIndices.enumerated().map { (model.featureIDs[$0.element],$0.offset) })
            local=value.features.map { selectedIDs[$0.id] ?? -1 }
            metadata=value;totals=Array(repeating: 0,count: n)
        },onEntry: { row,_,count in totals[row]=try vivoOmicsSum(totals[row],count) })
        guard let metadata else { throw VivoOmicsError.invalid("query metadata absent") }
        let projection = try VivoFrozenPCAProjection(centers: centers, loadings: reduction.loadings, target: model.plan.reduction.normalizationTarget)
        var scores=Array(repeating: projection.initial,count: metadata.cells.count),updates=0
        _=try VivoSingleCellH5AD.scanSnapshot(snapshot,plan: plan.mapping,limits: VivoH5ADPseudobulk.sourceLimits,onMetadata: {
            guard $0==metadata else { throw VivoOmicsError.invalid("query metadata changed between passes") }
        },onEntry: { row,j,count in
            let selected=local[j];if selected<0 { return }
            guard updates<=plan.maximumProjectionUpdates-d else { throw VivoOmicsError.limit("reference projection-update budget") }
            try scores[row].withUnsafeMutableBufferPointer { values in
                try projection.add(count: count, total: totals[row], selected: selected, to: values)
            };updates+=d
        })
        let labelIndex=Dictionary(uniqueKeysWithValues: model.classes.enumerated().map { ($0.element,$0.offset) })
        var predictions: [VivoSingleCellReferencePrediction]=[]
        for row in metadata.cells.indices {
            try Task.checkCancellation()
            let cell=metadata.cells[row],identity=VivoOmicsCellIdentity(sampleID: cell.sampleID,barcode: cell.barcode)
            if totals[row]==0 {
                operations-=d*reduction.cells.count
                predictions.append(.init(cell: identity,totalCounts: 0,scores: nil,neighborIndices: [],squaredDistances: [],votes: [],candidateLabel: nil));continue
            }
            guard scores[row].allSatisfy(\.isFinite) else { throw VivoOmicsError.invalid("query scores nonfinite") }
            let neighbors=try nearest(scores[row],reference: reduction.scores,neighbors: model.plan.neighbors)
            var votes=[Int](repeating: 0,count: model.classes.count)
            for neighbor in neighbors { votes[labelIndex[model.referenceLabels[neighbor.index]]!]+=1 }
            let winner=votes.firstIndex(of: votes.max()!)!
            predictions.append(.init(cell: identity,totalCounts: totals[row],scores: scores[row],neighborIndices: neighbors.map(\.index),
                squaredDistances: neighbors.map(\.squaredDistance),votes: votes,candidateLabel: model.classes[winner]))
        }
        let donors=Set(model.plan.mapping.samples.compactMap(\.donorID)).intersection(metadata.samples.compactMap(\.donorID)).sorted()
        return .init(method: model.method,referenceSource: model.source,classes: model.classes,cells: predictions,overlappingDonorIDs: donors,
            distanceOperations: operations,projectionUpdates: updates,hdf5Version: version,
            qualification: "Frozen training-only projection and uniform kNN votes; labels are candidates and votes are uncalibrated. Empty libraries remain unmapped. No novel-class rejection or biological identity guarantee. No query labels read or query fitting performed.")
    }
}
