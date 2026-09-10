import Foundation
import CryptoKit

/// Observation indices address exactly the source bytes named by `source`.
/// Selection changes aggregation only; the complete source remains in the bundle.
public struct VivoH5ADCellSelection: Codable, Sendable, Equatable {
    public let source: VivoFingerprint
    public let observationIndices: [Int]
    public let provenance: String
    public init(source: VivoFingerprint, observationIndices: [Int], provenance: String) {
        self.source=source; self.observationIndices=observationIndices; self.provenance=provenance
    }
    private enum CodingKeys: String, CodingKey { case source,observationIndices,provenance }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder,allowed: ["source","observationIndices","provenance"])
        let c=try decoder.container(keyedBy: CodingKeys.self)
        source=try c.decode(VivoFingerprint.self,forKey: .source)
        observationIndices=try c.decode([Int].self,forKey: .observationIndices)
        provenance=try c.decode(String.self,forKey: .provenance)
    }
    func validate() throws {
        guard !observationIndices.isEmpty, observationIndices.count<=VivoH5ADPseudobulk.sourceLimits.maximumCells,
              observationIndices.allSatisfy({ (0..<VivoH5ADPseudobulk.sourceLimits.maximumCells).contains($0) }),
              Set(observationIndices).count==observationIndices.count,
              !provenance.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty,
              provenance.utf8.count<=16_384 else { throw VivoOmicsError.invalid("streamed cell selection indices or provenance") }
    }
}

public struct VivoH5ADPseudobulkPlan: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let mapping: VivoH5ADImportPlan
    public let contrasts: [VivoOmicsExpressionContrast]
    public var programs: VivoSingleCellProgramOptions? = nil
    public var reduction: VivoH5ADReductionOptions? = nil
    public var cellSelection: VivoH5ADCellSelection? = nil
    public init(mapping: VivoH5ADImportPlan,contrasts: [VivoOmicsExpressionContrast] = [],reduction: VivoH5ADReductionOptions? = nil,programs: VivoSingleCellProgramOptions? = nil) {
        schemaVersion=1; self.mapping=mapping; self.contrasts=contrasts;self.reduction=reduction;self.programs=programs
    }
    private enum CodingKeys: String,CodingKey { case schemaVersion,mapping,contrasts,reduction,programs,cellSelection }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder,allowed: ["schemaVersion","mapping","contrasts","reduction","programs","cellSelection"])
        let c=try decoder.container(keyedBy: CodingKeys.self)
        schemaVersion=try c.decode(Int.self,forKey: .schemaVersion)
        mapping=try c.decode(VivoH5ADImportPlan.self,forKey: .mapping)
        programs=try c.decodeIfPresent(VivoSingleCellProgramOptions.self,forKey: .programs)
        reduction=try c.decodeIfPresent(VivoH5ADReductionOptions.self,forKey: .reduction)
        cellSelection=try c.decodeIfPresent(VivoH5ADCellSelection.self,forKey: .cellSelection)
        contrasts=try c.decodeIfPresent([VivoOmicsExpressionContrast].self,forKey: .contrasts) ?? []
    }
    public func validate() throws {
        guard schemaVersion == 1, contrasts.count <= 32, Set(contrasts.map(\.id)).count == contrasts.count else {
            throw VivoOmicsError.invalid("streamed pseudobulk plan schema or contrasts")
        }
        try programs?.validate()
        try reduction?.validate()
        try cellSelection?.validate()
        guard cellSelection == nil || (programs == nil && reduction == nil) else {
            throw VivoOmicsError.invalid("streamed cell selection currently supports aggregation and contrasts only")
        }
        for contrast in contrasts { try contrast.validate() }
    }
}
public struct VivoH5ADPseudobulkReport: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let method: String
    public let metadata: VivoSingleCellCountMetadata
    public let quality: [VivoCellQuality]
    public let pseudobulk: VivoPseudobulkCounts
    public let canonicalNonzeros: Int
    public let hdf5Version: String
    public let contrasts: [VivoOmicsExpressionResult]
    public var programs: VivoSingleCellProgramResult? = nil
    public var reduction: VivoSingleCellReductionResult? = nil
    public var reductionStorage: VivoH5ADReductionStorage? = nil
    /// Maps report-local rows (including group member indices) to original rows.
    public var sourceObservationIndices: [Int]? = nil
    public var sourceCellCount: Int? = nil
}
public struct VivoH5ADPseudobulkReceipt: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let source: VivoFingerprint
    public let plan: VivoFingerprint
    public let report: VivoFingerprint
    public let implementation: VivoFingerprint
}

/// The complete source file is retained. Optional source-bound cell selection
/// restricts raw aggregation without writing another cell-by-feature matrix.
public enum VivoH5ADPseudobulk {
    static var sourceLimits: VivoOmicsLimits {
        var limits=VivoOmicsLimits()
        limits.maximumCells=2_000_000; limits.maximumFeatures=100_000
        // Source entry count bounds scan work, not resident matrix allocation:
        // the reader keeps at most one sparse major segment and 65,536-entry
        // input slices. The separate aggregate allowance remains unchanged.
        limits.maximumNonzeros=4_000_000_000
        // Source snapshots use fixed 1 MiB I/O, independently of the retained
        // report/plan bounds below. Match the existing count-store file ceiling
        // so complete atlas files also enter streamed aggregation and PCA.
        limits.maximumInputBytes=64 * 1_024 * 1_024 * 1_024
        return limits
    }
    static let maximumAggregateNonzeros=5_000_000
    private struct Key: Hashable { let replicate: String; let condition: String; let group: String? }
    private final class Accumulator {
        let metadata: VivoSingleCellCountMetadata
        let groups: [VivoPseudobulkGroup]
        let groupForRow: [Int]
        var sums: [[Int: UInt64]]
        let quality: VivoSingleCellQualityAccumulator
        var nonzeros: Int { quality.nonzeros }
        var aggregateNonzeros=0
        init(_ metadata: VivoSingleCellCountMetadata) throws {
            self.metadata=metadata
            let samples=Dictionary(uniqueKeysWithValues: metadata.samples.map { ($0.id,$0) })
            var members: [Key:[Int]]=[:]
            for (i,cell) in metadata.cells.enumerated() {
                let sample=samples[cell.sampleID]!
                members[.init(replicate: sample.biologicalReplicateID,condition: sample.condition,group: cell.group),default: []].append(i)
            }
            let keys=members.keys.sorted {
                if $0.replicate != $1.replicate { return $0.replicate < $1.replicate }
                if $0.condition != $1.condition { return $0.condition < $1.condition }
                if $0.group == nil { return $1.group != nil }; if $1.group == nil { return false }
                return $0.group! < $1.group!
            }
            var groups: [VivoPseudobulkGroup]=[], assignment=[Int](repeating: 0,count: metadata.cells.count)
            for (i,key) in keys.enumerated() {
                let rows=members[key]!, first=samples[metadata.cells[rows[0]].sampleID]!
                for row in rows { assignment[row]=i }
                let ids=Set(rows.map { metadata.cells[$0].sampleID }).sorted()
                groups.append(.init(biologicalReplicateID: key.replicate,donorID: first.donorID,condition: key.condition,
                    organism: first.organism,cellGroup: key.group,sampleIDs: ids,
                    batchIDs: Set(ids.map { samples[$0]!.batchID }).sorted(),sourceCellIndices: rows))
            }
            self.groups=groups; groupForRow=assignment; sums=Array(repeating: [:],count: groups.count)
            quality = VivoSingleCellQualityAccumulator(metadata)
        }
        func add(row: Int,feature: Int,count: UInt64) throws {
            try quality.add(row: row, feature: feature, count: count)
            let group=groupForRow[row]
            if sums[group][feature] == nil {
                guard aggregateNonzeros < maximumAggregateNonzeros else { throw VivoOmicsError.limit("streamed aggregate exceeds five million nonzeros") }
                aggregateNonzeros+=1
            }
            sums[group][feature]=try vivoOmicsSum(sums[group][feature] ?? 0,count)
        }
        func finish(version: String,contrasts: [VivoOmicsExpressionContrast]) throws -> VivoH5ADPseudobulkReport {
            let quality = self.quality.finish()
            var offsets=[0],columns: [Int]=[],counts: [UInt64]=[]
            for group in sums {
                for feature in group.keys.sorted() { columns.append(feature); counts.append(group[feature]!) }
                offsets.append(counts.count)
            }
            let matrix=VivoSparseCounts(cellCount: groups.count,featureCount: metadata.features.count,rowOffsets: offsets,featureIndices: columns,counts: counts)
            var limits=sourceLimits; limits.maximumNonzeros=maximumAggregateNonzeros
            try matrix.validate(limits: limits)
            let bulk=VivoPseudobulkCounts(method: "raw-sum-by-replicate-condition-group-v1",countUnit: metadata.countUnit,
                groups: groups,featureIDs: metadata.features.map(\.id),matrix: matrix)
            let results=try contrasts.map { try VivoPseudobulkDifferentialExpression.evaluate(metadata: metadata,bulk: bulk,contrast: $0) }
            return .init(schemaVersion: 1,method: "H5AD-bounded-slices-all-source-cells-pseudobulk-v1",metadata: metadata,
                quality: quality,pseudobulk: bulk,canonicalNonzeros: nonzeros,hdf5Version: version,contrasts: results)
        }
    }
    static func evaluateSnapshot(_ url: URL,plan: VivoH5ADPseudobulkPlan) throws -> VivoH5ADPseudobulkReport {
        try plan.validate()
        if let selection=plan.cellSelection {
            guard try fingerprint(url)==selection.source else { throw VivoOmicsError.invalid("streamed cell selection source changed") }
        }
        var accumulator: Accumulator?
        var rowMap: [Int]?, sourceCellCount: Int?
        let version=try VivoSingleCellH5AD.scanSnapshot(url,plan: plan.mapping,limits: sourceLimits,onMetadata: { metadata in
            if let selection=plan.cellSelection {
                let indices=selection.observationIndices
                guard indices.allSatisfy({ $0<metadata.cells.count }) else {
                    throw VivoOmicsError.invalid("streamed cell selection index outside source axis")
                }
                sourceCellCount=metadata.cells.count
                var map=[Int](repeating: -1,count: metadata.cells.count)
                for (local,original) in indices.enumerated() { map[original]=local }
                rowMap=map
                let selected=VivoSingleCellCountMetadata(id: metadata.id,evidence: metadata.evidence,
                    sourceDescription: metadata.sourceDescription,countUnit: metadata.countUnit,
                    samples: metadata.samples,features: metadata.features,cells: indices.map { metadata.cells[$0] })
                try selected.validate(limits: sourceLimits)
                accumulator=try Accumulator(selected)
            } else {
                accumulator=try Accumulator(metadata)
            }
        },onEntry: { row,feature,count in
            guard let accumulator else { throw VivoOmicsError.invalid("stream has no metadata") }
            let local=rowMap?[row] ?? row
            if local>=0 { try accumulator.add(row: local,feature: feature,count: count) }
        })
        guard let accumulator else { throw VivoOmicsError.invalid("stream has no metadata") }
        var report=try accumulator.finish(version: version,contrasts: plan.contrasts)
        report.sourceObservationIndices=plan.cellSelection?.observationIndices
        report.sourceCellCount=sourceCellCount
        if let options=plan.programs {
            report.programs=try VivoSingleCellPrograms.run(snapshot: url,mapping: plan.mapping,metadata: report.metadata,quality: report.quality,normalizationTarget: plan.reduction?.normalizationTarget ?? 10_000,options: options)
        }
        if let options=plan.reduction {
            let (reduction,storage)=try VivoH5ADReduction.run(snapshot: url,mapping: plan.mapping,metadata: report.metadata,quality: report.quality,options: options)
            report.reduction=reduction;report.reductionStorage=storage
        }
        return report
    }
    /// Hash/copy in 1 MiB blocks, independent of matrix size. The copied bytes,
    /// not a subsequently reread live path, are the authority for computation.
    static func fingerprint(_ source: URL,copyTo destination: URL? = nil) throws -> VivoFingerprint {
        try VivoOmicsFileSnapshot.fingerprint(source, copyTo: destination, maximumBytes: sourceLimits.maximumInputBytes)
    }
    public static func publish(source: URL,plan: VivoH5ADPseudobulkPlan,implementation: VivoFingerprint,to destination: URL) throws -> VivoH5ADPseudobulkReceipt {
        try plan.validate()
        let planBytes=try VivoCanonicalJSON.encode(plan)
        guard planBytes.count <= 2_097_152 else { throw VivoOmicsError.limit("streamed plan exceeds 2 MiB") }
        guard !FileManager.default.fileExists(atPath: destination.path) else { throw VivoOmicsError.invalid("output already exists") }
        let staging=destination.deletingLastPathComponent().appendingPathComponent(".numivivo-stream-"+UUID().uuidString)
        try FileManager.default.createDirectory(at: staging,withIntermediateDirectories: false,attributes: [.posixPermissions: 0o700])
        defer { try? FileManager.default.removeItem(at: staging) }
        let snapshot=staging.appendingPathComponent("original.h5ad")
        let sourceID=try fingerprint(source,copyTo: snapshot)
        let report=try evaluateSnapshot(snapshot,plan: plan)
        let reportBytes=try VivoCanonicalJSON.encode(report)
        guard reportBytes.count <= 536_870_912 else { throw VivoOmicsError.limit("streamed report exceeds 512 MiB") }
        let receipt=try VivoH5ADPseudobulkReceipt(schemaVersion: 1,source: sourceID,plan: VivoCanonicalJSON.fingerprint(planBytes),
            report: VivoCanonicalJSON.fingerprint(reportBytes),implementation: implementation)
        try planBytes.write(to: staging.appendingPathComponent("plan.json"),options: .withoutOverwriting)
        try reportBytes.write(to: staging.appendingPathComponent("report.json"),options: .withoutOverwriting)
        try VivoCanonicalJSON.encode(receipt).write(to: staging.appendingPathComponent("receipt.json"),options: .withoutOverwriting)
        try Task.checkCancellation()
        try FileManager.default.moveItem(at: staging,to: destination)
        return receipt
    }
    public static func verify(_ directory: URL,implementation: VivoFingerprint) throws -> VivoH5ADPseudobulkReport {
        let receipt=try VivoCanonicalJSON.decode(VivoH5ADPseudobulkReceipt.self,from: VivoSingleCellCampaignIO.readDocument(directory.appendingPathComponent("receipt.json"),maximumBytes: 65_536))
        let planBytes=try VivoSingleCellCampaignIO.readDocument(directory.appendingPathComponent("plan.json"),maximumBytes: 2_097_152)
        let reportBytes=try VivoSingleCellCampaignIO.readDocument(directory.appendingPathComponent("report.json"),maximumBytes: 536_870_912)
        guard receipt.schemaVersion == 1, receipt.implementation == implementation,
              try VivoCanonicalJSON.fingerprint(planBytes) == receipt.plan, try VivoCanonicalJSON.fingerprint(reportBytes) == receipt.report else {
            throw VivoOmicsError.invalid("streamed bundle source, plan, report or implementation changed")
        }
        let temporary=FileManager.default.temporaryDirectory.appendingPathComponent("numivivo-verify-"+UUID().uuidString)
        try FileManager.default.createDirectory(at: temporary,withIntermediateDirectories: false,attributes: [.posixPermissions: 0o700])
        defer { try? FileManager.default.removeItem(at: temporary) }
        let snapshot=temporary.appendingPathComponent("original.h5ad")
        guard try fingerprint(directory.appendingPathComponent("original.h5ad"),copyTo: snapshot) == receipt.source else {
            throw VivoOmicsError.invalid("streamed bundle source changed")
        }
        let plan=try VivoCanonicalJSON.decode(VivoH5ADPseudobulkPlan.self,from: planBytes)
        let reconstructed=try evaluateSnapshot(snapshot,plan: plan)
        guard try VivoCanonicalJSON.encode(reconstructed) == reportBytes else { throw VivoOmicsError.invalid("streamed result does not reconstruct") }
        return reconstructed
    }
}
