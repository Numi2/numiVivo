import Foundation
import CryptoKit

public struct VivoH5ADPseudobulkPlan: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let mapping: VivoH5ADImportPlan
    public let contrasts: [VivoOmicsExpressionContrast]
    public var reduction: VivoH5ADReductionOptions? = nil
    public init(mapping: VivoH5ADImportPlan,contrasts: [VivoOmicsExpressionContrast] = [],reduction: VivoH5ADReductionOptions? = nil) {
        schemaVersion=1; self.mapping=mapping; self.contrasts=contrasts;self.reduction=reduction
    }
    private enum CodingKeys: String,CodingKey { case schemaVersion,mapping,contrasts,reduction }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder,allowed: ["schemaVersion","mapping","contrasts","reduction"])
        let c=try decoder.container(keyedBy: CodingKeys.self)
        schemaVersion=try c.decode(Int.self,forKey: .schemaVersion)
        mapping=try c.decode(VivoH5ADImportPlan.self,forKey: .mapping)
        reduction=try c.decodeIfPresent(VivoH5ADReductionOptions.self,forKey: .reduction)
        contrasts=try c.decodeIfPresent([VivoOmicsExpressionContrast].self,forKey: .contrasts) ?? []
    }
    public func validate() throws {
        guard schemaVersion == 1, contrasts.count <= 32, Set(contrasts.map(\.id)).count == contrasts.count else {
            throw VivoOmicsError.invalid("streamed pseudobulk plan schema or contrasts")
        }
        try reduction?.validate()
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
    public var reduction: VivoSingleCellReductionResult? = nil
    public var reductionStorage: VivoH5ADReductionStorage? = nil
}
public struct VivoH5ADPseudobulkReceipt: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let source: VivoFingerprint
    public let plan: VivoFingerprint
    public let report: VivoFingerprint
    public let implementation: VivoFingerprint
}

/// All source cells are retained. This is an exact raw aggregation route;
/// cell filtering and normalized per-cell matrix output are separate operations.
public enum VivoH5ADPseudobulk {
    static var sourceLimits: VivoOmicsLimits {
        var limits=VivoOmicsLimits()
        limits.maximumCells=1_000_000; limits.maximumFeatures=100_000
        limits.maximumNonzeros=100_000_000; limits.maximumInputBytes=1_073_741_824
        return limits
    }
    static let maximumAggregateNonzeros=5_000_000
    private struct Key: Hashable { let replicate: String; let condition: String; let group: String? }
    private final class Accumulator {
        let metadata: VivoSingleCellCountMetadata
        let groups: [VivoPseudobulkGroup]
        let groupForRow: [Int]
        var sums: [[Int: UInt64]]
        var totals: [UInt64], mitochondrial: [UInt64], detected: [Int]
        var nonzeros=0, aggregateNonzeros=0
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
            totals=Array(repeating: 0,count: metadata.cells.count)
            mitochondrial=totals; detected=Array(repeating: 0,count: metadata.cells.count)
        }
        func add(row: Int,feature: Int,count: UInt64) throws {
            totals[row]=try vivoOmicsSum(totals[row],count); detected[row]+=1; nonzeros+=1
            if metadata.features[feature].mitochondrial { mitochondrial[row]=try vivoOmicsSum(mitochondrial[row],count) }
            let group=groupForRow[row]
            if sums[group][feature] == nil {
                guard aggregateNonzeros < maximumAggregateNonzeros else { throw VivoOmicsError.limit("streamed aggregate exceeds five million nonzeros") }
                aggregateNonzeros+=1
            }
            sums[group][feature]=try vivoOmicsSum(sums[group][feature] ?? 0,count)
        }
        func finish(version: String,contrasts: [VivoOmicsExpressionContrast]) throws -> VivoH5ADPseudobulkReport {
            let mitoFeatures=metadata.features.filter(\.mitochondrial).count
            let quality=metadata.cells.indices.map { i in
                VivoCellQuality(sampleID: metadata.cells[i].sampleID,barcode: metadata.cells[i].barcode,totalCounts: totals[i],
                    detectedFeatures: detected[i],mitochondrialCounts: mitochondrial[i],mitochondrialFeatureCount: mitoFeatures,
                    mitochondrialFraction: totals[i] == 0 || mitoFeatures == 0 ? nil : Double(mitochondrial[i])/Double(totals[i]))
            }
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
        var accumulator: Accumulator?
        let version=try VivoSingleCellH5AD.scanSnapshot(url,plan: plan.mapping,limits: sourceLimits,onMetadata: {
            accumulator=try Accumulator($0)
        },onEntry: { row,feature,count in
            guard let accumulator else { throw VivoOmicsError.invalid("stream has no metadata") }
            try accumulator.add(row: row,feature: feature,count: count)
        })
        guard let accumulator else { throw VivoOmicsError.invalid("stream has no metadata") }
        var report=try accumulator.finish(version: version,contrasts: plan.contrasts)
        if let options=plan.reduction {
            let (reduction,storage)=try VivoH5ADReduction.run(snapshot: url,mapping: plan.mapping,metadata: report.metadata,quality: report.quality,options: options)
            report.reduction=reduction;report.reductionStorage=storage
        }
        return report
    }
    /// Hash/copy in 1 MiB blocks, independent of matrix size. The copied bytes,
    /// not a subsequently reread live path, are the authority for computation.
    static func fingerprint(_ source: URL,copyTo destination: URL? = nil) throws -> VivoFingerprint {
        let attrs=try FileManager.default.attributesOfItem(atPath: source.path)
        guard attrs[.type] as? FileAttributeType == .typeRegular else { throw VivoOmicsError.invalid("source must be a regular file") }
        let input=try FileHandle(forReadingFrom: source); defer { try? input.close() }
        var output: FileHandle?
        if let destination {
            guard FileManager.default.createFile(atPath: destination.path,contents: nil,attributes: [.posixPermissions: 0o600]) else { throw VivoOmicsError.invalid("cannot create source snapshot") }
            output=try FileHandle(forWritingTo: destination)
        }
        defer { try? output?.close() }
        var hasher=SHA256(), count=0
        while let block=try input.read(upToCount: 1_048_576), !block.isEmpty {
            try Task.checkCancellation()
            guard block.count <= sourceLimits.maximumInputBytes-count else { throw VivoOmicsError.limit("H5AD snapshot exceeds 1 GiB") }
            count+=block.count; hasher.update(data: block); try output?.write(contentsOf: block)
        }
        try output?.synchronize()
        return try VivoFingerprint(bytes: Array(hasher.finalize()))
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
