import Foundation

/// A standalone source mapping and supplied program definitions. No pseudobulk,
/// dimensional reduction, cell selection or inferred annotation is performed.
public struct VivoH5ADProgramPlan: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let mapping: VivoH5ADImportPlan
    public let programs: VivoSingleCellProgramOptions
    public let normalizationTarget: Double
    public let matchFeatureNames: Bool
    public init(mapping: VivoH5ADImportPlan, programs: VivoSingleCellProgramOptions, normalizationTarget: Double = 10_000, matchFeatureNames: Bool = false) {
        schemaVersion=1;self.mapping=mapping;self.programs=programs;self.normalizationTarget=normalizationTarget;self.matchFeatureNames=matchFeatureNames
    }
    private enum CodingKeys: String,CodingKey { case schemaVersion,mapping,programs,normalizationTarget,matchFeatureNames }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder,allowed: ["schemaVersion","mapping","programs","normalizationTarget","matchFeatureNames"])
        let c=try decoder.container(keyedBy: CodingKeys.self)
        schemaVersion=try c.decode(Int.self,forKey: .schemaVersion)
        mapping=try c.decode(VivoH5ADImportPlan.self,forKey: .mapping)
        programs=try c.decode(VivoSingleCellProgramOptions.self,forKey: .programs)
        matchFeatureNames=try c.decodeIfPresent(Bool.self,forKey: .matchFeatureNames) ?? false
        normalizationTarget=try c.decodeIfPresent(Double.self,forKey: .normalizationTarget) ?? 10_000
    }
    public func validate() throws {
        guard schemaVersion==1,normalizationTarget.isFinite,normalizationTarget>0 else {
            throw VivoOmicsError.invalid("standalone program schema or normalization")
        }
        try programs.validate()
    }
}
public struct VivoH5ADProgramModel: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let method: String
    public let normalizationTarget: Double
    public let cells: Int
    public let canonicalNonzeros: Int
    public let hdf5Version: String
    public let programs: [VivoSingleCellResolvedProgram]
    public let updates: Int
    public let sourcePasses: Int
    public let emptyCellScoreEncoding: String
    public let qualification: String
}
public struct VivoH5ADProgramReceipt: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let matrixFormat: String
    public let integerFormat: String
    public let source: VivoFingerprint
    public let plan: VivoFingerprint
    public let metadata: VivoFingerprint
    public let model: VivoFingerprint
    public let scores: VivoFingerprint
    public let detectedMembers: VivoFingerprint
    public let totalCounts: VivoFingerprint
    public let implementation: VivoFingerprint
}

public enum VivoH5ADPrograms {
    static let matrixFormat="complete-row-major-u32-row-u32-program-f64-le/v1"
    static let integerFormat="complete-row-major-u32-row-u32-column-u64-le/v1"
    private static let maximumMatrixBytes=20_000_000*16
    private static func staging(_ parent: URL) throws -> URL {
        let url=parent.appendingPathComponent(".numivivo-programs-"+UUID().uuidString)
        try FileManager.default.createDirectory(at: url,withIntermediateDirectories: false,attributes: [.posixPermissions: 0o700])
        return url
    }
    private static func writeJSON<T: Encodable>(_ value: T,name: String,root: URL,maximum: Int) throws -> VivoFingerprint {
        let bytes=try VivoCanonicalJSON.encode(value)
        guard bytes.count<=maximum else { throw VivoOmicsError.limit("program artifact bytes") }
        try bytes.write(to: root.appendingPathComponent(name),options: .withoutOverwriting)
        return try VivoCanonicalJSON.fingerprint(bytes)
    }
    /// All coordinates are explicit. A zero total is the authoritative missing
    /// score mask; its finite binary score placeholder is +0, never measured zero.
    static func writeScores(_ a: VivoSingleCellProgramAccumulator,root: URL) throws -> (VivoFingerprint,VivoFingerprint) {
        let scores=try VivoCountRecordWriter(root.appendingPathComponent("scores.bin"))
        let detected=try VivoCountRecordWriter(root.appendingPathComponent("detected-members.bin"))
        for row in a.cells.indices {
            try Task.checkCancellation()
            for p in a.programs.indices {
                let i=row*a.programs.count+p,value=a.values[i],count=a.detected[i]
                guard value.isFinite,count>=0,count<=a.programs[p].featureIndices.count,
                      a.hasLibrary[row] || (value==0 && count==0) else { throw VivoOmicsError.invalid("program output value or availability") }
                try scores.append(row: row,feature: p,bits: value.bitPattern)
                try detected.append(row: row,feature: p,bits: UInt64(count))
            }
        }
        return (try scores.finish(),try detected.finish())
    }
    public static func publish(source: URL,plan: VivoH5ADProgramPlan,implementation: VivoFingerprint,to destination: URL) throws -> VivoH5ADProgramReceipt {
        try plan.validate();try VivoH5ADCountStore.requireNew(destination)
        let temp=try staging(destination.deletingLastPathComponent());defer { try? FileManager.default.removeItem(at: temp) }
        let planHash=try writeJSON(plan,name: "plan.json",root: temp,maximum: 2_097_152)
        let snapshot=temp.appendingPathComponent("original.h5ad")
        let sourceHash=try VivoH5ADPseudobulk.fingerprint(source,copyTo: snapshot)
        var qc: VivoSingleCellQualityAccumulator?
        let version=try VivoSingleCellH5AD.scanSnapshot(snapshot,plan: plan.mapping,limits: VivoH5ADPseudobulk.sourceLimits,onMetadata: {
            qc=VivoSingleCellQualityAccumulator($0)
        },onEntry: { row,feature,count in
            guard let qc else { throw VivoOmicsError.invalid("program source metadata missing") }
            try qc.add(row: row,feature: feature,count: count)
        })
        guard let qc else { throw VivoOmicsError.invalid("program source metadata missing") }
        let metadata=qc.metadata,quality=qc.finish()
        let a=try VivoSingleCellPrograms.accumulate(snapshot: snapshot,mapping: plan.mapping,metadata: metadata,quality: quality,
                                                   normalizationTarget: plan.normalizationTarget,options: plan.programs,featureMatch: plan.matchFeatureNames ? .featureName : .featureID)
        let model=VivoH5ADProgramModel(schemaVersion: 1,method: "signed-l1-mean-log-normalized-expression-v1",
            normalizationTarget: plan.normalizationTarget,cells: metadata.cells.count,canonicalNonzeros: qc.nonzeros,
            hdf5Version: version,programs: a.programs,updates: a.updates,sourcePasses: 2,
            emptyCellScoreEncoding: "total-counts.bin column0 == 0 means unavailable score; scores.bin contains +0 placeholder; detected-members.bin contains 0",
            qualification: "Supplied expression summaries with explicit exact ID/name mapping in plan and original feature indices in model, not calibrated pathway activity or authoritative labels. Two sparse source scans; axes, cell QC and flat cell-by-program arrays are resident; binary publication uses bounded buffers. No cell-by-gene allocation or biological qualification implied.")
        let metadataHash=try writeJSON(metadata,name: "metadata.json",root: temp,maximum: 536_870_912)
        let modelHash=try writeJSON(model,name: "model.json",root: temp,maximum: 67_108_864)
        let (scores,detected)=try writeScores(a,root: temp)
        let totals=try VivoCountRecordWriter(temp.appendingPathComponent("total-counts.bin"))
        for (row,q) in quality.enumerated() { try totals.append(row: row,feature: 0,bits: q.totalCounts) }
        let receipt=VivoH5ADProgramReceipt(schemaVersion: 1,matrixFormat: matrixFormat,integerFormat: integerFormat,source: sourceHash,
            plan: planHash,metadata: metadataHash,model: modelHash,scores: scores,detectedMembers: detected,totalCounts: try totals.finish(),implementation: implementation)
        _=try writeJSON(receipt,name: "receipt.json",root: temp,maximum: 65_536)
        try Task.checkCancellation();try FileManager.default.moveItem(at: temp,to: destination)
        return receipt
    }
    public static func verify(_ directory: URL,implementation: VivoFingerprint) throws -> VivoH5ADProgramReceipt {
        let bytes=try VivoSingleCellCampaignIO.readDocument(directory.appendingPathComponent("receipt.json"),maximumBytes: 65_536)
        let receipt=try VivoCanonicalJSON.decode(VivoH5ADProgramReceipt.self,from: bytes)
        guard receipt.schemaVersion==1,receipt.matrixFormat==matrixFormat,receipt.integerFormat==integerFormat,
              receipt.implementation==implementation,try VivoCanonicalJSON.encode(receipt)==bytes else { throw VivoOmicsError.invalid("program receipt") }
        let planBytes=try VivoSingleCellCampaignIO.readDocument(directory.appendingPathComponent("plan.json"),maximumBytes: 2_097_152)
        let plan=try VivoCanonicalJSON.decode(VivoH5ADProgramPlan.self,from: planBytes)
        guard try VivoCanonicalJSON.encode(plan)==planBytes,try VivoCanonicalJSON.fingerprint(planBytes)==receipt.plan else { throw VivoOmicsError.invalid("program plan") }
        for (name,hash,maximum) in [("metadata.json",receipt.metadata,536_870_912),("model.json",receipt.model,67_108_864),
            ("scores.bin",receipt.scores,maximumMatrixBytes),("detected-members.bin",receipt.detectedMembers,maximumMatrixBytes),
            ("total-counts.bin",receipt.totalCounts,maximumMatrixBytes)] {
            guard try VivoOmicsFileSnapshot.fingerprint(directory.appendingPathComponent(name),maximumBytes: maximum)==hash else { throw VivoOmicsError.invalid("program artifact fingerprint differs") }
        }
        let temp=try staging(FileManager.default.temporaryDirectory);defer { try? FileManager.default.removeItem(at: temp) }
        let rebuilt=try publish(source: directory.appendingPathComponent("original.h5ad"),plan: plan,implementation: implementation,to: temp.appendingPathComponent("rebuilt"))
        guard rebuilt==receipt else { throw VivoOmicsError.invalid("program source reconstruction differs") }
        return receipt
    }
}
