import Foundation
import CryptoKit
#if canImport(Darwin)
import Darwin
#else
import Glibc
#endif

public struct VivoH5ADReductionOptions: Codable, Sendable, Equatable {
    public var pca: VivoSingleCellReductionOptions = .init()
    public var normalizationTarget: Double = 10_000
    public var maximumCacheBytes: Int = 2_000_000_000
    public var maximumEntryVisits: Int = 2_000_000_000
    public init() {}
    private enum CodingKeys: String,CodingKey { case pca,normalizationTarget,maximumCacheBytes,maximumEntryVisits }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder,allowed: ["pca","normalizationTarget","maximumCacheBytes","maximumEntryVisits"])
        let c=try decoder.container(keyedBy: CodingKeys.self)
        pca=try c.decodeIfPresent(VivoSingleCellReductionOptions.self,forKey: .pca) ?? .init()
        normalizationTarget=try c.decodeIfPresent(Double.self,forKey: .normalizationTarget) ?? 10_000
        maximumCacheBytes=try c.decodeIfPresent(Int.self,forKey: .maximumCacheBytes) ?? 2_000_000_000
        maximumEntryVisits=try c.decodeIfPresent(Int.self,forKey: .maximumEntryVisits) ?? 2_000_000_000
    }
    public func validate() throws {
        try pca.validate()
        guard normalizationTarget.isFinite,normalizationTarget>0,
              (16...2_000_000_000).contains(maximumCacheBytes), (1...20_000_000_000).contains(maximumEntryVisits) else {
            throw VivoOmicsError.invalid("streamed reduction options")
        }
    }
}
public struct VivoH5ADReductionStorage: Codable, Sendable, Equatable {
    public let method: String
    public let options: VivoH5ADReductionOptions
    public let sourcePasses: Int
    public let selectedEntries: Int
    public let cacheBytes: Int
    public let cacheFingerprint: VivoFingerprint
    public let entryVisits: Int
    public let qualification: String
}

/// Private scratch representation, reconstructed from the archived H5AD. Each
/// 16-byte little-endian record is UInt32 row, UInt32 selected-column, Float64
/// log-normalized value. It is not a new public count format or a second engine.
final class VivoMappedReductionEntries {
    private let records: VivoWindowedCountRecords
    private let rowCount: Int
    private let columnCount: Int
    let byteCount: Int
    let count: Int
    private(set) var visits = 0
    init(url: URL, expectedEntries: Int, rows: Int, columns: Int) throws {
        guard expectedEntries > 0, expectedEntries <= 125_000_000, rows > 0, columns > 0 else {
            throw VivoOmicsError.invalid("reduction cache dimensions")
        }
        records = try VivoWindowedCountRecords(url, entries: expectedEntries)
        rowCount = rows; columnCount = columns
        byteCount = expectedEntries * 16; count = expectedEntries
        // The owner retains an immutable private scratch file. Validate the
        // complete stream once, using the same bounded window as later passes.
        for i in 0..<count {
            let record = try records.record(i)
            let value = Double(bitPattern: record.bits)
            guard record.row < rows, record.feature < columns, value.isFinite, value > 0 else {
                throw VivoOmicsError.invalid("reduction cache record")
            }
        }
    }
    func project(_ vector: [Double], shift: Double, rows: Int) throws -> [Double] {
        guard vector.count == columnCount, rows == rowCount else {
            throw VivoOmicsError.invalid("reduction projection dimensions")
        }
        var result = [Double](repeating: -shift, count: rows)
        for i in 0..<count {
            let record = try records.record(i)
            result[record.row] += Double(bitPattern: record.bits) * vector[record.feature]
        }
        visits += count; return result
    }
    func transpose(_ vector: [Double], initial: [Double]) throws -> [Double] {
        guard vector.count == rowCount, initial.count == columnCount else {
            throw VivoOmicsError.invalid("reduction transpose dimensions")
        }
        var result = initial
        for i in 0..<count {
            let record = try records.record(i)
            result[record.feature] += Double(bitPattern: record.bits) * vector[record.row]
        }
        visits += count; return result
    }
}

enum VivoH5ADReduction {
    static func run(snapshot: URL,mapping: VivoH5ADImportPlan,metadata: VivoSingleCellCountMetadata,
                    quality: [VivoCellQuality],options: VivoH5ADReductionOptions) throws -> (VivoSingleCellReductionResult,VivoH5ADReductionStorage) {
        try options.validate()
        let n=metadata.cells.count,m=metadata.features.count
        guard n>options.pca.components, quality.count==n else { throw VivoOmicsError.invalid("streamed PCA cell count") }
        var seen=[Int](repeating: 0,count: m),means=[Double](repeating: 0,count: m),m2=means,logMeans=means,logM2=means
        func logValue(row: Int,count: UInt64) throws -> Double {
            guard quality[row].totalCounts>0 else { throw VivoOmicsError.invalid("nonzero entry in zero-library cell") }
            let value=log1p((Double(count)/Double(quality[row].totalCounts))*options.normalizationTarget)
            guard value.isFinite,value>0 else { throw VivoOmicsError.invalid("streamed normalization is nonfinite or underflowed") }
            return value
        }
        // Pass two after QC: every source feature contributes to normalization,
        // and implicit zero cells enter both the HVG and log-space moments.
        _ = try VivoSingleCellH5AD.scanSnapshot(snapshot,plan: mapping,limits: VivoH5ADPseudobulk.sourceLimits,onMetadata: {
            guard $0==metadata else { throw VivoOmicsError.invalid("streamed reduction metadata changed") }
        },onEntry: { row,j,count in
            let log=try logValue(row: row,count: count),value=expm1(log)
            seen[j]+=1
            let delta=value-means[j];means[j]+=delta/Double(seen[j]);m2[j]+=delta*(value-means[j])
            let logDelta=log-logMeans[j];logMeans[j]+=logDelta/Double(seen[j]);logM2[j]+=logDelta*(log-logMeans[j])
        })
        let selection=try VivoSingleCellReduction.selectFeatures(featureIDs: metadata.features.map(\.id),cells: n,seen: seen,nonzeroMeans: means,m2: m2,options: options.pca)
        let selected=selection.selected,width=selected.count
        var local=[Int](repeating: -1,count: m)
        for (j,source) in selected.enumerated() { local[source]=j }
        let centers=selected.map { logMeans[$0]*Double(seen[$0])/Double(n) }
        let totalVariance=selected.reduce(0.0) { sum,j in
            sum+(logM2[j]+logMeans[j]*logMeans[j]*Double(seen[j])*Double(n-seen[j])/Double(n))/Double(n-1)
        }
        let expectedEntries=selected.reduce(0) { $0+seen[$1] }
        guard expectedEntries>0,expectedEntries<=options.maximumCacheBytes/16 else { throw VivoOmicsError.limit("streamed PCA cache-byte budget") }
        let traversals=2*min(options.pca.maximumBasis,width)+3*options.pca.components
        let work=expectedEntries.multipliedReportingOverflow(by: traversals)
        guard !work.overflow,work.partialValue<=options.maximumEntryVisits else { throw VivoOmicsError.limit("streamed PCA entry-visit budget") }
        let url=snapshot.deletingLastPathComponent().appendingPathComponent(".reduction-"+UUID().uuidString+".bin")
        let fd=open(url.path,O_WRONLY|O_CREAT|O_EXCL|O_NOFOLLOW,mode_t(0o600))
        guard fd>=0 else { throw VivoOmicsError.invalid("cannot create reduction cache") }
        let output=FileHandle(fileDescriptor: fd,closeOnDealloc: true)
        defer { try? output.close(); try? FileManager.default.removeItem(at: url) }
        var buffer=Data(),hasher=SHA256(),written=0
        buffer.reserveCapacity(1_048_576)
        func append<T>(_ value: T) { var value=value; withUnsafeBytes(of: &value) { buffer.append(contentsOf: $0) } }
        func flush() throws { if !buffer.isEmpty { try output.write(contentsOf: buffer);hasher.update(data: buffer);buffer.removeAll(keepingCapacity: true) } }
        // Pass three writes only HVG entries. CSR and CSC sources retain their
        // canonical traversal order; axes in each record make both usable.
        _ = try VivoSingleCellH5AD.scanSnapshot(snapshot,plan: mapping,limits: VivoH5ADPseudobulk.sourceLimits,onMetadata: {
            guard $0==metadata else { throw VivoOmicsError.invalid("streamed reduction metadata changed") }
        },onEntry: { row,j,count in
            let column=local[j]
            if column>=0 {
                guard written<expectedEntries else { throw VivoOmicsError.invalid("streamed selected entry count changed") }
                append(UInt32(row).littleEndian);append(UInt32(column).littleEndian)
                append(try logValue(row: row,count: count).bitPattern.littleEndian);written+=1
                if buffer.count>=1_048_576 { try flush() }
            }
        })
        try flush();try output.synchronize();try output.close()
        guard written==expectedEntries else { throw VivoOmicsError.invalid("streamed selected entry count changed") }
        let fingerprint=try VivoFingerprint(bytes: Array(hasher.finalize()))
        let cache=try VivoMappedReductionEntries(url: url,expectedEntries: written,rows: n,columns: width)
        let result=try VivoSingleCellReduction.fit(cells: metadata.cells.map { .init(sampleID: $0.sampleID,barcode: $0.barcode) },statistics: selection.statistics,
            selected: selected,centers: centers,totalVariance: totalVariance,options: options.pca,
            project: { try cache.project($0,shift: $1,rows: n) },transpose: { try cache.transpose($0,initial: $1) })
        guard cache.visits==work.partialValue else { throw VivoOmicsError.invalid("streamed PCA entry-visit accounting") }
        return (result,.init(method: "three-source-passes-windowed-selected-COO-v2",options: options,sourcePasses: 3,
            selectedEntries: written,cacheBytes: written*16,cacheFingerprint: fingerprint,entryVisits: cache.visits,
            qualification: "HDF5 slices, bounded write buffer and 16 MiB POSIX mapping windows for selected expression; metadata, moments, Krylov basis, scores and report remain resident. No million-cell or GPU performance qualification."))
    }
}
