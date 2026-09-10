import Foundation
import CryptoKit
import NumiVivoCore
#if canImport(Darwin)
import Darwin
#else
import Glibc
#endif

/// Persistent exact counts. Each source-major record is little-endian UInt32
/// row, UInt32 feature, UInt64 count. Source traversal order is preserved; this
/// is not a CSR index and does not promise identical bytes across source layouts.
public struct VivoH5ADCountStoreReceipt: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let format: String
    public let entries: Int
    public let source: VivoFingerprint
    public let plan: VivoFingerprint
    public let metadata: VivoFingerprint
    public let quality: VivoFingerprint
    public let counts: VivoFingerprint
    public let implementation: VivoFingerprint
}
public struct VivoCountStoreQuality: Codable, Sendable, Equatable {
    public var rowTotals: [UInt64]
    public var rowNonzeros: [UInt64]
    public var featureTotals: [UInt64]
    public var featureNonzeros: [UInt64]
}
public struct VivoCountStoreNormalizationReceipt: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let format: String
    public let target: Double
    public let entries: Int
    public let input: VivoFingerprint
    public let metadata: VivoFingerprint
    public let values: VivoFingerprint
    public let implementation: VivoFingerprint
}

/// At most 1 MiB of record payload is buffered, independently of entry count.
final class VivoCountRecordWriter {
    private let output: FileHandle?
    private var buffer: [UInt64] = []
    private var hasher = SHA256()
    private(set) var entries = 0
    init(_ url: URL?) throws {
        if let url {
            let fd = open(url.path, O_WRONLY | O_CREAT | O_EXCL | O_NOFOLLOW, 0o600)
            guard fd >= 0 else { throw VivoOmicsError.invalid("cannot create count record file") }
            output = FileHandle(fileDescriptor: fd, closeOnDealloc: true)
        } else { output = nil }
        buffer.reserveCapacity(131_072)
    }
    func append(row: Int, feature: Int, bits: UInt64) throws {
        guard row >= 0, row <= Int(UInt32.max), feature >= 0, feature <= Int(UInt32.max),
              entries < VivoH5ADCountStore.maximumEntries else { throw VivoOmicsError.limit("count record coordinates or entries") }
        buffer.append((UInt64(row) | (UInt64(feature) << 32)).littleEndian)
        buffer.append(bits.littleEndian); entries += 1
        if buffer.count == 131_072 { try flush() }
    }
    private func flush() throws {
        guard !buffer.isEmpty else { return }
        try Task.checkCancellation()
        try buffer.withUnsafeBytes {
            let data = Data($0); hasher.update(data: data); try output?.write(contentsOf: data)
        }
        buffer.removeAll(keepingCapacity: true)
    }
    func finish() throws -> VivoFingerprint {
        try flush(); try output?.synchronize(); try output?.close()
        return try VivoFingerprint(bytes: Array(hasher.finalize()))
    }
}

/// Shared storage bounds for PCA artifacts and their downstream consumers.
enum VivoPCAStorageLimits {
    static let maximumRows = Int(NVIVO_OMICS_PCA_MAXIMUM_ROWS)
    static let maximumQualityBytes = 536_870_912
    static let maximumColumns = Int(NVIVO_OMICS_PCA_MAXIMUM_COLUMNS)
    static let maximumBytes = maximumRows * maximumColumns * 16
    static let maximumSymmetricGraphEntries = maximumRows * 2 * (Int(NVIVO_OMICS_HNSW_MAXIMUM_NEIGHBORS) - 1)
    static let maximumGraphBytes = maximumSymmetricGraphEntries * 16
}

/// Maps only a 16 MiB window. The caller owns an immutable private snapshot;
/// unmapping each previous window bounds the mapped working set even for large
/// files. Metadata and row/feature statistics are separate resident structures.
final class VivoWindowedCountRecords {
    static let windowBytes = 16 * 1_024 * 1_024
    private let fd: Int32
    let count: Int
    private var address: UnsafeMutableRawPointer?
    private var offset = -1
    private var length = 0
    init(_ url: URL, entries: Int) throws {
        guard entries >= 0, entries <= VivoH5ADCountStore.maximumEntries else { throw VivoOmicsError.limit("mapped count entries") }
        fd = open(url.path, O_RDONLY | O_NOFOLLOW | O_NONBLOCK)
        guard fd >= 0 else { throw VivoOmicsError.invalid("cannot open count records") }
        var info = stat()
        guard fstat(fd, &info) == 0, (info.st_mode & mode_t(S_IFMT)) == mode_t(S_IFREG), info.st_size == entries * 16 else {
            _ = close(fd); throw VivoOmicsError.invalid("count record size or type")
        }
        count = entries
    }
    deinit { if let address { _ = munmap(address, length) }; _ = close(fd) }
    func record(_ index: Int) throws -> (row: Int, feature: Int, bits: UInt64) {
        guard index >= 0, index < count else { throw VivoOmicsError.invalid("count record index") }
        let byte = index * 16, next = (byte / Self.windowBytes) * Self.windowBytes
        if next != offset {
            try Task.checkCancellation()
            if let address { _ = munmap(address, length); self.address = nil }
            length = min(Self.windowBytes, count * 16 - next)
            let mapped = mmap(nil, length, PROT_READ, MAP_PRIVATE, fd, off_t(next))
            guard mapped != MAP_FAILED, let mapped else { throw VivoOmicsError.limit("cannot map count window") }
            address = mapped; offset = next
        }
        let p = address!, local = byte - offset
        let packed = UInt64(littleEndian: p.load(fromByteOffset: local, as: UInt64.self))
        return (Int(packed & 0xffff_ffff), Int(packed >> 32), UInt64(littleEndian: p.load(fromByteOffset: local + 8, as: UInt64.self)))
    }
}

public enum VivoH5ADCountStore {
    static let maximumEntries = 2_000_000_000
    static let maximumFileBytes = 64 * 1_024 * 1_024 * 1_024
    static let format = "source-major-coo-u32-u32-u64-le/v1"
    static var limits: VivoOmicsLimits {
        var v = VivoOmicsLimits(); v.maximumCells = 2_000_000; v.maximumFeatures = 200_000
        v.maximumNonzeros = maximumEntries; v.maximumInputBytes = maximumFileBytes; return v
    }
    static func staging(_ parent: URL) throws -> URL {
        let url = parent.appendingPathComponent(".numivivo-count-store-" + UUID().uuidString)
        try FileManager.default.createDirectory(at: url, withIntermediateDirectories: false, attributes: [.posixPermissions: 0o700]); return url
    }
    static func requireNew(_ url: URL) throws {
        guard url.isFileURL, (try? FileManager.default.attributesOfItem(atPath: url.path)) == nil else { throw VivoOmicsError.invalid("count store destination exists or is not local") }
    }
    static func fingerprint(_ source: URL, copyTo destination: URL? = nil) throws -> VivoFingerprint {
        try VivoOmicsFileSnapshot.fingerprint(source, copyTo: destination, maximumBytes: maximumFileBytes)
    }
    static func read(_ root: URL, _ name: String, maximum: Int) throws -> Data {
        try VivoSingleCellCampaignIO.readDocument(root.appendingPathComponent(name), maximumBytes: maximum)
    }
    private static func build(_ source: URL, plan: VivoH5ADImportPlan, output: URL?) throws -> (VivoSingleCellCountMetadata, VivoCountStoreQuality, Int, VivoFingerprint) {
        let writer = try VivoCountRecordWriter(output)
        var metadata: VivoSingleCellCountMetadata?
        var quality = VivoCountStoreQuality(rowTotals: [], rowNonzeros: [], featureTotals: [], featureNonzeros: [])
        _ = try VivoSingleCellH5AD.scanSnapshot(source, plan: plan, limits: limits, onMetadata: {
            metadata = $0
            quality = .init(rowTotals: .init(repeating: 0, count: $0.cells.count), rowNonzeros: .init(repeating: 0, count: $0.cells.count),
                featureTotals: .init(repeating: 0, count: $0.features.count), featureNonzeros: .init(repeating: 0, count: $0.features.count))
        }, onEntry: { row, column, value in
            quality.rowTotals[row] = try vivoOmicsSum(quality.rowTotals[row], value)
            quality.featureTotals[column] = try vivoOmicsSum(quality.featureTotals[column], value)
            quality.rowNonzeros[row] += 1; quality.featureNonzeros[column] += 1
            try writer.append(row: row, feature: column, bits: value)
        })
        guard let metadata else { throw VivoOmicsError.invalid("missing count metadata") }
        return (metadata, quality, writer.entries, try writer.finish())
    }
    public static func publish(source: URL, plan: VivoH5ADImportPlan, implementation: VivoFingerprint, to destination: URL) throws -> VivoH5ADCountStoreReceipt {
        try requireNew(destination)
        let temp = try staging(destination.deletingLastPathComponent()); defer { try? FileManager.default.removeItem(at: temp) }
        let planData = try VivoCanonicalJSON.encode(plan)
        guard planData.count <= 2_097_152 else { throw VivoOmicsError.limit("count store plan bytes") }
        let sourceHash = try fingerprint(source, copyTo: temp.appendingPathComponent("original.h5ad"))
        let (metadata, quality, entries, counts) = try build(temp.appendingPathComponent("original.h5ad"), plan: plan, output: temp.appendingPathComponent("counts.bin"))
        let m = try VivoCanonicalJSON.encode(metadata), q = try VivoCanonicalJSON.encode(quality)
        guard m.count <= 536_870_912, q.count <= 268_435_456 else { throw VivoOmicsError.limit("count store metadata or quality bytes") }
        for (name, data) in [("plan.json", planData), ("metadata.json", m), ("quality.json", q)] { try data.write(to: temp.appendingPathComponent(name), options: .withoutOverwriting) }
        let receipt = try VivoH5ADCountStoreReceipt(schemaVersion: 1, format: format, entries: entries, source: sourceHash,
            plan: VivoCanonicalJSON.fingerprint(planData), metadata: VivoCanonicalJSON.fingerprint(m), quality: VivoCanonicalJSON.fingerprint(q), counts: counts, implementation: implementation)
        try VivoCanonicalJSON.encode(receipt).write(to: temp.appendingPathComponent("receipt.json"), options: .withoutOverwriting)
        try Task.checkCancellation(); try FileManager.default.moveItem(at: temp, to: destination); return receipt
    }
    private static func verified<T>(_ root: URL, implementation: VivoFingerprint,
        body: (URL, VivoH5ADCountStoreReceipt, Data, VivoCountStoreQuality) throws -> T) throws -> T {
        let receiptData = try read(root, "receipt.json", maximum: 65_536)
        let receipt = try VivoCanonicalJSON.decode(VivoH5ADCountStoreReceipt.self, from: receiptData)
        guard receipt.schemaVersion == 1, receipt.format == format, receipt.entries >= 0, receipt.entries <= maximumEntries,
              receipt.implementation == implementation, try VivoCanonicalJSON.encode(receipt) == receiptData else { throw VivoOmicsError.invalid("count store receipt") }
        let temp = try staging(FileManager.default.temporaryDirectory); defer { try? FileManager.default.removeItem(at: temp) }
        let plan = try read(root, "plan.json", maximum: 2_097_152), metadata = try read(root, "metadata.json", maximum: 536_870_912), quality = try read(root, "quality.json", maximum: 268_435_456)
        guard try VivoCanonicalJSON.fingerprint(plan) == receipt.plan, try VivoCanonicalJSON.fingerprint(metadata) == receipt.metadata,
              try VivoCanonicalJSON.fingerprint(quality) == receipt.quality,
              try fingerprint(root.appendingPathComponent("original.h5ad"), copyTo: temp.appendingPathComponent("original.h5ad")) == receipt.source,
              try fingerprint(root.appendingPathComponent("counts.bin"), copyTo: temp.appendingPathComponent("counts.bin")) == receipt.counts else { throw VivoOmicsError.invalid("count store artifact fingerprint") }
        let mapping = try VivoCanonicalJSON.decode(VivoH5ADImportPlan.self, from: plan)
        let (rebuilt, qc, count, hash) = try build(temp.appendingPathComponent("original.h5ad"), plan: mapping, output: nil)
        guard count == receipt.entries, hash == receipt.counts, try VivoCanonicalJSON.encode(rebuilt) == metadata,
              try VivoCanonicalJSON.encode(qc) == quality else { throw VivoOmicsError.invalid("count store source reconstruction differs") }
        return try body(temp, receipt, metadata, qc)
    }
    public static func verify(_ directory: URL, implementation: VivoFingerprint) throws -> VivoH5ADCountStoreReceipt {
        try verified(directory, implementation: implementation) { _, receipt, _, _ in receipt }
    }
    public static func normalize(_ directory: URL, target: Double, implementation: VivoFingerprint, to destination: URL) throws -> VivoCountStoreNormalizationReceipt {
        guard target.isFinite, target > 0, target <= 1_000_000_000 else { throw VivoOmicsError.invalid("count store normalization target") }
        try requireNew(destination)
        return try verified(directory, implementation: implementation) { snapshot, receipt, metadata, quality in
            let temp = try staging(destination.deletingLastPathComponent()); defer { try? FileManager.default.removeItem(at: temp) }
            let records = try VivoWindowedCountRecords(snapshot.appendingPathComponent("counts.bin"), entries: receipt.entries)
            let writer = try VivoCountRecordWriter(temp.appendingPathComponent("values.bin"))
            for i in 0..<records.count {
                let record = try records.record(i)
                guard record.row < quality.rowTotals.count, record.feature < quality.featureTotals.count, record.bits > 0,
                      quality.rowTotals[record.row] > 0 else { throw VivoOmicsError.invalid("invalid mapped count coordinate") }
                let value = log1p(Double(record.bits) * (target / Double(quality.rowTotals[record.row])))
                guard value.isFinite, value > 0 else { throw VivoOmicsError.invalid("nonfinite normalized count") }
                try writer.append(row: record.row, feature: record.feature, bits: value.bitPattern)
            }
            let result = try VivoCountStoreNormalizationReceipt(schemaVersion: 1, format: "source-major-coo-u32-u32-f64-log1p-le/v1",
                target: target, entries: records.count, input: VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(receipt)),
                metadata: receipt.metadata, values: writer.finish(), implementation: implementation)
            try metadata.write(to: temp.appendingPathComponent("metadata.json"), options: .withoutOverwriting)
            try VivoCanonicalJSON.encode(receipt).write(to: temp.appendingPathComponent("input-receipt.json"), options: .withoutOverwriting)
            try VivoCanonicalJSON.encode(result).write(to: temp.appendingPathComponent("receipt.json"), options: .withoutOverwriting)
            try Task.checkCancellation(); try FileManager.default.moveItem(at: temp, to: destination); return result
        }
    }
    public static func verifyNormalization(_ normalized: URL, store: URL, implementation: VivoFingerprint) throws -> VivoCountStoreNormalizationReceipt {
        let bytes = try read(normalized, "receipt.json", maximum: 65_536)
        let receipt = try VivoCanonicalJSON.decode(VivoCountStoreNormalizationReceipt.self, from: bytes)
        guard receipt.schemaVersion == 1, receipt.implementation == implementation, try VivoCanonicalJSON.encode(receipt) == bytes else { throw VivoOmicsError.invalid("normalized count receipt") }
        let temp = try staging(FileManager.default.temporaryDirectory); defer { try? FileManager.default.removeItem(at: temp) }
        let rebuilt = try normalize(store, target: receipt.target, implementation: implementation, to: temp.appendingPathComponent("rebuilt"))
        guard rebuilt == receipt, try fingerprint(normalized.appendingPathComponent("values.bin")) == receipt.values,
              try VivoCanonicalJSON.fingerprint(read(normalized, "metadata.json", maximum: 536_870_912)) == receipt.metadata,
              try VivoCanonicalJSON.fingerprint(read(normalized, "input-receipt.json", maximum: 65_536)) == receipt.input else { throw VivoOmicsError.invalid("normalized count reconstruction differs") }
        return receipt
    }
}
