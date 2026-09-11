import Foundation
import CryptoKit

public struct VivoFilePseudobulkGroup: Codable, Sendable, Equatable {
    public let biologicalReplicateID: String
    public let donorID: String?
    public let condition: String
    public let organism: String
    public let cellGroup: String?
    public let sampleIDs: [String]
    public let batchIDs: [String]
    public let sourceCellCount: Int
}
public struct VivoFileCountStreamReport: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let method: String
    public let cellCount: Int
    public let canonicalNonzeros: Int
    public let countUnit: VivoOmicsCountUnit
    public let featureIDs: [String]
    public let groups: [VivoFilePseudobulkGroup]
    public let matrix: VivoSparseCounts
}
public struct VivoFileCountStreamReceipt: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let encoding: String
    public let qualityEncoding: String
    public let streamBytes: Int
    public let stream: VivoFingerprint
    public let axis: VivoFingerprint
    public let quality: VivoFingerprint
    public let report: VivoFingerprint
    public let implementation: VivoFingerprint
}

/// Row identities, QC and memberships stay on disk; small dictionaries and
/// at most five million aggregate coordinates remain resident. No raw matrix is retained.
public enum VivoFileCountStream {
    public static let qualityEncoding = "total-u64-mito-u64-nnz-u32-group-u32-le-v1"
    static let maximumReportBytes = 268_435_456
    private struct Key: Hashable {
        let replicate: String; let condition: String; let group: String?
        static func ordered(_ a: Key, _ b: Key) -> Bool {
            if a.replicate != b.replicate { return a.replicate < b.replicate }
            if a.condition != b.condition { return a.condition < b.condition }
            if a.group == nil { return b.group != nil }; if b.group == nil { return false }
            return a.group! < b.group!
        }
    }
    private struct Membership { var samples: Set<Int> = []; var cells = 0 }
    private static func key(_ row: VivoCellAxisRow, _ axis: VivoFileCellAxis) -> Key {
        let sample = axis.header.metadata.samples[row.sampleIndex]
        return .init(replicate: sample.biologicalReplicateID, condition: sample.condition, group: row.cell.group)
    }
    static func evaluate(axis: VivoFileCellAxis, qualityURL: URL, read: () throws -> Data) throws -> (VivoFileCountStreamReport, VivoFingerprint, Int) {
        var memberships: [Key: Membership] = [:], associations = 0
        for i in 0..<axis.header.cellCount {
            try vivoAxisPool {
                try Task.checkCancellation()
                let row = try axis.row(i), key = key(row, axis)
                if memberships[key] == nil {
                    guard memberships.count < 100_000 else { throw VivoOmicsError.limit("file count stream exceeds 100000 aggregates") }
                    memberships[key] = Membership()
                }
                if memberships[key]!.samples.insert(row.sampleIndex).inserted {
                    associations += 1
                    guard associations <= 1_000_000 else { throw VivoOmicsError.limit("file count stream sample/group associations") }
                }
                memberships[key]!.cells += 1
            }
        }
        let keys = memberships.keys.sorted(by: Key.ordered)
        let groupIndex = Dictionary(uniqueKeysWithValues: keys.enumerated().map { ($0.element, $0.offset) })
        let groups: [VivoFilePseudobulkGroup] = keys.map { key in
            let member = memberships[key]!, samples = member.samples.map { axis.header.metadata.samples[$0] }.sorted { $0.id < $1.id }
            let first = samples[0]
            return .init(biologicalReplicateID: key.replicate, donorID: first.donorID, condition: key.condition,
                         organism: first.organism, cellGroup: key.group, sampleIDs: samples.map(\.id),
                         batchIDs: Set(samples.map(\.batchID)).sorted(), sourceCellCount: member.cells)
        }
        memberships.removeAll()
        var sums: [[Int: UInt64]] = Array(repeating: [:], count: groups.count), aggregateNonzeros = 0
        let writer = try VivoAxisWriter(qualityURL)
        var nextRow = 0, current: VivoCellAxisRow?, currentGroup = 0, total: UInt64 = 0, mitochondrial: UInt64 = 0
        var detected = 0, previousFeature = -1, entries = 0, buffer = Data(), digest = SHA256()
        func startRow() throws {
            let row = try axis.row(nextRow); current = row; currentGroup = groupIndex[key(row, axis)]!
        }
        func finishRow() throws {
            if current == nil { try startRow() }
            let row = current!
            guard detected == row.nonzeros, row.totalCounts.map({ $0 == total }) ?? true else {
                throw VivoOmicsError.invalid("file count stream differs from declared row cardinality or total")
            }
            var bytes = Data(); bytes.vivoAppendLE(total); bytes.vivoAppendLE(mitochondrial)
            bytes.vivoAppendLE(UInt32(detected)); bytes.vivoAppendLE(UInt32(currentGroup)); try writer.append(bytes)
            nextRow += 1; current = nil; total = 0; mitochondrial = 0; detected = 0; previousFeature = -1
        }
        while true {
            let ended = try vivoAxisPool {
                try Task.checkCancellation()
                let chunk = try read()
                guard chunk.count <= VivoCountStreamPseudobulk.maximumChunkBytes else { throw VivoOmicsError.limit("count read exceeds 1 MiB") }
                if chunk.isEmpty { return true }
                guard chunk.count <= axis.header.nonzeros * 16 - entries * 16 - buffer.count else { throw VivoOmicsError.invalid("file count stream exceeds declared records") }
                digest.update(data: chunk); buffer.append(chunk)
                let complete = buffer.count / 16 * 16
                try buffer.withUnsafeBytes { bytes in
                    for offset in stride(from: 0, to: complete, by: 16) {
                        let packed = UInt64(littleEndian: bytes.loadUnaligned(fromByteOffset: offset, as: UInt64.self))
                        let row = Int(packed & 0xffff_ffff), feature = Int(packed >> 32)
                        let count = UInt64(littleEndian: bytes.loadUnaligned(fromByteOffset: offset + 8, as: UInt64.self))
                        guard row >= nextRow, row < axis.header.cellCount, feature < axis.header.metadata.features.count, count > 0 else {
                            throw VivoOmicsError.invalid("file count stream row, feature or positive-count requirement")
                        }
                        while nextRow < row { try finishRow() }
                        if current == nil { try startRow() }
                        guard feature > previousFeature else { throw VivoOmicsError.invalid("file count stream requires sorted unique columns") }
                        total = try vivoOmicsSum(total, count)
                        if axis.header.metadata.features[feature].mitochondrial { mitochondrial = try vivoOmicsSum(mitochondrial, count) }
                        if sums[currentGroup][feature] == nil {
                            guard aggregateNonzeros < VivoH5ADPseudobulk.maximumAggregateNonzeros else { throw VivoOmicsError.limit("file count stream aggregate nonzero limit") }
                            aggregateNonzeros += 1
                        }
                        sums[currentGroup][feature] = try vivoOmicsSum(sums[currentGroup][feature] ?? 0, count)
                        detected += 1; entries += 1; previousFeature = feature
                    }
                }
                buffer = Data(buffer.suffix(buffer.count - complete)); return false
            }
            if ended { break }
        }
        guard buffer.isEmpty, entries == axis.header.nonzeros else { throw VivoOmicsError.invalid("file count stream is truncated") }
        while nextRow < axis.header.cellCount { try vivoAxisPool { try finishRow() } }
        try writer.finish()
        var offsets = [0], columns: [Int] = [], counts: [UInt64] = []
        for group in sums {
            for feature in group.keys.sorted() { columns.append(feature); counts.append(group[feature]!) }
            offsets.append(counts.count)
        }
        let matrix = VivoSparseCounts(cellCount: groups.count, featureCount: axis.header.metadata.features.count,
                                      rowOffsets: offsets, featureIndices: columns, counts: counts)
        var limits = VivoH5ADPseudobulk.sourceLimits; limits.maximumNonzeros = VivoH5ADPseudobulk.maximumAggregateNonzeros
        try matrix.validate(limits: limits)
        let report = VivoFileCountStreamReport(schemaVersion: 1, method: "file-cell-axis-count-stream-pseudobulk-v1",
            cellCount: axis.header.cellCount, canonicalNonzeros: entries, countUnit: axis.header.metadata.countUnit,
            featureIDs: axis.header.metadata.features.map(\.id), groups: groups, matrix: matrix)
        return (report, try VivoFingerprint(bytes: Array(digest.finalize())), entries * 16)
    }
    public static func publish(axis directory: URL, input: FileHandle, implementation: VivoFingerprint,
                               to destination: URL) throws -> VivoFileCountStreamReceipt {
        guard !FileManager.default.fileExists(atPath: destination.path) else { throw VivoOmicsError.invalid("file count output exists") }
        let axis = try VivoFileCellAxis.open(directory, implementation: implementation)
        let staging = destination.deletingLastPathComponent().appendingPathComponent(".numivivo-file-count-" + UUID().uuidString)
        try FileManager.default.createDirectory(at: staging, withIntermediateDirectories: false, attributes: [.posixPermissions: 0o700])
        defer { try? FileManager.default.removeItem(at: staging) }
        let (report, stream, bytes) = try evaluate(axis: axis, qualityURL: staging.appendingPathComponent("quality.bin")) {
            try input.read(upToCount: VivoCountStreamPseudobulk.maximumChunkBytes) ?? Data()
        }
        let raw = try VivoCanonicalJSON.encode(report)
        guard raw.count <= maximumReportBytes else { throw VivoOmicsError.limit("file count aggregate report size") }
        try raw.write(to: staging.appendingPathComponent("report.json"), options: .withoutOverwriting)
        try axis.copy(to: staging.appendingPathComponent("axis", isDirectory: true))
        let receipt = try VivoFileCountStreamReceipt(schemaVersion: 1, encoding: VivoCountStreamPseudobulk.encoding,
            qualityEncoding: qualityEncoding, streamBytes: bytes, stream: stream,
            axis: VivoOmicsFileSnapshot.fingerprint(axis.snapshot.appendingPathComponent("receipt.json"), maximumBytes: 65_536),
            quality: VivoOmicsFileSnapshot.fingerprint(staging.appendingPathComponent("quality.bin"), maximumBytes: axis.header.cellCount * 24),
            report: VivoCanonicalJSON.fingerprint(raw), implementation: implementation)
        try VivoCanonicalJSON.encode(receipt).write(to: staging.appendingPathComponent("receipt.json"), options: .withoutOverwriting)
        try Task.checkCancellation(); try FileManager.default.moveItem(at: staging, to: destination); return receipt
    }
    public static func verify(_ directory: URL, input: FileHandle, implementation: VivoFingerprint) throws -> VivoFileCountStreamReceipt {
        let saved = try VivoFileCountSnapshot.open(directory, implementation: implementation)
        let scratch = FileManager.default.temporaryDirectory.appendingPathComponent("numivivo-file-count-replay-" + UUID().uuidString)
        try FileManager.default.createDirectory(at: scratch, withIntermediateDirectories: false)
        defer { try? FileManager.default.removeItem(at: scratch) }
        let quality = scratch.appendingPathComponent("quality.bin")
        let (report, stream, bytes) = try evaluate(axis: saved.axis, qualityURL: quality) {
            try input.read(upToCount: VivoCountStreamPseudobulk.maximumChunkBytes) ?? Data()
        }
        guard stream == saved.receipt.stream, bytes == saved.receipt.streamBytes,
              try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(report)) == saved.receipt.report,
              try VivoOmicsFileSnapshot.fingerprint(quality, maximumBytes: saved.axis.header.cellCount * 24) == saved.receipt.quality else {
            throw VivoOmicsError.invalid("file count report or QC does not reconstruct")
        }
        return saved.receipt
    }
}

/// Snapshot reads validate hashes and shapes. Full numerical reconstruction
/// additionally requires VivoFileCountStream.verify with the original stream.
public final class VivoFileCountSnapshot {
    public let axis: VivoFileCellAxis
    public let report: VivoFileCountStreamReport
    public let receipt: VivoFileCountStreamReceipt
    private let root: URL
    private let qualityFile: VivoAxisBinaryFile
    private let mitochondrialFeatures: Int
    private init(root: URL, axis: VivoFileCellAxis, report: VivoFileCountStreamReport, receipt: VivoFileCountStreamReceipt) throws {
        self.root = root; self.axis = axis; self.report = report; self.receipt = receipt
        qualityFile = try VivoAxisBinaryFile(root.appendingPathComponent("quality.bin"), bytes: axis.header.cellCount * 24)
        mitochondrialFeatures = axis.header.metadata.features.filter(\.mitochondrial).count
    }
    deinit { try? FileManager.default.removeItem(at: root) }
    func receiptFingerprint() throws -> VivoFingerprint {
        try VivoOmicsFileSnapshot.fingerprint(root.appendingPathComponent("receipt.json"), maximumBytes: 65_536)
    }
    /// Copy the owned snapshots, never reopen mutable source paths during analysis publication.
    func copy(to destination: URL) throws {
        guard !FileManager.default.fileExists(atPath: destination.path) else { throw VivoOmicsError.invalid("file count copy destination exists") }
        try FileManager.default.createDirectory(at: destination, withIntermediateDirectories: false, attributes: [.posixPermissions: 0o700])
        var accepted = false
        defer { if !accepted { try? FileManager.default.removeItem(at: destination) } }
        for (name, limit) in [("receipt.json", 65_536), ("report.json", VivoFileCountStream.maximumReportBytes), ("quality.bin", axis.header.cellCount * 24)] {
            _ = try VivoOmicsFileSnapshot.fingerprint(root.appendingPathComponent(name), copyTo: destination.appendingPathComponent(name), maximumBytes: limit)
        }
        try axis.copy(to: destination.appendingPathComponent("axis")); accepted = true
    }
    public func quality(_ row: Int) throws -> (quality: VivoCellQuality, group: Int) {
        let cell = try axis.row(row), bytes = try qualityFile.read(row * 24, 24)
        let total = bytes.vivoLE(UInt64.self, at: 0), mito = bytes.vivoLE(UInt64.self, at: 8)
        let nnz = Int(bytes.vivoLE(UInt32.self, at: 16)), group = Int(bytes.vivoLE(UInt32.self, at: 20))
        guard group < report.groups.count, nnz == cell.nonzeros, (total == 0) == (nnz == 0),
              UInt64(nnz) <= total, mito <= total, mitochondrialFeatures > 0 || mito == 0,
              cell.totalCounts.map({ $0 == total }) ?? true else { throw VivoOmicsError.invalid("file quality row") }
        let sample = axis.header.metadata.samples[cell.sampleIndex], info = report.groups[group]
        guard info.biologicalReplicateID == sample.biologicalReplicateID, info.condition == sample.condition,
              info.cellGroup == cell.cell.group, info.donorID == sample.donorID, info.organism == sample.organism,
              info.sampleIDs.contains(sample.id), info.batchIDs.contains(sample.batchID) else { throw VivoOmicsError.invalid("file quality membership") }
        return (.init(sampleID: cell.cell.sampleID, barcode: cell.cell.barcode, totalCounts: total, detectedFeatures: nnz,
                      mitochondrialCounts: mito, mitochondrialFeatureCount: mitochondrialFeatures,
                      mitochondrialFraction: total == 0 || mitochondrialFeatures == 0 ? nil : Double(mito) / Double(total)), group)
    }
    public static func open(_ directory: URL, implementation: VivoFingerprint) throws -> VivoFileCountSnapshot {
        let root = FileManager.default.temporaryDirectory.appendingPathComponent("numivivo-file-count-read-" + UUID().uuidString)
        try FileManager.default.createDirectory(at: root, withIntermediateDirectories: false, attributes: [.posixPermissions: 0o700])
        var accepted = false
        defer { if !accepted { try? FileManager.default.removeItem(at: root) } }
        func copy(_ name: String, _ limit: Int) throws -> VivoFingerprint {
            try VivoOmicsFileSnapshot.fingerprint(directory.appendingPathComponent(name), copyTo: root.appendingPathComponent(name), maximumBytes: limit)
        }
        _ = try copy("receipt.json", 65_536)
        let receipt = try VivoCanonicalJSON.decode(VivoFileCountStreamReceipt.self, from: Data(contentsOf: root.appendingPathComponent("receipt.json")))
        guard receipt.schemaVersion == 1, receipt.encoding == VivoCountStreamPseudobulk.encoding,
              receipt.qualityEncoding == VivoFileCountStream.qualityEncoding, receipt.implementation == implementation else { throw VivoOmicsError.invalid("file count receipt") }
        let axis = try VivoFileCellAxis.open(directory.appendingPathComponent("axis"), implementation: implementation)
        guard try VivoOmicsFileSnapshot.fingerprint(axis.snapshot.appendingPathComponent("receipt.json"), maximumBytes: 65_536) == receipt.axis,
              try copy("report.json", VivoFileCountStream.maximumReportBytes) == receipt.report,
              try copy("quality.bin", axis.header.cellCount * 24) == receipt.quality else { throw VivoOmicsError.invalid("file count payload fingerprint") }
        let report = try VivoCanonicalJSON.decode(VivoFileCountStreamReport.self, from: Data(contentsOf: root.appendingPathComponent("report.json")))
        guard report.schemaVersion == 1, report.method == "file-cell-axis-count-stream-pseudobulk-v1",
              report.cellCount == axis.header.cellCount, report.canonicalNonzeros == axis.header.nonzeros,
              receipt.streamBytes == axis.header.nonzeros * 16, report.featureIDs == axis.header.metadata.features.map(\.id),
              report.countUnit == axis.header.metadata.countUnit, report.groups.count <= 100_000,
              report.matrix.cellCount == report.groups.count, report.matrix.featureCount == report.featureIDs.count,
              report.groups.allSatisfy({ $0.sourceCellCount > 0 && $0.sourceCellCount <= axis.header.cellCount }),
              report.groups.reduce(0, { $0 + $1.sourceCellCount }) == axis.header.cellCount else { throw VivoOmicsError.invalid("file count report axes") }
        var limits = VivoH5ADPseudobulk.sourceLimits; limits.maximumNonzeros = VivoH5ADPseudobulk.maximumAggregateNonzeros
        try report.matrix.validate(limits: limits)
        let result = try VivoFileCountSnapshot(root: root, axis: axis, report: report, receipt: receipt); accepted = true; return result
    }
}
