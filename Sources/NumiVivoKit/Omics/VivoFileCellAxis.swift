import Foundation
import CryptoKit
#if canImport(Darwin)
import Darwin
#else
import Glibc
#endif

/// Small resident dictionaries. The cell axis itself is never decoded into an array.
public struct VivoCellAxisHeader: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let metadata: VivoSingleCellCountMetadata
    public let cellCount: Int
    public let nonzeros: Int
    public let hasRowTotals: Bool
    public let annotations: [String]
    public let sourceDeclaration: String
    public init(metadata: VivoSingleCellCountMetadata, cellCount: Int, nonzeros: Int,
                hasRowTotals: Bool, annotations: [String], sourceDeclaration: String) {
        schemaVersion = 1; self.metadata = metadata; self.cellCount = cellCount; self.nonzeros = nonzeros
        self.hasRowTotals = hasRowTotals; self.annotations = annotations; self.sourceDeclaration = sourceDeclaration
    }
    private enum CodingKeys: String, CodingKey { case schemaVersion, metadata, cellCount, nonzeros, hasRowTotals, annotations, sourceDeclaration }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["schemaVersion", "metadata", "cellCount", "nonzeros", "hasRowTotals", "annotations", "sourceDeclaration"])
        let c = try decoder.container(keyedBy: CodingKeys.self)
        schemaVersion = try c.decode(Int.self, forKey: .schemaVersion)
        metadata = try c.decode(VivoSingleCellCountMetadata.self, forKey: .metadata)
        cellCount = try c.decode(Int.self, forKey: .cellCount); nonzeros = try c.decode(Int.self, forKey: .nonzeros)
        hasRowTotals = try c.decode(Bool.self, forKey: .hasRowTotals)
        annotations = try c.decode([String].self, forKey: .annotations)
        sourceDeclaration = try c.decode(String.self, forKey: .sourceDeclaration)
    }
    public func validate() throws {
        var limits = VivoH5ADPseudobulk.sourceLimits; limits.maximumCells = 20_000_000
        try metadata.validate(limits: limits)
        guard schemaVersion == 1, metadata.cells.isEmpty, (0...20_000_000).contains(cellCount),
              (0...40_000_000_000).contains(nonzeros), metadata.samples.count <= 100_000,
              annotations.count <= 100_000, annotations.allSatisfy(vivoOmicsID),
              Set(annotations).count == annotations.count,
              !sourceDeclaration.isEmpty, sourceDeclaration.utf8.count <= 16_384 else {
            throw VivoOmicsError.invalid("file cell axis header, bounds or dictionaries")
        }
    }
}

/// One JSONL observation. Optional totals refer to the retained matrix, not historical QC.
public struct VivoCellAxisInputRow: Codable, Sendable, Equatable {
    public let barcode: String
    public let sampleID: String
    public let group: String?
    public let nonzeros: Int
    public let totalCounts: UInt64?
    public init(cell: VivoOmicsCell, nonzeros: Int, totalCounts: UInt64?) {
        barcode = cell.barcode; sampleID = cell.sampleID; group = cell.group
        self.nonzeros = nonzeros; self.totalCounts = totalCounts
    }
    private enum CodingKeys: String, CodingKey { case barcode, sampleID, group, nonzeros, totalCounts }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["barcode", "sampleID", "group", "nonzeros", "totalCounts"])
        let c = try decoder.container(keyedBy: CodingKeys.self)
        barcode = try c.decode(String.self, forKey: .barcode); sampleID = try c.decode(String.self, forKey: .sampleID)
        group = try c.decodeIfPresent(String.self, forKey: .group); nonzeros = try c.decode(Int.self, forKey: .nonzeros)
        totalCounts = try c.decodeIfPresent(UInt64.self, forKey: .totalCounts)
    }
}
public struct VivoCellAxisRow: Sendable, Equatable {
    public let cell: VivoOmicsCell
    public let sampleIndex: Int
    public let annotationIndex: Int?
    public let nonzeros: Int
    public let totalCounts: UInt64?
}
public struct VivoCellAxisReceipt: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let encoding: String
    public let header: VivoFingerprint
    public let rows: VivoFingerprint
    public let strings: VivoFingerprint
    public let stringBytes: Int
    public let inputJSONL: VivoFingerprint
    public let implementation: VivoFingerprint
}

func vivoAxisPool<T>(_ body: () throws -> T) rethrows -> T {
    #if canImport(ObjectiveC)
    return try autoreleasepool(invoking: body)
    #else
    return try body()
    #endif
}
extension Data {
    mutating func vivoAppendLE<T: FixedWidthInteger>(_ value: T) {
        var bits = value.littleEndian
        Swift.withUnsafeBytes(of: &bits) { append(contentsOf: $0) }
    }
    func vivoLE<T: FixedWidthInteger>(_ type: T.Type, at offset: Int) -> T {
        withUnsafeBytes { T(littleEndian: $0.loadUnaligned(fromByteOffset: offset, as: T.self)) }
    }
}
final class VivoAxisWriter {
    let file: FileHandle
    var buffer = Data()
    var closed = false
    init(_ path: URL) throws {
        guard FileManager.default.createFile(atPath: path.path, contents: nil, attributes: [.posixPermissions: 0o600]) else {
            throw VivoOmicsError.invalid("cannot create axis output")
        }
        file = try FileHandle(forWritingTo: path)
    }
    deinit { if !closed { try? file.close() } }
    func append(_ data: Data) throws {
        buffer.append(data)
        if buffer.count >= 65_536 { try file.write(contentsOf: buffer); buffer.removeAll(keepingCapacity: true) }
    }
    func finish() throws {
        try file.write(contentsOf: buffer); buffer.removeAll(); try file.synchronize(); try file.close(); closed = true
    }
}
/// Descriptor reads stay bounded and cannot follow a replaced input symlink.
final class VivoAxisBinaryFile {
    let fd: Int32
    let bytes: Int
    init(_ path: URL, bytes: Int, scratch: Bool = false) throws {
        guard bytes >= 0 else { throw VivoOmicsError.invalid("negative axis file extent") }
        let descriptor = open(path.path, (scratch ? O_RDWR | O_CREAT | O_EXCL : O_RDONLY) | O_NOFOLLOW | O_NONBLOCK, 0o600)
        guard descriptor >= 0 else { throw VivoOmicsError.invalid("cannot open axis file") }
        var accepted = false
        defer { if !accepted { _ = close(descriptor) } }
        if scratch { guard ftruncate(descriptor, off_t(bytes)) == 0 else { throw VivoOmicsError.invalid("axis index allocation") } }
        var info = stat()
        guard fstat(descriptor, &info) == 0, info.st_mode & mode_t(S_IFMT) == mode_t(S_IFREG), info.st_size == bytes else {
            throw VivoOmicsError.invalid("axis file type or extent")
        }
        fd = descriptor; self.bytes = bytes; accepted = true
    }
    deinit { _ = close(fd) }
    func read(_ offset: Int, _ count: Int) throws -> Data {
        guard offset >= 0, count >= 0, offset <= bytes, count <= bytes - offset else { throw VivoOmicsError.invalid("axis read bounds") }
        var data = Data(count: count)
        try data.withUnsafeMutableBytes { raw in
            var done = 0
            while done < count {
                let n = pread(fd, raw.baseAddress!.advanced(by: done), count - done, off_t(offset + done))
                if n < 0 && errno == EINTR { continue }
                guard n > 0 else { throw VivoOmicsError.invalid("axis read failed") }; done += n
            }
        }
        return data
    }
    func writeSlot(_ index: Int, _ value: UInt64) throws {
        guard index >= 0, index < bytes / 8 else { throw VivoOmicsError.invalid("axis index bounds") }
        var bits = value.littleEndian
        try Swift.withUnsafeBytes(of: &bits) { raw in
            var done = 0
            while done < 8 {
                let n = pwrite(fd, raw.baseAddress!.advanced(by: done), 8 - done, off_t(index * 8 + done))
                if n < 0 && errno == EINTR { continue }
                guard n > 0 else { throw VivoOmicsError.invalid("axis index write failed") }; done += n
            }
        }
    }
}

/// Owns a private immutable snapshot. Identity uniqueness uses a transient disk
/// hash table with exact Swift string comparison, including Unicode equivalence.
public final class VivoFileCellAxis {
    public static let recordBytes = 40
    public static let encoding = "string-offset-u64-barcode-length-u32-sample-length-u32-group-length-u32-sample-u32-annotation-u32-nnz-u32-total-u64-le-v1"
    public static let maximumHeaderBytes = 67_108_864
    public let header: VivoCellAxisHeader
    public let receipt: VivoCellAxisReceipt
    let snapshot: URL
    private let rows: VivoAxisBinaryFile
    private let strings: VivoAxisBinaryFile
    private init(snapshot: URL, header: VivoCellAxisHeader, receipt: VivoCellAxisReceipt) throws {
        self.snapshot = snapshot; self.header = header; self.receipt = receipt
        rows = try VivoAxisBinaryFile(snapshot.appendingPathComponent("rows.bin"), bytes: header.cellCount * Self.recordBytes)
        strings = try VivoAxisBinaryFile(snapshot.appendingPathComponent("strings.bin"), bytes: receipt.stringBytes)
    }
    deinit { try? FileManager.default.removeItem(at: snapshot) }
    public func row(_ index: Int) throws -> VivoCellAxisRow {
        guard (0..<header.cellCount).contains(index) else { throw VivoOmicsError.invalid("cell row outside axis") }
        let data = try rows.read(index * Self.recordBytes, Self.recordBytes)
        let offset = data.vivoLE(UInt64.self, at: 0), barcodeLength = Int(data.vivoLE(UInt32.self, at: 8))
        let sampleLength = Int(data.vivoLE(UInt32.self, at: 12)), groupLength = Int(data.vivoLE(UInt32.self, at: 16))
        let sample = Int(data.vivoLE(UInt32.self, at: 20)), annotation = data.vivoLE(UInt32.self, at: 24)
        let nnz = Int(data.vivoLE(UInt32.self, at: 28)), total = data.vivoLE(UInt64.self, at: 32)
        guard offset <= UInt64(receipt.stringBytes), (1...1024).contains(barcodeLength), (1...1024).contains(sampleLength),
              (0...1024).contains(groupLength), sample < header.metadata.samples.count,
              annotation == UInt32.max || Int(annotation) < header.annotations.count,
              (annotation == UInt32.max) == (groupLength == 0), nnz <= header.metadata.features.count,
              header.hasRowTotals ? (UInt64(nnz) <= total && (nnz == 0) == (total == 0)) : total == 0 else {
            throw VivoOmicsError.invalid("invalid file cell row extent")
        }
        let bytes = try strings.read(Int(offset), barcodeLength + sampleLength + groupLength)
        func string(_ start: Int, _ length: Int) throws -> String {
            guard let value = String(data: bytes.subdata(in: start..<(start + length)), encoding: .utf8), vivoOmicsID(value) else {
                throw VivoOmicsError.invalid("invalid cell identity string")
            }
            return value
        }
        let barcode = try string(0, barcodeLength), sampleID = try string(barcodeLength, sampleLength)
        let group = groupLength == 0 ? nil : try string(barcodeLength + sampleLength, groupLength)
        let label = annotation == UInt32.max ? nil : Int(annotation)
        guard sampleID == header.metadata.samples[sample].id, group == label.map({ header.annotations[$0] }) else {
            throw VivoOmicsError.invalid("cell identity dictionary reference")
        }
        return .init(cell: .init(barcode: barcode, sampleID: sampleID, group: group), sampleIndex: sample,
                     annotationIndex: label, nonzeros: nnz, totalCounts: header.hasRowTotals ? total : nil)
    }
    private func validateRows() throws {
        var capacity = 1
        while capacity < max(2, header.cellCount * 2) { capacity *= 2 }
        let path = snapshot.appendingPathComponent(".identity-index")
        defer { try? FileManager.default.removeItem(at: path) }
        let index = try VivoAxisBinaryFile(path, bytes: capacity * 8, scratch: true)
        var stringOffset: UInt64 = 0, entries = 0
        for i in 0..<header.cellCount {
            try vivoAxisPool {
                try Task.checkCancellation()
                let record = try row(i), raw = try rows.read(i * Self.recordBytes, Self.recordBytes)
                guard raw.vivoLE(UInt64.self, at: 0) == stringOffset,
                      record.nonzeros <= header.nonzeros - entries else { throw VivoOmicsError.invalid("axis offsets or record total") }
                stringOffset += UInt64(record.cell.barcode.utf8.count + record.cell.sampleID.utf8.count + (record.cell.group?.utf8.count ?? 0)); entries += record.nonzeros
                var hash = Hasher(); hash.combine(record.sampleIndex); hash.combine(record.cell.barcode)
                var slot = Int(UInt(bitPattern: hash.finalize()) & UInt(capacity - 1)), inserted = false
                for _ in 0..<capacity {
                    let previous = try index.read(slot * 8, 8).vivoLE(UInt64.self, at: 0)
                    if previous == 0 { try index.writeSlot(slot, UInt64(i + 1)); inserted = true; break }
                    guard previous <= UInt64(i) else { throw VivoOmicsError.invalid("axis uniqueness index") }
                    let other = try row(Int(previous - 1))
                    guard other.sampleIndex != record.sampleIndex || other.cell.barcode != record.cell.barcode else {
                        throw VivoOmicsError.invalid("duplicate sample/barcode identity")
                    }
                    slot = (slot + 1) & (capacity - 1)
                }
                guard inserted else { throw VivoOmicsError.invalid("axis identity index exhausted") }
            }
        }
        guard stringOffset == receipt.stringBytes, entries == header.nonzeros else { throw VivoOmicsError.invalid("axis terminal extents") }
    }
    public static func open(_ directory: URL, implementation: VivoFingerprint) throws -> VivoFileCellAxis {
        let snapshot = FileManager.default.temporaryDirectory.appendingPathComponent("numivivo-cell-axis-" + UUID().uuidString)
        try FileManager.default.createDirectory(at: snapshot, withIntermediateDirectories: false, attributes: [.posixPermissions: 0o700])
        var accepted = false
        defer { if !accepted { try? FileManager.default.removeItem(at: snapshot) } }
        func copy(_ name: String, _ limit: Int) throws -> VivoFingerprint {
            try VivoOmicsFileSnapshot.fingerprint(directory.appendingPathComponent(name), copyTo: snapshot.appendingPathComponent(name), maximumBytes: limit)
        }
        _ = try copy("receipt.json", 65_536)
        let receipt = try VivoCanonicalJSON.decode(VivoCellAxisReceipt.self, from: Data(contentsOf: snapshot.appendingPathComponent("receipt.json")))
        guard receipt.schemaVersion == 1, receipt.encoding == encoding, receipt.implementation == implementation,
              receipt.stringBytes >= 0, receipt.stringBytes <= 20_000_000 * 3072 else { throw VivoOmicsError.invalid("axis receipt") }
        guard try copy("header.json", maximumHeaderBytes) == receipt.header else { throw VivoOmicsError.invalid("axis header fingerprint") }
        let header = try VivoCanonicalJSON.decode(VivoCellAxisHeader.self, from: Data(contentsOf: snapshot.appendingPathComponent("header.json")))
        try header.validate()
        guard try copy("rows.bin", header.cellCount * Self.recordBytes) == receipt.rows,
              try copy("strings.bin", receipt.stringBytes) == receipt.strings else { throw VivoOmicsError.invalid("axis payload fingerprint") }
        let axis = try VivoFileCellAxis(snapshot: snapshot, header: header, receipt: receipt)
        try axis.validateRows(); accepted = true; return axis
    }
    func copy(to destination: URL) throws {
        try FileManager.default.createDirectory(at: destination, withIntermediateDirectories: false, attributes: [.posixPermissions: 0o700])
        for name in ["header.json", "receipt.json", "rows.bin", "strings.bin"] {
            _ = try VivoOmicsFileSnapshot.fingerprint(snapshot.appendingPathComponent(name), copyTo: destination.appendingPathComponent(name), maximumBytes: 64 * 1_024 * 1_024 * 1_024)
        }
    }
    public static func publish(header: VivoCellAxisHeader, input: FileHandle, implementation: VivoFingerprint,
                               to destination: URL) throws -> VivoCellAxisReceipt {
        try header.validate()
        let headerBytes = try VivoCanonicalJSON.encode(header)
        guard headerBytes.count <= maximumHeaderBytes, !FileManager.default.fileExists(atPath: destination.path) else { throw VivoOmicsError.invalid("axis header size or existing output") }
        let staging = destination.deletingLastPathComponent().appendingPathComponent(".numivivo-axis-" + UUID().uuidString)
        try FileManager.default.createDirectory(at: staging, withIntermediateDirectories: false, attributes: [.posixPermissions: 0o700])
        defer { try? FileManager.default.removeItem(at: staging) }
        let rows = try VivoAxisWriter(staging.appendingPathComponent("rows.bin")), strings = try VivoAxisWriter(staging.appendingPathComponent("strings.bin"))
        let samples = Dictionary(uniqueKeysWithValues: header.metadata.samples.enumerated().map { ($0.element.id, $0.offset) })
        let annotations = Dictionary(uniqueKeysWithValues: header.annotations.enumerated().map { ($0.element, $0.offset) })
        var pending = Data(), hash = SHA256(), count = 0, offset = 0
        while true {
            let ended = try vivoAxisPool {
                try Task.checkCancellation()
                let chunk = try input.read(upToCount: 1_048_576) ?? Data()
                if chunk.isEmpty { return true }
                hash.update(data: chunk); pending.append(chunk)
                var start = pending.startIndex
                while let end = pending[start...].firstIndex(of: 10) {
                    guard end - start <= 16_384, end > start, count < header.cellCount else { throw VivoOmicsError.invalid("axis JSONL line or cell bound") }
                    let row = try VivoCanonicalJSON.decode(VivoCellAxisInputRow.self, from: Data(pending[start..<end]))
                    guard vivoOmicsID(row.barcode), vivoOmicsID(row.sampleID), row.group.map(vivoOmicsID) ?? true, let sample = samples[row.sampleID],
                          row.group.map({ annotations[$0] != nil }) ?? true,
                          row.nonzeros >= 0, row.nonzeros <= header.metadata.features.count,
                          (row.totalCounts != nil) == header.hasRowTotals else { throw VivoOmicsError.invalid("axis JSONL identity or counts") }
                    let barcode = Data(row.barcode.utf8), sampleName = Data(row.sampleID.utf8), groupName = Data((row.group ?? "").utf8)
                    var text = barcode; text.append(sampleName); text.append(groupName)
                    var record = Data(); record.vivoAppendLE(UInt64(offset))
                    record.vivoAppendLE(UInt32(barcode.count)); record.vivoAppendLE(UInt32(sampleName.count)); record.vivoAppendLE(UInt32(groupName.count))
                    record.vivoAppendLE(UInt32(sample)); record.vivoAppendLE(row.group.map { UInt32(annotations[$0]!) } ?? UInt32.max)
                    record.vivoAppendLE(UInt32(row.nonzeros)); record.vivoAppendLE(row.totalCounts ?? 0)
                    try rows.append(record); try strings.append(text); offset += text.count; count += 1
                    start = end + 1
                }
                pending = Data(pending[start...])
                guard pending.count <= 16_384 else { throw VivoOmicsError.invalid("axis JSONL line exceeds bound") }
                return false
            }
            if ended { break }
        }
        guard pending.isEmpty, count == header.cellCount else { throw VivoOmicsError.invalid("axis JSONL requires complete newline-terminated rows") }
        try rows.finish(); try strings.finish()
        try headerBytes.write(to: staging.appendingPathComponent("header.json"), options: .withoutOverwriting)
        let receipt = try VivoCellAxisReceipt(schemaVersion: 1, encoding: encoding, header: VivoCanonicalJSON.fingerprint(headerBytes),
            rows: VivoOmicsFileSnapshot.fingerprint(staging.appendingPathComponent("rows.bin"), maximumBytes: header.cellCount * Self.recordBytes),
            strings: VivoOmicsFileSnapshot.fingerprint(staging.appendingPathComponent("strings.bin"), maximumBytes: offset), stringBytes: offset,
            inputJSONL: VivoFingerprint(bytes: Array(hash.finalize())), implementation: implementation)
        try VivoCanonicalJSON.encode(receipt).write(to: staging.appendingPathComponent("receipt.json"), options: .withoutOverwriting)
        _ = try open(staging, implementation: implementation)
        try Task.checkCancellation(); try FileManager.default.moveItem(at: staging, to: destination)
        return receipt
    }
}
