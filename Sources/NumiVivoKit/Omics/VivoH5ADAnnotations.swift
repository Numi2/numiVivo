import Foundation

public struct VivoH5ADAnnotationEdit: Codable, Sendable, Equatable {
    public enum Mode: String, Codable, Sendable { case add, replace }
    public let path: String
    public let mode: Mode
    public let value: VivoH5ADElement
    public init(path: String, mode: Mode = .add, value: VivoH5ADElement) { self.path = path; self.mode = mode; self.value = value }
    private enum CodingKeys: String, CodingKey { case path, mode, value }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["path", "mode", "value"])
        let c = try decoder.container(keyedBy: CodingKeys.self)
        path = try c.decode(String.self, forKey: .path)
        mode = try c.decode(Mode.self, forKey: .mode)
        value = try c.decode(VivoH5ADElement.self, forKey: .value)
    }
}

/// Describes derived fields in the unchanged source cell/feature order. Source
/// SHA-256 prevents applying results to a different ordering or assay snapshot.
public struct VivoH5ADAnnotationPlan: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let source: VivoFingerprint
    public let provenance: String
    public let edits: [VivoH5ADAnnotationEdit]
    public init(source: VivoFingerprint, provenance: String, edits: [VivoH5ADAnnotationEdit]) {
        schemaVersion = 1; self.source = source; self.provenance = provenance; self.edits = edits
    }
    private enum CodingKeys: String, CodingKey { case schemaVersion, source, provenance, edits }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["schemaVersion", "source", "provenance", "edits"])
        let c = try decoder.container(keyedBy: CodingKeys.self)
        schemaVersion = try c.decode(Int.self, forKey: .schemaVersion)
        source = try c.decode(VivoFingerprint.self, forKey: .source)
        provenance = try c.decode(String.self, forKey: .provenance)
        edits = try c.decode([VivoH5ADAnnotationEdit].self, forKey: .edits)
    }
    func validate() throws {
        guard schemaVersion == 1, !provenance.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty,
              provenance.utf8.count <= 16_384, (1...256).contains(edits.count), Set(edits.map(\.path)).count == edits.count else {
            throw VivoOmicsError.invalid("H5AD annotation schema, provenance or duplicate edits")
        }
        var remaining = 2_000_000, bytes = 64 * 1_024 * 1_024
        for edit in edits {
            let parts = edit.path.split(separator: "/", omittingEmptySubsequences: false).map(String.init)
            guard parts.count == 2, ["obs", "var", "obsm", "varm", "obsp", "varp", "layers", "uns"].contains(parts[0]),
                  edit.path != "uns/numivivo_edits" else { throw VivoOmicsError.invalid("unsupported annotation path") }
            try vivoH5ADComponent(parts[1]); try edit.value.validate(remaining: &remaining, stringBytes: &bytes)
        }
    }
}

public struct VivoH5ADAnnotationReceipt: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let source: VivoFingerprint
    public let plan: VivoFingerprint
    public let output: VivoFingerprint
    public let implementation: VivoFingerprint
    public let hdf5Version: String
}

extension VivoHDF5 {
    func exists(_ parent: ID, _ path: String) throws -> Bool {
        let fn: @convention(c) (ID, UnsafePointer<CChar>, ID) -> Int32 = try symbol("H5Lexists")
        let result = fn(parent, path, 0); try check(result, "query link"); return result > 0
    }
    func unlink(_ parent: ID, _ path: String) throws {
        let fn: @convention(c) (ID, UnsafePointer<CChar>, ID) -> Int32 = try symbol("H5Ldelete")
        try check(fn(parent, path, 0), "remove replaced link")
    }
    /// Copy before changing a group's attributes or children. Hard-link aliases
    /// elsewhere in the source continue to refer to the original unedited data.
    func detachGroup(_ file: ID, _ path: String) throws {
        let current = try object(file, path); defer { close(current, "H5Oclose") }
        let temporary = ".numivivo-edit-" + UUID().uuidString
        let detached = try group(file, temporary); defer { close(detached, "H5Gclose") }
        // Edits replace immediate child links; untouched datasets and nested
        // groups need no physical copy. Copy attributes onto a new group so
        // aliases to the original group retain their original metadata.
        try projectionAttributes(current, detached)
        let link: @convention(c) (ID, UnsafePointer<CChar>, ID, UnsafePointer<CChar>, ID, ID) -> Int32 = try symbol("H5Lcreate_hard")
        for child in try projectionNames(current) {
            try check(link(current, child, detached, child, 0, 0), "retain unedited annotation child")
        }
        try unlink(file, path)
        let move: @convention(c) (ID, UnsafePointer<CChar>, ID, UnsafePointer<CChar>, ID, ID) -> Int32 = try symbol("H5Lmove")
        try check(move(file, temporary, file, path, 0, 0), "install edited group")
    }
    func replaceStringAttribute(_ parent: ID, _ name: String, _ strings: [String]) throws {
        let delete: @convention(c) (ID, UnsafePointer<CChar>) -> Int32 = try symbol("H5Adelete")
        try check(delete(parent, name), "replace dataframe attribute")
        try writeStrings(parent, name, strings, attribute: true)
    }
}

extension VivoSingleCellH5AD {
    /// Publishes a NEW annotated AnnData file. X, raw and untouched fields retain
    /// their stored dtypes and values; neither a dense cell-gene projection nor
    /// a JSON translation of the original object is constructed.
    public static func annotate(_ sourceURL: URL, plan: VivoH5ADAnnotationPlan, implementation: VivoFingerprint,
                                to destination: URL) throws -> VivoH5ADAnnotationReceipt {
        try plan.validate()
        let outputLimit = 32 * 1_024 * 1_024 * 1_024
        let planBytes = try VivoCanonicalJSON.encode(plan), planID = try VivoCanonicalJSON.fingerprint(planBytes)
        let directory = FileManager.default.temporaryDirectory.appendingPathComponent("numivivo-h5ad-annotation-" + UUID().uuidString)
        try FileManager.default.createDirectory(at: directory, withIntermediateDirectories: false, attributes: [.posixPermissions: 0o700])
        defer { try? FileManager.default.removeItem(at: directory) }
        let staging = directory.appendingPathComponent("annotated.h5ad")
        // Hash the exact private copy HDF5 will edit, without retaining the
        // complete file in Data. Axis/value limits remain independently bounded.
        guard try VivoOmicsFileSnapshot.fingerprint(sourceURL, copyTo: staging, maximumBytes: 16 * 1_024 * 1_024 * 1_024) == plan.source else {
            throw VivoOmicsError.invalid("annotation source fingerprint mismatch")
        }
        let version = try VivoHDF5.lock.withLock {
            let h = try VivoHDF5()
            h.projectionMaximumOutputBytes = outputLimit
            let file = try h.file(staging.path, writable: true)
            defer { h.close(file, "H5Fclose") }
            guard try h.text(file, "encoding-type") == "anndata", try h.text(file, "encoding-version") == "0.1.0" else {
                throw VivoOmicsError.invalid("unsupported AnnData root encoding")
            }
            try h.validateAnnotationStorage(file)
            struct Frame { let count: Int; let index: String; var columns: [String] }
            func frame(_ path: String) throws -> Frame {
                let g = try h.object(file, path); defer { h.close(g, "H5Oclose") }
                guard try h.text(g, "encoding-type") == "dataframe", try h.text(g, "encoding-version") == "0.2.0" else {
                    throw VivoOmicsError.invalid("unsupported annotation dataframe")
                }
                let name = try h.text(g, "_index"); try vivoH5ADComponent(name)
                let index = try h.object(g, name); defer { h.close(index, "H5Oclose") }
                let encoding = try h.text(index, "encoding-type"), version = try h.text(index, "encoding-version")
                let dataPath: String
                if ["string-array", "array"].contains(encoding), version == "0.2.0" { dataPath = name }
                else if ["nullable-string-array", "nullable-integer", "nullable-boolean"].contains(encoding), version == "0.1.0" { dataPath = name + "/values" }
                else if encoding == "categorical", version == "0.2.0" { dataPath = name + "/codes" }
                else { throw VivoOmicsError.invalid("unsupported dataframe index encoding") }
                let d = try h.dataset(g, dataPath); defer { h.close(d, "H5Dclose") }
                let shape = try h.shape(d, attribute: false)
                guard shape.count == 1, shape[0] <= 2_000_000 else { throw VivoOmicsError.limit("annotation dataframe index shape") }
                if ["nullable-string-array", "nullable-integer", "nullable-boolean"].contains(encoding) {
                    let mask = try h.dataset(index, "mask"); defer { h.close(mask, "H5Dclose") }
                    guard try h.shape(mask, attribute: false) == shape, try h.mask(mask, maximum: 2_000_000).allSatisfy({ !$0 }) else {
                        throw VivoOmicsError.invalid("missing or inconsistent dataframe index")
                    }
                } else if encoding == "categorical" {
                    let categories = try h.dataset(index, "categories"); defer { h.close(categories, "H5Dclose") }
                    let labels = try h.strings(categories, maximum: 100_000)
                    guard try h.categoricalCodes(d, maximum: 2_000_000).allSatisfy({ $0 >= 0 && $0 < labels.count }) else {
                        throw VivoOmicsError.invalid("missing or out-of-range categorical index")
                    }
                }
                let a = try h.attribute(g, "column-order"); defer { h.close(a, "H5Aclose") }
                guard try h.shape(a, attribute: true).count == 1 else { throw VivoOmicsError.invalid("dataframe column-order must be a vector") }
                let columns = try h.strings(a, attribute: true, maximum: 10_000)
                guard Set(columns).count == columns.count, !columns.contains(name) else { throw VivoOmicsError.invalid("dataframe column-order is invalid") }
                return Frame(count: Int(shape[0]), index: name, columns: columns)
            }
            var frames = try ["obs": frame("obs"), "var": frame("var")]
            let cells = frames["obs"]!.count, features = frames["var"]!.count
            let parents = Set(plan.edits.map { String($0.path.split(separator: "/")[0]) }).union(["uns"])
            for parent in parents.sorted() {
                if try h.exists(file, parent) {
                    let g = try h.object(file, parent); defer { h.close(g, "H5Oclose") }
                    if parent != "obs" && parent != "var" {
                        guard try h.text(g, "encoding-type") == "dict", try h.text(g, "encoding-version") == "0.1.0" else {
                            throw VivoOmicsError.invalid("annotation parent is not an AnnData mapping")
                        }
                    }
                    try h.detachGroup(file, parent)
                } else {
                    let g = try h.group(file, parent); defer { h.close(g, "H5Gclose") }; try h.encoding(g, "dict", "0.1.0")
                }
            }
            for edit in plan.edits {
                try Task.checkCancellation()
                let parts = edit.path.split(separator: "/").map(String.init), parent = parts[0], name = parts[1]
                let shape = edit.value.shape
                if let frame = frames[parent] {
                    guard name != frame.index, name != "_index", shape == [frame.count] else { throw VivoOmicsError.invalid("annotation column changes index or has wrong length") }
                    guard frame.columns.contains(name) == (edit.mode == .replace) else { throw VivoOmicsError.invalid("column add/replace mode disagrees with dataframe") }
                } else if parent == "obsm" || parent == "varm" {
                    guard let shape, shape.count == 2, shape[0] == (parent == "obsm" ? cells : features),
                          edit.value.isSparse || shape[1] <= 256 else { throw VivoOmicsError.invalid("embedding shape or dense component allowance") }
                    switch edit.value { case .float64, .int64, .uint64, .boolean, .csrFloat64, .csrUInt64: break
                    default: throw VivoOmicsError.invalid("embedding must be numeric") }
                } else if parent == "obsp" || parent == "varp" || parent == "layers" {
                    let expected = parent == "layers" ? [cells, features] : Array(repeating: parent == "obsp" ? cells : features, count: 2)
                    guard edit.value.isSparse, shape == expected else { throw VivoOmicsError.invalid("graph/layer must be sparse with matching axes") }
                }
                let g = try h.object(file, parent); defer { h.close(g, "H5Oclose") }
                let exists = try h.exists(g, name)
                guard exists == (edit.mode == .replace) else { throw VivoOmicsError.invalid("annotation add/replace mode disagrees with existing field") }
                if exists { try h.unlink(g, name) }
                try h.writeElement(g, name, edit.value)
                if edit.mode == .add, frames[parent] != nil { frames[parent]!.columns.append(name) }
            }
            for parent in ["obs", "var"] where parents.contains(parent) {
                let g = try h.object(file, parent); defer { h.close(g, "H5Oclose") }
                try h.replaceStringAttribute(g, "column-order", frames[parent]!.columns)
            }
            // Each immutable edit record binds the complete plan and unchanged
            // source order. Repeated edits append history instead of erasing it.
            let journal = "uns/numivivo_edits"
            if try h.exists(file, journal) {
                let g = try h.object(file, journal); defer { h.close(g, "H5Oclose") }
                guard try h.text(g, "encoding-type") == "dict", try h.text(g, "encoding-version") == "0.1.0" else {
                    throw VivoOmicsError.invalid("reserved edit history is not a mapping")
                }
                try h.detachGroup(file, journal)
            } else {
                let g = try h.group(file, journal); defer { h.close(g, "H5Gclose") }; try h.encoding(g, "dict", "0.1.0")
            }
            let g = try h.object(file, journal); defer { h.close(g, "H5Oclose") }
            let runtime = try h.version()
            try h.writeElement(g, planID.hex, .dictionary(values: [
                "source_sha256": .string(shape: [], values: [plan.source.hex]),
                "plan_sha256": .string(shape: [], values: [planID.hex]),
                "implementation_sha256": .string(shape: [], values: [implementation.hex]),
                "hdf5_version": .string(shape: [], values: [runtime]),
                "plan_json": .string(shape: [], values: [String(decoding: planBytes, as: UTF8.self)])
            ]))
            let flush: @convention(c) (Int64, Int32) -> Int32 = try h.symbol("H5Fflush")
            try h.check(flush(file, 1), "flush annotation output")
            return runtime
        }
        try Task.checkCancellation()
        let output = try VivoOmicsFileSnapshot.fingerprint(staging, maximumBytes: outputLimit)
        let receipt = try VivoH5ADAnnotationReceipt(schemaVersion: 1, source: plan.source, plan: planID,
            output: output, implementation: implementation, hdf5Version: version)
        guard destination.isFileURL else { throw VivoOmicsError.invalid("local H5AD output required") }
        let files = try VivoRootedFileStore(rootURL: destination.deletingLastPathComponent(), createIfNeeded: false)
        guard try files.writeFile(from: staging, relative: destination.lastPathComponent, maximumBytes: outputLimit, immutable: true) else {
            throw VivoOmicsError.invalid("H5AD output already exists")
        }
        return receipt
    }
}

/// Validate storage dependencies without reading matrix values. A new file cannot
/// preserve relative external paths, virtual sources or opaque object references
/// by copying bytes and replacing groups, so such inputs must fail explicitly.
private final class VivoH5ADStorageCheck {
    let h: VivoHDF5
    var error: Error?
    var visited = 0
    init(_ h: VivoHDF5) { self.h = h }
    func noReferences(_ type: Int64) throws {
        let detect: @convention(c) (Int64, Int32) -> Int32 = try h.symbol("H5Tdetect_class")
        guard detect(type, 7) == 0 else { throw VivoOmicsError.invalid("H5AD object references cannot be preserved by annotation edits") }
    }
    func attributes(_ object: Int64) throws {
        typealias Callback = @convention(c) (Int64, UnsafePointer<CChar>?, UnsafeRawPointer?, UnsafeMutableRawPointer?) -> Int32
        let iterate: @convention(c) (Int64, Int32, Int32, UnsafeMutablePointer<UInt64>, Callback, UnsafeMutableRawPointer?) -> Int32 = try h.symbol("H5Aiterate2")
        let callback: Callback = { object, name, _, pointer in
            guard let name, let pointer else { return -1 }
            let context = Unmanaged<VivoH5ADStorageCheck>.fromOpaque(pointer).takeUnretainedValue()
            do {
                context.visited += 1
                guard context.visited <= 100_000 else { throw VivoOmicsError.limit("H5AD storage metadata allowance") }
                let a = try context.h.attribute(object, String(cString: name)); defer { context.h.close(a, "H5Aclose") }
                let type = try context.h.type(a, attribute: true); defer { context.h.close(type, "H5Tclose") }
                try context.noReferences(type)
                return 0
            } catch { context.error = error; return 1 }
        }
        var index: UInt64 = 0
        try h.check(iterate(object, 0, 2, &index, callback, Unmanaged.passUnretained(self).toOpaque()), "validate attribute storage")
        if let error { throw error }
    }
}

extension VivoHDF5 {
    func validateAnnotationStorage(_ file: ID) throws {
        let context = VivoH5ADStorageCheck(self)
        try context.attributes(file)
        typealias Callback = @convention(c) (Int64, UnsafePointer<CChar>?, UnsafeRawPointer?, UnsafeMutableRawPointer?) -> Int32
        let visit: @convention(c) (ID, Int32, Int32, Callback, UnsafeMutableRawPointer?) -> Int32 = try symbol("H5Lvisit2")
        let callback: Callback = { group, name, info, pointer in
            guard let name, let info, let pointer else { return -1 }
            let context = Unmanaged<VivoH5ADStorageCheck>.fromOpaque(pointer).takeUnretainedValue()
            do {
                // H5L_info2_t begins with the public C enum H5L_type_t.
                guard info.load(as: Int32.self) == 0 else { throw VivoOmicsError.invalid("annotation input must use self-contained hard links") }
                context.visited += 1
                guard context.visited <= 100_000 else { throw VivoOmicsError.limit("H5AD storage metadata allowance") }
                let h = context.h, path = String(cString: name), object = try h.object(group, path)
                defer { h.close(object, "H5Oclose") }
                let kind: @convention(c) (ID) -> Int32 = try h.symbol("H5Iget_type")
                if kind(object) == 5 {
                    // The dataset opener checks external storage and VDS layout.
                    let d = try h.dataset(group, path); defer { h.close(d, "H5Dclose") }
                    let type = try h.type(d, attribute: false); defer { h.close(type, "H5Tclose") }
                    try context.noReferences(type)
                } else if kind(object) != 2 { throw VivoOmicsError.invalid("unsupported HDF5 object in annotation input") }
                try context.attributes(object)
                return 0
            } catch { context.error = error; return 1 }
        }
        try check(visit(file, 0, 2, callback, Unmanaged.passUnretained(context).toOpaque()), "validate source storage")
        if let error = context.error { throw error }
    }
}
