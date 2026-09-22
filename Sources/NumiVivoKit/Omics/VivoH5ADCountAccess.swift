import Foundation

/// AnnData's root metadata changed after the dataframe-0.1 writer used by
/// early `write_h5ad` releases.  The legacy dialect is intentionally narrow:
/// it is selected only after the root and both frames have been admitted by
/// `VivoSingleCellH5AD`, never by a per-column heuristic.
enum VivoH5ADDialect: Sendable, Equatable {
    case encodedAnnData010
    case legacyDataframe010
}

/// Shared AnnData frame decoding for native H5AD and H5MU paths.
struct VivoH5ADFrameReader {
    let h: VivoHDF5
    let file: Int64
    let dialect: VivoH5ADDialect
    init(h: VivoHDF5, file: Int64, dialect: VivoH5ADDialect = .encodedAnnData010) {
        self.h = h
        self.file = file
        self.dialect = dialect
    }
    func component(_ value: String) throws -> String {
        guard !value.isEmpty, !value.contains("/"), value != ".", value != "..", !value.contains("\0") else {
            throw VivoOmicsError.invalid("invalid AnnData component")
        }; return value
    }
    private func dataframeVersion(_ frame: String) throws -> String {
        let group = try h.object(file, frame); defer { h.close(group, "H5Oclose") }
        guard try h.text(group, "encoding-type") == "dataframe" else {
            throw VivoOmicsError.invalid("unsupported AnnData dataframe")
        }
        let version = try h.text(group, "encoding-version")
        switch dialect {
        case .encodedAnnData010:
            guard version == "0.2.0" else { throw VivoOmicsError.invalid("unsupported AnnData dataframe") }
        case .legacyDataframe010:
            guard version == "0.1.0" else { throw VivoOmicsError.invalid("mixed legacy AnnData dataframe") }
        }
        return version
    }
    func requireFrame(_ frame: String) throws { _ = try dataframeVersion(frame) }
    private func legacyStringColumn(_ path: String, maximum: Int) throws -> [String?] {
        let object = try h.dataset(file, path); defer { h.close(object, "H5Dclose") }
        let shape = try h.shape(object, attribute: false)
        guard shape.count == 1, shape[0] <= UInt64(maximum) else {
            throw VivoOmicsError.limit("legacy identity column shape")
        }
        let hasEncoding = try h.legacyHasAttribute(object, "encoding-type")
        let hasVersion = try h.legacyHasAttribute(object, "encoding-version")
        guard hasEncoding == hasVersion, !hasEncoding else {
            throw VivoOmicsError.invalid("legacy identity column encoding")
        }
        if try h.legacyHasAttribute(object, "categories") {
            let categories = try h.legacyObjectReference(object, "categories")
            defer { h.close(categories, "H5Oclose") }
            guard try h.legacyObjectKind(categories) == 5,
                  !(try h.legacyHasAttribute(categories, "encoding-type")),
                  !(try h.legacyHasAttribute(categories, "encoding-version")),
                  try h.shape(categories, attribute: false).count == 1 else {
                throw VivoOmicsError.invalid("legacy categorical categories")
            }
            let labels = try h.strings(categories, maximum: maximum)
            guard Set(labels).count == labels.count else {
                throw VivoOmicsError.invalid("legacy categorical labels must be unique")
            }
            let codes = try h.categoricalCodes(object, maximum: maximum)
            guard codes.count == Int(shape[0]) else { throw VivoOmicsError.invalid("legacy categorical code shape") }
            return try codes.map {
                if $0 == -1 { return nil }
                guard $0 >= 0, $0 < labels.count else { throw VivoOmicsError.invalid("legacy categorical code out of range") }
                return labels[Int($0)]
            }
        }
        return try h.strings(object, maximum: maximum)
    }
    func column(_ frame: String, _ name: String, maximum: Int) throws -> [String?] {
        _ = try dataframeVersion(frame)
        if dialect == .legacyDataframe010 { return try legacyStringColumn(frame + "/" + (try component(name)), maximum: maximum) }
        let path = frame + "/" + (try component(name)), object = try h.object(file, path)
        defer { h.close(object, "H5Oclose") }
        let encoding = try h.text(object, "encoding-type"), version = try h.text(object, "encoding-version")
        if encoding == "categorical", version == "0.2.0" {
            let categories = try h.dataset(file, path + "/categories"), codes = try h.dataset(file, path + "/codes")
            defer { h.close(categories, "H5Dclose"); h.close(codes, "H5Dclose") }
            guard try h.shape(categories, attribute: false).count == 1, try h.shape(codes, attribute: false).count == 1 else {
                throw VivoOmicsError.invalid("categorical arrays must be vectors")
            }
            let labels = try h.strings(categories, maximum: maximum)
            // Missing categorical values stay missing; required design identities reject them below.
            return try h.categoricalCodes(codes, maximum: maximum).map {
                if $0 == -1 { return nil }
                guard $0 >= 0, $0 < labels.count else { throw VivoOmicsError.invalid("categorical code out of range") }
                return labels[Int($0)]
            }
        }
        if encoding == "nullable-string-array", version == "0.1.0" {
            let values = try h.dataset(file, path + "/values"), mask = try h.dataset(file, path + "/mask")
            defer { h.close(values, "H5Dclose"); h.close(mask, "H5Dclose") }
            guard try h.shape(values, attribute: false).count == 1, try h.shape(mask, attribute: false).count == 1 else {
                throw VivoOmicsError.invalid("nullable arrays must be vectors")
            }
            let strings = try h.strings(values, maximum: maximum), missing = try h.mask(mask, maximum: maximum)
            guard strings.count == missing.count else { throw VivoOmicsError.invalid("nullable mask shape differs") }
            return strings.indices.map { missing[$0] ? nil : strings[$0] }
        }
        guard encoding == "string-array", version == "0.2.0" else { throw VivoOmicsError.invalid("selected identity column must be string-array or categorical") }
        let d = try h.dataset(file, path); defer { h.close(d, "H5Dclose") }
        guard try h.shape(d, attribute: false).count == 1 else { throw VivoOmicsError.invalid("identity column must be a vector") }
        return try h.strings(d, maximum: maximum)
    }
    func required(_ values: [String?]) throws -> [String] {
        try values.map { guard let value = $0 else { throw VivoOmicsError.invalid("missing selected design identity") }; return value }
    }
    /// Validate axis shape independently of identity values. Explicit mapping
    /// columns must not require coercion of an unused dataframe index.
    func indexLength(_ frame: String, maximum: Int) throws -> Int {
        let group = try h.object(file, frame); defer { h.close(group, "H5Oclose") }
        _ = try dataframeVersion(frame)
        let name = try component(h.text(group, "_index")), object = try h.object(group, name)
        defer { h.close(object, "H5Oclose") }
        if dialect == .legacyDataframe010 {
            guard !(try h.legacyHasAttribute(object, "encoding-type")),
                  !(try h.legacyHasAttribute(object, "encoding-version")),
                  !(try h.legacyHasAttribute(object, "categories")),
                  try h.legacyObjectKind(object) == 5 else {
                throw VivoOmicsError.invalid("legacy dataframe index encoding")
            }
            let type = try h.type(object, attribute: false); defer { h.close(type, "H5Tclose") }
            guard try h.legacyTypeKind(type) == 3 else { throw VivoOmicsError.invalid("legacy dataframe index must be strings") }
            let shape = try h.shape(object, attribute: false)
            guard shape.count == 1, shape[0] <= UInt64(maximum) else { throw VivoOmicsError.limit("dataframe index shape") }
            return Int(shape[0])
        }
        let kind = try h.text(object, "encoding-type"), version = try h.text(object, "encoding-version")
        let path: String
        if ["array", "string-array"].contains(kind), version == "0.2.0" { path = name }
        else if ["nullable-integer", "nullable-boolean", "nullable-string-array"].contains(kind), version == "0.1.0" { path = name + "/values" }
        else if kind == "categorical", version == "0.2.0" { path = name + "/codes" }
        else { throw VivoOmicsError.invalid("unsupported dataframe index encoding") }
        let data = try h.dataset(group, path); defer { h.close(data, "H5Dclose") }
        let shape = try h.shape(data, attribute: false)
        guard shape.count == 1, shape[0] <= maximum else { throw VivoOmicsError.limit("dataframe index shape") }
        if kind.hasPrefix("nullable-") {
            let mask = try h.dataset(object, "mask"); defer { h.close(mask, "H5Dclose") }
            guard try h.shape(mask, attribute: false) == shape else { throw VivoOmicsError.invalid("nullable index mask shape differs") }
        }
        return Int(shape[0])
    }
    func index(_ frame: String, maximum: Int) throws -> [String] {
        let group = try h.object(file, frame); defer { h.close(group, "H5Oclose") }
        _ = try dataframeVersion(frame)
        return try required(column(frame, h.text(group, "_index"), maximum: maximum))
    }
}

/// The caller owns HDF5 serialization and an immutable file snapshot.
enum VivoH5ADCountReader {
    static func scan(_ h: VivoHDF5, file: Int64, path: String, rows: Int, features: Int, limits: VivoOmicsLimits,
                     dialect: VivoH5ADDialect = .encodedAnnData010, onEntry: (Int, Int, UInt64) throws -> Void) throws {
        let matrix = try h.object(file, path); defer { h.close(matrix, "H5Oclose") }
        let hasEncoding = try h.legacyHasAttribute(matrix, "encoding-type")
        let hasVersion = try h.legacyHasAttribute(matrix, "encoding-version")
        guard hasEncoding == hasVersion else { throw VivoOmicsError.invalid("partial H5AD matrix encoding") }
        let encoding = hasEncoding ? try h.text(matrix, "encoding-type") : ""
        if encoding == "array" || (!hasEncoding && dialect == .legacyDataframe010) {
        if hasEncoding { guard try h.text(matrix, "encoding-version") == "0.2.0" else { throw VivoOmicsError.invalid("unsupported dense encoding") } }
        let d = try h.dataset(file, path); defer { h.close(d, "H5Dclose") }
        guard try h.shape(d, attribute: false) == [UInt64(rows), UInt64(features)] else {
            throw VivoOmicsError.invalid("dense matrix shape differs from obs/var")
        }
        var retained = 0
        for row in 0..<rows {
            try Task.checkCancellation()
            let values = try h.integers(d, maximum: limits.maximumFeatures, allowFloat: true, row: row)
            for (feature, count) in values.enumerated() where count > 0 {
                guard retained < limits.maximumNonzeros else { throw VivoOmicsError.limit("dense input exceeds sparse nonzero allowance") }
                try onEntry(row,feature,count); retained += 1
            }
        }
    } else {
        guard dialect == .encodedAnnData010,
              ["csr_matrix", "csc_matrix"].contains(encoding), try h.text(matrix, "encoding-version") == "0.1.0" else {
            throw VivoOmicsError.invalid("count import requires CSR, CSC or dense numeric array")
        }
        let shapeAttribute = try h.attribute(matrix, "shape"); defer { h.close(shapeAttribute, "H5Aclose") }
        guard try h.shape(shapeAttribute, attribute: true) == [2] else { throw VivoOmicsError.invalid("sparse shape must have two dimensions") }
        let shape = try h.integers(shapeAttribute, attribute: true, maximum: 2)
        guard shape == [UInt64(rows), UInt64(features)] else { throw VivoOmicsError.invalid("matrix shape differs from obs/var") }
        func integers(_ name: String, maximum: Int, counts: Bool = false) throws -> [UInt64] {
            let d = try h.dataset(file, path + "/" + name); defer { h.close(d, "H5Dclose") }
            guard try h.shape(d, attribute: false).count == 1 else { throw VivoOmicsError.invalid("sparse arrays must be vectors") }
            return try h.integers(d, maximum: maximum, allowFloat: counts)
        }
        let isCSR = encoding == "csr_matrix"
        let major = isCSR ? rows : features
        let minor = isCSR ? features : rows
        let offsets = try integers("indptr", maximum: major + 1)
        let indexDataset = try h.dataset(file,path+"/indices")
        defer { h.close(indexDataset,"H5Dclose") }
        let countDataset = try h.dataset(file,path+"/data")
        defer { h.close(countDataset,"H5Dclose") }
        let dimensions = try h.shape(countDataset,attribute: false)
        guard dimensions.count == 1, dimensions[0] <= UInt64(limits.maximumNonzeros),
              try h.shape(indexDataset,attribute: false) == dimensions,
              offsets.count == major+1, offsets.first == 0, offsets.last == dimensions[0] else {
            throw VivoOmicsError.invalid("malformed or oversized sparse array lengths")
        }
        for i in 0..<major {
            try Task.checkCancellation()
            guard offsets[i] <= offsets[i+1], offsets[i+1] <= dimensions[0] else { throw VivoOmicsError.invalid("malformed sparse offsets") }
            var canonical: [Int: UInt64] = [:]
            var cursor = Int(offsets[i])
            let end = Int(offsets[i+1])
            // Typical cell-major inputs already have strictly increasing feature
            // indices and fit in one bounded read. Emit those validated entries
            // directly; retain the canonicalization path for unsorted/duplicate
            // coordinates and major segments crossing the slice boundary.
            if cursor == end { continue }
            if end - cursor <= 65_536 {
                let indices = try h.integers(indexDataset,maximum: 65_536,range: cursor..<end)
                let counts = try h.integers(countDataset,maximum: 65_536,allowFloat: true,range: cursor..<end)
                var increasing = true
                for k in indices.indices {
                    guard indices[k] < minor else { throw VivoOmicsError.invalid("sparse index out of range") }
                    if k > 0, indices[k] <= indices[k-1] { increasing = false }
                }
                if increasing {
                    for k in indices.indices where counts[k] > 0 {
                        let index = Int(indices[k])
                        try onEntry(isCSR ? i : index,isCSR ? index : i,counts[k])
                    }
                    continue
                }
                for k in indices.indices where counts[k] > 0 {
                    let index = Int(indices[k])
                    canonical[index] = try vivoOmicsSum(canonical[index] ?? 0,counts[k])
                }
                cursor = end
            }
            while cursor < end {
                try Task.checkCancellation()
                let next = min(end,cursor+65_536)
                let indices = try h.integers(indexDataset,maximum: 65_536,range: cursor..<next)
                let counts = try h.integers(countDataset,maximum: 65_536,allowFloat: true,range: cursor..<next)
                for k in indices.indices {
                    guard indices[k] < minor else { throw VivoOmicsError.invalid("sparse index out of range") }
                    if counts[k] > 0 {
                        let index = Int(indices[k])
                        canonical[index] = try vivoOmicsSum(canonical[index] ?? 0,counts[k])
                    }
                }
                cursor = next
            }
            for j in canonical.keys.sorted() {
                try onEntry(isCSR ? i : j,isCSR ? j : i,canonical[j]!)
            }
        }
    }
    }
}
