import Foundation

/// Typed values for NEW or explicitly replaced AnnData fields. Original fields
/// remain in HDF5 without decoding or narrowing their precision/encoding.
public indirect enum VivoH5ADElement: Codable, Sendable, Equatable {
    case float64(shape: [Int], values: [Double])
    case int64(shape: [Int], values: [Int64])
    case uint64(shape: [Int], values: [UInt64])
    case boolean(shape: [Int], values: [Bool])
    case string(shape: [Int], values: [String])
    case categorical(codes: [Int64], categories: [String], ordered: Bool)
    /// Missing floating results use the conventional NumPy NaN representation.
    case nullableFloat64(values: [Double?])
    case nullableString(values: [String?])
    case nullableInt64(values: [Int64?])
    case nullableBoolean(values: [Bool?])
    case csrFloat64(shape: [Int], indptr: [Int], indices: [Int], values: [Double])
    case csrUInt64(shape: [Int], indptr: [Int], indices: [Int], values: [UInt64])
    case dictionary(values: [String: VivoH5ADElement])

    var shape: [Int]? {
        switch self {
        case .float64(let shape, _), .int64(let shape, _), .uint64(let shape, _), .boolean(let shape, _), .string(let shape, _),
             .csrFloat64(let shape, _, _, _), .csrUInt64(let shape, _, _, _): return shape
        case .categorical(let values, _, _): return [values.count]
        case .nullableFloat64(let values): return [values.count]
        case .nullableString(let values): return [values.count]
        case .nullableInt64(let values): return [values.count]
        case .nullableBoolean(let values): return [values.count]
        case .dictionary: return nil
        }
    }
    var isSparse: Bool {
        switch self { case .csrFloat64, .csrUInt64: true; default: false }
    }
    func validate(depth: Int = 0, remaining: inout Int, stringBytes: inout Int) throws {
        guard depth < 16, remaining > 0 else { throw VivoOmicsError.limit("H5AD element depth or element allowance") }
        remaining -= 1
        func array(_ shape: [Int], _ n: Int) throws {
            guard shape.count <= 8, n <= remaining else { throw VivoOmicsError.limit("H5AD array rank or element allowance") }
            var size = 1
            for dimension in shape {
                guard dimension >= 0, dimension <= 2_000_000, dimension == 0 || size <= 2_000_000 / dimension else {
                    throw VivoOmicsError.limit("H5AD array shape")
                }
                size *= dimension
            }
            guard size == n else { throw VivoOmicsError.invalid("H5AD shape and values disagree") }
            remaining -= n
        }
        func strings(_ values: [String]) throws {
            for value in values {
                guard !value.contains("\0"), value.utf8.count <= 16_384, value.utf8.count <= stringBytes else {
                    throw VivoOmicsError.limit("H5AD string allowance or embedded NUL")
                }
                stringBytes -= value.utf8.count
            }
        }
        func sparse(_ shape: [Int], _ indptr: [Int], _ indices: [Int], _ n: Int) throws {
            guard shape.count == 2, shape.allSatisfy({ (0...100_000).contains($0) }),
                  indptr.count == shape[0] + 1, indptr.first == 0, indptr.last == n, indices.count == n else {
                throw VivoOmicsError.invalid("H5AD CSR shape or array lengths")
            }
            guard n <= remaining / 2, indptr.count <= remaining - 2 * n else { throw VivoOmicsError.limit("H5AD CSR allowance") }
            remaining -= 2 * n + indptr.count
            for row in 0..<shape[0] {
                guard indptr[row] >= 0, indptr[row] <= indptr[row + 1], indptr[row + 1] <= n else { throw VivoOmicsError.invalid("H5AD CSR offsets") }
                var previous = -1
                for k in indptr[row]..<indptr[row + 1] {
                    guard indices[k] > previous, indices[k] < shape[1] else { throw VivoOmicsError.invalid("H5AD CSR indices must be sorted, unique and in range") }
                    previous = indices[k]
                }
            }
        }
        switch self {
        case .float64(let shape, let values):
            try array(shape, values.count)
            guard values.allSatisfy(\.isFinite) else { throw VivoOmicsError.invalid("nonfinite annotation value") }
        case .int64(let shape, let values): try array(shape, values.count)
        case .uint64(let shape, let values): try array(shape, values.count)
        case .boolean(let shape, let values): try array(shape, values.count)
        case .string(let shape, let values): try array(shape, values.count); try strings(values)
        case .categorical(let codes, let categories, _):
            try array([codes.count], codes.count); try array([categories.count], categories.count); try strings(categories)
            guard Set(categories).count == categories.count, codes.allSatisfy({ $0 >= -1 && $0 < categories.count }) else {
                throw VivoOmicsError.invalid("H5AD categorical code or duplicate category")
            }
        case .nullableFloat64(let values):
            try array([values.count], values.count)
            guard values.compactMap({ $0 }).allSatisfy(\.isFinite) else { throw VivoOmicsError.invalid("nonfinite supplied annotation value") }
        case .nullableString(let values): try array([values.count], values.count); try strings(values.compactMap { $0 })
        case .nullableInt64(let values): try array([values.count], values.count)
        case .nullableBoolean(let values): try array([values.count], values.count)
        case .csrFloat64(let shape, let indptr, let indices, let values):
            try sparse(shape, indptr, indices, values.count)
            guard values.allSatisfy(\.isFinite) else { throw VivoOmicsError.invalid("nonfinite sparse annotation value") }
        case .csrUInt64(let shape, let indptr, let indices, let values): try sparse(shape, indptr, indices, values.count)
        case .dictionary(let values):
            for key in values.keys.sorted() { try vivoH5ADComponent(key); try strings([key]); try values[key]!.validate(depth: depth + 1, remaining: &remaining, stringBytes: &stringBytes) }
        }
    }
}

func vivoH5ADComponent(_ value: String) throws {
    guard !value.isEmpty, value.utf8.count <= 1_024, !value.contains("/"), !value.contains("\0"), value != ".", value != ".." else {
        throw VivoOmicsError.invalid("invalid H5AD component")
    }
}

extension VivoHDF5 {
    func writeBooleans(_ object: ID, _ name: String, shape: [Int], values: [Bool], attribute: Bool = false) throws {
        let create: @convention(c) (ID) -> ID = try symbol("H5Tenum_create")
        let t = try id(create(native("NATIVE_UCHAR")), "create boolean type"); defer { close(t, "H5Tclose") }
        let insert: @convention(c) (ID, UnsafePointer<CChar>, UnsafeRawPointer) -> Int32 = try symbol("H5Tenum_insert")
        var zero: UInt8 = 0, one: UInt8 = 1
        try check(insert(t, "FALSE", &zero), "false enum"); try check(insert(t, "TRUE", &one), "true enum")
        let bytes: [UInt8] = values.map { $0 ? 1 : 0 }
        try bytes.withUnsafeBytes { try write(object, name, type: t, dimensions: shape.map(UInt64.init), attribute: attribute, buffer: $0.baseAddress) }
    }
    func writeElement(_ parent: ID, _ name: String, _ element: VivoH5ADElement) throws {
        func array<T>(_ values: [T], _ shape: [Int], _ type: String) throws {
            try values.withUnsafeBytes { try write(parent, name, type: native(type), dimensions: shape.map(UInt64.init), attribute: false, buffer: $0.baseAddress) }
        }
        var kind = "array", version = "0.2.0"
        switch element {
        case .float64(let shape, let values): try array(values, shape, "NATIVE_DOUBLE"); kind = shape.isEmpty ? "numeric-scalar" : "array"
        case .int64(let shape, let values): try array(values, shape, "NATIVE_LLONG"); kind = shape.isEmpty ? "numeric-scalar" : "array"
        case .uint64(let shape, let values): try array(values, shape, "NATIVE_ULLONG"); kind = shape.isEmpty ? "numeric-scalar" : "array"
        case .boolean(let shape, let values): try writeBooleans(parent, name, shape: shape, values: values); kind = shape.isEmpty ? "numeric-scalar" : "array"
        case .string(let shape, let values):
            try writeStrings(parent, name, values, scalar: shape.isEmpty, dimensions: shape.map(UInt64.init))
            kind = shape.isEmpty ? "string" : "string-array"
        case .categorical(let codes, let categories, let ordered):
            let g = try group(parent, name); defer { close(g, "H5Gclose") }
            try writeElement(g, "codes", .int64(shape: [codes.count], values: codes))
            try writeElement(g, "categories", .string(shape: [categories.count], values: categories))
            try writeBooleans(g, "ordered", shape: [], values: [ordered], attribute: true)
            kind = "categorical"
        case .nullableFloat64(let values):
            try array(values.map { $0 ?? .nan }, [values.count], "NATIVE_DOUBLE")
        case .nullableString(let values):
            let g = try group(parent, name); defer { close(g, "H5Gclose") }
            try writeElement(g, "values", .string(shape: [values.count], values: values.map { $0 ?? "" }))
            try writeElement(g, "mask", .boolean(shape: [values.count], values: values.map { $0 == nil }))
            kind = "nullable-string-array"; version = "0.1.0"
        case .nullableInt64(let values):
            let g = try group(parent, name); defer { close(g, "H5Gclose") }
            try writeElement(g, "values", .int64(shape: [values.count], values: values.map { $0 ?? 0 }))
            try writeElement(g, "mask", .boolean(shape: [values.count], values: values.map { $0 == nil }))
            kind = "nullable-integer"; version = "0.1.0"
        case .nullableBoolean(let values):
            let g = try group(parent, name); defer { close(g, "H5Gclose") }
            try writeElement(g, "values", .boolean(shape: [values.count], values: values.map { $0 ?? false }))
            try writeElement(g, "mask", .boolean(shape: [values.count], values: values.map { $0 == nil }))
            kind = "nullable-boolean"; version = "0.1.0"
        case .csrFloat64(let shape, let indptr, let indices, let values):
            let g = try group(parent, name); defer { close(g, "H5Gclose") }
            try writeIntegers(g, "shape", shape.map(UInt64.init), attribute: true)
            try writeIndices(g, "indptr", indptr); try writeIndices(g, "indices", indices)
            try writeElement(g, "data", .float64(shape: [values.count], values: values))
            kind = "csr_matrix"; version = "0.1.0"
        case .csrUInt64(let shape, let indptr, let indices, let values):
            let g = try group(parent, name); defer { close(g, "H5Gclose") }
            try writeIntegers(g, "shape", shape.map(UInt64.init), attribute: true)
            try writeIndices(g, "indptr", indptr); try writeIndices(g, "indices", indices)
            try writeIntegers(g, "data", values)
            kind = "csr_matrix"; version = "0.1.0"
        case .dictionary(let values):
            let g = try group(parent, name); defer { close(g, "H5Gclose") }
            for key in values.keys.sorted() { try writeElement(g, key, values[key]!) }
            kind = "dict"; version = "0.1.0"
        }
        let object = try object(parent, name); defer { close(object, "H5Oclose") }
        try encoding(object, kind, version)
    }
}
