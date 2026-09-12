import Foundation
import NumiVivoCore
#if canImport(Darwin)
import Darwin
#else
import Glibc
#endif

/// HDF5 is loaded only for H5AD operations. No Python interpreter, subprocess,
/// or mandatory HDF5 dependency is added to the molecular/iOS runtime.
/// All calls are serialized because a system HDF5 build may not be thread safe.
final class VivoHDF5 {
    /// Per-instance projection budget; each projection opens its own locked handle.
    var projectionMaximumOutputBytes = 1_073_741_824
    static let lock = NSLock()
    let library: UnsafeMutableRawPointer
    typealias ID = Int64
    init() throws {
        let candidates = [ProcessInfo.processInfo.environment["NUMIVIVO_HDF5_LIBRARY"],
            "/opt/homebrew/opt/hdf5/lib/libhdf5.dylib", "/usr/local/opt/hdf5/lib/libhdf5.dylib",
            "libhdf5.dylib", "libhdf5_serial.so", "libhdf5.so"].compactMap { $0 }
        guard let handle = candidates.lazy.compactMap({ dlopen($0, RTLD_NOW | RTLD_LOCAL) }).first else {
            throw VivoOmicsError.invalid("H5AD requires native HDF5; install hdf5 or set NUMIVIVO_HDF5_LIBRARY")
        }
        library = handle
        do {
            let open: @convention(c) () -> Int32 = try symbol("H5open")
            try check(open(), "initialize HDF5")
            try registerLZFDecoder()
        } catch { dlclose(handle); throw error }
    }
    deinit { dlclose(library) }
    /// Decode h5py's LZF chunks without loading Python or an external plugin.
    /// Existing runtime filters take precedence. New output remains uncompressed
    /// unless the caller explicitly configures an available encoder.
    private func registerLZFDecoder() throws {
        let available: @convention(c) (Int32) -> Int32 = try symbol("H5Zfilter_avail")
        if available(32000)>0 { return }
        let register: @convention(c) (UnsafeRawPointer) -> Int32 = try symbol("H5Zregister")
        var descriptor=NVivoHDF5FilterDescriptor()
        descriptor.version=1;descriptor.id=32000;descriptor.encoder_present=0;descriptor.decoder_present=1
        descriptor.filter = { flags,n,parameters,bytes,capacity,buffer in
            guard flags & 0x0100 != 0,let capacity,let buffer,let input=buffer.pointee,
                  bytes>0,bytes<=capacity.pointee else { return 0 }
            // Older h5py releases omit the uncompressed-size hint. Variable-
            // length storage can also exceed a supplied hint; grow only when
            // the decoder specifically reports insufficient output capacity.
            let maximum=64*1_024*1_024
            var size=(n>=3 && parameters != nil && parameters![2]>0) ? Int(parameters![2]) : capacity.pointee
            guard size>0,size<=maximum else { return 0 }
            while true {
                guard let output=malloc(size) else { return 0 }
                let count=VivoHDF5.decodeLZF(input: input.assumingMemoryBound(to: UInt8.self),count: bytes,
                    output: output.assumingMemoryBound(to: UInt8.self),capacity: size)
                if count>0 { free(input);buffer.pointee=output;capacity.pointee=size;return count }
                free(output)
                guard count == -1,size<maximum else { return 0 }
                size=min(maximum,size*2)
            }
        }
        try check(withUnsafePointer(to: &descriptor) { register(UnsafeRawPointer($0)) },"register native LZF decoder")
    }
    /// LZF literal/back-reference stream. Offsets are checked before pointer
    /// access; overlapping references copy forward, as required by the format.
    /// Zero means malformed input; -1 means insufficient output capacity.
    static func decodeLZF(input: UnsafePointer<UInt8>,count: Int,output: UnsafeMutablePointer<UInt8>,capacity: Int) -> Int {
        guard count>0,capacity>0 else { return 0 }
        var read=0,written=0
        while read<count {
            let control=Int(input[read]);read+=1
            if control<32 {
                let length=control+1
                guard length<=count-read else { return 0 }
                guard length<=capacity-written else { return -1 }
                for _ in 0..<length { output[written]=input[read];written+=1;read+=1 }
            } else {
                var length=control>>5
                if length==7 { guard read<count else { return 0 };length+=Int(input[read]);read+=1 }
                guard read<count else { return 0 }
                let distance=((control & 31)<<8)+Int(input[read])+1;read+=1;length+=2
                guard distance<=written else { return 0 }
                guard length<=capacity-written else { return -1 }
                for _ in 0..<length { output[written]=output[written-distance];written+=1 }
            }
        }
        return written
    }
    func symbol<T>(_ name: String) throws -> T {
        guard let address = dlsym(library, name) else { throw VivoOmicsError.invalid("HDF5 lacks \(name)") }
        return unsafeBitCast(address, to: T.self)
    }
    func check(_ status: Int32, _ operation: String) throws {
        guard status >= 0 else { throw VivoOmicsError.invalid("HDF5 \(operation) failed") }
    }
    func id(_ value: ID, _ operation: String) throws -> ID {
        guard value >= 0 else { throw VivoOmicsError.invalid("HDF5 \(operation) failed") }; return value
    }
    func native(_ name: String) throws -> ID {
        guard let address = dlsym(library, "H5T_" + name + "_g") else { throw VivoOmicsError.invalid("HDF5 missing datatype \(name)") }
        return address.load(as: ID.self)
    }
    func close(_ value: ID, _ name: String) {
        if let function: @convention(c) (ID) -> Int32 = try? symbol(name) { _ = function(value) }
    }
    func file(_ path: String, create: Bool = false, writable: Bool = false) throws -> ID {
        if create {
            let fn: @convention(c) (UnsafePointer<CChar>, UInt32, ID, ID) -> ID = try symbol("H5Fcreate")
            return try id(fn(path, 4, 0, 0), "create exclusive file")
        }
        let fn: @convention(c) (UnsafePointer<CChar>, UInt32, ID) -> ID = try symbol("H5Fopen")
        return try id(fn(path, writable ? 1 : 0, 0), "open file")
    }
    func object(_ file: ID, _ path: String) throws -> ID {
        let fn: @convention(c) (ID, UnsafePointer<CChar>, ID) -> ID = try symbol("H5Oopen")
        let value = try id(fn(file, path, 0), "open \(path)")
        do { try sameFile(file, value); return value } catch { close(value, "H5Oclose"); throw error }
    }
    func dataset(_ file: ID, _ path: String) throws -> ID {
        let fn: @convention(c) (ID, UnsafePointer<CChar>, ID) -> ID = try symbol("H5Dopen2")
        let value = try id(fn(file, path, 0), "open dataset \(path)")
        do {
            try sameFile(file, value)
            let get: @convention(c) (ID) -> ID = try symbol("H5Dget_create_plist")
            let properties = try id(get(value), "dataset properties"); defer { close(properties, "H5Pclose") }
            let external: @convention(c) (ID) -> Int32 = try symbol("H5Pget_external_count")
            let layout: @convention(c) (ID) -> Int32 = try symbol("H5Pget_layout")
            guard external(properties) == 0, (0...2).contains(layout(properties)) else {
                throw VivoOmicsError.invalid("external or virtual datasets are not self-contained H5AD inputs")
            }
            return value
        } catch { close(value, "H5Dclose"); throw error }
    }
    func sameFile(_ first: ID, _ second: ID) throws {
        let name: @convention(c) (ID, UnsafeMutablePointer<CChar>?, Int) -> Int = try symbol("H5Fget_name")
        func path(_ object: ID) throws -> String {
            let length = name(object, nil, 0)
            guard length >= 0, length <= 8192 else { throw VivoOmicsError.invalid("HDF5 file identity") }
            var buffer = [CChar](repeating: 0, count: length + 1)
            guard name(object, &buffer, buffer.count) == length else { throw VivoOmicsError.invalid("HDF5 file identity changed") }
            return String(decoding: buffer.prefix(length).map { UInt8(bitPattern: $0) }, as: UTF8.self)
        }
        guard try path(first) == path(second) else { throw VivoOmicsError.invalid("external HDF5 links cannot supply imported data") }
    }
    func attribute(_ object: ID, _ name: String) throws -> ID {
        let fn: @convention(c) (ID, UnsafePointer<CChar>, ID) -> ID = try symbol("H5Aopen")
        return try id(fn(object, name, 0), "open attribute \(name)")
    }
    func type(_ value: ID, attribute: Bool) throws -> ID {
        let fn: @convention(c) (ID) -> ID = try symbol(attribute ? "H5Aget_type" : "H5Dget_type")
        return try id(fn(value), "get datatype")
    }
    func shape(_ value: ID, attribute: Bool) throws -> [UInt64] {
        let get: @convention(c) (ID) -> ID = try symbol(attribute ? "H5Aget_space" : "H5Dget_space")
        let space = try id(get(value), "get dataspace"); defer { close(space, "H5Sclose") }
        let rank: @convention(c) (ID) -> Int32 = try symbol("H5Sget_simple_extent_ndims")
        let n = rank(space)
        guard n >= 0, n <= 2 else { throw VivoOmicsError.invalid("expected scalar, vector or matrix") }
        var dims = [UInt64](repeating: 0, count: Int(n))
        let extent: @convention(c) (ID, UnsafeMutablePointer<UInt64>?, UnsafeMutablePointer<UInt64>?) -> Int32 = try symbol("H5Sget_simple_extent_dims")
        try check(extent(space, &dims, nil), "get dimensions")
        return dims
    }
    func read(_ value: ID, type: ID, attribute: Bool, into buffer: UnsafeMutableRawPointer?, row: Int? = nil, range: Range<Int>? = nil) throws {
        if attribute {
            let fn: @convention(c) (ID, ID, UnsafeMutableRawPointer?) -> Int32 = try symbol("H5Aread")
            try check(fn(value, type, buffer), "read attribute")
        } else {
            let fn: @convention(c) (ID, ID, ID, ID, ID, UnsafeMutableRawPointer?) -> Int32 = try symbol("H5Dread")
            if let range {
                let dimensions = try shape(value,attribute: false)
                guard row == nil, dimensions.count == 1, range.lowerBound >= 0, UInt64(range.upperBound) <= dimensions[0] else {
                    throw VivoOmicsError.invalid("vector slice selection")
                }
                if range.isEmpty { return }
                let get: @convention(c) (ID) -> ID = try symbol("H5Dget_space")
                let fileSpace = try id(get(value),"vector dataspace"); defer { close(fileSpace,"H5Sclose") }
                let memorySpace = try space([UInt64(range.count)]); defer { close(memorySpace,"H5Sclose") }
                let select: @convention(c) (ID,Int32,UnsafePointer<UInt64>,UnsafePointer<UInt64>?,UnsafePointer<UInt64>,UnsafePointer<UInt64>?) -> Int32 = try symbol("H5Sselect_hyperslab")
                try check(select(fileSpace,0,[UInt64(range.lowerBound)],nil,[UInt64(range.count)],nil),"select vector slice")
                try check(fn(value,type,memorySpace,fileSpace,0,buffer),"read vector slice")
            } else if let row {
                let dimensions = try shape(value, attribute: false)
                guard dimensions.count == 2, row >= 0, UInt64(row) < dimensions[0] else { throw VivoOmicsError.invalid("dense row selection") }
                let get: @convention(c) (ID) -> ID = try symbol("H5Dget_space")
                let fileSpace = try id(get(value), "dense dataspace"); defer { close(fileSpace, "H5Sclose") }
                let memorySpace = try space([1, dimensions[1]]); defer { close(memorySpace, "H5Sclose") }
                let select: @convention(c) (ID, Int32, UnsafePointer<UInt64>, UnsafePointer<UInt64>?, UnsafePointer<UInt64>, UnsafePointer<UInt64>?) -> Int32 = try symbol("H5Sselect_hyperslab")
                try check(select(fileSpace, 0, [UInt64(row), 0], nil, [1, dimensions[1]], nil), "select dense row")
                try check(fn(value, type, memorySpace, fileSpace, 0, buffer), "read dense row")
            } else { try check(fn(value, type, 0, 0, 0, buffer), "read dataset") }
        }
    }
    func count(_ shape: [UInt64], maximum: Int) throws -> Int {
        var result = 1
        for dimension in shape {
            guard dimension <= UInt64(maximum), let d = Int(exactly: dimension), d == 0 || result <= maximum / d else {
                throw VivoOmicsError.limit("HDF5 array exceeds admission limit")
            }
            result *= Int(dimension)
        }
        guard result <= maximum else { throw VivoOmicsError.limit("HDF5 array exceeds admission limit") }
        return result
    }
    func strings(_ value: ID, attribute: Bool = false, maximum: Int) throws -> [String] {
        let n = try count(shape(value, attribute: attribute), maximum: maximum)
        // AnnData/h5py encode an empty dataframe column-order as an empty
        // numeric attribute. There are no values to coerce in this case.
        if n == 0 { return [] }
        let t = try type(value, attribute: attribute); defer { close(t, "H5Tclose") }
        let kind: @convention(c) (ID) -> Int32 = try symbol("H5Tget_class")
        guard kind(t) == 3 else { throw VivoOmicsError.invalid("expected HDF5 strings") }
        let variable: @convention(c) (ID) -> Int32 = try symbol("H5Tis_variable_str")
        if variable(t) > 0 {
            var pointers = [UnsafeMutablePointer<CChar>?](repeating: nil, count: n)
            let free: @convention(c) (UnsafeMutableRawPointer?) -> Int32 = try symbol("H5free_memory")
            defer { for pointer in pointers { if let pointer { _ = free(pointer) } } }
            try pointers.withUnsafeMutableBytes { try read(value, type: t, attribute: attribute, into: $0.baseAddress) }
            return try pointers.map {
                guard let p = $0, let string = String(validatingCString: p), string.utf8.count <= 16_384 else {
                    throw VivoOmicsError.invalid("invalid or oversized UTF-8 string")
                }; return string
            }
        }
        let size: @convention(c) (ID) -> Int = try symbol("H5Tget_size")
        let width = size(t)
        guard width > 0, width <= 16_384, n <= 64 * 1_024 * 1_024 / width else { throw VivoOmicsError.limit("HDF5 string storage") }
        var bytes = [UInt8](repeating: 0, count: n * width)
        try bytes.withUnsafeMutableBytes { try read(value, type: t, attribute: attribute, into: $0.baseAddress) }
        return try (0..<n).map { i in
            let slice = bytes[(i * width)..<((i + 1) * width)].prefix { $0 != 0 }
            guard let string = String(bytes: slice, encoding: .utf8) else { throw VivoOmicsError.invalid("invalid UTF-8") }
            return string
        }
    }
    func text(_ object: ID, _ name: String) throws -> String {
        let a = try attribute(object, name); defer { close(a, "H5Aclose") }
        guard try shape(a, attribute: true).isEmpty else { throw VivoOmicsError.invalid("expected scalar attribute \(name)") }
        return try strings(a, attribute: true, maximum: 1)[0]
    }
    func integers(_ value: ID, attribute: Bool = false, maximum: Int, allowFloat: Bool = false, row: Int? = nil, range: Range<Int>? = nil) throws -> [UInt64] {
        let dimensions = try shape(value, attribute: attribute)
        if row != nil && (attribute || dimensions.count != 2) { throw VivoOmicsError.invalid("dense row requires a matrix dataset") }
        if let range {
            guard !attribute, row == nil, dimensions.count == 1, range.lowerBound >= 0,
                  UInt64(range.upperBound) <= dimensions[0] else { throw VivoOmicsError.invalid("integer vector slice") }
        }
        let selected = range.map { [UInt64($0.count)] } ?? (row == nil ? dimensions : [dimensions[1]])
        let n = try count(selected,maximum: maximum)
        let t = try type(value, attribute: attribute); defer { close(t, "H5Tclose") }
        let kind: @convention(c) (ID) -> Int32 = try symbol("H5Tget_class")
        let size: @convention(c) (ID) -> Int = try symbol("H5Tget_size")
        if kind(t) == 1 && allowFloat {
            guard size(t) <= 8 else { throw VivoOmicsError.invalid("unsupported floating count precision") }
            var values = [Double](repeating: 0, count: n)
            try values.withUnsafeMutableBytes { try read(value, type: native("NATIVE_DOUBLE"), attribute: attribute, into: $0.baseAddress, row: row, range: range) }
            return try values.map {
                guard $0.isFinite, $0 >= 0, $0.rounded(.towardZero) == $0, $0 <= 9_007_199_254_740_992,
                      let count = UInt64(exactly: $0) else { throw VivoOmicsError.invalid("selected count matrix contains fractional, negative, nonfinite or inexact floating counts") }
                return count
            }
        }
        guard kind(t) == 0, size(t) <= 8 else { throw VivoOmicsError.invalid("expected integer array") }
        let sign: @convention(c) (ID) -> Int32 = try symbol("H5Tget_sign")
        if sign(t) == 1 {
            var values = [Int64](repeating: 0, count: n)
            try values.withUnsafeMutableBytes { try read(value, type: native("NATIVE_LLONG"), attribute: attribute, into: $0.baseAddress, row: row, range: range) }
            return try values.map { guard $0 >= 0 else { throw VivoOmicsError.invalid("negative integer") }; return UInt64($0) }
        }
        var values = [UInt64](repeating: 0, count: n)
        try values.withUnsafeMutableBytes { try read(value, type: native("NATIVE_ULLONG"), attribute: attribute, into: $0.baseAddress, row: row, range: range) }
        return values
    }
}

extension VivoHDF5 {
    func group(_ file: ID, _ path: String) throws -> ID {
        let fn: @convention(c) (ID, UnsafePointer<CChar>, ID, ID, ID) -> ID = try symbol("H5Gcreate2")
        return try id(fn(file, path, 0, 0, 0), "create group")
    }
    func space(_ dimensions: [UInt64]) throws -> ID {
        if dimensions.isEmpty {
            let fn: @convention(c) (Int32) -> ID = try symbol("H5Screate")
            return try id(fn(0), "create scalar space")
        }
        let fn: @convention(c) (Int32, UnsafePointer<UInt64>, UnsafePointer<UInt64>?) -> ID = try symbol("H5Screate_simple")
        return try id(fn(Int32(dimensions.count), dimensions, nil), "create space")
    }
    func write(_ object: ID, _ name: String, type: ID, dimensions: [UInt64], attribute: Bool, buffer: UnsafeRawPointer?) throws {
        let s = try space(dimensions); defer { close(s, "H5Sclose") }
        if attribute {
            let create: @convention(c) (ID, UnsafePointer<CChar>, ID, ID, ID, ID) -> ID = try symbol("H5Acreate2")
            let a = try id(create(object, name, type, s, 0, 0), "create attribute"); defer { close(a, "H5Aclose") }
            let fn: @convention(c) (ID, ID, UnsafeRawPointer?) -> Int32 = try symbol("H5Awrite")
            try check(fn(a, type, buffer), "write attribute")
        } else {
            let create: @convention(c) (ID, UnsafePointer<CChar>, ID, ID, ID, ID, ID) -> ID = try symbol("H5Dcreate2")
            let d = try id(create(object, name, type, s, 0, 0, 0), "create dataset"); defer { close(d, "H5Dclose") }
            let fn: @convention(c) (ID, ID, ID, ID, ID, UnsafeRawPointer?) -> Int32 = try symbol("H5Dwrite")
            try check(fn(d, type, 0, 0, 0, buffer), "write dataset")
        }
    }
    func writeStrings(_ object: ID, _ name: String, _ values: [String], attribute: Bool = false, scalar: Bool = false, dimensions: [UInt64]? = nil) throws {
        guard !scalar || values.count == 1, values.allSatisfy({ !$0.contains("\0") }) else { throw VivoOmicsError.invalid("invalid string write") }
        let copy: @convention(c) (ID) -> ID = try symbol("H5Tcopy")
        let t = try id(copy(native("C_S1")), "copy string type"); defer { close(t, "H5Tclose") }
        let size: @convention(c) (ID, Int) -> Int32 = try symbol("H5Tset_size")
        let charset: @convention(c) (ID, Int32) -> Int32 = try symbol("H5Tset_cset")
        try check(size(t, -1), "set variable string"); try check(charset(t, 1), "set UTF-8")
        let pointers = values.map { strdup($0) }
        defer { pointers.forEach { free($0) } }
        guard pointers.allSatisfy({ $0 != nil }) else { throw VivoOmicsError.limit("string allocation failed") }
        try pointers.withUnsafeBytes { try write(object, name, type: t, dimensions: dimensions ?? (scalar ? [] : [UInt64(values.count)]), attribute: attribute, buffer: $0.baseAddress) }
    }
    func encoding(_ object: ID, _ type: String, _ version: String) throws {
        try writeStrings(object, "encoding-type", [type], attribute: true, scalar: true)
        try writeStrings(object, "encoding-version", [version], attribute: true, scalar: true)
    }
    func writeIntegers(_ object: ID, _ name: String, _ values: [UInt64], attribute: Bool = false) throws {
        try values.withUnsafeBytes { try write(object, name, type: native("NATIVE_ULLONG"), dimensions: [UInt64(values.count)], attribute: attribute, buffer: $0.baseAddress) }
    }
}

extension VivoHDF5 {
    func categoricalCodes(_ value: ID, maximum: Int) throws -> [Int64] {
        let t = try type(value, attribute: false); defer { close(t, "H5Tclose") }
        let kind: @convention(c) (ID) -> Int32 = try symbol("H5Tget_class")
        let sign: @convention(c) (ID) -> Int32 = try symbol("H5Tget_sign")
        let size: @convention(c) (ID) -> Int = try symbol("H5Tget_size")
        guard kind(t) == 0, size(t) <= 8 else { throw VivoOmicsError.invalid("categorical codes must be integers") }
        if sign(t) == 0 {
            return try integers(value, maximum: maximum).map {
                guard let code = Int64(exactly: $0) else { throw VivoOmicsError.invalid("categorical code exceeds signed index range") }
                return code
            }
        }
        var values = [Int64](repeating: 0, count: try count(shape(value, attribute: false), maximum: maximum))
        try values.withUnsafeMutableBytes { try read(value, type: native("NATIVE_LLONG"), attribute: false, into: $0.baseAddress) }
        return values
    }
    func mask(_ value: ID, maximum: Int) throws -> [Bool] {
        let t = try type(value, attribute: false); defer { close(t, "H5Tclose") }
        let kind: @convention(c) (ID) -> Int32 = try symbol("H5Tget_class")
        let size: @convention(c) (ID) -> Int = try symbol("H5Tget_size")
        guard kind(t) == 8, size(t) == 1 else { throw VivoOmicsError.invalid("nullable mask must be boolean") }
        var values = [UInt8](repeating: 0, count: try count(shape(value, attribute: false), maximum: maximum))
        try values.withUnsafeMutableBytes { try read(value, type: t, attribute: false, into: $0.baseAddress) }
        guard values.allSatisfy({ $0 <= 1 }) else { throw VivoOmicsError.invalid("invalid boolean mask") }
        return values.map { $0 == 1 }
    }
}

extension VivoHDF5 {
    func writeIndices(_ object: ID, _ name: String, _ values: [Int]) throws {
        let signed = values.map(Int64.init)
        try signed.withUnsafeBytes { try write(object, name, type: native("NATIVE_LLONG"), dimensions: [UInt64(values.count)], attribute: false, buffer: $0.baseAddress) }
    }
}

extension VivoHDF5 {
    func version() throws -> String {
        let fn: @convention(c) (UnsafeMutablePointer<UInt32>, UnsafeMutablePointer<UInt32>, UnsafeMutablePointer<UInt32>) -> Int32 = try symbol("H5get_libversion")
        var major: UInt32 = 0, minor: UInt32 = 0, release: UInt32 = 0
        try check(fn(&major, &minor, &release), "read library version")
        return "\(major).\(minor).\(release)"
    }
}
