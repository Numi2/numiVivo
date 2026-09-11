import Foundation
#if canImport(Darwin)
import Darwin
#else
import Glibc
#endif

extension VivoHDF5 {
    func legacyHasAttribute(_ object: ID,_ name: String) throws -> Bool {
        let fn: @convention(c) (ID,UnsafePointer<CChar>) -> Int32 = try symbol("H5Aexists")
        let value=fn(object,name);try check(value,"legacy attribute presence");return value>0
    }
    func legacyObjectKind(_ object: ID) throws -> Int32 {
        let fn: @convention(c) (ID) -> Int32 = try symbol("H5Iget_type");return fn(object)
    }
    func legacyTypeKind(_ type: ID) throws -> Int32 {
        let fn: @convention(c) (ID) -> Int32 = try symbol("H5Tget_class");return fn(type)
    }
    func legacyCompound(_ object: ID) throws -> Bool {
        guard try legacyObjectKind(object)==5 else { return false }
        let t=try type(object,attribute: false);defer { close(t,"H5Tclose") };return try legacyTypeKind(t)==6
    }
    func legacyMembers(_ dataset: ID) throws -> [String] {
        let t=try type(dataset,attribute: false);defer { close(t,"H5Tclose") }
        guard try legacyTypeKind(t)==6 else { throw VivoOmicsError.invalid("legacy aligned object must be compound") }
        let count: @convention(c) (ID) -> Int32 = try symbol("H5Tget_nmembers")
        let name: @convention(c) (ID,UInt32) -> UnsafeMutablePointer<CChar>? = try symbol("H5Tget_member_name")
        let free: @convention(c) (UnsafeMutableRawPointer?) -> Int32 = try symbol("H5free_memory")
        let n=count(t);guard (1...100_000).contains(n) else { throw VivoOmicsError.limit("legacy compound members") }
        return try (0..<n).map { i in
            guard let pointer=name(t,UInt32(i)) else { throw VivoOmicsError.invalid("legacy member name") }
            defer { _=free(pointer) };let value=String(cString: pointer);try vivoH5ADComponent(value);return value
        }
    }
    func legacyMemberType(_ dataset: ID,_ index: Int) throws -> ID {
        let t=try type(dataset,attribute: false);defer { close(t,"H5Tclose") }
        let fn: @convention(c) (ID,UInt32) -> ID = try symbol("H5Tget_member_type")
        return try id(fn(t,UInt32(index)),"legacy member datatype")
    }
    func legacyMemoryType(_ name: String,_ member: ID) throws -> ID {
        let size: @convention(c) (ID) -> Int = try symbol("H5Tget_size")
        let create: @convention(c) (Int32,Int) -> ID = try symbol("H5Tcreate")
        let insert: @convention(c) (ID,UnsafePointer<CChar>,Int,ID) -> Int32 = try symbol("H5Tinsert")
        let t=try id(create(6,size(member)),"legacy member memory type")
        do { try check(insert(t,name,0,member),"select legacy member");return t }
        catch { close(t,"H5Tclose");throw error }
    }
    func legacyArrayShape(_ type: ID) throws -> [UInt64] {
        guard try legacyTypeKind(type)==10 else { return [] }
        let rank: @convention(c) (ID) -> Int32 = try symbol("H5Tget_array_ndims")
        let dims: @convention(c) (ID,UnsafeMutablePointer<UInt64>) -> Int32 = try symbol("H5Tget_array_dims2")
        let n=rank(type);guard (1...7).contains(n) else { throw VivoOmicsError.limit("legacy member array rank") }
        var shape=[UInt64](repeating: 0,count: Int(n));try check(dims(type,&shape),"legacy member array dimensions");return shape
    }
    func legacyBaseType(_ type: ID) throws -> ID {
        let fn: @convention(c) (ID) -> ID = try symbol(try legacyTypeKind(type)==10 ? "H5Tget_super" : "H5Tcopy")
        return try id(fn(type),"legacy component datatype")
    }
    /// Read one compound member without converting its numeric representation,
    /// then write its scalar or fixed-array components to a modern dataset.
    func legacyTransfer(_ source: ID,_ output: ID,memoryType: ID,outputType: ID,rows: [UInt64],components: Int,outputStart: [UInt64],outputCount: [UInt64],scalarSource: Bool = false,reserveVariableBytes: (Int) throws -> Void) throws {
        guard !rows.isEmpty else { return }
        let n=rows.count;guard n<=65_536,components>0,n<=65_536/components else { throw VivoOmicsError.limit("legacy transfer values") }
        let get: @convention(c) (ID) -> ID = try symbol("H5Dget_space")
        let input=try id(get(source),"legacy input space");defer { close(input,"H5Sclose") }
        let out=try id(get(output),"legacy output space");defer { close(out,"H5Sclose") }
        if !scalarSource {
            let select: @convention(c) (ID,Int32,Int,UnsafePointer<UInt64>) -> Int32 = try symbol("H5Sselect_elements")
            try check(select(input,0,n,rows),"legacy selected rows")
        } else { guard n==1 else { throw VivoOmicsError.invalid("legacy scalar transfer") } }
        let select: @convention(c) (ID,Int32,UnsafePointer<UInt64>,UnsafePointer<UInt64>?,UnsafePointer<UInt64>,UnsafePointer<UInt64>?) -> Int32 = try symbol("H5Sselect_hyperslab")
        try check(select(out,0,outputStart,nil,outputCount,nil),"legacy output slab")
        let readSpace=try space([UInt64(n)]);defer { close(readSpace,"H5Sclose") }
        let writeSpace=try space([UInt64(n*components)]);defer { close(writeSpace,"H5Sclose") }
        let size: @convention(c) (ID) -> Int = try symbol("H5Tget_size")
        let width=size(memoryType);guard width>0,width<=65_536,n<=64*1_024*1_024/width else { throw VivoOmicsError.limit("legacy compound buffer") }
        let detect: @convention(c) (ID,Int32) -> Int32 = try symbol("H5Tdetect_class")
        let transfer=try scalarSource ? legacyScalarTransferProperties() : 0
        defer { if transfer != 0 { close(transfer,"H5Pclose") } }
        if scalarSource { try reserveVariableBytes(16_385) }
        if !scalarSource && (detect(memoryType,3)>0 || detect(memoryType,9)>0) {
            let needed: @convention(c) (ID,ID,ID,UnsafeMutablePointer<UInt64>) -> Int32 = try symbol("H5Dvlen_get_buf_size")
            var bytes: UInt64=0;try check(needed(source,memoryType,input,&bytes),"legacy variable buffer")
            guard bytes<=64*1_024*1_024 else { throw VivoOmicsError.limit("legacy variable buffer") };try reserveVariableBytes(Int(bytes))
        }
        let buffer=UnsafeMutableRawPointer.allocate(byteCount: n*width,alignment: 16);defer { buffer.deallocate() }
        let read: @convention(c) (ID,ID,ID,ID,ID,UnsafeMutableRawPointer?) -> Int32 = try symbol("H5Dread")
        try check(read(source,memoryType,readSpace,input,transfer,buffer),"read legacy member")
        let reclaim: @convention(c) (ID,ID,ID,UnsafeMutableRawPointer?) -> Int32 = try symbol("H5Dvlen_reclaim")
        defer { _=reclaim(memoryType,readSpace,transfer,buffer) }
        let write: @convention(c) (ID,ID,ID,ID,ID,UnsafeRawPointer?) -> Int32 = try symbol("H5Dwrite")
        try check(write(output,outputType,writeSpace,out,0,buffer),"write legacy member");try projectionStorageLimit(output)
    }
    /// HDF5 2.2.0's variable-buffer-size query crashes on scalar strings.
    /// A bounded transfer allocator admits those strings without that query.
    func legacyFalseMask(_ parent: ID,count: Int) throws {
        let create: @convention(c) (ID) -> ID = try symbol("H5Tenum_create")
        let t=try id(create(native("NATIVE_UCHAR")),"legacy mask type");defer { close(t,"H5Tclose") }
        let insert: @convention(c) (ID,UnsafePointer<CChar>,UnsafeRawPointer) -> Int32 = try symbol("H5Tenum_insert")
        var zero: UInt8=0,one: UInt8=1
        try check(insert(t,"FALSE",&zero),"legacy false mask");try check(insert(t,"TRUE",&one),"legacy true mask")
        let d=try projectionDataset(parent,"mask",type: t,shape: [UInt64(count)]);defer { close(d,"H5Dclose") }
        try encoding(d,"array","0.2.0")
        if count>0 {
            let values=[UInt8](repeating: 0,count: count)
            let write: @convention(c) (ID,ID,ID,ID,ID,UnsafeRawPointer?) -> Int32 = try symbol("H5Dwrite")
            try values.withUnsafeBytes { try check(write(d,t,0,0,0,$0.baseAddress),"legacy string mask") }
        }
        try projectionStorageLimit(d)
    }
    func legacyScalarTransferProperties() throws -> ID {
        guard let address=dlsym(library,"H5P_CLS_DATASET_XFER_ID_g") else { throw VivoOmicsError.invalid("legacy transfer property class") }
        let create: @convention(c) (ID) -> ID = try symbol("H5Pcreate")
        let p=try id(create(address.load(as: ID.self)),"legacy scalar transfer properties")
        typealias Allocate = @convention(c) (Int,UnsafeMutableRawPointer?) -> UnsafeMutableRawPointer?
        typealias Free = @convention(c) (UnsafeMutableRawPointer?,UnsafeMutableRawPointer?) -> Void
        let allocate: Allocate = { n,_ in n>0 && n<=16_385 ? malloc(n) : nil }
        let release: Free = { pointer,_ in free(pointer) }
        let set: @convention(c) (ID,Allocate,UnsafeMutableRawPointer?,Free,UnsafeMutableRawPointer?) -> Int32 = try symbol("H5Pset_vlen_mem_manager")
        do { try check(set(p,allocate,nil,release,nil),"bounded legacy string allocator");return p }
        catch { close(p,"H5Pclose");throw error }
    }
    func legacyScalarString(_ source: ID,type: ID) throws -> String {
        let variable: @convention(c) (ID) -> Int32 = try symbol("H5Tis_variable_str")
        if variable(type)<=0 { return try strings(source,maximum: 1)[0] }
        let p=try legacyScalarTransferProperties();defer { close(p,"H5Pclose") }
        var pointer: UnsafeMutablePointer<CChar>?=nil
        defer { free(pointer) }
        let read: @convention(c) (ID,ID,ID,ID,ID,UnsafeMutableRawPointer?) -> Int32 = try symbol("H5Dread")
        try check(read(source,type,0,0,p,&pointer),"bounded legacy scalar string")
        guard let pointer,let value=String(validatingCString: pointer),value.utf8.count<=16_384 else { throw VivoOmicsError.invalid("legacy scalar UTF-8 string") };return value
    }
    func legacyStringReadBound(_ source: ID,type: ID) throws {
        let get: @convention(c) (ID) -> ID = try symbol("H5Dget_space")
        let s=try id(get(source),"legacy category space");defer { close(s,"H5Sclose") }
        let needed: @convention(c) (ID,ID,ID,UnsafeMutablePointer<UInt64>) -> Int32 = try symbol("H5Dvlen_get_buf_size")
        var bytes: UInt64=0;try check(needed(source,type,s,&bytes),"legacy category string storage")
        guard bytes<=64*1_024*1_024 else { throw VivoOmicsError.limit("legacy category string buffer") }
    }
    func legacyCodes(_ source: ID,member: String,range: Range<Int>) throws -> [Int64] {
        let memory=try legacyMemoryType(member,native("NATIVE_LLONG"));defer { close(memory,"H5Tclose") }
        let get: @convention(c) (ID) -> ID = try symbol("H5Dget_space")
        let input=try id(get(source),"legacy codes space");defer { close(input,"H5Sclose") }
        let select: @convention(c) (ID,Int32,UnsafePointer<UInt64>,UnsafePointer<UInt64>?,UnsafePointer<UInt64>,UnsafePointer<UInt64>?) -> Int32 = try symbol("H5Sselect_hyperslab")
        try check(select(input,0,[UInt64(range.lowerBound)],nil,[UInt64(range.count)],nil),"legacy code rows")
        let output=try space([UInt64(range.count)]);defer { close(output,"H5Sclose") }
        var values=[Int64](repeating: 0,count: range.count)
        let read: @convention(c) (ID,ID,ID,ID,ID,UnsafeMutableRawPointer?) -> Int32 = try symbol("H5Dread")
        try values.withUnsafeMutableBytes { try check(read(source,memory,output,input,0,$0.baseAddress),"legacy code domain") };return values
    }
}
