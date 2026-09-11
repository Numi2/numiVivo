import Foundation
#if canImport(Darwin)
import Darwin
#else
import Glibc
#endif

private final class VivoH5ADNames {
    var names: [String] = []
    var error: Error?
}

/// Raw HDF5 transfers preserve stored datatypes (including float payloads,
/// integer widths, enums and strings). Numeric matrices never pass through Double.
extension VivoHDF5 {
    func projectionStorageLimit(_ object: ID,reserving bytes: Int = 0) throws {
        let get: @convention(c) (ID) -> ID = try symbol("H5Iget_file_id")
        let file=try id(get(object),"projection storage file");defer { close(file,"H5Fclose") }
        let size: @convention(c) (ID,UnsafeMutablePointer<UInt64>) -> Int32 = try symbol("H5Fget_filesize")
        var current: UInt64=0;try check(size(file,&current),"projection storage size")
        guard bytes>=0,current<=UInt64(projectionMaximumOutputBytes),UInt64(bytes)<=UInt64(projectionMaximumOutputBytes)-current else { throw VivoOmicsError.limit("projected H5AD exceeds output byte allowance") }
    }
    func projectionNames(_ object: ID, attributes: Bool = false) throws -> [String] {
        let state=VivoH5ADNames()
        typealias Callback = @convention(c) (ID, UnsafePointer<CChar>?, UnsafeRawPointer?, UnsafeMutableRawPointer?) -> Int32
        let callback: Callback = { _, name, _, pointer in
            guard let name,let pointer else { return -1 }
            let state=Unmanaged<VivoH5ADNames>.fromOpaque(pointer).takeUnretainedValue()
            do {
                let text=String(cString: name);try vivoH5ADComponent(text)
                guard state.names.count<100_000 else { throw VivoOmicsError.limit("projection metadata entries") }
                state.names.append(text);return 0
            } catch { state.error=error;return 1 }
        }
        let iterate: @convention(c) (ID,Int32,Int32,UnsafeMutablePointer<UInt64>,Callback,UnsafeMutableRawPointer?) -> Int32 = try symbol(attributes ? "H5Aiterate2" : "H5Literate2")
        var index: UInt64=0
        try check(iterate(object,0,0,&index,callback,Unmanaged.passUnretained(state).toOpaque()),"enumerate projection metadata")
        if let error=state.error { throw error }
        return state.names.sorted()
    }
    func projectionShape(_ value: ID, attribute: Bool = false) throws -> [UInt64] {
        let get: @convention(c) (ID) -> ID = try symbol(attribute ? "H5Aget_space" : "H5Dget_space")
        let space=try id(get(value),"projection dataspace");defer { close(space,"H5Sclose") }
        let spaceKind: @convention(c) (ID) -> Int32 = try symbol("H5Sget_simple_extent_type")
        guard spaceKind(space) != 2 else { throw VivoOmicsError.invalid("null HDF5 dataspace cannot be projected") }
        let rank: @convention(c) (ID) -> Int32 = try symbol("H5Sget_simple_extent_ndims")
        let n=rank(space);guard (0...8).contains(n) else { throw VivoOmicsError.invalid("projection supports ranks zero through eight") }
        var shape=[UInt64](repeating: 0,count: Int(n))
        let dims: @convention(c) (ID,UnsafeMutablePointer<UInt64>?,UnsafeMutablePointer<UInt64>?) -> Int32 = try symbol("H5Sget_simple_extent_dims")
        try check(dims(space,&shape,nil),"projection dimensions");return shape
    }
    func projectionProperties(_ kind: String) throws -> ID {
        guard let address=dlsym(library,"H5P_CLS_"+kind+"_CREATE_ID_g") else { throw VivoOmicsError.invalid("HDF5 projection property class") }
        let create: @convention(c) (ID) -> ID = try symbol("H5Pcreate")
        let p=try id(create(address.load(as: ID.self)),"projection creation properties")
        do {
            let times: @convention(c) (ID,UInt32) -> Int32 = try symbol("H5Pset_obj_track_times")
            try check(times(p,0),"disable projection timestamps");return p
        } catch { close(p,"H5Pclose");throw error }
    }
    func projectionFile(_ path: String) throws -> ID {
        let p=try projectionProperties("FILE");defer { close(p,"H5Pclose") }
        let create: @convention(c) (UnsafePointer<CChar>,UInt32,ID,ID) -> ID = try symbol("H5Fcreate")
        return try id(create(path,4,p,0),"create projection file")
    }
    func projectionGroup(_ parent: ID,_ name: String) throws -> ID {
        let p=try projectionProperties("GROUP");defer { close(p,"H5Pclose") }
        let create: @convention(c) (ID,UnsafePointer<CChar>,ID,ID,ID) -> ID = try symbol("H5Gcreate2")
        return try id(create(parent,name,0,p,0),"create projection group")
    }
    func projectionDataset(_ parent: ID,_ name: String,type: ID,shape: [UInt64]) throws -> ID {
        let p=try projectionProperties("DATASET");defer { close(p,"H5Pclose") }
        let s=try space(shape);defer { close(s,"H5Sclose") }
        let create: @convention(c) (ID,UnsafePointer<CChar>,ID,ID,ID,ID,ID) -> ID = try symbol("H5Dcreate2")
        return try id(create(parent,name,type,s,0,p,0),"create projected dataset")
    }
    func projectionCopy(_ source: ID,_ sourceName: String,_ output: ID,_ name: String) throws {
        let copy: @convention(c) (ID,UnsafePointer<CChar>,ID,UnsafePointer<CChar>,ID,ID) -> Int32 = try symbol("H5Ocopy")
        try check(copy(source,sourceName,output,name,0,0),"copy unaligned H5AD object")
        try projectionStorageLimit(output)
    }
    func projectionAttributes(_ source: ID,_ output: ID,excluding: Set<String> = []) throws {
        let size: @convention(c) (ID) -> Int = try symbol("H5Tget_size")
        let reclaim: @convention(c) (ID,ID,ID,UnsafeMutableRawPointer?) -> Int32 = try symbol("H5Dvlen_reclaim")
        for name in try projectionNames(source,attributes: true) where !excluding.contains(name) {
            let a=try attribute(source,name);defer { close(a,"H5Aclose") }
            let t=try type(a,attribute: true);defer { close(t,"H5Tclose") }
            let dims=try projectionShape(a,attribute: true),n=try count(dims,maximum: 2_000_000),width=size(t)
            guard width>0,width<=65_536,n<=64*1_024*1_024/width else { throw VivoOmicsError.limit("projection attribute buffer") }
            let buffer=UnsafeMutableRawPointer.allocate(byteCount: max(1,n*width),alignment: 16)
            defer { buffer.deallocate() }
            let s=try space(dims);defer { close(s,"H5Sclose") }
            try read(a,type: t,attribute: true,into: buffer)
            defer { _=reclaim(t,s,0,buffer) }
            try write(output,name,type: t,dimensions: dims,attribute: true,buffer: buffer)
            try projectionStorageLimit(output)
        }
    }
    /// Gather selected source coordinates into one contiguous output interval.
    /// The datatype is identical in memory and on disk: no numerical conversion.
    func projectionTransfer(_ source: ID,_ output: ID,type: ID,coordinates: [UInt64],count n: Int,outputStart: [UInt64],outputCount: [UInt64],reserveVariableBytes: ((Int) throws -> Void)? = nil) throws {
        guard n>0 else { return }
        guard n<=65_536 else { throw VivoOmicsError.limit("projection transfer entries") }
        let get: @convention(c) (ID) -> ID = try symbol("H5Dget_space")
        let inputSpace=try id(get(source),"projection input space");defer { close(inputSpace,"H5Sclose") }
        let outputSpace=try id(get(output),"projection output space");defer { close(outputSpace,"H5Sclose") }
        let memory=try space([UInt64(n)]);defer { close(memory,"H5Sclose") }
        let points: @convention(c) (ID,Int32,Int,UnsafePointer<UInt64>) -> Int32 = try symbol("H5Sselect_elements")
        try check(points(inputSpace,0,n,coordinates),"select projection input points")
        let slab: @convention(c) (ID,Int32,UnsafePointer<UInt64>,UnsafePointer<UInt64>?,UnsafePointer<UInt64>,UnsafePointer<UInt64>?) -> Int32 = try symbol("H5Sselect_hyperslab")
        try check(slab(outputSpace,0,outputStart,nil,outputCount,nil),"select projection output interval")
        let widthFn: @convention(c) (ID) -> Int = try symbol("H5Tget_size")
        let width=widthFn(type)
        guard width>0,width<=65_536,n<=64*1_024*1_024/width else { throw VivoOmicsError.limit("projection data buffer") }
        let detect: @convention(c) (ID,Int32) -> Int32 = try symbol("H5Tdetect_class")
        if detect(type,3)>0 || detect(type,9)>0 {
            let needed: @convention(c) (ID,ID,ID,UnsafeMutablePointer<UInt64>) -> Int32 = try symbol("H5Dvlen_get_buf_size")
            var bytes: UInt64=0;try check(needed(source,type,inputSpace,&bytes),"projection variable storage size")
            guard bytes<=64*1_024*1_024 else { throw VivoOmicsError.limit("projection variable data buffer") }
            try reserveVariableBytes?(Int(bytes))
        }
        let buffer=UnsafeMutableRawPointer.allocate(byteCount: n*width,alignment: 16);defer { buffer.deallocate() }
        let read: @convention(c) (ID,ID,ID,ID,ID,UnsafeMutableRawPointer?) -> Int32 = try symbol("H5Dread")
        try check(read(source,type,memory,inputSpace,0,buffer),"read projection block")
        let reclaim: @convention(c) (ID,ID,ID,UnsafeMutableRawPointer?) -> Int32 = try symbol("H5Dvlen_reclaim")
        defer { _=reclaim(type,memory,0,buffer) }
        let write: @convention(c) (ID,ID,ID,ID,ID,UnsafeRawPointer?) -> Int32 = try symbol("H5Dwrite")
        try check(write(output,type,memory,outputSpace,0,buffer),"write projection block")
        try projectionStorageLimit(output)
    }
    func projectionIntegers(_ output: ID,values: [Int64],offset: Int = 0) throws {
        guard !values.isEmpty else { return }
        let get: @convention(c) (ID) -> ID = try symbol("H5Dget_space")
        let s=try id(get(output),"index output space");defer { close(s,"H5Sclose") }
        let m=try space([UInt64(values.count)]);defer { close(m,"H5Sclose") }
        let select: @convention(c) (ID,Int32,UnsafePointer<UInt64>,UnsafePointer<UInt64>?,UnsafePointer<UInt64>,UnsafePointer<UInt64>?) -> Int32 = try symbol("H5Sselect_hyperslab")
        try check(select(s,0,[UInt64(offset)],nil,[UInt64(values.count)],nil),"select projected indices")
        let write: @convention(c) (ID,ID,ID,ID,ID,UnsafeRawPointer?) -> Int32 = try symbol("H5Dwrite")
        try values.withUnsafeBytes { try check(write(output,native("NATIVE_LLONG"),m,s,0,$0.baseAddress),"write projected indices") }
    }
}
