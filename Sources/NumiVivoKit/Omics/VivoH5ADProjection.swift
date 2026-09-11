import Foundation

public struct VivoH5ADProjectionPlan: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let source: VivoFingerprint
    public let provenance: String
    /// nil keeps the complete axis; [] selects an empty axis. Indices address
    /// the fingerprinted source order. Repeated indices repeat aligned values.
    public let observationIndices: [Int]?
    public let featureIndices: [Int]?
    public let maximumElementVisits: Int
    /// nil preserves the historical 1 GiB output limit and encoded plan.
    public let maximumOutputBytes: Int?
    public init(source: VivoFingerprint,provenance: String,observationIndices: [Int]? = nil,featureIndices: [Int]? = nil,maximumElementVisits: Int = 500_000_000, maximumOutputBytes: Int? = nil) {
        self.maximumOutputBytes=maximumOutputBytes
        schemaVersion=1;self.source=source;self.provenance=provenance;self.observationIndices=observationIndices
        self.featureIndices=featureIndices;self.maximumElementVisits=maximumElementVisits
    }
    private enum CodingKeys: String,CodingKey { case schemaVersion,source,provenance,observationIndices,featureIndices,maximumElementVisits,maximumOutputBytes }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder,allowed: ["schemaVersion","source","provenance","observationIndices","featureIndices","maximumElementVisits","maximumOutputBytes"])
        let c=try decoder.container(keyedBy: CodingKeys.self)
        schemaVersion=try c.decode(Int.self,forKey: .schemaVersion);source=try c.decode(VivoFingerprint.self,forKey: .source)
        provenance=try c.decode(String.self,forKey: .provenance)
        observationIndices=try c.decodeIfPresent([Int].self,forKey: .observationIndices)
        featureIndices=try c.decodeIfPresent([Int].self,forKey: .featureIndices)
        maximumOutputBytes=try c.decodeIfPresent(Int.self,forKey: .maximumOutputBytes)
        maximumElementVisits=try c.decodeIfPresent(Int.self,forKey: .maximumElementVisits) ?? 500_000_000
    }
    func validate() throws {
        guard schemaVersion==1,!provenance.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty,provenance.utf8.count<=16_384,
              (1...2_000_000_000).contains(maximumElementVisits) else { throw VivoOmicsError.invalid("projection schema, provenance or work allowance") }
        if let maximumOutputBytes {
            guard (1...8_589_934_592).contains(maximumOutputBytes) else { throw VivoOmicsError.invalid("projection output byte allowance") }
        }
        for indices in [observationIndices,featureIndices].compactMap({ $0 }) {
            guard indices.count<=2_000_000,indices.allSatisfy({ (0..<2_000_000).contains($0) }) else {
                throw VivoOmicsError.invalid("projection indices must be nonnegative and bounded")
            }
        }
    }
}
public struct VivoH5ADProjectedField: Codable,Sendable,Equatable {
    public let path: String
    public let encoding: String
    public let sourceShape: [Int]
    public let outputShape: [Int]
}
public struct VivoH5ADProjectionReport: Codable,Sendable,Equatable {
    public let method: String
    public let sourceShape: [Int]
    public let outputShape: [Int]
    public let elementVisits: Int
    public let fields: [VivoH5ADProjectedField]
    public let hdf5Version: String
    public let qualification: String
}
public struct VivoH5ADProjectionReceipt: Codable,Sendable,Equatable {
    public let schemaVersion: Int
    public let source: VivoFingerprint
    public let plan: VivoFingerprint
    public let output: VivoFingerprint
    public let report: VivoFingerprint
    public let implementation: VivoFingerprint
}

private final class VivoH5ADProjector {
    let h: VivoHDF5
    let legacy: Bool
    let root: Int64
    var movedLegacyCategories=Set<String>()
    var remaining: Int
    var remainingStorage: Int
    var fields: [VivoH5ADProjectedField]=[]
    init(_ h: VivoHDF5,limit: Int,root: Int64,legacy: Bool) { self.h=h;remaining=limit;remainingStorage=h.projectionMaximumOutputBytes;self.root=root;self.legacy=legacy }
    func consume(_ count: Int) throws {
        try Task.checkCancellation()
        guard count>=0,count<=remaining else { throw VivoOmicsError.limit("H5AD projection element visits") };remaining-=count
    }
    func reserveStorage(_ bytes: Int,output: Int64) throws {
        guard bytes>=0,bytes<=remainingStorage else { throw VivoOmicsError.limit("projection output allocation exceeds byte allowance") }
        try h.projectionStorageLimit(output,reserving: bytes);remainingStorage-=bytes
    }
    func typeWidth(_ type: Int64) throws -> Int {
        let size: @convention(c) (Int64) -> Int = try h.symbol("H5Tget_size")
        let width=size(type);guard (1...65_536).contains(width) else { throw VivoOmicsError.limit("projection datatype width") };return width
    }
    func encoding(_ object: Int64,_ kind: String,_ version: String) throws {
        if legacy,try !h.legacyHasAttribute(object,"encoding-type"),try !h.legacyHasAttribute(object,"encoding-version") {
            if kind=="anndata" || kind=="dict" || kind=="raw" {
                guard try [1,2].contains(h.legacyObjectKind(object)) else { throw VivoOmicsError.invalid("legacy group encoding") };return
            }
            if kind=="csr_matrix" || kind=="csc_matrix" {
                guard try h.text(object,"h5sparse_format")==String(kind.prefix(3)) else { throw VivoOmicsError.invalid("legacy sparse format") };return
            }
            if kind=="array" || kind=="string-array" { return }
        }
        guard try h.text(object,"encoding-type")==kind,try h.text(object,"encoding-version")==version else { throw VivoOmicsError.invalid("unsupported projection encoding: "+kind) }
    }
    func frameLength(_ frame: Int64) throws -> Int {
        if legacy,try h.legacyCompound(frame) {
            let dims=try h.projectionShape(frame);guard dims.count==1,dims[0]<=2_000_000 else { throw VivoOmicsError.limit("legacy dataframe axis") };return Int(dims[0])
        }
        try encoding(frame,"dataframe","0.2.0")
        let name=try h.text(frame,"_index");try vivoH5ADComponent(name)
        let obj=try h.object(frame,name);defer { h.close(obj,"H5Oclose") }
        let kind=try h.text(obj,"encoding-type"),path: String
        switch kind {
        case "string-array":path=name
        case "nullable-string-array":path=name+"/values"
        case "categorical":path=name+"/codes"
        default:throw VivoOmicsError.invalid("unsupported projection dataframe index")
        }
        let d=try h.dataset(frame,path);defer { h.close(d,"H5Dclose") }
        let shape=try h.projectionShape(d)
        guard shape.count==1,shape[0]<=2_000_000 else { throw VivoOmicsError.limit("projection dataframe axis") }
        return Int(shape[0])
    }
    func indices(_ selection: [Int]?,length: Int) throws -> [Int] {
        let value=selection ?? Array(0..<length)
        guard value.allSatisfy({ $0<length }) else { throw VivoOmicsError.invalid("projection index outside source axis") }
        return value
    }
    func rawArray(_ source: Int64,_ destination: Int64,_ name: String,selections: [[Int]?],expected: [Int?],path: String,kind: String) throws {
        let dims=try h.projectionShape(source)
        guard dims.count>=selections.count,dims.count>=expected.count,!dims.isEmpty,dims.allSatisfy({ $0<=2_000_000 }) else { throw VivoOmicsError.invalid("projected array rank or axis bound: "+path) }
        let shape=dims.map(Int.init)
        if ["categorical","nullable-integer","nullable-boolean","nullable-string-array"].contains(kind) {
            guard shape.count==1 else { throw VivoOmicsError.invalid("encoded column components must be vectors") }
        }
        if path=="X" || path=="raw/X" || path.hasPrefix("layers/") {
            guard shape.count==2 else { throw VivoOmicsError.invalid("expression matrix must have two axes") }
        }
        for (i,n) in expected.enumerated() { if let n { guard shape[i]==n else { throw VivoOmicsError.invalid("misaligned array: "+path) } } }
        var axes=[[Int]?](repeating: nil,count: dims.count)
        for (i,s) in selections.enumerated() { if let s { axes[i]=try indices(s,length: shape[i]) } }
        let outputShape=shape.indices.map { axes[$0]?.count ?? shape[$0] }
        let count=outputShape.contains(0) ? 0 : try h.count(outputShape.map(UInt64.init),maximum: remaining)
        try consume(count)
        let type=try h.type(source,attribute: false);defer { h.close(type,"H5Tclose") }
        try reserveStorage(count*typeWidth(type),output: destination)
        let output=try h.projectionDataset(destination,name,type: type,shape: outputShape.map(UInt64.init));defer { h.close(output,"H5Dclose") }
        try h.projectionAttributes(source,output)
        if legacy,try !h.legacyHasAttribute(source,"encoding-type") { try h.encoding(output,kind,"0.2.0") }
        fields.append(.init(path: path,encoding: kind,sourceShape: shape,outputShape: outputShape))
        guard count>0 else { return }
        let rank=shape.count,last=rank-1,width=outputShape[last],outer=count/width
        for row in 0..<outer {
            try Task.checkCancellation()
            var prefix=[Int](repeating: 0,count: rank),q=row
            if rank>1 { for j in stride(from: last-1,through: 0,by: -1) { prefix[j]=q%outputShape[j];q/=outputShape[j] } }
            for start in stride(from: 0,to: width,by: 65_536) {
                let length=min(65_536,width-start)
                var coordinates=[UInt64]();coordinates.reserveCapacity(length*rank)
                for k in start..<(start+length) {
                    for j in 0..<rank { let index=j==last ? k : prefix[j];coordinates.append(UInt64(axes[j]?[index] ?? index)) }
                }
                var origin=prefix.map(UInt64.init);origin[last]=UInt64(start)
                var extent=[UInt64](repeating: 1,count: rank);extent[last]=UInt64(length)
                try h.projectionTransfer(source,output,type: type,coordinates: coordinates,count: length,outputStart: origin,outputCount: extent,
                    reserveVariableBytes: { try self.reserveStorage($0,output: output) })
            }
        }
    }
    func frame(_ source: Int64,_ destination: Int64,_ name: String,selection: [Int]?,length: Int,path: String) throws {
        if legacy,try h.legacyCompound(source) { try legacyFrame(source,destination,name,selection: selection,length: length,path: path);return }
        guard try frameLength(source)==length else { throw VivoOmicsError.invalid("misaligned dataframe: "+path) }
        let index=try h.text(source,"_index"),attribute=try h.attribute(source,"column-order");defer { h.close(attribute,"H5Aclose") }
        guard try h.projectionShape(attribute,attribute: true).count==1 else { throw VivoOmicsError.invalid("dataframe column-order must be a vector") }
        let columns=try h.strings(attribute,attribute: true,maximum: 100_000)
        guard Set(columns).count==columns.count,!columns.contains(index),Set(try h.projectionNames(source))==Set(columns+[index]) else { throw VivoOmicsError.invalid("dataframe columns disagree: "+path) }
        let group=try h.projectionGroup(destination,name);defer { h.close(group,"H5Gclose") }
        try h.projectionAttributes(source,group)
        for column in [index]+columns {
            try element(source,column,group,selections: [selection],expected: [length],path: path+"/"+column,frameAllowed: false)
        }
    }
    func sparse(_ source: Int64,_ destination: Int64,_ name: String,selections: [[Int]?],expected: [Int?],path: String,kind: String) throws {
        try encoding(source,kind,"0.1.0")
        guard Set(try h.projectionNames(source))==["data","indices","indptr"] else { throw VivoOmicsError.invalid("unknown sparse components") }
        let legacySparse=try legacy && !h.legacyHasAttribute(source,"encoding-type")
        let attr=try h.attribute(source,legacySparse ? "h5sparse_shape" : "shape");defer { h.close(attr,"H5Aclose") }
        let dims=try h.integers(attr,attribute: true,maximum: 2)
        guard dims.count==2,dims.allSatisfy({ $0<=2_000_000 }),selections.count<=2 else { throw VivoOmicsError.invalid("sparse projection shape") }
        let shape=dims.map(Int.init)
        for (i,n) in expected.enumerated() { guard i<2 else { throw VivoOmicsError.invalid("sparse rank") };if let n { guard shape[i]==n else { throw VivoOmicsError.invalid("misaligned sparse field: "+path) } } }
        let rows=try indices(selections.first ?? nil,length: shape[0]),cols=try indices(selections.count>1 ? selections[1] : nil,length: shape[1])
        let csr=kind=="csr_matrix",major=csr ? rows : cols,minor=csr ? cols : rows
        let majorLength=csr ? shape[0] : shape[1],minorLength=csr ? shape[1] : shape[0]
        var selectedMajor=[Bool](repeating: false,count: majorLength)
        for k in major { selectedMajor[k]=true }
        // Each source minor index owns an ordered chain of output positions.
        // This preserves sparse entry order while expanding repeated selections,
        // without a dense matrix or one heap allocation per source index.
        var minorHead=[Int](repeating: -1,count: minorLength)
        var minorNext=[Int](repeating: -1,count: minor.count)
        var minorMultiplicity=[Int](repeating: 0,count: minorLength)
        for i in minor.indices.reversed() {
            let k=minor[i];minorNext[i]=minorHead[k];minorHead[k]=i;minorMultiplicity[k]+=1
        }
        let data=try h.dataset(source,"data");defer { h.close(data,"H5Dclose") }
        let ii=try h.dataset(source,"indices");defer { h.close(ii,"H5Dclose") }
        let pp=try h.dataset(source,"indptr");defer { h.close(pp,"H5Dclose") }
        let ds=try h.shape(data,attribute: false),is_=try h.shape(ii,attribute: false)
        guard ds.count==1,ds==is_,try h.shape(pp,attribute: false)==[UInt64(majorLength+1)] else { throw VivoOmicsError.invalid("sparse component dimensions") }
        // Entry count bounds scan work, not resident payload. Check before Int conversion.
        guard ds[0]<=UInt64(remaining) else { throw VivoOmicsError.limit("H5AD projection element visits") }
        let ptr=try h.integers(pp,maximum: 2_000_001),nnz=Int(ds[0])
        guard ptr.first==0,ptr.last==ds[0],zip(ptr,ptr.dropFirst()).allSatisfy({ $0<=$1 }) else { throw VivoOmicsError.invalid("sparse offset domain") }
        try consume(nnz)
        var counts=[Int](repeating: 0,count: majorLength),sourceMajor=0
        for start in stride(from: 0,to: nnz,by: 65_536) {
            try Task.checkCancellation()
            let values=try h.integers(ii,maximum: 65_536,range: start..<min(start+65_536,nnz))
            for (j,index) in values.enumerated() {
                guard index<UInt64(minorLength) else { throw VivoOmicsError.invalid("sparse index outside axis") }
                while sourceMajor<majorLength && UInt64(start+j)>=ptr[sourceMajor+1] { sourceMajor+=1 }
                guard sourceMajor<majorLength else { throw VivoOmicsError.invalid("sparse offset coverage") }
                if selectedMajor[sourceMajor] { counts[sourceMajor]+=minorMultiplicity[Int(index)] }
            }
        }
        let scanned=major.reduce(0) { $0+Int(ptr[$1+1]-ptr[$1]) }
        try consume(scanned)
        var outputPtr: [Int64]=[0],total=0
        for row in major {
            // Bound expansion before summing or allocating. Both selected axes
            // can repeat, so an unchecked output-size product can overflow Int.
            guard counts[row]<=remaining-total else { throw VivoOmicsError.limit("H5AD projection element visits") }
            total+=counts[row];outputPtr.append(Int64(total))
        }
        try consume(total)
        let group=try h.projectionGroup(destination,name);defer { h.close(group,"H5Gclose") }
        try h.projectionAttributes(source,group,excluding: legacySparse ? ["h5sparse_shape","h5sparse_format"] : ["shape"])
        if legacySparse { try h.encoding(group,kind,"0.1.0") }
        try h.writeIntegers(group,"shape",[UInt64(rows.count),UInt64(cols.count)],attribute: true)
        let type=try h.type(data,attribute: false);defer { h.close(type,"H5Tclose") }
        try reserveStorage(total*(typeWidth(type)+8)+outputPtr.count*8,output: group)
        let od=try h.projectionDataset(group,"data",type: type,shape: [UInt64(total)]);defer { h.close(od,"H5Dclose") }
        let oi=try h.projectionDataset(group,"indices",type: h.native("NATIVE_LLONG"),shape: [UInt64(total)]);defer { h.close(oi,"H5Dclose") }
        let op=try h.projectionDataset(group,"indptr",type: h.native("NATIVE_LLONG"),shape: [UInt64(outputPtr.count)]);defer { h.close(op,"H5Dclose") }
        try h.projectionAttributes(data,od);try h.projectionAttributes(ii,oi);try h.projectionAttributes(pp,op)
        try h.projectionIntegers(op,values: outputPtr)
        var cursor=0
        for row in major {
            for start in stride(from: Int(ptr[row]),to: Int(ptr[row+1]),by: 65_536) {
                try Task.checkCancellation()
                let values=try h.integers(ii,maximum: 65_536,range: start..<min(start+65_536,Int(ptr[row+1])))
                var positions: [UInt64]=[],mapped: [Int64]=[]
                func flush() throws {
                    guard !positions.isEmpty else { return }
                    try Task.checkCancellation()
                    try h.projectionTransfer(data,od,type: type,coordinates: positions,count: positions.count,outputStart: [UInt64(cursor)],outputCount: [UInt64(positions.count)],
                        reserveVariableBytes: { try self.reserveStorage($0,output: od) })
                    try h.projectionIntegers(oi,values: mapped,offset: cursor);cursor+=mapped.count
                    positions.removeAll(keepingCapacity: true);mapped.removeAll(keepingCapacity: true)
                }
                for (j,index) in values.enumerated() {
                    var position=minorHead[Int(index)]
                    while position>=0 {
                        positions.append(UInt64(start+j));mapped.append(Int64(position))
                        if positions.count==65_536 { try flush() }
                        position=minorNext[position]
                    }
                }
                try flush()
            }
        }
        guard cursor==total else { throw VivoOmicsError.invalid("projection sparse pass disagreement") }
        fields.append(.init(path: path,encoding: kind,sourceShape: shape,outputShape: [rows.count,cols.count]))
    }
    func element(_ source: Int64,_ name: String,_ destination: Int64,selections: [[Int]?],expected: [Int?],path: String,frameAllowed: Bool = true,outputName: String? = nil) throws {
        guard fields.count<100_000 else { throw VivoOmicsError.limit("projection field count") }
        let object=try h.object(source,name);defer { h.close(object,"H5Oclose") }
        let kind=try elementKind(object)
        switch kind {
        case "array","string-array":
            try encoding(object,kind,"0.2.0")
            let dataset=try h.dataset(source,name);defer { h.close(dataset,"H5Dclose") }
            if !frameAllowed,selections.count==1 {
                guard try h.projectionShape(dataset).count==1 else { throw VivoOmicsError.invalid("dataframe columns must be vectors") }
            }
            try rawArray(dataset,destination,outputName ?? name,selections: selections,expected: expected,path: path,kind: kind)
        case "csr_matrix","csc_matrix":try sparse(object,destination,outputName ?? name,selections: selections,expected: expected,path: path,kind: kind)
        case "dataframe" where frameAllowed:
            guard selections.count==1,expected.count==1,let n=expected[0] else { throw VivoOmicsError.invalid("dataframe projection context") }
            try frame(object,destination,outputName ?? name,selection: selections[0],length: n,path: path)
        case "categorical","nullable-integer","nullable-boolean","nullable-string-array":
            guard selections.count==1,expected.count==1 else { throw VivoOmicsError.invalid("column projection context") }
            try encoding(object,kind,kind=="categorical" ? "0.2.0" : "0.1.0")
            let children=try h.projectionNames(object)
            guard Set(children)==(kind=="categorical" ? Set(["codes","categories"]) : Set(["values","mask"])) else { throw VivoOmicsError.invalid("unknown encoded column components") }
            let group=try h.projectionGroup(destination,outputName ?? name);defer { h.close(group,"H5Gclose") }
            try h.projectionAttributes(object,group)
            for child in children {
                if child=="categories" { try h.projectionCopy(object,child,group,child);continue }
                let dataset=try h.dataset(object,child);defer { h.close(dataset,"H5Dclose") }
                try rawArray(dataset,group,child,selections: selections,expected: expected,path: path+"/"+child,kind: kind)
            }
        default:throw VivoOmicsError.invalid("unsupported aligned encoding at "+path+": "+kind)
        }
    }
    func mapping(_ source: Int64,_ destination: Int64,_ name: String,selections: [[Int]?],expected: [Int?],path: String,outputName: String? = nil) throws {
        guard try h.exists(source,name) else { return }
        let input=try h.object(source,name);defer { h.close(input,"H5Oclose") }
        if legacy,try h.legacyCompound(input) { try legacyMapping(input,destination,outputName ?? name,selections: selections,expected: expected,path: path);return }
        try encoding(input,"dict","0.1.0")
        let output=try h.projectionGroup(destination,outputName ?? name);defer { h.close(output,"H5Gclose") };try h.projectionAttributes(input,output)
        if legacy,try !h.legacyHasAttribute(input,"encoding-type") { try h.encoding(output,"dict","0.1.0") }
        for key in try h.projectionNames(input) { try element(input,key,output,selections: selections,expected: expected,path: path+"/"+key) }
    }
}

public enum VivoH5ADProjection {
    static func evaluate(_ source: URL,plan: VivoH5ADProjectionPlan,to output: URL) throws -> VivoH5ADProjectionReport {
        try plan.validate()
        return try VivoHDF5.lock.withLock {
            let h=try VivoHDF5(),input=try h.file(source.path);defer { h.close(input,"H5Fclose") }
            h.projectionMaximumOutputBytes=plan.maximumOutputBytes ?? 1_073_741_824
            let legacy=try !h.legacyHasAttribute(input,"encoding-type") && !h.legacyHasAttribute(input,"encoding-version")
            let engine=VivoH5ADProjector(h,limit: plan.maximumElementVisits,root: input,legacy: legacy)
            try engine.encoding(input,"anndata","0.1.0");try h.validateAnnotationStorage(input)
            let children=Set(try h.projectionNames(input))
            guard children.isSubset(of: Set(["X","obs","var","layers","obsm","varm","obsp","varp","raw","uns"]+(legacy ? ["raw.X","raw.var","raw.varm"] : []))) else { throw VivoOmicsError.invalid("unknown AnnData root fields cannot be aligned") }
            let obs=try h.object(input,"obs");defer { h.close(obs,"H5Oclose") }
            let vars=try h.object(input,"var");defer { h.close(vars,"H5Oclose") }
            let n=try engine.frameLength(obs),p=try engine.frameLength(vars)
            let oi=try engine.indices(plan.observationIndices,length: n),vi=try engine.indices(plan.featureIndices,length: p)
            let out=try h.projectionFile(output.path);defer { h.close(out,"H5Fclose") }
            try h.projectionAttributes(input,out)
            if legacy { try h.encoding(out,"anndata","0.1.0") }
            try engine.frame(obs,out,"obs",selection: oi,length: n,path: "obs")
            try engine.frame(vars,out,"var",selection: vi,length: p,path: "var")
            if children.contains("X") { try engine.element(input,"X",out,selections: [oi,vi],expected: [n,p],path: "X",frameAllowed: false) }
            try engine.mapping(input,out,"layers",selections: [oi,vi],expected: [n,p],path: "layers")
            try engine.mapping(input,out,"obsm",selections: [oi],expected: [n],path: "obsm")
            try engine.mapping(input,out,"varm",selections: [vi],expected: [p],path: "varm")
            try engine.mapping(input,out,"obsp",selections: [oi,oi],expected: [n,n],path: "obsp")
            try engine.mapping(input,out,"varp",selections: [vi,vi],expected: [p,p],path: "varp")
            let dotted=children.contains("raw.X") || children.contains("raw.var") || children.contains("raw.varm")
            guard !(dotted && children.contains("raw")) else { throw VivoOmicsError.invalid("mixed legacy raw representations") }
            if children.contains("raw") || dotted {
                let raw=try h.object(input,dotted ? "." : "raw");defer { h.close(raw,"H5Oclose") }
                let varName=dotted ? "raw.var" : "var",xName=dotted ? "raw.X" : "X",varmName=dotted ? "raw.varm" : "varm"
                if !dotted {
                    try engine.encoding(raw,"raw","0.1.0")
                    guard Set(try h.projectionNames(raw)).isSubset(of: ["X","var","varm"]) else { throw VivoOmicsError.invalid("unknown raw fields") }
                }
                let rg=try h.projectionGroup(out,"raw");defer { h.close(rg,"H5Gclose") }
                if !dotted { try h.projectionAttributes(raw,rg) }
                if try dotted || (legacy && !h.legacyHasAttribute(raw,"encoding-type")) { try h.encoding(rg,"raw","0.1.0") }
                let rv=try h.object(raw,varName);defer { h.close(rv,"H5Oclose") };let rawP=try engine.frameLength(rv)
                try engine.frame(rv,rg,"var",selection: nil,length: rawP,path: "raw/var")
                if try h.exists(raw,xName) { try engine.element(raw,xName,rg,selections: [oi,nil],expected: [n,rawP],path: "raw/X",frameAllowed: false,outputName: "X") }
                try engine.mapping(raw,rg,varmName,selections: [nil],expected: [rawP],path: "raw/varm",outputName: "varm")
            }
            if children.contains("uns") {
                try h.projectionCopy(input,"uns",out,"uns")
                let uns=try h.object(out,"uns");defer { h.close(uns,"H5Oclose") }
                // AnnData accepts legacy dicts through its fallback reader;
                // downstream native annotation requires a versioned parent.
                // Preserve every child and add only the missing dict tag.
                if legacy,try h.legacyObjectKind(uns)==2,try !h.legacyHasAttribute(uns,"encoding-type"),try !h.legacyHasAttribute(uns,"encoding-version") {
                    try h.encoding(uns,"dict","0.1.0")
                }
                let remove: @convention(c) (Int64,UnsafePointer<CChar>,Int64) -> Int32 = try h.symbol("H5Ldelete")
                for name in engine.movedLegacyCategories.sorted() { try h.check(remove(uns,name,0),"relocate legacy categories") }
            }
            let repeated=Set(oi).count != oi.count || Set(vi).count != vi.count
            return try .init(method: legacy ? "native-H5AD-legacy-axis-projection-v2" : repeated ? "native-H5AD-repeated-axis-projection-v1" : "native-H5AD-axis-projection-v1",sourceShape: [n,p],outputShape: [oi.count,vi.count],elementVisits: plan.maximumElementVisits-engine.remaining,
                fields: engine.fields,hdf5Version: h.version(),qualification: (legacy ? "Legacy compound annotations and embeddings split without numeric conversion; category definitions relocated from uns and legacy uns mapping versioned; original source retained. Cell/feature selection and reordering;" : repeated ? "Repeated cell/feature selection and reordering;" : "Unique cell/feature selection and reordering;")+" raw retains its feature axis. Stored values/datatypes and categories retained; sparse structural indices become int64. Unstructured data copied without inferred axis semantics. No biological validation claim.")
        }
    }
    public static func publish(source: URL,plan: VivoH5ADProjectionPlan,implementation: VivoFingerprint,to destination: URL) throws -> VivoH5ADProjectionReceipt {
        try plan.validate();let bytes=try VivoCanonicalJSON.encode(plan)
        guard bytes.count<=64*1_024*1_024 else { throw VivoOmicsError.limit("projection plan bytes") }
        guard !FileManager.default.fileExists(atPath: destination.path) else { throw VivoOmicsError.invalid("projection output exists") }
        let staging=destination.deletingLastPathComponent().appendingPathComponent(".numivivo-project-"+UUID().uuidString)
        try FileManager.default.createDirectory(at: staging,withIntermediateDirectories: false,attributes: [.posixPermissions: 0o700])
        defer { try? FileManager.default.removeItem(at: staging) }
        let snapshot=staging.appendingPathComponent("original.h5ad"),output=staging.appendingPathComponent("projected.h5ad")
        guard try VivoH5ADPseudobulk.fingerprint(source,copyTo: snapshot)==plan.source else { throw VivoOmicsError.invalid("projection source fingerprint mismatch") }
        let report=try evaluate(snapshot,plan: plan,to: output),reportBytes=try VivoCanonicalJSON.encode(report)
        guard reportBytes.count<=64*1_024*1_024 else { throw VivoOmicsError.limit("projection report exceeds replay allowance") }
        let receipt=try VivoH5ADProjectionReceipt(schemaVersion: 1,source: plan.source,plan: VivoCanonicalJSON.fingerprint(bytes),output: VivoH5ADPseudobulk.fingerprint(output),
            report: VivoCanonicalJSON.fingerprint(reportBytes),implementation: implementation)
        try bytes.write(to: staging.appendingPathComponent("plan.json"),options: .withoutOverwriting)
        try reportBytes.write(to: staging.appendingPathComponent("report.json"),options: .withoutOverwriting)
        try VivoCanonicalJSON.encode(receipt).write(to: staging.appendingPathComponent("receipt.json"),options: .withoutOverwriting)
        try Task.checkCancellation();try FileManager.default.moveItem(at: staging,to: destination);return receipt
    }
    public static func verify(_ directory: URL,implementation: VivoFingerprint) throws -> VivoH5ADProjectionReport {
        let receipt=try VivoCanonicalJSON.decode(VivoH5ADProjectionReceipt.self,from: VivoSingleCellCampaignIO.readDocument(directory.appendingPathComponent("receipt.json"),maximumBytes: 65_536))
        let bytes=try VivoSingleCellCampaignIO.readDocument(directory.appendingPathComponent("plan.json"),maximumBytes: 64*1_024*1_024)
        let reportBytes=try VivoSingleCellCampaignIO.readDocument(directory.appendingPathComponent("report.json"),maximumBytes: 64*1_024*1_024)
        guard receipt.schemaVersion==1,receipt.implementation==implementation,try VivoCanonicalJSON.fingerprint(bytes)==receipt.plan,
              try VivoCanonicalJSON.fingerprint(reportBytes)==receipt.report,try VivoH5ADPseudobulk.fingerprint(directory.appendingPathComponent("projected.h5ad"))==receipt.output else {
            throw VivoOmicsError.invalid("projection output, plan, report or implementation changed")
        }
        let temporary=FileManager.default.temporaryDirectory.appendingPathComponent("numivivo-project-verify-"+UUID().uuidString)
        try FileManager.default.createDirectory(at: temporary,withIntermediateDirectories: false,attributes: [.posixPermissions: 0o700]);defer { try? FileManager.default.removeItem(at: temporary) }
        let snapshot=temporary.appendingPathComponent("original.h5ad"),output=temporary.appendingPathComponent("projected.h5ad")
        guard try VivoH5ADPseudobulk.fingerprint(directory.appendingPathComponent("original.h5ad"),copyTo: snapshot)==receipt.source else { throw VivoOmicsError.invalid("projection original changed") }
        let plan=try VivoCanonicalJSON.decode(VivoH5ADProjectionPlan.self,from: bytes)
        guard plan.source==receipt.source else { throw VivoOmicsError.invalid("projection plan source differs") }
        let report=try evaluate(snapshot,plan: plan,to: output)
        guard try VivoCanonicalJSON.encode(report)==reportBytes,try VivoH5ADPseudobulk.fingerprint(output)==receipt.output else { throw VivoOmicsError.invalid("projection does not reconstruct") }
        return report
    }
}

private extension VivoH5ADProjector {
    func legacyScalarKind(_ type: Int64) throws -> String {
        switch try h.legacyTypeKind(type) {
        case 0,1,8:return "array"
        case 3:return "string-array"
        default:throw VivoOmicsError.invalid("unsupported legacy scalar datatype")
        }
    }
    func elementKind(_ object: Int64) throws -> String {
        if try h.legacyHasAttribute(object,"encoding-type") { return try h.text(object,"encoding-type") }
        guard legacy,try !h.legacyHasAttribute(object,"encoding-version") else { throw VivoOmicsError.invalid("missing aligned encoding") }
        if try h.legacyObjectKind(object)==2 {
            guard try h.legacyHasAttribute(object,"h5sparse_format") else { throw VivoOmicsError.invalid("unsupported legacy aligned group") }
            let format=try h.text(object,"h5sparse_format")
            guard ["csr","csc"].contains(format) else { throw VivoOmicsError.invalid("unsupported legacy sparse format") };return format+"_matrix"
        }
        guard try h.legacyObjectKind(object)==5 else { throw VivoOmicsError.invalid("legacy aligned object kind") }
        let t=try h.type(object,attribute: false);defer { h.close(t,"H5Tclose") };return try legacyScalarKind(t)
    }
    func legacyMember(_ source: Int64,_ destination: Int64,_ name: String,member: String,type: Int64,rows: [Int],length: Int,trailing: [UInt64],path: String) throws {
        guard fields.count<100_000,trailing.allSatisfy({ $0>0 && $0<=2_000_000 }) else { throw VivoOmicsError.limit("legacy field shape") }
        let base=try h.legacyBaseType(type);defer { h.close(base,"H5Tclose") }
        let kind=try legacyScalarKind(base),components=try h.count(trailing,maximum: 65_536)
        guard rows.count<=remaining/components else { throw VivoOmicsError.limit("legacy aligned element visits") }
        let count=rows.count*components;try consume(count);try reserveStorage(count*typeWidth(base),output: destination)
        let shape=[UInt64(rows.count)]+trailing
        let out=try h.projectionDataset(destination,name,type: base,shape: shape);defer { h.close(out,"H5Dclose") }
        try h.encoding(out,kind,"0.2.0")
        let memory=try h.legacyMemoryType(member,type);defer { h.close(memory,"H5Tclose") }
        guard try typeWidth(memory)==components*typeWidth(base) else { throw VivoOmicsError.invalid("legacy member storage width") }
        for start in stride(from: 0,to: rows.count,by: 65_536/components) {
            let end=min(rows.count,start+65_536/components)
            try h.legacyTransfer(source,out,memoryType: memory,outputType: base,rows: rows[start..<end].map(UInt64.init),components: components,
                outputStart: [UInt64(start)]+trailing.map { _ in 0 },outputCount: [UInt64(end-start)]+trailing,
                reserveVariableBytes: { try self.reserveStorage($0,output: out) })
        }
        fields.append(.init(path: path,encoding: "legacy-compound/"+kind,sourceShape: [length]+trailing.map(Int.init),outputShape: shape.map(Int.init)))
    }
    func legacyCategory(_ source: Int64,_ destination: Int64,column: String,type: Int64,rows: [Int],length: Int,path: String) throws -> Bool {
        guard ["obs","var"].contains(path),try h.exists(root,"uns") else { return false }
        let obs=try h.object(root,"obs");defer { h.close(obs,"H5Oclose") }
        guard try h.legacyCompound(obs) else { return false }
        let uns=try h.object(root,"uns");defer { h.close(uns,"H5Oclose") }
        let matches=try h.projectionNames(uns).filter { $0.hasSuffix("_categories") && $0.replacingOccurrences(of: "_categories",with: "")==column }
        guard matches.count<=1 else { throw VivoOmicsError.invalid("ambiguous legacy category definitions") }
        guard let key=matches.first else { return false }
        guard try h.legacyTypeKind(type)==0,try typeWidth(type)<=8 else { throw VivoOmicsError.invalid("legacy category codes must be bounded integers") }
        let cats=try h.dataset(uns,key);defer { h.close(cats,"H5Dclose") }
        let shape=try h.projectionShape(cats)
        guard shape.count<=1,shape.allSatisfy({ $0<=100_000 }) else { throw VivoOmicsError.invalid("legacy categories must be bounded scalar or vector") }
        let count=shape.isEmpty ? 1 : Int(shape[0]);guard count<=100_000 else { throw VivoOmicsError.limit("legacy categories") }
        var tooLarge=false,tooSmall=false
        for start in stride(from: 0,to: length,by: 65_536) {
            let range=start..<min(length,start+65_536);try consume(range.count)
            let codes=try h.legacyCodes(source,member: column,range: range)
            tooLarge = tooLarge || codes.contains { $0>=Int64(count) };tooSmall = tooSmall || codes.contains { $0 < -1 }
        }
        // AnnData checks every original code before slicing or moving categories.
        if tooLarge { return false }
        guard !tooSmall else { throw VivoOmicsError.invalid("legacy categorical code below missing sentinel") }
        let ct=try h.type(cats,attribute: false);defer { h.close(ct,"H5Tclose") }
        guard try h.legacyTypeKind(ct)==3 else { throw VivoOmicsError.invalid("legacy category migration currently requires string labels") }
        if !shape.isEmpty { try h.legacyStringReadBound(cats,type: ct) }
        let labels=try shape.isEmpty ? [h.legacyScalarString(cats,type: ct)] : h.strings(cats,maximum: 100_000)
        guard labels.count==count,Set(labels).count==count else { throw VivoOmicsError.invalid("duplicate legacy category labels") }
        let group=try h.projectionGroup(destination,column);defer { h.close(group,"H5Gclose") }
        try h.encoding(group,"categorical","0.2.0");try h.writeBooleans(group,"ordered",shape: [],values: [false],attribute: true)
        try legacyMember(source,group,"codes",member: column,type: type,rows: rows,length: length,trailing: [],path: path+"/"+column+"/codes")
        try consume(count);try reserveStorage(count*typeWidth(ct),output: group)
        let output=try h.projectionDataset(group,"categories",type: ct,shape: [UInt64(count)]);defer { h.close(output,"H5Dclose") }
        try h.projectionAttributes(cats,output,excluding: ["encoding-type","encoding-version"]);try h.encoding(output,"string-array","0.2.0")
        for start in stride(from: 0,to: count,by: 65_536) {
            let end=min(count,start+65_536)
            try h.legacyTransfer(cats,output,memoryType: ct,outputType: ct,rows: (start..<end).map(UInt64.init),components: 1,
                outputStart: [UInt64(start)],outputCount: [UInt64(end-start)],scalarSource: shape.isEmpty,
                reserveVariableBytes: { try self.reserveStorage($0,output: output) })
        }
        movedLegacyCategories.insert(key);return true
    }
    func legacyFrame(_ source: Int64,_ destination: Int64,_ name: String,selection: [Int]?,length: Int,path: String) throws {
        guard try frameLength(source)==length else { throw VivoOmicsError.invalid("legacy dataframe alignment") }
        guard try Set(h.projectionNames(source,attributes: true)).isDisjoint(with: ["encoding-type","encoding-version","_index","column-order"]) else { throw VivoOmicsError.invalid("contradictory legacy dataframe metadata") }
        let names=try h.legacyMembers(source),rows=try indices(selection,length: length)
        guard Set(names).count==names.count else { throw VivoOmicsError.invalid("duplicate legacy dataframe fields") }
        let group=try h.projectionGroup(destination,name);defer { h.close(group,"H5Gclose") }
        try h.projectionAttributes(source,group,excluding: ["encoding-type","encoding-version","_index","column-order"])
        try h.encoding(group,"dataframe","0.2.0");try h.writeStrings(group,"_index",[names[0]],attribute: true,scalar: true)
        try h.writeStrings(group,"column-order",Array(names.dropFirst()),attribute: true)
        for (i,column) in names.enumerated() {
            let type=try h.legacyMemberType(source,i);defer { h.close(type,"H5Tclose") }
            guard try h.legacyArrayShape(type).isEmpty else { throw VivoOmicsError.invalid("legacy dataframe field must be scalar") }
            if i==0 { guard try h.legacyTypeKind(type)==3 else { throw VivoOmicsError.invalid("legacy dataframe index must be string") } }
            if i>0,try legacyCategory(source,group,column: column,type: type,rows: rows,length: length,path: path) { continue }
            if try h.legacyTypeKind(type)==3 {
                // Preserve string semantics on empty axes as well as nonempty
                // pandas 3 frames, without changing the original string bytes.
                let text=try h.projectionGroup(group,column);defer { h.close(text,"H5Gclose") }
                try h.encoding(text,"nullable-string-array","0.1.0")
                try h.writeStrings(text,"na-value",["NaN"],attribute: true,scalar: true)
                try legacyMember(source,text,"values",member: column,type: type,rows: rows,length: length,trailing: [],path: path+"/"+column+"/values")
                try consume(rows.count);try reserveStorage(rows.count,output: text)
                try h.legacyFalseMask(text,count: rows.count)
            } else {
                try legacyMember(source,group,column,member: column,type: type,rows: rows,length: length,trailing: [],path: path+"/"+column)
            }
        }
    }
    func legacyMapping(_ source: Int64,_ destination: Int64,_ name: String,selections: [[Int]?],expected: [Int?],path: String) throws {
        guard ["obsm","varm","raw/varm"].contains(path),selections.count==1,expected.count==1,let length=expected[0],try h.projectionShape(source)==[UInt64(length)] else { throw VivoOmicsError.invalid("legacy compound mapping alignment") }
        guard try Set(h.projectionNames(source,attributes: true)).isDisjoint(with: ["encoding-type","encoding-version"]) else { throw VivoOmicsError.invalid("contradictory legacy mapping metadata") }
        let rows=try indices(selections[0],length: length),names=try h.legacyMembers(source)
        let output=try h.projectionGroup(destination,name);defer { h.close(output,"H5Gclose") }
        try h.projectionAttributes(source,output,excluding: ["encoding-type","encoding-version"]);try h.encoding(output,"dict","0.1.0")
        for (i,member) in names.enumerated() {
            let type=try h.legacyMemberType(source,i);defer { h.close(type,"H5Tclose") }
            let trailing=try h.legacyArrayShape(type)
            guard !trailing.isEmpty else { throw VivoOmicsError.invalid("legacy embedding member must be a fixed array") }
            try legacyMember(source,output,member,member: member,type: type,rows: rows,length: length,trailing: trailing,path: path+"/"+member)
        }
    }
}
