import Foundation
#if canImport(Metal)
@preconcurrency import Metal
#endif

/// Exhaustive candidate search with explicit FP32 distance arithmetic and bounded
/// score/distance tiles. Graph assembly and published distances remain FP64.
enum VivoMetalPCANeighbors {
    static let method = "metal-FP32-windowed-PCA-knn-v1"
    static let qualification = "Exhaustive candidate search with Metal FP32 squared distances, widened before FP64 square root and fuzzy graph assembly. Row/index ties are deterministic within this backend; FP32 may alter membership, ordering and weights versus FP64. Bounded score and distance tiles; input metadata and graph state retain separate resident limits. No cross-device bitwise, biological-preservation or million-cell qualification."
    static let source = """
    #include <metal_stdlib>
    using namespace metal;
    kernel void distances(device const float *q [[buffer(0)]], device const float *c [[buffer(1)]],
        device float *out [[buffer(2)]], constant uint4 &a [[buffer(3)]], constant uint &origin [[buffer(4)]],
        uint2 i [[thread_position_in_grid]]) {
        if(i.x>=a.x || i.y>=a.y) return;
        float sum=0;
        if(origin+i.y != a.w+i.x) {
            for(uint d=0;d<a.z;d++) { float delta=q[i.y*a.z+d]-c[i.x*a.z+d]; sum=sum+delta*delta; }
        }
        out[i.y*a.x+i.x]=sum;
    }
    """
    static func stream(source url: URL, rows n: Int, dimensions d: Int, options: VivoSingleCellNeighborOptions,
                       execution: VivoPCANeighborExecution, sink: (Int,[Int],[Double]) throws -> Void) throws {
        try options.validate(); try execution.validate(); try Task.checkCancellation()
        guard execution.backend == .metalFP32, n >= options.neighbors, n <= VivoPCAStorageLimits.maximumRows,
              (1...64).contains(d), n*(n-1)/2 <= options.maximumDistancePairs else { throw VivoOmicsError.limit("Metal PCA neighbor axes or pair budget") }
        #if canImport(Metal)
        guard let device=MTLCreateSystemDefaultDevice(), device.hasUnifiedMemory, device.supportsFamily(.apple7),
              !device.name.lowercased().contains("paravirtual"), let queue=device.makeCommandQueue() else {
            throw VivoOmicsError.invalid("Metal PCA neighbors require a physical Apple GPU")
        }
        let compile=MTLCompileOptions(); compile.fastMathEnabled=false
        let library=try device.makeLibrary(source: source, options: compile)
        guard let function=library.makeFunction(name:"distances") else { throw VivoOmicsError.invalid("Metal PCA neighbor kernel") }
        let pipeline=try device.makeComputePipelineState(function:function)
        let qr=execution.queryBlockRows, cr=execution.candidateBlockRows
        guard let qbuffer=device.makeBuffer(length:qr*d*4,options:.storageModeShared),
              let cbuffer=device.makeBuffer(length:cr*d*4,options:.storageModeShared),
              let obuffer=device.makeBuffer(length:qr*cr*4,options:.storageModeShared) else { throw VivoOmicsError.limit("Metal PCA neighbor tile allocation") }
        let qp=qbuffer.contents().assumingMemoryBound(to:Float.self), cp=cbuffer.contents().assumingMemoryBound(to:Float.self)
        let distances=obuffer.contents().assumingMemoryBound(to:Float.self)
        let reader=try VivoPCAScoreReader(url,rows:n,columns:d), k=options.neighbors
        func copy(_ values:[Double], to pointer:UnsafeMutablePointer<Float>) throws {
            for i in values.indices { let value=Float(values[i]); guard value.isFinite else { throw VivoOmicsError.invalid("PCA score exceeds FP32 range") }; pointer[i]=value }
        }
        for first in stride(from:0,to:n,by:qr) {
            try Task.checkCancellation(); let count=min(qr,n-first)
            try copy(reader.readRows(first..<first+count),to:qp)
            var heaps=Array(repeating:[VivoSingleCellNeighbors.Neighbor](),count:count)
            for start in stride(from:0,to:n,by:cr) {
                try Task.checkCancellation(); let candidates=min(cr,n-start)
                try copy(reader.readRows(start..<start+candidates),to:cp)
                guard let command=queue.makeCommandBuffer(), let encoder=command.makeComputeCommandEncoder() else { throw VivoOmicsError.invalid("Metal PCA command allocation") }
                encoder.setComputePipelineState(pipeline); encoder.setBuffer(qbuffer,offset:0,index:0); encoder.setBuffer(cbuffer,offset:0,index:1); encoder.setBuffer(obuffer,offset:0,index:2)
                var axes=SIMD4<UInt32>(UInt32(candidates),UInt32(count),UInt32(d),UInt32(start)), origin=UInt32(first)
                encoder.setBytes(&axes,length:16,index:3);encoder.setBytes(&origin,length:4,index:4)
                let width=min(pipeline.threadExecutionWidth,pipeline.maxTotalThreadsPerThreadgroup)
                encoder.dispatchThreads(MTLSize(width:candidates,height:count,depth:1),threadsPerThreadgroup:MTLSize(width:width,height:1,depth:1))
                encoder.endEncoding();command.commit();command.waitUntilCompleted()
                guard command.status == .completed else { throw VivoOmicsError.invalid("Metal PCA distance command failed") }
                try Task.checkCancellation()
                for local in 0..<count {
                    for other in 0..<candidates where start+other != first+local {
                        let value=distances[local*candidates+other]
                        guard value.isFinite, value>=0 else { throw VivoOmicsError.invalid("nonfinite Metal PCA distance") }
                        VivoSingleCellNeighbors.retain(.init(index:start+other,squaredDistance:Double(value)),in:&heaps[local],capacity:k-1)
                    }
                }
            }
            for local in 0..<count {
                try Task.checkCancellation()
                let sorted=heaps[local].sorted(by:{$0.precedes($1)})
                try sink(first+local,[first+local]+sorted.map(\.index),[0]+sorted.map{sqrt($0.squaredDistance)})
            }
        }
        #else
        throw VivoOmicsError.invalid("Metal PCA neighbors unavailable on this platform")
        #endif
    }
}
