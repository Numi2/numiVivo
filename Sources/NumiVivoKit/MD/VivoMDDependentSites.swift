import Foundation
@preconcurrency import Metal
import NumiVivoShaders

private struct VivoMDDependentSiteRecord {
    var identity: SIMD4<UInt32>
    var parents: SIMD4<UInt32>
    var originWeights: SIMD4<Float>
    var xWeights: SIMD4<Float>
    var yWeights: SIMD4<Float>
    var local: SIMD4<Float>
}
private struct VivoMDDependentSitePass { var particleCount,siteCount,depth,reserved: UInt32 }

/// GPU implementation of VivoVirtualSiteGraph. Used for every site when any
/// generalized definition is present; legacy-only systems retain their exact
/// established numerical profile. No CPU coordinate readback is introduced.
final class VivoMDDependentSites: @unchecked Sendable {
    let graph: VivoVirtualSiteGraph
    private let records,jacobians,offsets,incidences: MTLBuffer
    private let construction,velocities,forces,zero: NumiVivoPipeline
    static func make(device: MTLDevice,catalog: NumiVivoPipelineCatalog,system: VivoClassicalSystem) async throws -> VivoMDDependentSites {
        let graph = try system.resolvedVirtualSiteGraph()
        guard MemoryLayout<VivoMDDependentSiteRecord>.stride == 96 else { throw VivoMDRuntimeError.metal("dependent site ABI mismatch") }
        func float4(_ x: [Double]) throws -> SIMD4<Float> {
            guard x.count <= 4,x.allSatisfy({ $0.isFinite && Float($0).isFinite && ($0 == 0 || Float($0) != 0) }) else {
                throw VivoMDRuntimeError.unsupported(["dependent-site parameters not representable in FP32"])
            }
            var output = SIMD4<Float>.zero
            for i in x.indices { output[i] = Float(x[i]) };return output
        }
        var source: [VivoMDDependentSiteRecord] = [],buckets = [[SIMD2<UInt32>]](repeating: [],count: graph.particleCount)
        for (index,site) in graph.sites.enumerated() {
            var parent = SIMD4<UInt32>.zero,o = SIMD4<Float>.zero,x = o,y = o,local = o,kind: UInt32 = 0
            for i in site.parentParticles.indices {
                parent[i] = site.parentParticles[i];buckets[Int(parent[i])].append(.init(UInt32(index),UInt32(i)))
            }
            switch site.rule {
            case .linear(let weights): o = try float4(weights)
            case .outOfPlane(let a,let b,let c): kind = 1;local = try float4([a,b,c])
            case .localCoordinates(let ow,let xw,let yw,let p):
                kind = 2;o = try float4(ow);x = try float4(xw);y = try float4(yw);local = try float4([p.x,p.y,p.z])
            }
            source.append(.init(identity: .init(site.siteParticle,UInt32(site.parentParticles.count),UInt32(graph.depths[index]),kind),
                parents: parent,originWeights: o,xWeights: x,yWeights: y,local: local))
        }
        func buffer<T>(_ data: [T],_ label: String) throws -> MTLBuffer {
            let bytes = data.count.multipliedReportingOverflow(by: MemoryLayout<T>.stride)
            guard !bytes.overflow,bytes.partialValue <= device.maxBufferLength,
                  let buffer = device.makeBuffer(length: max(16,bytes.partialValue),options: .storageModeShared) else {
                throw VivoMDRuntimeError.metal("dependent-site \(label) allocation")
            }
            if !data.isEmpty { data.withUnsafeBytes { buffer.contents().copyMemory(from: $0.baseAddress!,byteCount: $0.count) } }
            buffer.label = "NumiVivo.DependentSites."+label;return buffer
        }
        var offset: [UInt32] = [0],incidence: [SIMD2<UInt32>] = []
        for bucket in buckets {
            incidence += bucket
            guard let count = UInt32(exactly: incidence.count) else { throw VivoMDRuntimeError.metal("dependent-site incidence width") }
            offset.append(count)
        }
        guard graph.sites.count <= device.maxBufferLength/192,
              let derivatives = device.makeBuffer(length: max(16,graph.sites.count*192),options: .storageModePrivate) else {
            throw VivoMDRuntimeError.metal("dependent-site Jacobian allocation")
        }
        return try await .init(graph: graph,records: buffer(source,"records"),jacobians: derivatives,
            offsets: buffer(offset,"offsets"),incidences: buffer(incidence,"incidences"),
            construction: catalog.pipeline(.mdConstructDependentSites),velocities: catalog.pipeline(.mdDependentSiteVelocities),
            forces: catalog.pipeline(.mdDependentSiteForces),zero: catalog.pipeline(.mdZeroDependentSiteForces))
    }
    private init(graph: VivoVirtualSiteGraph,records: MTLBuffer,jacobians: MTLBuffer,offsets: MTLBuffer,incidences: MTLBuffer,
                 construction: NumiVivoPipeline,velocities: NumiVivoPipeline,forces: NumiVivoPipeline,zero: NumiVivoPipeline) {
        self.graph = graph;self.records = records;self.jacobians = jacobians;self.offsets = offsets;self.incidences = incidences
        self.construction = construction;self.velocities = velocities;self.forces = forces;self.zero = zero
    }
    func construct(_ buffer: MTLCommandBuffer,positions: MTLBuffer,status: MTLBuffer,command: VivoMDMetalCommand) throws {
        guard !graph.sites.isEmpty else { return }
        for depth in 0...graph.maximumDepth {
            try dispatch(construction,buffer: buffer,buffers: [positions,records,jacobians,status],depth: depth,command: command,particles: false)
        }
    }
    func normalize(_ buffer: MTLCommandBuffer,positions: MTLBuffer,velocity: MTLBuffer,status: MTLBuffer,command: VivoMDMetalCommand) throws {
        try construct(buffer,positions: positions,status: status,command: command)
        guard !graph.sites.isEmpty else { return }
        for depth in 0...graph.maximumDepth {
            try dispatch(velocities,buffer: buffer,buffers: [velocity,records,jacobians,status],depth: depth,particles: false)
        }
    }
    func redistribute(_ buffer: MTLCommandBuffer,forceEnergy: MTLBuffer,status: MTLBuffer) throws {
        guard !graph.sites.isEmpty else { return }
        for depth in stride(from: graph.maximumDepth,through: 0,by: -1) {
            try dispatch(forces,buffer: buffer,buffers: [forceEnergy,records,jacobians,offsets,incidences,status],depth: depth,particles: true)
        }
        try dispatch(zero,buffer: buffer,buffers: [forceEnergy,records],depth: 0,particles: false)
    }
    private func dispatch(_ pipeline: NumiVivoPipeline,buffer: MTLCommandBuffer,buffers: [MTLBuffer],depth: Int,
                          command: VivoMDMetalCommand? = nil,particles: Bool) throws {
        guard let encoder = buffer.makeComputeCommandEncoder() else { throw VivoMDRuntimeError.metal("dependent-site encoder") }
        defer { encoder.endEncoding() }
        encoder.setComputePipelineState(pipeline.state)
        for (index,resource) in buffers.enumerated() { encoder.setBuffer(resource,offset: 0,index: index) }
        var index = buffers.count
        if var command { encoder.setBytes(&command,length: MemoryLayout<VivoMDMetalCommand>.stride,index: index);index += 1 }
        var pass = VivoMDDependentSitePass(particleCount: UInt32(graph.particleCount),siteCount: UInt32(graph.sites.count),depth: UInt32(depth),reserved: 0)
        encoder.setBytes(&pass,length: MemoryLayout<VivoMDDependentSitePass>.stride,index: index)
        let count = particles ? graph.particleCount : graph.sites.count
        encoder.dispatchThreads(pipeline.gridSize(for: count),threadsPerThreadgroup: pipeline.threadgroupSize(for: count))
    }
}
