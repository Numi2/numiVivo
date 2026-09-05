import Foundation
@preconcurrency import Metal
import NumiVivoShaders

public struct VivoMDBarostatMoveCertificate: Codable, Sendable, Equatable {
    public var stepIndex: UInt64
    public var attempted: Bool
    public var accepted: Bool
    public var volumeBeforeNM3: Double
    public var volumeAfterNM3: Double
    public var logVolumeDelta: Double
    public var logAcceptanceRatio: Double
    /// Non-nil for an inadmissible geometric proposal, not a Metropolis draw.
    public var rejectionReason: String? = nil
}

/// Molecule membership includes every explicit bonded term, holonomic constraint
/// and virtual-site parent relation. BFS order gives an unambiguous unwrapping
/// tree; non-tree edges are checked for periodic winding on the GPU.
public struct VivoMDBarostatPlan: Sendable, Equatable {
    public var componentCount: UInt32
    public var componentOffsets: [UInt32]
    public var componentParticles: [UInt32]
    public var componentIndexByParticle: [UInt32]
    public var parentByParticle: [UInt32]
    public var edgeOffsets: [UInt32]
    public var edges: [SIMD2<UInt32>]

    public static func make(system: VivoClassicalSystem) throws -> Self {
        try VivoClassicalSystemValidator.validate(system)
        let n = system.particles.count
        var pairKeys = Set<UInt64>()
        func connect(_ a: UInt32, _ b: UInt32) throws {
            guard a != b, Int(a) < n, Int(b) < n,
                  system.particles[Int(a)].massDa > 0, system.particles[Int(b)].massDa > 0 else {
                throw VivoMDRuntimeError.unsupported(["NPT bonded/constraint connectivity must join massive particles"])
            }
            pairKeys.insert(UInt64(min(a,b)) << 32 | UInt64(max(a,b)))
        }
        for q in system.bonds { try connect(q.a,q.b) }
        for q in system.constraints { try connect(q.a,q.b) }
        for q in system.angles { try connect(q.a,q.b); try connect(q.b,q.c) }
        for q in system.torsions { try connect(q.a,q.b); try connect(q.b,q.c); try connect(q.c,q.d) }
        for site in system.linearVirtualSites ?? [] {
            guard let first = site.parentParticles.first else { throw VivoMDRuntimeError.metal("empty virtual-site parents") }
            for parent in site.parentParticles.dropFirst() { try connect(first,parent) }
        }
        guard pairKeys.count <= Int(UInt32.max) else { throw VivoMDRuntimeError.metal("NPT edge count exceeds UInt32") }
        let pairs = pairKeys.sorted().map { SIMD2<UInt32>(UInt32($0 >> 32),UInt32(truncatingIfNeeded:$0)) }
        var adjacency = [[UInt32]](repeating: [], count: n)
        for q in pairs { adjacency[Int(q.x)].append(q.y); adjacency[Int(q.y)].append(q.x) }
        var components = [UInt32](repeating: .max, count: n)
        var parents = [UInt32](repeating: .max, count: n)
        var members: [UInt32] = [], offsets: [UInt32] = [0]
        var count: UInt32 = 0
        for root in 0..<n where system.particles[root].massDa > 0 && components[root] == .max {
            var queue: [UInt32] = [UInt32(root)], cursor = 0
            components[root] = count
            while cursor < queue.count {
                let current = queue[cursor]; cursor += 1
                for neighbor in adjacency[Int(current)].sorted() where components[Int(neighbor)] == .max {
                    components[Int(neighbor)] = count
                    parents[Int(neighbor)] = current
                    queue.append(neighbor)
                }
            }
            members.append(contentsOf: queue)
            offsets.append(UInt32(members.count)); count += 1
        }
        for site in system.linearVirtualSites ?? [] {
            let component = components[Int(site.parentParticles[0])]
            guard component != .max, site.parentParticles.allSatisfy({ components[Int($0)] == component }) else {
                throw VivoMDRuntimeError.metal("virtual-site parents have inconsistent molecular ownership")
            }
            components[Int(site.siteParticle)] = component
        }
        guard count > 0, !components.contains(.max) else {
            throw VivoMDRuntimeError.unsupported(["NPT contains an unresolved massless component"])
        }
        var perComponent = [[SIMD2<UInt32>]](repeating: [], count: Int(count))
        for pair in pairs { perComponent[Int(components[Int(pair.x)])].append(pair) }
        var edges: [SIMD2<UInt32>] = [], edgeOffsets: [UInt32] = [0]
        for group in perComponent { edges.append(contentsOf: group); edgeOffsets.append(UInt32(edges.count)) }
        return .init(componentCount: count, componentOffsets: offsets, componentParticles: members,
                     componentIndexByParticle: components, parentByParticle: parents,
                     edgeOffsets: edgeOffsets, edges: edges)
    }
}

private struct VivoMDBarostatMetalCommand {
    var particleCount: UInt32
    var componentCount: UInt32
    var scale: Float
    var toleranceNM: Float
    var cellA: SIMD4<Float>
    var cellB: SIMD4<Float>
    var cellC: SIMD4<Float>
    var reciprocalA: SIMD4<Float>
    var reciprocalB: SIMD4<Float>
    var reciprocalC: SIMD4<Float>
}

final class VivoMDBarostatEngine: @unchecked Sendable {
    let plan: VivoMDBarostatPlan
    private let offsets, members, componentIndices, parents, edgeOffsets, edges, centers, unwrapped: MTLBuffer
    private let centerPipeline, scalePipeline: NumiVivoPipeline

    static func make(device: MTLDevice, catalog: NumiVivoPipelineCatalog,
                     system: VivoClassicalSystem) async throws -> VivoMDBarostatEngine {
        let plan = try VivoMDBarostatPlan.make(system: system)
        func immutable<T>(_ values: [T], _ label: String) throws -> MTLBuffer {
            let length = max(values.count * MemoryLayout<T>.stride,16)
            guard length <= device.maxBufferLength, let buffer = device.makeBuffer(length: length,options:.storageModeShared) else {
                throw VivoMDRuntimeError.metal("NPT \(label) allocation failed")
            }
            buffer.label = "NumiVivo.NPT." + label
            if !values.isEmpty { values.withUnsafeBytes { buffer.contents().copyMemory(from:$0.baseAddress!,byteCount:$0.count) } }
            return buffer
        }
        func scratch(_ count: Int, _ label: String) throws -> MTLBuffer {
            guard count <= device.maxBufferLength / 16, let buffer = device.makeBuffer(length:max(16,count*16),options:.storageModePrivate) else {
                throw VivoMDRuntimeError.metal("NPT \(label) allocation failed")
            }
            buffer.label = "NumiVivo.NPT." + label
            return buffer
        }
        let centers = try scratch(Int(plan.componentCount),"centers")
        let unwrapped = try scratch(system.particles.count,"unwrapped")
        let centerPipeline = try await catalog.pipeline(.mdBarostatCenters)
        let scalePipeline = try await catalog.pipeline(.mdBarostatScale)
        return try .init(plan:plan,offsets:immutable(plan.componentOffsets,"offsets"),
                         members:immutable(plan.componentParticles,"members"),
                         componentIndices:immutable(plan.componentIndexByParticle,"componentIndices"),
                         parents:immutable(plan.parentByParticle,"parents"),
                         edgeOffsets:immutable(plan.edgeOffsets,"edgeOffsets"),edges:immutable(plan.edges,"edges"),
                         centers:centers,unwrapped:unwrapped,centerPipeline:centerPipeline,scalePipeline:scalePipeline)
    }
    private init(plan:VivoMDBarostatPlan,offsets:MTLBuffer,members:MTLBuffer,componentIndices:MTLBuffer,
                 parents:MTLBuffer,edgeOffsets:MTLBuffer,edges:MTLBuffer,centers:MTLBuffer,unwrapped:MTLBuffer,
                 centerPipeline:NumiVivoPipeline,scalePipeline:NumiVivoPipeline) {
        self.plan=plan;self.offsets=offsets;self.members=members;self.componentIndices=componentIndices
        self.parents=parents;self.edgeOffsets=edgeOffsets;self.edges=edges;self.centers=centers;self.unwrapped=unwrapped
        self.centerPipeline=centerPipeline;self.scalePipeline=scalePipeline
    }
    func encodeProposal(commandBuffer:MTLCommandBuffer,positions:MTLBuffer,dynamics:MTLBuffer,
                        cell:VivoPeriodicCell,scale:Double,status:MTLBuffer) throws {
        guard scale.isFinite,scale>0,Float(scale).isFinite,Float(scale)>0 else { throw VivoMDRuntimeError.metal("invalid NPT scale") }
        let d=cell.a.dot(cell.b.cross(cell.c))
        guard d.isFinite,abs(d)>1e-18 else { throw VivoMDRuntimeError.metal("singular NPT cell") }
        func f(_ v:VivoVector3D)->SIMD4<Float>{ .init(Float(v.x),Float(v.y),Float(v.z),0) }
        var abi=VivoMDBarostatMetalCommand(particleCount:UInt32(plan.componentIndexByParticle.count),
            componentCount:plan.componentCount,scale:Float(scale),toleranceNM:1e-4,
            cellA:f(cell.a),cellB:f(cell.b),cellC:f(cell.c),
            reciprocalA:f(cell.b.cross(cell.c)/d),reciprocalB:f(cell.c.cross(cell.a)/d),reciprocalC:f(cell.a.cross(cell.b)/d))
        guard MemoryLayout<VivoMDBarostatMetalCommand>.stride==112 else { throw VivoMDRuntimeError.metal("NPT ABI mismatch") }
        try encode(centerPipeline,command:commandBuffer,buffers:[positions,dynamics,offsets,members,parents,edgeOffsets,edges,unwrapped,centers,status],
                   abi:&abi,count:Int(plan.componentCount))
        try encode(scalePipeline,command:commandBuffer,buffers:[positions,dynamics,componentIndices,centers,unwrapped,status],
                   abi:&abi,count:Int(abi.particleCount))
    }
    private func encode(_ pipeline:NumiVivoPipeline,command:MTLCommandBuffer,buffers:[MTLBuffer],
                        abi:inout VivoMDBarostatMetalCommand,count:Int) throws {
        guard let encoder=command.makeComputeCommandEncoder() else { throw VivoMDRuntimeError.metal("NPT encoder unavailable") }
        defer { encoder.endEncoding() }
        encoder.setComputePipelineState(pipeline.state)
        for (i,buffer) in buffers.enumerated(){encoder.setBuffer(buffer,offset:0,index:i)}
        encoder.setBytes(&abi,length:MemoryLayout<VivoMDBarostatMetalCommand>.stride,index:buffers.count)
        encoder.dispatchThreads(pipeline.gridSize(for:count),threadsPerThreadgroup:pipeline.threadgroupSize(for:count))
    }
}
