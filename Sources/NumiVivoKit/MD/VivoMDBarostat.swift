import Foundation
import Metal
import NumiVivoShaders

public struct VivoMDBarostatMoveCertificate: Codable, Sendable, Equatable {
    public var stepIndex: UInt64
    public var attempted: Bool
    public var accepted: Bool
    public var volumeBeforeNM3: Double
    public var volumeAfterNM3: Double
    public var logVolumeDelta: Double
    public var logAcceptanceRatio: Double
}

public struct VivoMDBarostatPlan: Sendable, Equatable {
    public var componentCount: UInt32
    public var componentOffsets: [UInt32]
    public var componentParticles: [UInt32]
    public var componentIndexByParticle: [UInt32]

    public static func make(system: VivoClassicalSystem) throws -> Self {
        try VivoClassicalSystemValidator.validate(system)
        let n = system.particles.count
        guard n > 0, n <= Int(UInt32.max) else {
            throw VivoArtifactValidationError.invalid("barostat particle count exceeds UInt32")
        }
        var parent = Array(0..<n)
        func root(_ start: Int, _ parent: inout [Int]) -> Int {
            var x = start
            while parent[x] != x { parent[x] = parent[parent[x]]; x = parent[x] }
            return x
        }
        func unite(_ a: Int, _ b: Int, _ parent: inout [Int]) {
            let ra = root(a, &parent), rb = root(b, &parent)
            if ra == rb { return }
            if ra < rb { parent[rb] = ra } else { parent[ra] = rb }
        }
        for bond in system.bonds {
            let a = Int(bond.a), b = Int(bond.b)
            if system.particles[a].massDa > 0, system.particles[b].massDa > 0 {
                unite(a, b, &parent)
            }
        }
        for i in 0..<n where system.particles[i].massDa > 0 { _ = root(i, &parent) }

        var virtualParent: [UInt32: UInt32] = [:]
        for site in system.linearVirtualSites ?? [] {
            guard let first = site.parentParticles.first else {
                throw VivoArtifactValidationError.invalid("virtual site has no parent for barostat component assignment")
            }
            virtualParent[site.siteParticle] = first
        }
        var representative = [Int](repeating: -1, count: n)
        for i in 0..<n {
            if system.particles[i].massDa > 0 {
                representative[i] = root(i, &parent)
            } else if system.particles[i].role == .virtualSite,
                      let p = virtualParent[UInt32(i)] {
                representative[i] = root(Int(p), &parent)
            } else {
                representative[i] = i
            }
        }
        let roots = Array(Set(representative)).sorted()
        guard roots.count <= Int(UInt32.max) else {
            throw VivoArtifactValidationError.invalid("barostat component count exceeds UInt32")
        }
        let componentForRoot = Dictionary(uniqueKeysWithValues: roots.enumerated().map { ($0.element, UInt32($0.offset)) })
        var byComponent = [[UInt32]](repeating: [], count: roots.count)
        var indexByParticle = [UInt32](repeating: 0, count: n)
        for i in 0..<n {
            guard let component = componentForRoot[representative[i]] else {
                throw VivoArtifactValidationError.invalid("barostat component construction failed")
            }
            indexByParticle[i] = component
            byComponent[Int(component)].append(UInt32(i))
        }
        var offsets: [UInt32] = [0]
        var members: [UInt32] = []
        for component in byComponent {
            members.append(contentsOf: component)
            guard let next = UInt32(exactly: members.count) else {
                throw VivoArtifactValidationError.invalid("barostat component membership exceeds UInt32")
            }
            offsets.append(next)
        }
        return .init(componentCount: UInt32(roots.count), componentOffsets: offsets,
                     componentParticles: members, componentIndexByParticle: indexByParticle)
    }
}

private struct VivoMDBarostatMetalCommand {
    var particleCount: UInt32
    var componentCount: UInt32
    var scale: Float
    var reserved: Float = 0
    var cellA: SIMD4<Float>
    var cellB: SIMD4<Float>
    var cellC: SIMD4<Float>
    var reciprocalA: SIMD4<Float>
    var reciprocalB: SIMD4<Float>
    var reciprocalC: SIMD4<Float>
}

final class VivoMDBarostatEngine: @unchecked Sendable {
    let plan: VivoMDBarostatPlan
    private let componentOffsets: MTLBuffer
    private let componentParticles: MTLBuffer
    private let componentIndexByParticle: MTLBuffer
    private let centers: MTLBuffer
    private let centerPipeline: NumiVivoPipeline
    private let scalePipeline: NumiVivoPipeline

    static func make(device: MTLDevice,
                     catalog: NumiVivoPipelineCatalog,
                     system: VivoClassicalSystem) async throws -> VivoMDBarostatEngine {
        let plan = try VivoMDBarostatPlan.make(system: system)
        func immutable<T>(_ values: [T], _ label: String) throws -> MTLBuffer {
            let bytes = max(values.count * MemoryLayout<T>.stride, 16)
            guard bytes <= device.maxBufferLength else {
                throw VivoMDRuntimeError.metal("barostat \(label) exceeds Metal buffer limit")
            }
            if values.isEmpty {
                guard let b = device.makeBuffer(length: bytes, options: .storageModePrivate) else {
                    throw VivoMDRuntimeError.metal("barostat \(label) allocation failed")
                }
                b.label = label; return b
            }
            let buffer = values.withUnsafeBytes { raw -> MTLBuffer? in
                guard let base = raw.baseAddress else { return nil }
                return device.makeBuffer(bytes: base, length: raw.count, options: .storageModeShared)
            }
            guard let buffer else { throw VivoMDRuntimeError.metal("barostat \(label) upload failed") }
            buffer.label = label
            return buffer
        }
        let centerBytes = max(Int(plan.componentCount) * MemoryLayout<SIMD4<Float>>.stride, 16)
        guard centerBytes <= device.maxBufferLength,
              let centers = device.makeBuffer(length: centerBytes, options: .storageModePrivate) else {
            throw VivoMDRuntimeError.metal("barostat center buffer allocation failed")
        }
        centers.label = "NumiVivo.NPT.componentCenters"
        return try await .init(
            plan: plan,
            componentOffsets: immutable(plan.componentOffsets, "NumiVivo.NPT.componentOffsets"),
            componentParticles: immutable(plan.componentParticles, "NumiVivo.NPT.componentParticles"),
            componentIndexByParticle: immutable(plan.componentIndexByParticle, "NumiVivo.NPT.componentIndexByParticle"),
            centers: centers,
            centerPipeline: catalog.pipeline(.mdBarostatCenters),
            scalePipeline: catalog.pipeline(.mdBarostatScale)
        )
    }

    private init(plan: VivoMDBarostatPlan,
                 componentOffsets: MTLBuffer,
                 componentParticles: MTLBuffer,
                 componentIndexByParticle: MTLBuffer,
                 centers: MTLBuffer,
                 centerPipeline: NumiVivoPipeline,
                 scalePipeline: NumiVivoPipeline) {
        self.plan = plan; self.componentOffsets = componentOffsets
        self.componentParticles = componentParticles
        self.componentIndexByParticle = componentIndexByParticle
        self.centers = centers; self.centerPipeline = centerPipeline
        self.scalePipeline = scalePipeline
    }

    func encodeProposal(commandBuffer: MTLCommandBuffer,
                        positions: MTLBuffer,
                        dynamics: MTLBuffer,
                        cell: VivoPeriodicCell,
                        scale: Double) throws {
        guard scale.isFinite, scale > 0, scale <= Double(Float.greatestFiniteMagnitude) else {
            throw VivoMDRuntimeError.metal("barostat scale is invalid")
        }
        let d = cell.a.dot(cell.b.cross(cell.c))
        guard d.isFinite, abs(d) > 1e-18 else { throw VivoMDRuntimeError.metal("barostat cell is singular") }
        let reciprocal = (cell.b.cross(cell.c)/d, cell.c.cross(cell.a)/d, cell.a.cross(cell.b)/d)
        func f4(_ v: VivoVector3D) -> SIMD4<Float> { .init(Float(v.x),Float(v.y),Float(v.z),0) }
        var abi = VivoMDBarostatMetalCommand(particleCount: UInt32(plan.componentIndexByParticle.count),
                                             componentCount: plan.componentCount,
                                             scale: Float(scale),
                                             cellA: f4(cell.a), cellB: f4(cell.b), cellC: f4(cell.c),
                                             reciprocalA: f4(reciprocal.0), reciprocalB: f4(reciprocal.1), reciprocalC: f4(reciprocal.2))
        guard MemoryLayout<VivoMDBarostatMetalCommand>.stride == 112 else {
            throw VivoMDRuntimeError.metal("Swift/Metal barostat command ABI mismatch")
        }
        try encode(centerPipeline, commandBuffer: commandBuffer,
                   buffers: [positions,dynamics,componentOffsets,componentParticles,centers],
                   abi: &abi, elements: Int(plan.componentCount))
        try encode(scalePipeline, commandBuffer: commandBuffer,
                   buffers: [positions,componentIndexByParticle,centers],
                   abi: &abi, elements: plan.componentIndexByParticle.count)
    }

    private func encode(_ pipeline: NumiVivoPipeline,
                        commandBuffer: MTLCommandBuffer,
                        buffers: [MTLBuffer],
                        abi: inout VivoMDBarostatMetalCommand,
                        elements: Int) throws {
        guard let encoder = commandBuffer.makeComputeCommandEncoder() else {
            throw VivoMDRuntimeError.metal("barostat compute encoder unavailable")
        }
        encoder.setComputePipelineState(pipeline.state)
        for (index,buffer) in buffers.enumerated(){encoder.setBuffer(buffer,offset:0,index:index)}
        encoder.setBytes(&abi,length:MemoryLayout<VivoMDBarostatMetalCommand>.stride,index:buffers.count)
        encoder.dispatchThreads(pipeline.gridSize(for:elements),threadsPerThreadgroup:pipeline.threadgroupSize(for:elements))
        encoder.endEncoding()
    }
}
