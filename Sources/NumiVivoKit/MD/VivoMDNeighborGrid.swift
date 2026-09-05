import Foundation
import Metal
import NumiVivoShaders

public struct VivoMDNeighborGridPlan: Sendable, Equatable, Codable {
    public var dimensions: SIMD3<UInt32>
    public var cellCount: UInt32
    public var cellCapacity: UInt32
    public var neighborCapacity: UInt32
    public var neighborRadiusNM: Double

    public static func make(particleCount: UInt32,
                            cell: VivoPeriodicCell,
                            neighborRadiusNM: Double,
                            neighborCapacity: UInt32) throws -> Self {
        guard particleCount > 0, cell.isValid,
              neighborRadiusNM.isFinite, neighborRadiusNM > 0,
              neighborCapacity > 0 else {
            throw VivoMDRuntimeError.metal("invalid spatial neighbor-grid request")
        }
        let volume = cell.a.dot(cell.b.cross(cell.c))
        let heights = [abs(volume) / cell.b.cross(cell.c).norm,
                       abs(volume) / cell.c.cross(cell.a).norm,
                       abs(volume) / cell.a.cross(cell.b).norm]
        guard heights.allSatisfy({ $0.isFinite && $0 > 2 * neighborRadiusNM }) else {
            throw VivoMDRuntimeError.unsupported([
                "periodic cell is too narrow for the requested Verlet radius"
            ])
        }
        let raw = heights.map { max(2, Int(floor($0 / neighborRadiusNM))) }
        guard raw.allSatisfy({ $0 <= Int(UInt32.max) }) else {
            throw VivoMDRuntimeError.metal("neighbor-grid dimension exceeds UInt32")
        }
        let dims = SIMD3<UInt32>(UInt32(raw[0]), UInt32(raw[1]), UInt32(raw[2]))
        let xy = UInt64(dims.x) * UInt64(dims.y)
        let cells64 = xy * UInt64(dims.z)
        guard cells64 > 0, cells64 <= UInt64(UInt32.max) else {
            throw VivoMDRuntimeError.metal("neighbor-grid cell count exceeds UInt32")
        }
        let cells = UInt32(cells64)
        let mean = max(1, (UInt64(particleCount) + cells64 - 1) / cells64)
        // Deliberately generous bounded occupancy. Overflow is certified on GPU,
        // so this affects memory/performance but never silently affects physics.
        let proposed = max(UInt64(64), mean * 8)
        let capacity = UInt32(min(UInt64(particleCount), min(proposed, UInt64(UInt32.max))))
        return .init(dimensions: dims, cellCount: cells,
                     cellCapacity: max(capacity, 1),
                     neighborCapacity: neighborCapacity,
                     neighborRadiusNM: neighborRadiusNM)
    }
}

struct VivoMDNeighborGridCommand {
    var particleCount: UInt32
    var cellsX: UInt32
    var cellsY: UInt32
    var cellsZ: UInt32
    var cellCount: UInt32
    var cellCapacity: UInt32
    var neighborCapacity: UInt32
    var periodic: UInt32
    var neighborRadiusNM: Float
    var reserved0: Float = 0
    var reserved1: Float = 0
    var reserved2: Float = 0
    var cellA: SIMD4<Float>
    var cellB: SIMD4<Float>
    var cellC: SIMD4<Float>
    var reciprocalA: SIMD4<Float>
    var reciprocalB: SIMD4<Float>
    var reciprocalC: SIMD4<Float>
}

/// GPU spatial binning backend for the existing bounded Verlet-list buffers.
/// The force kernel remains unchanged: this module only replaces list construction.
final class VivoMDNeighborGrid: @unchecked Sendable {
    let plan: VivoMDNeighborGridPlan
    let cellCounts: MTLBuffer
    let cellParticles: MTLBuffer
    private let clearPipeline: NumiVivoPipeline
    private let binPipeline: NumiVivoPipeline
    private let buildPipeline: NumiVivoPipeline
    private let command: VivoMDNeighborGridCommand

    static func make(device: MTLDevice,
                     catalog: NumiVivoPipelineCatalog,
                     particleCount: UInt32,
                     cell: VivoPeriodicCell,
                     neighborRadiusNM: Double,
                     neighborCapacity: UInt32) async throws -> VivoMDNeighborGrid {
        let plan = try VivoMDNeighborGridPlan.make(particleCount: particleCount,
                                                   cell: cell,
                                                   neighborRadiusNM: neighborRadiusNM,
                                                   neighborCapacity: neighborCapacity)
        let countBytes = Int(plan.cellCount) * MemoryLayout<UInt32>.stride
        let slots64 = UInt64(plan.cellCount) * UInt64(plan.cellCapacity)
        guard slots64 <= UInt64(Int.max / MemoryLayout<UInt32>.stride) else {
            throw VivoMDRuntimeError.metal("neighbor-grid slot count exceeds host address space")
        }
        let particleBytes = Int(slots64) * MemoryLayout<UInt32>.stride
        guard countBytes <= device.maxBufferLength, particleBytes <= device.maxBufferLength,
              let counts = device.makeBuffer(length: max(countBytes, 16), options: .storageModePrivate),
              let particles = device.makeBuffer(length: max(particleBytes, 16), options: .storageModePrivate) else {
            throw VivoMDRuntimeError.metal("neighbor-grid allocation exceeds Metal buffer limits")
        }
        counts.label = "NumiVivo.MD.gridCellCounts"
        particles.label = "NumiVivo.MD.gridCellParticles"
        let reciprocal = try reciprocalCell(cell)
        func f4(_ value: VivoVector3D) -> SIMD4<Float> {
            .init(Float(value.x), Float(value.y), Float(value.z), 0)
        }
        let abi = VivoMDNeighborGridCommand(
            particleCount: particleCount,
            cellsX: plan.dimensions.x, cellsY: plan.dimensions.y, cellsZ: plan.dimensions.z,
            cellCount: plan.cellCount, cellCapacity: plan.cellCapacity,
            neighborCapacity: plan.neighborCapacity, periodic: 1,
            neighborRadiusNM: Float(plan.neighborRadiusNM),
            cellA: f4(cell.a), cellB: f4(cell.b), cellC: f4(cell.c),
            reciprocalA: f4(reciprocal.0), reciprocalB: f4(reciprocal.1),
            reciprocalC: f4(reciprocal.2)
        )
        guard MemoryLayout<VivoMDNeighborGridCommand>.stride == 144 else {
            throw VivoMDRuntimeError.metal("Swift spatial-grid ABI stride is not 144 bytes")
        }
        return try await .init(plan: plan, cellCounts: counts, cellParticles: particles,
                               command: abi,
                               clearPipeline: catalog.pipeline(.mdGridClear),
                               binPipeline: catalog.pipeline(.mdGridBin),
                               buildPipeline: catalog.pipeline(.mdGridBuildNeighbors))
    }

    private init(plan: VivoMDNeighborGridPlan,
                 cellCounts: MTLBuffer,
                 cellParticles: MTLBuffer,
                 command: VivoMDNeighborGridCommand,
                 clearPipeline: NumiVivoPipeline,
                 binPipeline: NumiVivoPipeline,
                 buildPipeline: NumiVivoPipeline) {
        self.plan = plan; self.cellCounts = cellCounts; self.cellParticles = cellParticles
        self.command = command; self.clearPipeline = clearPipeline
        self.binPipeline = binPipeline; self.buildPipeline = buildPipeline
    }

    func encode(commandBuffer: MTLCommandBuffer,
                positions: MTLBuffer,
                neighborCounts: MTLBuffer,
                neighborIndices: MTLBuffer,
                referencePositions: MTLBuffer,
                status: MTLBuffer) throws {
        var abi = command
        try dispatch(clearPipeline, commandBuffer: commandBuffer,
                     buffers: [cellCounts, status], abi: &abi,
                     elements: Int(plan.cellCount), label: "NumiVivo.MD.grid.clear")
        try dispatch(binPipeline, commandBuffer: commandBuffer,
                     buffers: [positions, cellCounts, cellParticles, status], abi: &abi,
                     elements: Int(abi.particleCount), label: "NumiVivo.MD.grid.bin")
        try dispatch(buildPipeline, commandBuffer: commandBuffer,
                     buffers: [positions, cellCounts, cellParticles,
                               neighborCounts, neighborIndices, referencePositions, status],
                     abi: &abi, elements: Int(abi.particleCount),
                     label: "NumiVivo.MD.grid.neighbors")
    }

    private func dispatch(_ pipeline: NumiVivoPipeline,
                          commandBuffer: MTLCommandBuffer,
                          buffers: [MTLBuffer],
                          abi: inout VivoMDNeighborGridCommand,
                          elements: Int,
                          label: String) throws {
        guard let encoder = commandBuffer.makeComputeCommandEncoder() else {
            throw VivoMDRuntimeError.metal("cannot create \(label) encoder")
        }
        encoder.label = label
        encoder.setComputePipelineState(pipeline.state)
        for (index, buffer) in buffers.enumerated() { encoder.setBuffer(buffer, offset: 0, index: index) }
        encoder.setBytes(&abi, length: MemoryLayout<VivoMDNeighborGridCommand>.stride, index: buffers.count)
        encoder.dispatchThreads(pipeline.gridSize(for: elements),
                                threadsPerThreadgroup: pipeline.threadgroupSize(for: elements))
        encoder.endEncoding()
    }

    private static func reciprocalCell(_ cell: VivoPeriodicCell) throws
    -> (VivoVector3D, VivoVector3D, VivoVector3D) {
        let determinant = cell.a.dot(cell.b.cross(cell.c))
        guard determinant.isFinite, abs(determinant) > 1e-18 else {
            throw VivoMDRuntimeError.metal("neighbor-grid periodic cell is singular")
        }
        return (cell.b.cross(cell.c) / determinant,
                cell.c.cross(cell.a) / determinant,
                cell.a.cross(cell.b) / determinant)
    }
}
