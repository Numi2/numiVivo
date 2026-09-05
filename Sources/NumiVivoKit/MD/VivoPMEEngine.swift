import Foundation
import Metal
import NumiVivoShaders

private struct VivoPMEMetalCommand {
    var particleCount: UInt32
    var gridX: UInt32
    var gridY: UInt32
    var gridZ: UInt32
    var gridPointCount: UInt32
    var axis: UInt32
    var stage: UInt32
    var inverse: UInt32
    var betaPerNM: Float
    var volumeNM3: Float
    var coulombPrefactor: Float
    var inverseGridCount: Float
    var reciprocalA: SIMD4<Float>
    var reciprocalB: SIMD4<Float>
    var reciprocalC: SIMD4<Float>
}

/// Persistent GPU-resident particle-mesh Ewald reciprocal backend.
/// Two complex grids are ping-ponged through a radix-2 3-D FFT. No host
/// readback occurs in the force path.
final class VivoPMEEngine: @unchecked Sendable {
    let plan: VivoPMEPlan
    private let gridA: MTLBuffer
    private let gridB: MTLBuffer
    private let pipelines: [NumiVivoKernel: NumiVivoPipeline]
    private let baseCommand: VivoPMEMetalCommand

    static func make(device: MTLDevice,
                     catalog: NumiVivoPipelineCatalog,
                     particleCount: UInt32,
                     cell: VivoPeriodicCell,
                     configuration: VivoMDConfiguration) async throws -> VivoPMEEngine {
        guard configuration.electrostatics == .pme else {
            throw VivoMDRuntimeError.metal("PME engine requested for non-PME configuration")
        }
        let plan = try VivoPMEPlan.make(cell: cell,
                                        cutoffNM: configuration.cutoffNM,
                                        tolerance: configuration.resolvedPMETolerance,
                                        targetGridSpacingNM: configuration.resolvedPMEGridSpacingNM)
        guard plan.gridPointCount <= UInt64(UInt32.max),
              plan.cellVolumeNM3 <= Double(Float.greatestFiniteMagnitude),
              plan.ewaldBetaPerNM <= Double(Float.greatestFiniteMagnitude) else {
            throw VivoMDRuntimeError.metal("PME plan exceeds current UInt32/FP32 execution ABI")
        }
        let pointCount = UInt32(plan.gridPointCount)
        let byteCount64 = plan.gridPointCount * UInt64(MemoryLayout<SIMD2<Float>>.stride)
        guard byteCount64 <= UInt64(device.maxBufferLength), byteCount64 <= UInt64(Int.max),
              let a = device.makeBuffer(length: Int(byteCount64), options: .storageModePrivate),
              let b = device.makeBuffer(length: Int(byteCount64), options: .storageModePrivate) else {
            throw VivoMDRuntimeError.metal("PME complex-grid allocation exceeds Metal limits")
        }
        a.label = "NumiVivo.PME.gridA"
        b.label = "NumiVivo.PME.gridB"
        func f4(_ v: VivoVector3D) -> SIMD4<Float> {
            .init(Float(v.x), Float(v.y), Float(v.z), 0)
        }
        let command = VivoPMEMetalCommand(
            particleCount: particleCount,
            gridX: plan.gridX, gridY: plan.gridY, gridZ: plan.gridZ,
            gridPointCount: pointCount, axis: 0, stage: 0, inverse: 0,
            betaPerNM: Float(plan.ewaldBetaPerNM),
            volumeNM3: Float(plan.cellVolumeNM3),
            coulombPrefactor: Float(138.935456 / configuration.relativeDielectric),
            inverseGridCount: 1 / Float(pointCount),
            reciprocalA: f4(plan.reciprocalA),
            reciprocalB: f4(plan.reciprocalB),
            reciprocalC: f4(plan.reciprocalC)
        )
        guard MemoryLayout<VivoPMEMetalCommand>.stride == 96 else {
            throw VivoMDRuntimeError.metal("Swift/Metal PME command ABI mismatch")
        }
        let names: [NumiVivoKernel] = [
            .mdPMEClearGrid, .mdPMESpread, .mdPMEBitReverse, .mdPMEFFTStage,
            .mdPMEInfluence, .mdPMEScaleInverse, .mdPMEGather
        ]
        var pipelines: [NumiVivoKernel: NumiVivoPipeline] = [:]
        for name in names { pipelines[name] = try await catalog.pipeline(name) }
        return .init(plan: plan, gridA: a, gridB: b,
                     pipelines: pipelines, baseCommand: command)
    }

    private init(plan: VivoPMEPlan, gridA: MTLBuffer, gridB: MTLBuffer,
                 pipelines: [NumiVivoKernel: NumiVivoPipeline],
                 baseCommand: VivoPMEMetalCommand) {
        self.plan = plan; self.gridA = gridA; self.gridB = gridB
        self.pipelines = pipelines; self.baseCommand = baseCommand
    }

    /// Adds reciprocal-space electrostatic force and energy into an already
    /// populated per-particle force/energy buffer.
    func encodeReciprocal(commandBuffer: MTLCommandBuffer,
                          positions: MTLBuffer,
                          dynamics: MTLBuffer,
                          forceEnergy: MTLBuffer,
                          status: MTLBuffer) throws {
        var command = baseCommand
        try dispatch(.mdPMEClearGrid, commandBuffer: commandBuffer,
                     buffers: [gridA], command: &command,
                     elements: Int(command.gridPointCount))
        try dispatch(.mdPMESpread, commandBuffer: commandBuffer,
                     buffers: [positions, dynamics, gridA, status], command: &command,
                     elements: Int(command.particleCount))

        var current = gridA
        var scratch = gridB
        try fft3D(commandBuffer: commandBuffer, current: &current, scratch: &scratch,
                  inverse: false, command: &command)
        try dispatch(.mdPMEInfluence, commandBuffer: commandBuffer,
                     buffers: [current, scratch], command: &command,
                     elements: Int(command.gridPointCount))
        swap(&current, &scratch)
        try fft3D(commandBuffer: commandBuffer, current: &current, scratch: &scratch,
                  inverse: true, command: &command)
        try dispatch(.mdPMEScaleInverse, commandBuffer: commandBuffer,
                     buffers: [current], command: &command,
                     elements: Int(command.gridPointCount))
        try dispatch(.mdPMEGather, commandBuffer: commandBuffer,
                     buffers: [positions, dynamics, current, forceEnergy, status],
                     command: &command, elements: Int(command.particleCount))
    }

    private func fft3D(commandBuffer: MTLCommandBuffer,
                       current: inout MTLBuffer,
                       scratch: inout MTLBuffer,
                       inverse: Bool,
                       command: inout VivoPMEMetalCommand) throws {
        for axis in UInt32(0)...UInt32(2) {
            command.axis = axis
            command.inverse = inverse ? 1 : 0
            command.stage = 0
            try dispatch(.mdPMEBitReverse, commandBuffer: commandBuffer,
                         buffers: [current, scratch], command: &command,
                         elements: Int(command.gridPointCount))
            swap(&current, &scratch)
            let length = axis == 0 ? command.gridX : (axis == 1 ? command.gridY : command.gridZ)
            let stages = exactLog2(length)
            for stage in UInt32(0)..<stages {
                command.stage = stage
                try dispatch(.mdPMEFFTStage, commandBuffer: commandBuffer,
                             buffers: [current, scratch], command: &command,
                             elements: Int(command.gridPointCount / 2))
                swap(&current, &scratch)
            }
        }
    }

    private func dispatch(_ kernel: NumiVivoKernel,
                          commandBuffer: MTLCommandBuffer,
                          buffers: [MTLBuffer],
                          command: inout VivoPMEMetalCommand,
                          elements: Int) throws {
        guard let pipeline = pipelines[kernel],
              let encoder = commandBuffer.makeComputeCommandEncoder() else {
            throw VivoMDRuntimeError.metal("PME encoder unavailable for \(kernel.rawValue)")
        }
        encoder.label = kernel.rawValue
        encoder.setComputePipelineState(pipeline.state)
        for (index, buffer) in buffers.enumerated() { encoder.setBuffer(buffer, offset: 0, index: index) }
        encoder.setBytes(&command, length: MemoryLayout<VivoPMEMetalCommand>.stride,
                         index: buffers.count)
        encoder.dispatchThreads(pipeline.gridSize(for: elements),
                                threadsPerThreadgroup: pipeline.threadgroupSize(for: elements))
        encoder.endEncoding()
    }

    private func exactLog2(_ value: UInt32) -> UInt32 {
        precondition(value > 0 && (value & (value - 1)) == 0)
        return UInt32(31 - value.leadingZeroBitCount)
    }
}
