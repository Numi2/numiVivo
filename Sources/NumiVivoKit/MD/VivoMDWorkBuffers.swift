import Foundation
import Metal

/// Reused for scalar reductions, constrained descent and drift corrections.
/// Scratch allocation is outside the microstep and is covered by one local
/// working-set preflight. This is not a replacement for a cross-runtime budget.
final class VivoMDWorkBuffers: @unchecked Sendable {
    let referencePosition: MTLBuffer
    let directionA: MTLBuffer
    let directionB: MTLBuffer
    let unitDynamics: MTLBuffer
    let termsA: MTLBuffer
    let termsB: MTLBuffer
    let scalar: MTLBuffer

    init(device: MTLDevice, particles: [VivoClassicalParticle]) throws {
        let count = particles.count
        guard count > 0, count <= Int(UInt32.max), count <= device.maxBufferLength / 16 else {
            throw VivoMDRuntimeError.metal("MD workspace shape exceeds device limits")
        }
        let vectorBytes = max(count * 16, 16), termBytes = max(count * 8, 16)
        let required = UInt64(vectorBytes) * 4 + UInt64(termBytes) * 2 + 16
        let budget = UInt64(Double(device.recommendedMaxWorkingSetSize) * 0.8)
        guard required <= budget, UInt64(device.currentAllocatedSize) <= budget - required else {
            throw VivoMDRuntimeError.metal("MD workspace exceeds remaining working-set headroom")
        }
        func buffer(_ bytes: Int, _ label: String, shared: Bool = false) throws -> MTLBuffer {
            guard let result = device.makeBuffer(length: bytes, options: shared ? .storageModeShared : .storageModePrivate) else {
                throw VivoMDRuntimeError.metal("cannot allocate \(label)")
            }
            result.label = "NumiVivo.MD." + label
            return result
        }
        referencePosition = try buffer(vectorBytes, "unconstrainedPosition")
        directionA = try buffer(vectorBytes, "projectedDirectionA")
        directionB = try buffer(vectorBytes, "projectedDirectionB")
        unitDynamics = try buffer(vectorBytes, "unitMetric", shared: true)
        let pointer = unitDynamics.contents().assumingMemoryBound(to: SIMD4<Float>.self)
        for (i, particle) in particles.enumerated() {
            pointer[i] = .init(0, particle.massDa > 0 ? 1 : 0, 0, 0)
        }
        termsA = try buffer(termBytes, "reductionA")
        termsB = try buffer(termBytes, "reductionB")
        scalar = try buffer(16, "reductionResult", shared: true)
    }
}
