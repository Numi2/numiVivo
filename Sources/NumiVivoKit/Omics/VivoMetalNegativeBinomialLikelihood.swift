import Foundation
#if canImport(Metal)
@preconcurrency import Metal
#endif

/// Bounded FP32 Metal evaluation of the mean-dependent part of an NB2 log
/// likelihood. Count-only log-gamma terms stay on the exact CPU owner, so this
/// operator is an explicit objective accelerator rather than a second
/// statistical authority.
final class VivoMetalNegativeBinomialLikelihood {
    static let method = "metal-FP32-batched-NB-objective-v1"
    static let maximumCount = 16_777_216
    static let maximumObservations = 4_000_000
    static let batchObservations = 262_144
    private let counts: [UInt32]
    static let source = """
    #include <metal_stdlib>
    using namespace metal;
    kernel void vivo_nb_variable_objective(
        device const uint *counts [[buffer(0)]],
        device const float *means [[buffer(1)]],
        device float *terms [[buffer(2)]],
        constant float &dispersion [[buffer(3)]],
        constant uint &count [[buffer(4)]],
        uint i [[thread_position_in_grid]]) {
        if (i >= count) return;
        float y = float(counts[i]);
        float mu = means[i];
        float r = 1.0f / dispersion;
        terms[i] = y * log(mu) - (y + r) * log(1.0f + dispersion * mu);
    }
    """

    #if canImport(Metal)
    private let queue: MTLCommandQueue
    private let pipeline: MTLComputePipelineState
    private let countBuffer: MTLBuffer
    private let meanBuffer: MTLBuffer
    private let termBuffer: MTLBuffer
    private let dispersion: Float
    #endif

    init(counts: [UInt64], dispersion: Double) throws {
        guard !counts.isEmpty, counts.count <= Self.maximumObservations,
              counts.allSatisfy({ $0 <= UInt64(Self.maximumCount) }),
              dispersion.isFinite, (1e-8...100).contains(dispersion),
              Float(dispersion).isFinite else {
            throw VivoOmicsStatisticsError.invalid("Metal NB objective count or dispersion domain")
        }
        self.counts = counts.map(UInt32.init)
        #if canImport(Metal)
        try Task.checkCancellation()
        guard let device = MTLCreateSystemDefaultDevice(), device.hasUnifiedMemory,
              device.supportsFamily(.apple7), !device.name.lowercased().contains("paravirtual"),
              let queue = device.makeCommandQueue() else {
            throw VivoOmicsStatisticsError.invalid("Metal NB objective requires a physical Apple GPU")
        }
        let options = MTLCompileOptions()
        options.fastMathEnabled = false
        let library = try device.makeLibrary(source: Self.source, options: options)
        guard let function = library.makeFunction(name: "vivo_nb_variable_objective") else {
            throw VivoOmicsStatisticsError.invalid("missing Metal NB objective kernel")
        }
        let pipeline = try device.makeComputePipelineState(function: function)
        guard let countBuffer = device.makeBuffer(length: Self.batchObservations * MemoryLayout<UInt32>.stride,
                                                   options: .storageModeShared),
              let meanBuffer = device.makeBuffer(length: Self.batchObservations * MemoryLayout<Float>.stride,
                                                  options: .storageModeShared),
              let termBuffer = device.makeBuffer(length: Self.batchObservations * MemoryLayout<Float>.stride,
                                                  options: .storageModeShared) else {
            throw VivoOmicsStatisticsError.invalid("Metal NB objective buffers")
        }
        self.queue = queue
        self.pipeline = pipeline
        self.countBuffer = countBuffer
        self.meanBuffer = meanBuffer
        self.termBuffer = termBuffer
        self.dispersion = Float(dispersion)
        #else
        throw VivoOmicsStatisticsError.invalid("Metal NB objective is unavailable")
        #endif
    }

    /// Evaluates only y*log(mu) - (y + 1/a)*log(1+a*mu), which is the part
    /// that changes during coefficient line search. Results are widened and
    /// compensated on the CPU after each bounded dispatch.
    func evaluate(means: [Double]) throws -> Double {
        guard means.count == counts.count else {
            throw VivoOmicsStatisticsError.invalid("Metal NB objective mean dimensions")
        }
        #if canImport(Metal)
        try Task.checkCancellation()
        let fp32Means = try means.map { value -> Float in
            guard value.isFinite, value > 0, value >= Double(Float.leastNormalMagnitude),
                  value <= Double(Float.greatestFiniteMagnitude) else {
                throw VivoOmicsStatisticsError.invalid("Metal NB objective mean outside FP32 domain")
            }
            let converted = Float(value)
            guard converted.isFinite, converted > 0,
                  (converted * dispersion).isFinite else {
                throw VivoOmicsStatisticsError.invalid("Metal NB objective mean or product outside FP32 domain")
            }
            return converted
        }
        let countPointer = countBuffer.contents().bindMemory(to: UInt32.self, capacity: Self.batchObservations)
        let meanPointer = meanBuffer.contents().bindMemory(to: Float.self, capacity: Self.batchObservations)
        let termPointer = termBuffer.contents().bindMemory(to: Float.self, capacity: Self.batchObservations)
        var total = 0.0
        var compensation = 0.0
        for start in stride(from: 0, to: counts.count, by: Self.batchObservations) {
            try Task.checkCancellation()
            let count = min(Self.batchObservations, counts.count - start)
            counts[start..<start + count].withUnsafeBufferPointer { values in
                countPointer.update(from: values.baseAddress!, count: count)
            }
            fp32Means[start..<start + count].withUnsafeBufferPointer { values in
                meanPointer.update(from: values.baseAddress!, count: count)
            }
            guard let command = queue.makeCommandBuffer(), let encoder = command.makeComputeCommandEncoder() else {
                throw VivoOmicsStatisticsError.invalid("Metal NB objective command")
            }
            encoder.setComputePipelineState(pipeline)
            encoder.setBuffer(countBuffer, offset: 0, index: 0)
            encoder.setBuffer(meanBuffer, offset: 0, index: 1)
            encoder.setBuffer(termBuffer, offset: 0, index: 2)
            var a = dispersion, n = UInt32(count)
            encoder.setBytes(&a, length: MemoryLayout<Float>.size, index: 3)
            encoder.setBytes(&n, length: MemoryLayout<UInt32>.size, index: 4)
            let width = min(pipeline.threadExecutionWidth, pipeline.maxTotalThreadsPerThreadgroup)
            encoder.dispatchThreads(MTLSize(width: count, height: 1, depth: 1),
                                    threadsPerThreadgroup: MTLSize(width: width, height: 1, depth: 1))
            encoder.endEncoding()
            command.commit()
            command.waitUntilCompleted()
            guard command.status == .completed else {
                throw VivoOmicsStatisticsError.invalid("Metal NB objective failed: " +
                    (command.error?.localizedDescription ?? "unknown"))
            }
            for i in 0..<count {
                let term = Double(termPointer[i])
                guard term.isFinite else { throw VivoOmicsStatisticsError.invalid("nonfinite Metal NB objective term") }
                let corrected = term - compensation
                let next = total + corrected
                compensation = (next - total) - corrected
                total = next
            }
        }
        guard total.isFinite else { throw VivoOmicsStatisticsError.invalid("nonfinite Metal NB objective") }
        return total
        #else
        throw VivoOmicsStatisticsError.invalid("Metal NB objective is unavailable")
        #endif
    }
}
