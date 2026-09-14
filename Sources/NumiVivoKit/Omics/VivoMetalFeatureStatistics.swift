import Foundation
#if canImport(Metal)
@preconcurrency import Metal
#endif

/// FP32 Metal reduction for sparse, nonzero log-normalized expression.
///
/// The caller supplies a feature-major sparse stream. One GPU thread owns one
/// feature and reduces only that feature's nonzero entries; implicit zero cells
/// are accounted for by the CPU owner after the reduction. No cell-by-feature
/// dense matrix is allocated.
enum VivoMetalFeatureStatistics {
    static let method = "metal-FP32-sparse-feature-statistics-v1"
    static let maximumEntries = 16_777_216
    static let maximumFeatures = 200_000
    static let source = """
    #include <metal_stdlib>
    using namespace metal;
    kernel void vivo_sparse_feature_statistics(
        device const uint *offsets [[buffer(0)]],
        device const float *values [[buffer(1)]],
        device float *sums [[buffer(2)]],
        device float *squares [[buffer(3)]],
        uint feature [[thread_position_in_grid]]) {
        uint start = offsets[feature];
        uint end = offsets[feature + 1];
        float sum = 0.0f;
        float square = 0.0f;
        for (uint i = start; i < end; ++i) {
            float y = exp(values[i]);
            // exp(log1p(x)) - 1 loses tiny x in FP32. The normalized log
            // values are nonnegative, so the input itself is the stable
            // fallback when exp rounds to one.
            float x = y == 1.0f ? values[i] : y - 1.0f;
            sum += x;
            square += x * x;
        }
        sums[feature] = sum;
        squares[feature] = square;
    }
    """

    /// Returns nonzero counts, nonzero means and nonzero Welford m2 values.
    /// The sums are accumulated in FP32 on the GPU and widened before the
    /// implicit-zero correction. This is an opt-in numerical profile; it does
    /// not silently replace the exact CPU FP64 path.
    static func run(values: [Double], featureIndices: [Int], featureCount: Int) throws ->
        (seen: [Int], means: [Double], m2: [Double]) {
        guard featureCount > 0, featureCount <= maximumFeatures,
              values.count == featureIndices.count, !values.isEmpty,
              values.count <= maximumEntries else {
            throw VivoOmicsError.limit("Metal feature-statistics sparse dimensions")
        }
        #if canImport(Metal)
        try Task.checkCancellation()
        guard let device = MTLCreateSystemDefaultDevice(), device.hasUnifiedMemory,
              device.supportsFamily(.apple7), !device.name.lowercased().contains("paravirtual"),
              let queue = device.makeCommandQueue() else {
            throw VivoOmicsError.invalid("Metal feature statistics require a physical Apple GPU")
        }
        var counts = [Int](repeating: 0, count: featureCount)
        for (index, pair) in zip(featureIndices, values).enumerated() {
            if index % 65_536 == 0 { try Task.checkCancellation() }
            let (feature, value) = pair
            guard (0..<featureCount).contains(feature), value.isFinite, value >= 0,
                  value <= 80, exp(value).isFinite else {
                throw VivoOmicsError.invalid("Metal feature-statistics coordinate or value")
            }
            counts[feature] += 1
        }
        var offsets = [UInt32](repeating: 0, count: featureCount + 1)
        for feature in 0..<featureCount {
            let next = offsets[feature] + UInt32(counts[feature])
            guard next >= offsets[feature] else { throw VivoOmicsError.limit("Metal feature-statistics offset overflow") }
            offsets[feature + 1] = next
        }
        var grouped = [Float](repeating: 0, count: values.count)
        var cursor = offsets.dropLast().map { Int($0) }
        for (index, pair) in zip(featureIndices, values).enumerated() {
            if index % 65_536 == 0 { try Task.checkCancellation() }
            let (feature, value) = pair
            let position = cursor[feature]
            grouped[position] = Float(value)
            cursor[feature] += 1
        }
        let options = MTLCompileOptions()
        options.fastMathEnabled = false
        let library = try device.makeLibrary(source: source, options: options)
        guard let function = library.makeFunction(name: "vivo_sparse_feature_statistics") else {
            throw VivoOmicsError.invalid("missing Metal feature-statistics kernel")
        }
        let pipeline = try device.makeComputePipelineState(function: function)
        let offsetBytes = offsets.count * MemoryLayout<UInt32>.stride
        let valueBytes = grouped.count * MemoryLayout<Float>.stride
        let outputBytes = featureCount * MemoryLayout<Float>.stride
        guard let offsetBuffer = device.makeBuffer(bytes: offsets, length: offsetBytes, options: .storageModeShared),
              let valueBuffer = device.makeBuffer(bytes: grouped, length: valueBytes, options: .storageModeShared),
              let sumBuffer = device.makeBuffer(length: outputBytes, options: .storageModeShared),
              let squareBuffer = device.makeBuffer(length: outputBytes, options: .storageModeShared),
              let command = queue.makeCommandBuffer(),
              let encoder = command.makeComputeCommandEncoder() else {
            throw VivoOmicsError.limit("Metal feature-statistics buffers or command")
        }
        encoder.setComputePipelineState(pipeline)
        encoder.setBuffer(offsetBuffer, offset: 0, index: 0)
        encoder.setBuffer(valueBuffer, offset: 0, index: 1)
        encoder.setBuffer(sumBuffer, offset: 0, index: 2)
        encoder.setBuffer(squareBuffer, offset: 0, index: 3)
        let width = min(pipeline.threadExecutionWidth, pipeline.maxTotalThreadsPerThreadgroup)
        encoder.dispatchThreads(MTLSize(width: featureCount, height: 1, depth: 1),
                                threadsPerThreadgroup: MTLSize(width: width, height: 1, depth: 1))
        encoder.endEncoding()
        command.commit()
        command.waitUntilCompleted()
        guard command.status == .completed else {
            throw VivoOmicsError.invalid("Metal feature-statistics command failed: " + (command.error?.localizedDescription ?? "unknown"))
        }
        try Task.checkCancellation()
        let sums = sumBuffer.contents().bindMemory(to: Float.self, capacity: featureCount)
        let squares = squareBuffer.contents().bindMemory(to: Float.self, capacity: featureCount)
        var means = [Double](repeating: 0, count: featureCount)
        var m2 = [Double](repeating: 0, count: featureCount)
        for feature in 0..<featureCount where counts[feature] > 0 {
            let sum = Double(sums[feature]), square = Double(squares[feature]), count = Double(counts[feature])
            guard sum.isFinite, square.isFinite, sum >= 0, square >= 0 else {
                throw VivoOmicsError.invalid("nonfinite Metal feature-statistics output")
            }
            means[feature] = sum / count
            // Rounding can make the recovered centered sum very slightly
            // negative. Clamp only that FP32 reconstruction residue.
            let centered = square - sum * sum / count
            guard centered.isFinite, centered >= -max(1e-8, square * 1e-6) else {
                throw VivoOmicsError.invalid("Metal feature-statistics variance")
            }
            m2[feature] = max(0, centered)
        }
        return (counts, means, m2)
        #else
        throw VivoOmicsError.invalid("Metal feature statistics are unavailable")
        #endif
    }
}
