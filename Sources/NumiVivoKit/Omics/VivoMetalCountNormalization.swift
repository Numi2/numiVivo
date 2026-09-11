import Foundation
#if canImport(Metal)
@preconcurrency import Metal
#endif

/// Explicit FP32 accelerator for the already verified sparse count transform.
/// One owner/queue per synchronous call; buffers are reused only after completion.
/// The output container remains FP64, holding exactly widened FP32 results.
enum VivoMetalCountNormalization {
    static let batchEntries = 262_144
    static let source = """
    #include <metal_stdlib>
    using namespace metal;
    kernel void vivo_count_log1p(device const uint4 *records [[buffer(0)]],
        device const float *scales [[buffer(1)]], device float *output [[buffer(2)]],
        constant uint &count [[buffer(3)]], uint i [[thread_position_in_grid]]) {
        if (i >= count) return;
        uint4 r = records[i];
        ulong molecules = (ulong(r.w) << 32) | ulong(r.z);
        float x = float(molecules) * scales[r.x];
        // Compensate rounding of 1+x, including x below half an FP32 ulp at one.
        // Fast math is disabled; the result is an explicit FP32 approximation.
        float y = 1.0f + x;
        output[i] = y == 1.0f ? x : log(y) * (x / (y - 1.0f));
    }
    """

    static func run(records: VivoWindowedCountRecords, quality: VivoCountStoreQuality,
                    target: Double, writer: VivoCountRecordWriter) throws -> VivoCountStoreNormalizationExecution {
        #if canImport(Metal)
        guard target.isFinite, (1...1_000_000_000).contains(target),
              let device = MTLCreateSystemDefaultDevice(), device.hasUnifiedMemory,
              device.supportsFamily(.apple7), !device.name.lowercased().contains("paravirtual"),
              let queue = device.makeCommandQueue() else {
            throw VivoOmicsError.invalid("Metal count normalization requires a physical Apple GPU and target in 1...1e9")
        }
        let options = MTLCompileOptions()
        options.fastMathEnabled = false
        let library = try device.makeLibrary(source: source, options: options)
        guard let function = library.makeFunction(name: "vivo_count_log1p") else {
            throw VivoOmicsError.invalid("missing Metal count normalization kernel")
        }
        let pipeline = try device.makeComputePipelineState(function: function)
        let scales = quality.rowTotals.map { $0 == 0 ? Float(0) : Float(target / Double($0)) }
        guard zip(scales, quality.rowTotals).allSatisfy({ $0.0.isFinite && ($0.1 == 0 || $0.0 > 0) }),
              let input = device.makeBuffer(length: batchEntries * 16, options: .storageModeShared),
              let output = device.makeBuffer(length: batchEntries * 4, options: .storageModeShared),
              let scaleBuffer = device.makeBuffer(length: max(4, scales.count * 4), options: .storageModeShared) else {
            throw VivoOmicsError.limit("Metal normalization buffers or FP32 scale range")
        }
        if !scales.isEmpty {
            scales.withUnsafeBytes { scaleBuffer.contents().copyMemory(from: $0.baseAddress!, byteCount: $0.count) }
        }
        let inputWords = input.contents().bindMemory(to: UInt64.self, capacity: batchEntries * 2)
        let outputValues = output.contents().bindMemory(to: Float.self, capacity: batchEntries)
        for start in stride(from: 0, to: records.count, by: batchEntries) {
            try Task.checkCancellation()
            let count = min(batchEntries, records.count - start)
            for i in 0..<count {
                let record = try records.record(start + i)
                guard record.row < quality.rowTotals.count, record.feature < quality.featureTotals.count,
                      record.bits > 0, record.bits <= quality.rowTotals[record.row] else {
                    throw VivoOmicsError.invalid("invalid Metal mapped count coordinate or total")
                }
                inputWords[i * 2] = (UInt64(record.row) | (UInt64(record.feature) << 32)).littleEndian
                inputWords[i * 2 + 1] = record.bits.littleEndian
            }
            guard let command = queue.makeCommandBuffer(), let encoder = command.makeComputeCommandEncoder() else {
                throw VivoOmicsError.limit("Metal normalization command")
            }
            encoder.setComputePipelineState(pipeline)
            encoder.setBuffer(input, offset: 0, index: 0)
            encoder.setBuffer(scaleBuffer, offset: 0, index: 1)
            encoder.setBuffer(output, offset: 0, index: 2)
            var n = UInt32(count)
            encoder.setBytes(&n, length: MemoryLayout<UInt32>.size, index: 3)
            let width = min(256, pipeline.maxTotalThreadsPerThreadgroup)
            encoder.dispatchThreads(MTLSize(width: count, height: 1, depth: 1),
                                    threadsPerThreadgroup: MTLSize(width: width, height: 1, depth: 1))
            encoder.endEncoding()
            command.commit(); command.waitUntilCompleted()
            try Task.checkCancellation()
            guard command.status == .completed else {
                throw VivoOmicsError.invalid("Metal normalization failed: \(command.error?.localizedDescription ?? "unknown")")
            }
            for i in 0..<count {
                let value = Double(outputValues[i]), packed = UInt64(littleEndian: inputWords[i * 2])
                guard value.isFinite, value > 0 else { throw VivoOmicsError.invalid("invalid Metal normalized value") }
                try writer.append(row: Int(packed & 0xffff_ffff), feature: Int(packed >> 32), bits: value.bitPattern)
            }
        }
        return .init(backend: .metalFP32, numericalProfile: "fp64-row-scale-to-fp32-count-multiply-compensated-log1p/v1",
                     deviceName: device.name, registryID: device.registryID,
                     kernel: try VivoCanonicalJSON.fingerprint(Data(source.utf8)), batchEntries: batchEntries)
        #else
        throw VivoOmicsError.invalid("Metal count normalization is unavailable")
        #endif
    }
}
