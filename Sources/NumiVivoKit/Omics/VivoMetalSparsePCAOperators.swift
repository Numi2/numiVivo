import Foundation
#if canImport(Metal)
@preconcurrency import Metal
#endif

/// Sparse FP32 Metal operators for the resident Krylov PCA owner.
///
/// Projection uses the source CSR stream and transpose uses a feature-major
/// CSC stream. Both streams retain only nonzero entries; no cells-by-features
/// dense matrix or covariance matrix is materialized.
final class VivoMetalSparsePCAOperators {
    static let method = "metal-FP32-sparse-PCA-operators-v1"
    static let maximumEntries = 16_777_216
    static let source = """
    #include <metal_stdlib>
    using namespace metal;
    kernel void vivo_sparse_pca_project(
        device const uint *row_offsets [[buffer(0)]],
        device const uint *columns [[buffer(1)]],
        device const float *values [[buffer(2)]],
        device const float *vector [[buffer(3)]],
        device float *output [[buffer(4)]],
        constant float &shift [[buffer(5)]],
        uint row [[thread_position_in_grid]]) {
        uint start = row_offsets[row];
        uint end = row_offsets[row + 1];
        float sum = 0.0f;
        for (uint i = start; i < end; ++i) {
            sum += values[i] * vector[columns[i]];
        }
        output[row] = sum - shift;
    }
    kernel void vivo_sparse_pca_transpose(
        device const uint *column_offsets [[buffer(0)]],
        device const uint *rows [[buffer(1)]],
        device const float *values [[buffer(2)]],
        device const float *vector [[buffer(3)]],
        device const float *initial [[buffer(4)]],
        device float *output [[buffer(5)]],
        uint column [[thread_position_in_grid]]) {
        uint start = column_offsets[column];
        uint end = column_offsets[column + 1];
        float sum = initial[column];
        for (uint i = start; i < end; ++i) {
            sum += values[i] * vector[rows[i]];
        }
        output[column] = sum;
    }
    """

    #if canImport(Metal)
    private let queue: MTLCommandQueue
    private let projectPipeline: MTLComputePipelineState
    private let transposePipeline: MTLComputePipelineState
    private let rowCount: Int
    private let columnCount: Int
    private let vectorBuffer: MTLBuffer
    private let initialBuffer: MTLBuffer
    private let outputBuffer: MTLBuffer
    #endif

    init(rowOffsets: [Int], columns: [Int], values: [Double], rows: Int, columnCount: Int) throws {
        guard rows > 0, columnCount > 0, rowOffsets.count == rows + 1,
              rowOffsets.first == 0, rowOffsets.last == values.count,
              columns.count == values.count, values.count > 0,
              values.count <= Self.maximumEntries,
              rows <= Int(UInt32.max), columnCount <= Int(UInt32.max),
              rowOffsets.allSatisfy({ $0 >= 0 && $0 <= values.count }),
              zip(rowOffsets, rowOffsets.dropFirst()).allSatisfy({ $0 <= $1 }) else {
            throw VivoOmicsError.invalid("Metal sparse PCA CSR dimensions")
        }
        #if canImport(Metal)
        try Task.checkCancellation()
        guard let device = MTLCreateSystemDefaultDevice(), device.hasUnifiedMemory,
              device.supportsFamily(.apple7), !device.name.lowercased().contains("paravirtual"),
              let queue = device.makeCommandQueue() else {
            throw VivoOmicsError.invalid("Metal sparse PCA operators require a physical Apple GPU")
        }
        var rowOffsets32 = [UInt32](repeating: 0, count: rowOffsets.count)
        for i in rowOffsets.indices { rowOffsets32[i] = UInt32(rowOffsets[i]) }
        var columns32 = [UInt32](repeating: 0, count: columns.count)
        var values32 = [Float](repeating: 0, count: values.count)
        var columnCounts = [Int](repeating: 0, count: columnCount)
        for (index, pair) in zip(columns, values).enumerated() {
            if index % 65_536 == 0 { try Task.checkCancellation() }
            let (column, value) = pair
            guard (0..<columnCount).contains(column), value.isFinite, value >= 0,
                  value <= 80, Float(value).isFinite else {
                throw VivoOmicsError.invalid("Metal sparse PCA coordinate or value")
            }
            columns32[index] = UInt32(column)
            values32[index] = Float(value)
            columnCounts[column] += 1
        }
        var columnOffsets32 = [UInt32](repeating: 0, count: columnCount + 1)
        for column in 0..<columnCount {
            let next = columnOffsets32[column] + UInt32(columnCounts[column])
            guard next >= columnOffsets32[column] else { throw VivoOmicsError.limit("Metal sparse PCA CSC offset overflow") }
            columnOffsets32[column + 1] = next
        }
        var cscRows = [UInt32](repeating: 0, count: values.count)
        var cscValues = [Float](repeating: 0, count: values.count)
        var cursor = columnOffsets32.dropLast().map { Int($0) }
        for row in 0..<rows {
            if row % 4_096 == 0 { try Task.checkCancellation() }
            for index in rowOffsets[row]..<rowOffsets[row + 1] {
                let column = Int(columns32[index]), position = cursor[column]
                cscRows[position] = UInt32(row)
                cscValues[position] = values32[index]
                cursor[column] += 1
            }
        }
        let options = MTLCompileOptions()
        options.fastMathEnabled = false
        let library = try device.makeLibrary(source: Self.source, options: options)
        guard let projectFunction = library.makeFunction(name: "vivo_sparse_pca_project"),
              let transposeFunction = library.makeFunction(name: "vivo_sparse_pca_transpose") else {
            throw VivoOmicsError.invalid("missing Metal sparse PCA operator")
        }
        self.projectPipeline = try device.makeComputePipelineState(function: projectFunction)
        self.transposePipeline = try device.makeComputePipelineState(function: transposeFunction)
        let vectorBytes = max(rows, columnCount) * MemoryLayout<Float>.stride
        let outputBytes = vectorBytes
        guard let rowOffsetsBuffer = device.makeBuffer(bytes: rowOffsets32, length: rowOffsets32.count * 4, options: .storageModeShared),
              let columnsBuffer = device.makeBuffer(bytes: columns32, length: max(4, columns32.count * 4), options: .storageModeShared),
              let valuesBuffer = device.makeBuffer(bytes: values32, length: max(4, values32.count * 4), options: .storageModeShared),
              let columnOffsetsBuffer = device.makeBuffer(bytes: columnOffsets32, length: columnOffsets32.count * 4, options: .storageModeShared),
              let cscRowsBuffer = device.makeBuffer(bytes: cscRows, length: max(4, cscRows.count * 4), options: .storageModeShared),
              let cscValuesBuffer = device.makeBuffer(bytes: cscValues, length: max(4, cscValues.count * 4), options: .storageModeShared),
              let vectorBuffer = device.makeBuffer(length: vectorBytes, options: .storageModeShared),
              let initialBuffer = device.makeBuffer(length: max(4, columnCount * 4), options: .storageModeShared),
              let outputBuffer = device.makeBuffer(length: outputBytes, options: .storageModeShared) else {
            throw VivoOmicsError.limit("Metal sparse PCA buffers")
        }
        self.queue = queue
        self.rowCount = rows
        self.columnCount = columnCount
        self.vectorBuffer = vectorBuffer
        self.initialBuffer = initialBuffer
        self.outputBuffer = outputBuffer
        self.rowOffsetsBuffer = rowOffsetsBuffer
        self.columnsBuffer = columnsBuffer
        self.valuesBuffer = valuesBuffer
        self.columnOffsetsBuffer = columnOffsetsBuffer
        self.cscRowsBuffer = cscRowsBuffer
        self.cscValuesBuffer = cscValuesBuffer
        #else
        throw VivoOmicsError.invalid("Metal sparse PCA operators are unavailable")
        #endif
    }

    #if canImport(Metal)
    private let rowOffsetsBuffer: MTLBuffer
    private let columnsBuffer: MTLBuffer
    private let valuesBuffer: MTLBuffer
    private let columnOffsetsBuffer: MTLBuffer
    private let cscRowsBuffer: MTLBuffer
    private let cscValuesBuffer: MTLBuffer

    private func copy(_ values: [Double], to buffer: MTLBuffer, maximumMagnitude: Double = Double(Float.greatestFiniteMagnitude)) throws {
        guard values.allSatisfy({ $0.isFinite && abs($0) <= maximumMagnitude }) else {
            throw VivoOmicsError.invalid("Metal sparse PCA vector")
        }
        var fp32 = values.map(Float.init)
        guard fp32.allSatisfy(\.isFinite) else { throw VivoOmicsError.invalid("Metal sparse PCA FP32 vector") }
        fp32.withUnsafeBytes { bytes in
            buffer.contents().copyMemory(from: bytes.baseAddress!, byteCount: bytes.count)
        }
    }

    private func command(_ pipeline: MTLComputePipelineState) throws -> (MTLCommandBuffer, MTLComputeCommandEncoder) {
        guard let command = queue.makeCommandBuffer(), let encoder = command.makeComputeCommandEncoder() else {
            throw VivoOmicsError.limit("Metal sparse PCA command")
        }
        encoder.setComputePipelineState(pipeline)
        return (command, encoder)
    }

    func project(_ vector: [Double], shift: Double) throws -> [Double] {
        guard vector.count == columnCount, shift.isFinite, abs(shift) <= Double(Float.greatestFiniteMagnitude) else {
            throw VivoOmicsError.invalid("Metal sparse PCA projection dimensions")
        }
        try Task.checkCancellation(); try copy(vector, to: vectorBuffer)
        var scalar = Float(shift)
        let (command, encoder) = try command(projectPipeline)
        encoder.setBuffer(rowOffsetsBuffer, offset: 0, index: 0)
        encoder.setBuffer(columnsBuffer, offset: 0, index: 1)
        encoder.setBuffer(valuesBuffer, offset: 0, index: 2)
        encoder.setBuffer(vectorBuffer, offset: 0, index: 3)
        encoder.setBuffer(outputBuffer, offset: 0, index: 4)
        encoder.setBytes(&scalar, length: MemoryLayout<Float>.size, index: 5)
        let width = min(projectPipeline.threadExecutionWidth, projectPipeline.maxTotalThreadsPerThreadgroup)
        encoder.dispatchThreads(MTLSize(width: rowCount, height: 1, depth: 1),
                                threadsPerThreadgroup: MTLSize(width: width, height: 1, depth: 1))
        encoder.endEncoding(); command.commit(); command.waitUntilCompleted()
        guard command.status == .completed else {
            throw VivoOmicsError.invalid("Metal sparse PCA projection failed: " + (command.error?.localizedDescription ?? "unknown"))
        }
        try Task.checkCancellation()
        let output = outputBuffer.contents().bindMemory(to: Float.self, capacity: rowCount)
        var result = [Double](repeating: 0, count: rowCount)
        for row in 0..<rowCount {
            result[row] = Double(output[row])
            guard result[row].isFinite else { throw VivoOmicsError.invalid("nonfinite Metal sparse PCA projection") }
        }
        return result
    }

    func transpose(_ vector: [Double], initial: [Double]) throws -> [Double] {
        guard vector.count == rowCount, initial.count == columnCount else {
            throw VivoOmicsError.invalid("Metal sparse PCA transpose dimensions")
        }
        try Task.checkCancellation(); try copy(vector, to: vectorBuffer)
        try copy(initial, to: initialBuffer, maximumMagnitude: Double(Float.greatestFiniteMagnitude))
        let (command, encoder) = try command(transposePipeline)
        encoder.setBuffer(columnOffsetsBuffer, offset: 0, index: 0)
        encoder.setBuffer(cscRowsBuffer, offset: 0, index: 1)
        encoder.setBuffer(cscValuesBuffer, offset: 0, index: 2)
        encoder.setBuffer(vectorBuffer, offset: 0, index: 3)
        encoder.setBuffer(initialBuffer, offset: 0, index: 4)
        encoder.setBuffer(outputBuffer, offset: 0, index: 5)
        let width = min(transposePipeline.threadExecutionWidth, transposePipeline.maxTotalThreadsPerThreadgroup)
        encoder.dispatchThreads(MTLSize(width: columnCount, height: 1, depth: 1),
                                threadsPerThreadgroup: MTLSize(width: width, height: 1, depth: 1))
        encoder.endEncoding(); command.commit(); command.waitUntilCompleted()
        guard command.status == .completed else {
            throw VivoOmicsError.invalid("Metal sparse PCA transpose failed: " + (command.error?.localizedDescription ?? "unknown"))
        }
        try Task.checkCancellation()
        let output = outputBuffer.contents().bindMemory(to: Float.self, capacity: columnCount)
        var result = [Double](repeating: 0, count: columnCount)
        for column in 0..<columnCount {
            result[column] = Double(output[column])
            guard result[column].isFinite else { throw VivoOmicsError.invalid("nonfinite Metal sparse PCA transpose") }
        }
        return result
    }
    #else
    func project(_ vector: [Double], shift: Double) throws -> [Double] { throw VivoOmicsError.invalid("Metal sparse PCA operators are unavailable") }
    func transpose(_ vector: [Double], initial: [Double]) throws -> [Double] { throw VivoOmicsError.invalid("Metal sparse PCA operators are unavailable") }
    #endif
}
