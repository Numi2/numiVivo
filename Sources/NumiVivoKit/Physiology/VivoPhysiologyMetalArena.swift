import Foundation
@preconcurrency import Metal

final class VivoPhysiologyMetalArena: @unchecked Sendable {
    struct Capacities: Sendable {
        let stateElements: Int
        let preTransforms: Int
        let postTransforms: Int
        let publications: Int
    }

    let device: MTLDevice
    let heap: MTLHeap
    let capabilities: VivoMetalCapabilities
    let model: PreparedVivoPhysiologyModel
    let configuration: VivoPhysiologyRuntimeConfiguration
    let capacities: Capacities

    var currentState: MTLBuffer
    let baseState: MTLBuffer
    let stageState: MTLBuffer
    var candidateState: MTLBuffer
    let derivativeK1: MTLBuffer
    let incidenceOffsets: MTLBuffer
    let incidence: MTLBuffer
    let clearances: MTLBuffer
    let bounds: MTLBuffer

    let command: MTLBuffer
    let status: MTLBuffer
    let preTransforms: MTLBuffer
    let postTransforms: MTLBuffer
    let publicationRequests: MTLBuffer
    let publicationOutput: MTLBuffer
    let stateReadback: MTLBuffer

    private static let privateOptions: MTLResourceOptions = .storageModePrivate
    private static let sharedOptions: MTLResourceOptions = .storageModeShared

    init(
        device: MTLDevice,
        commandQueue: MTLCommandQueue,
        model: PreparedVivoPhysiologyModel,
        configuration: VivoPhysiologyRuntimeConfiguration
    ) throws {
        try VivoPhysiologyRuntimeCommandABI.validateMemoryLayout()
        try VivoPhysiologyRuntimeStatusABI.validateMemoryLayout()
        try configuration.validate(for: model)

        self.device = device
        self.model = model
        self.configuration = configuration
        self.capabilities = VivoMetalCapabilities(device: device)

        let stateElements = model.initialState.count
        self.capacities = Capacities(
            stateElements: stateElements,
            preTransforms: configuration.maximumTransformsPerStep,
            postTransforms: configuration.maximumTransformsPerStep,
            publications: configuration.maximumPublicationsPerStep
        )

        let stateBytes = try Self.checkedBytes(stateElements, MemoryLayout<Float>.stride, label: "physiology state")
        let offsetBytes = try Self.checkedBytes(model.incidenceOffsets.count, MemoryLayout<UInt32>.stride, label: "physiology incidence offsets")
        let incidenceBytes = try Self.checkedBytes(max(model.incidence.count, 1), MemoryLayout<VivoPhysiologyIncidenceABI>.stride, label: "physiology incidence")
        let clearanceBytes = try Self.checkedBytes(model.clearances.count, MemoryLayout<VivoPhysiologyClearanceABI>.stride, label: "physiology clearance")
        let boundsBytes = try Self.checkedBytes(Int(model.pairCount), MemoryLayout<SIMD2<Float>>.stride, label: "physiology bounds")

        let requirements: [(String, Int)] = [
            ("current-state", stateBytes),
            ("base-state", stateBytes),
            ("stage-state", stateBytes),
            ("candidate-state", stateBytes),
            ("derivative-k1", stateBytes),
            ("incidence-offsets", offsetBytes),
            ("incidence", incidenceBytes),
            ("clearances", clearanceBytes),
            ("bounds", boundsBytes)
        ]
        let requiredHeap = try Self.requiredHeapSize(
            device: device,
            requirements: requirements,
            options: Self.privateOptions
        )
        let heapBytes = try Self.applyHeadroom(requiredHeap, factor: configuration.privateHeapHeadroom)
        try capabilities.validate(programBytes: 0, privateHeapBytes: UInt64(heapBytes))

        let descriptor = MTLHeapDescriptor()
        descriptor.label = "NumiVivo.Physiology.PrivateArena"
        descriptor.storageMode = .private
        descriptor.cpuCacheMode = .defaultCache
        descriptor.hazardTrackingMode = .tracked
        descriptor.size = heapBytes
        guard let heap = device.makeHeap(descriptor: descriptor) else {
            throw VivoRuntimeError.allocationFailed("Metal could not allocate the physiology private heap")
        }
        self.heap = heap

        self.currentState = try Self.makePrivateBuffer(heap: heap, length: stateBytes, label: "NumiVivo.Physiology.State.Current")
        self.baseState = try Self.makePrivateBuffer(heap: heap, length: stateBytes, label: "NumiVivo.Physiology.State.Base")
        self.stageState = try Self.makePrivateBuffer(heap: heap, length: stateBytes, label: "NumiVivo.Physiology.State.Stage")
        self.candidateState = try Self.makePrivateBuffer(heap: heap, length: stateBytes, label: "NumiVivo.Physiology.State.Candidate")
        self.derivativeK1 = try Self.makePrivateBuffer(heap: heap, length: stateBytes, label: "NumiVivo.Physiology.State.K1")
        self.incidenceOffsets = try Self.makePrivateBuffer(heap: heap, length: offsetBytes, label: "NumiVivo.Physiology.IncidenceOffsets")
        self.incidence = try Self.makePrivateBuffer(heap: heap, length: incidenceBytes, label: "NumiVivo.Physiology.Incidence")
        self.clearances = try Self.makePrivateBuffer(heap: heap, length: clearanceBytes, label: "NumiVivo.Physiology.Clearances")
        self.bounds = try Self.makePrivateBuffer(heap: heap, length: boundsBytes, label: "NumiVivo.Physiology.Bounds")

        self.command = try Self.makeSharedBuffer(
            device: device,
            length: MemoryLayout<VivoPhysiologyRuntimeCommandABI>.stride,
            label: "NumiVivo.Physiology.Command"
        )
        self.status = try Self.makeSharedBuffer(
            device: device,
            length: MemoryLayout<VivoPhysiologyRuntimeStatusABI>.stride,
            label: "NumiVivo.Physiology.Status"
        )
        self.preTransforms = try Self.makeSharedBuffer(
            device: device,
            length: try Self.checkedBytes(configuration.maximumTransformsPerStep, MemoryLayout<VivoPhysiologyStateTransformABI>.stride, label: "pre transforms"),
            label: "NumiVivo.Physiology.PreTransforms"
        )
        self.postTransforms = try Self.makeSharedBuffer(
            device: device,
            length: try Self.checkedBytes(configuration.maximumTransformsPerStep, MemoryLayout<VivoPhysiologyStateTransformABI>.stride, label: "post transforms"),
            label: "NumiVivo.Physiology.PostTransforms"
        )
        self.publicationRequests = try Self.makeSharedBuffer(
            device: device,
            length: try Self.checkedBytes(configuration.maximumPublicationsPerStep, MemoryLayout<VivoPhysiologyPublicationRequestABI>.stride, label: "physiology publication requests"),
            label: "NumiVivo.Physiology.PublicationRequests"
        )
        self.publicationOutput = try Self.makeSharedBuffer(
            device: device,
            length: try Self.checkedBytes(configuration.maximumPublicationsPerStep, MemoryLayout<Float>.stride, label: "physiology publication output"),
            label: "NumiVivo.Physiology.PublicationOutput"
        )
        self.stateReadback = try Self.makeSharedBuffer(
            device: device,
            length: stateBytes,
            label: "NumiVivo.Physiology.State.Readback"
        )

        try Self.uploadInitialContents(
            commandQueue: commandQueue,
            arena: self,
            model: model
        )
    }

    func write(command value: VivoPhysiologyRuntimeCommandABI) {
        var copy = value
        withUnsafeBytes(of: &copy) { bytes in
            command.contents().copyMemory(from: bytes.baseAddress!, byteCount: bytes.count)
        }
    }

    func write(pre values: [VivoPhysiologyStateTransformABI]) throws {
        try writeTransforms(values, destination: preTransforms, capacity: capacities.preTransforms, label: "pre")
    }

    func write(post values: [VivoPhysiologyStateTransformABI]) throws {
        try writeTransforms(values, destination: postTransforms, capacity: capacities.postTransforms, label: "post")
    }

    func write(publications values: [VivoPhysiologyPublicationRequestABI]) throws {
        guard values.count <= capacities.publications else {
            throw VivoRuntimeError.invalidConfiguration("physiology publication count exceeds capacity")
        }
        guard !values.isEmpty else { return }
        values.withUnsafeBytes { bytes in
            publicationRequests.contents().copyMemory(from: bytes.baseAddress!, byteCount: bytes.count)
        }
    }

    func runtimeStatus() -> VivoPhysiologyRuntimeStatus {
        VivoPhysiologyRuntimeStatus(
            raw: status.contents().assumingMemoryBound(to: VivoPhysiologyRuntimeStatusABI.self).pointee
        )
    }

    func publicationValues(count: Int) throws -> [Float] {
        guard count >= 0, count <= capacities.publications else {
            throw VivoRuntimeError.invalidConfiguration("physiology publication count is outside capacity")
        }
        guard count > 0 else { return [] }
        let pointer = publicationOutput.contents().assumingMemoryBound(to: Float.self)
        return Array(UnsafeBufferPointer(start: pointer, count: count))
    }

    func readbackValues() -> [Float] {
        let pointer = stateReadback.contents().assumingMemoryBound(to: Float.self)
        return Array(UnsafeBufferPointer(start: pointer, count: capacities.stateElements))
    }

    func commit() {
        swap(&currentState, &candidateState)
    }

    func uploadCurrentState(
        _ values: [Float],
        commandQueue: MTLCommandQueue,
        label: String
    ) throws {
        guard values.count == capacities.stateElements,
              values.allSatisfy(\.isFinite) else {
            throw VivoRuntimeError.invalidConfiguration("restored physiology state has invalid size or non-finite values")
        }
        let staging = try Self.makeStagingBuffer(
            device: device,
            values: values,
            label: "\(label).Staging"
        )
        guard let commandBuffer = commandQueue.makeCommandBuffer(),
              let blit = commandBuffer.makeBlitCommandEncoder() else {
            throw VivoRuntimeError.commandQueueUnavailable
        }
        commandBuffer.label = label
        let byteCount = values.count * MemoryLayout<Float>.stride
        for destination in [currentState, baseState, stageState, candidateState] {
            blit.copy(from: staging, sourceOffset: 0, to: destination, destinationOffset: 0, size: byteCount)
        }
        blit.fill(buffer: derivativeK1, range: 0..<derivativeK1.length, value: 0)
        blit.endEncoding()
        commandBuffer.commit()
        commandBuffer.waitUntilCompleted()
        if let error = commandBuffer.error {
            throw VivoRuntimeError.commandFailed("\(label): \(error)")
        }
    }

    private func writeTransforms(
        _ values: [VivoPhysiologyStateTransformABI],
        destination: MTLBuffer,
        capacity: Int,
        label: String
    ) throws {
        guard values.count <= capacity else {
            throw VivoRuntimeError.invalidConfiguration("physiology \(label) transform count exceeds capacity")
        }
        guard !values.isEmpty else { return }
        values.withUnsafeBytes { bytes in
            destination.contents().copyMemory(from: bytes.baseAddress!, byteCount: bytes.count)
        }
    }

    private static func uploadInitialContents(
        commandQueue: MTLCommandQueue,
        arena: VivoPhysiologyMetalArena,
        model: PreparedVivoPhysiologyModel
    ) throws {
        let incidenceValues = model.incidence.isEmpty
            ? [VivoPhysiologyIncidenceABI(sourcePairIndex: 0, coefficientPerSecond: 0)]
            : model.incidence
        let boundsValues: [SIMD2<Float>] = try model.analytes.flatMap { analyte in
            guard analyte.minimum.isFinite,
                  analyte.maximum.isFinite,
                  abs(analyte.minimum) <= Double(Float.greatestFiniteMagnitude),
                  abs(analyte.maximum) <= Double(Float.greatestFiniteMagnitude) else {
                throw VivoRuntimeError.invalidConfiguration("physiology analyte bounds are not FP32 representable")
            }
            return [SIMD2<Float>](
                repeating: SIMD2(Float(analyte.minimum), Float(analyte.maximum)),
                count: model.compartments.count
            )
        }

        let stateStaging = try makeStagingBuffer(device: arena.device, values: model.initialState, label: "NumiVivo.Physiology.Stage.State")
        let offsetsStaging = try makeStagingBuffer(device: arena.device, values: model.incidenceOffsets, label: "NumiVivo.Physiology.Stage.Offsets")
        let incidenceStaging = try makeStagingBuffer(device: arena.device, values: incidenceValues, label: "NumiVivo.Physiology.Stage.Incidence")
        let clearanceStaging = try makeStagingBuffer(device: arena.device, values: model.clearances, label: "NumiVivo.Physiology.Stage.Clearance")
        let boundsStaging = try makeStagingBuffer(device: arena.device, values: boundsValues, label: "NumiVivo.Physiology.Stage.Bounds")

        guard let command = commandQueue.makeCommandBuffer(),
              let blit = command.makeBlitCommandEncoder() else {
            throw VivoRuntimeError.commandQueueUnavailable
        }
        command.label = "NumiVivo.Physiology.InitializeArena"
        let stateBytes = model.initialState.count * MemoryLayout<Float>.stride
        for destination in [arena.currentState, arena.baseState, arena.stageState, arena.candidateState] {
            blit.copy(from: stateStaging, sourceOffset: 0, to: destination, destinationOffset: 0, size: stateBytes)
        }
        blit.fill(buffer: arena.derivativeK1, range: 0..<arena.derivativeK1.length, value: 0)
        blit.copy(
            from: offsetsStaging,
            sourceOffset: 0,
            to: arena.incidenceOffsets,
            destinationOffset: 0,
            size: model.incidenceOffsets.count * MemoryLayout<UInt32>.stride
        )
        blit.copy(
            from: incidenceStaging,
            sourceOffset: 0,
            to: arena.incidence,
            destinationOffset: 0,
            size: incidenceValues.count * MemoryLayout<VivoPhysiologyIncidenceABI>.stride
        )
        blit.copy(
            from: clearanceStaging,
            sourceOffset: 0,
            to: arena.clearances,
            destinationOffset: 0,
            size: model.clearances.count * MemoryLayout<VivoPhysiologyClearanceABI>.stride
        )
        blit.copy(
            from: boundsStaging,
            sourceOffset: 0,
            to: arena.bounds,
            destinationOffset: 0,
            size: boundsValues.count * MemoryLayout<SIMD2<Float>>.stride
        )
        blit.endEncoding()
        command.commit()
        command.waitUntilCompleted()
        if let error = command.error {
            throw VivoRuntimeError.commandFailed("physiology arena initialization: \(error)")
        }
    }

    private static func makePrivateBuffer(heap: MTLHeap, length: Int, label: String) throws -> MTLBuffer {
        guard let buffer = heap.makeBuffer(length: max(length, 1), options: privateOptions) else {
            throw VivoRuntimeError.allocationFailed("physiology private heap could not allocate \(label)")
        }
        buffer.label = label
        return buffer
    }

    private static func makeSharedBuffer(device: MTLDevice, length: Int, label: String) throws -> MTLBuffer {
        guard let buffer = device.makeBuffer(length: max(length, 1), options: sharedOptions) else {
            throw VivoRuntimeError.allocationFailed("Metal could not allocate \(label)")
        }
        buffer.label = label
        memset(buffer.contents(), 0, buffer.length)
        return buffer
    }

    private static func makeStagingBuffer<T>(device: MTLDevice, values: [T], label: String) throws -> MTLBuffer {
        let byteCount = max(values.count * MemoryLayout<T>.stride, 1)
        let buffer: MTLBuffer? = values.withUnsafeBytes { bytes in
            if let base = bytes.baseAddress, !bytes.isEmpty {
                return device.makeBuffer(bytes: base, length: bytes.count, options: sharedOptions)
            }
            return device.makeBuffer(length: byteCount, options: sharedOptions)
        }
        guard let buffer else {
            throw VivoRuntimeError.allocationFailed("Metal could not allocate \(label)")
        }
        buffer.label = label
        return buffer
    }

    private static func requiredHeapSize(
        device: MTLDevice,
        requirements: [(String, Int)],
        options: MTLResourceOptions
    ) throws -> Int {
        var cursor = 0
        for (label, length) in requirements {
            let sizeAndAlign = device.heapBufferSizeAndAlign(length: max(length, 1), options: options)
            cursor = try align(cursor, to: sizeAndAlign.align, label: label)
            let result = cursor.addingReportingOverflow(sizeAndAlign.size)
            guard !result.overflow else {
                throw VivoRuntimeError.allocationFailed("physiology heap size overflow while adding \(label)")
            }
            cursor = result.partialValue
        }
        return try align(cursor, to: 4_096, label: "physiology heap")
    }

    private static func applyHeadroom(_ value: Int, factor: Double) throws -> Int {
        let expanded = Double(value) * factor
        guard expanded.isFinite, expanded <= Double(Int.max) else {
            throw VivoRuntimeError.allocationFailed("physiology heap headroom overflow")
        }
        return try align(Int(expanded.rounded(.up)), to: 4_096, label: "physiology heap headroom")
    }

    private static func align(_ value: Int, to alignment: Int, label: String) throws -> Int {
        guard alignment > 0 else {
            throw VivoRuntimeError.allocationFailed("invalid alignment for \(label)")
        }
        let remainder = value % alignment
        if remainder == 0 { return value }
        let result = value.addingReportingOverflow(alignment - remainder)
        guard !result.overflow else {
            throw VivoRuntimeError.allocationFailed("alignment overflow for \(label)")
        }
        return result.partialValue
    }

    private static func checkedBytes(_ count: Int, _ stride: Int, label: String) throws -> Int {
        let result = count.multipliedReportingOverflow(by: stride)
        guard !result.overflow else {
            throw VivoRuntimeError.allocationFailed("\(label) byte count overflow")
        }
        return max(result.partialValue, 1)
    }
}
