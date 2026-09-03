import Foundation
import Metal

final class VivoMetalArena: @unchecked Sendable {
    struct Capacities: Sendable {
        let stateElements: Int
        let temporalElements: Int
        let reactionEventElements: Int
        let couplingUpdates: Int
        let publicationRequests: Int
        let eventRecords: Int
    }

    let device: MTLDevice
    let heap: MTLHeap
    let capabilities: VivoMetalCapabilities
    let pack: VivoProgramPack
    let configuration: VivoRuntimeConfiguration
    let capacities: Capacities

    let programBuffer: MTLBuffer
    let parameterBuffer: MTLBuffer
    var currentState: MTLBuffer
    let baseState: MTLBuffer
    let stageState: MTLBuffer
    var candidateState: MTLBuffer
    let derivativeK1: MTLBuffer
    var temporalCurrent: MTLBuffer
    var temporalCandidate: MTLBuffer
    let reactionEvents: MTLBuffer
    let transport: MTLBuffer
    let velocity: MTLBuffer
    let volumeFraction: MTLBuffer

    let commandBuffer: MTLBuffer
    let statusBuffer: MTLBuffer
    let eventBuffer: MTLBuffer
    let couplingBuffer: MTLBuffer
    let publicationRequestBuffer: MTLBuffer
    let publicationOutputBuffer: MTLBuffer
    let stateReadbackBuffer: MTLBuffer

    private static let privateOptions: MTLResourceOptions = .storageModePrivate
    private static let sharedOptions: MTLResourceOptions = .storageModeShared

    init(
        device: MTLDevice,
        commandQueue: MTLCommandQueue,
        pack: VivoProgramPack,
        configuration: VivoRuntimeConfiguration,
        couplingCapacity: Int = 65_536,
        publicationCapacity: Int = 65_536
    ) throws {
        try VivoRuntimeCommandABI.validateMemoryLayout()
        try VivoRuntimeStatusABI.validateMemoryLayout()
        try configuration.validate(for: pack)

        guard MemoryLayout<VivoSpeciesTransportABI>.stride == 16,
              MemoryLayout<VivoCouplingUpdateABI>.stride == 16,
              MemoryLayout<VivoPublicationRequestABI>.stride == 16,
              MemoryLayout<VivoEventABI>.stride == 32 else {
            throw VivoRuntimeError.incompatibleDevice("Swift–Metal record layouts do not match the runtime ABI")
        }

        self.device = device
        self.pack = pack
        self.configuration = configuration
        self.capabilities = VivoMetalCapabilities(device: device)

        let speciesCount = Int(pack.runtimeContract.speciesCount)
        let parameterCount = Int(pack.runtimeContract.parameterCount)
        let reactionCount = Int(pack.runtimeContract.reactionCount)
        let laneCount = Int(configuration.laneCount)
        let environmentCount = Int(configuration.environmentCount)
        let temporalCount = Int(pack.runtimeContract.temporalStateCount)
        let stateElements = try Self.checkedProduct(speciesCount, laneCount, label: "state")
        let temporalElements = max(try Self.checkedProduct(temporalCount, laneCount, label: "temporal state"), 1)
        let reactionEvents = max(try Self.checkedProduct(reactionCount, laneCount, label: "reaction events"), 1)
        let parameterElements = max(try Self.checkedProduct(parameterCount, environmentCount, label: "parameters"), 1)
        let eventCapacity = Int(configuration.eventCapacity)
        let effectiveCouplingCapacity = max(couplingCapacity, laneCount, 1)
        let effectivePublicationCapacity = max(publicationCapacity, 1)
        self.capacities = Capacities(
            stateElements: stateElements,
            temporalElements: temporalElements,
            reactionEventElements: reactionEvents,
            couplingUpdates: effectiveCouplingCapacity,
            publicationRequests: effectivePublicationCapacity,
            eventRecords: eventCapacity
        )

        let stateBytes = try Self.checkedBytes(stateElements, MemoryLayout<Float>.stride, label: "state")
        let temporalBytes = try Self.checkedBytes(temporalElements, MemoryLayout<SIMD2<Float>>.stride, label: "temporal state")
        let reactionEventBytes = try Self.checkedBytes(reactionEvents, MemoryLayout<Int32>.stride, label: "reaction events")
        let parameterBytes = try Self.checkedBytes(parameterElements, MemoryLayout<Float>.stride, label: "parameters")
        let transportBytes = try Self.checkedBytes(max(speciesCount, 1), MemoryLayout<VivoSpeciesTransportABI>.stride, label: "transport")
        let velocityBytes = try Self.checkedBytes(max(laneCount, 1), MemoryLayout<SIMD4<Float>>.stride, label: "velocity")
        let volumeFractionBytes = try Self.checkedBytes(max(laneCount, 1), MemoryLayout<Float>.stride, label: "volume fractions")
        let programBytes = max(pack.data.count, 1)

        let privateRequirements: [(String, Int)] = [
            ("ProgramPack", programBytes),
            ("parameters", parameterBytes),
            ("current-state", stateBytes),
            ("base-state", stateBytes),
            ("stage-state", stateBytes),
            ("candidate-state", stateBytes),
            ("derivative-k1", stateBytes),
            ("temporal-current", temporalBytes),
            ("temporal-candidate", temporalBytes),
            ("reaction-events", reactionEventBytes),
            ("species-transport", transportBytes),
            ("velocity", velocityBytes),
            ("volume-fraction", volumeFractionBytes)
        ]
        let requiredHeapSize = try Self.requiredHeapSize(
            device: device,
            requirements: privateRequirements,
            options: Self.privateOptions
        )
        let requestedHeapSize = try Self.applyHeadroom(
            requiredHeapSize,
            factor: configuration.privateHeapHeadroom
        )
        try capabilities.validate(
            programBytes: UInt64(programBytes),
            privateHeapBytes: UInt64(requestedHeapSize)
        )

        let descriptor = MTLHeapDescriptor()
        descriptor.label = "NumiVivo.PrivateArena"
        descriptor.storageMode = .private
        descriptor.cpuCacheMode = .defaultCache
        descriptor.hazardTrackingMode = .tracked
        descriptor.size = requestedHeapSize
        guard let heap = device.makeHeap(descriptor: descriptor) else {
            throw VivoRuntimeError.allocationFailed("Metal could not allocate a \(requestedHeapSize)-byte private heap")
        }
        self.heap = heap

        self.programBuffer = try Self.makePrivateBuffer(heap: heap, length: programBytes, label: "NumiVivo.ProgramPack")
        self.parameterBuffer = try Self.makePrivateBuffer(heap: heap, length: parameterBytes, label: "NumiVivo.Parameters")
        self.currentState = try Self.makePrivateBuffer(heap: heap, length: stateBytes, label: "NumiVivo.State.Current")
        self.baseState = try Self.makePrivateBuffer(heap: heap, length: stateBytes, label: "NumiVivo.State.Base")
        self.stageState = try Self.makePrivateBuffer(heap: heap, length: stateBytes, label: "NumiVivo.State.Stage")
        self.candidateState = try Self.makePrivateBuffer(heap: heap, length: stateBytes, label: "NumiVivo.State.Candidate")
        self.derivativeK1 = try Self.makePrivateBuffer(heap: heap, length: stateBytes, label: "NumiVivo.State.K1")
        self.temporalCurrent = try Self.makePrivateBuffer(heap: heap, length: temporalBytes, label: "NumiVivo.Temporal.Current")
        self.temporalCandidate = try Self.makePrivateBuffer(heap: heap, length: temporalBytes, label: "NumiVivo.Temporal.Candidate")
        self.reactionEvents = try Self.makePrivateBuffer(heap: heap, length: reactionEventBytes, label: "NumiVivo.ReactionEvents")
        self.transport = try Self.makePrivateBuffer(heap: heap, length: transportBytes, label: "NumiVivo.Transport")
        self.velocity = try Self.makePrivateBuffer(heap: heap, length: velocityBytes, label: "NumiVivo.Velocity")
        self.volumeFraction = try Self.makePrivateBuffer(heap: heap, length: volumeFractionBytes, label: "NumiVivo.VolumeFraction")

        self.commandBuffer = try Self.makeSharedBuffer(device: device, length: MemoryLayout<VivoRuntimeCommandABI>.stride, label: "NumiVivo.Command")
        self.statusBuffer = try Self.makeSharedBuffer(device: device, length: MemoryLayout<VivoRuntimeStatusABI>.stride, label: "NumiVivo.Status")
        self.eventBuffer = try Self.makeSharedBuffer(
            device: device,
            length: try Self.checkedBytes(eventCapacity, MemoryLayout<VivoEventABI>.stride, label: "event log"),
            label: "NumiVivo.Events"
        )
        self.couplingBuffer = try Self.makeSharedBuffer(
            device: device,
            length: try Self.checkedBytes(effectiveCouplingCapacity, MemoryLayout<VivoCouplingUpdateABI>.stride, label: "coupling updates"),
            label: "NumiVivo.CouplingUpdates"
        )
        self.publicationRequestBuffer = try Self.makeSharedBuffer(
            device: device,
            length: try Self.checkedBytes(effectivePublicationCapacity, MemoryLayout<VivoPublicationRequestABI>.stride, label: "publication requests"),
            label: "NumiVivo.PublicationRequests"
        )
        self.publicationOutputBuffer = try Self.makeSharedBuffer(
            device: device,
            length: try Self.checkedBytes(effectivePublicationCapacity, MemoryLayout<Float>.stride, label: "publication output"),
            label: "NumiVivo.PublicationOutput"
        )
        self.stateReadbackBuffer = try Self.makeSharedBuffer(device: device, length: stateBytes, label: "NumiVivo.State.Readback")

        try Self.uploadInitialContents(
            commandQueue: commandQueue,
            pack: pack,
            configuration: configuration,
            arena: self
        )
    }

    func bindProgramSection(
        _ kind: VivoProgramPack.SectionKind,
        to encoder: MTLComputeCommandEncoder,
        index: Int
    ) throws {
        let section = try pack.section(kind)
        guard section.offset <= UInt64(Int.max) else {
            throw VivoRuntimeError.packError("section \(kind) offset exceeds Int.max")
        }
        encoder.setBuffer(programBuffer, offset: Int(section.offset), index: index)
    }

    func write(command: VivoRuntimeCommandABI) {
        var command = command
        withUnsafeBytes(of: &command) { bytes in
            commandBuffer.contents().copyMemory(from: bytes.baseAddress!, byteCount: bytes.count)
        }
    }

    func write(couplingUpdates: [VivoCouplingUpdateABI]) throws {
        guard couplingUpdates.count <= capacities.couplingUpdates else {
            throw VivoRuntimeError.invalidConfiguration(
                "coupling update count \(couplingUpdates.count) exceeds capacity \(capacities.couplingUpdates)"
            )
        }
        guard !couplingUpdates.isEmpty else { return }
        couplingUpdates.withUnsafeBytes { bytes in
            couplingBuffer.contents().copyMemory(from: bytes.baseAddress!, byteCount: bytes.count)
        }
    }

    func write(publicationRequests: [VivoPublicationRequestABI]) throws {
        guard publicationRequests.count <= capacities.publicationRequests else {
            throw VivoRuntimeError.invalidConfiguration(
                "publication request count \(publicationRequests.count) exceeds capacity \(capacities.publicationRequests)"
            )
        }
        guard !publicationRequests.isEmpty else { return }
        publicationRequests.withUnsafeBytes { bytes in
            publicationRequestBuffer.contents().copyMemory(from: bytes.baseAddress!, byteCount: bytes.count)
        }
    }

    func status() -> VivoRuntimeStatus {
        let raw = statusBuffer.contents().assumingMemoryBound(to: VivoRuntimeStatusABI.self).pointee
        return VivoRuntimeStatus(raw: raw)
    }

    func events(status: VivoRuntimeStatus) -> [VivoEvent] {
        let count = min(Int(status.eventCount), capacities.eventRecords)
        guard count > 0 else { return [] }
        let pointer = eventBuffer.contents().assumingMemoryBound(to: VivoEventABI.self)
        return UnsafeBufferPointer(start: pointer, count: count).map(VivoEvent.init)
    }

    func publicationValues(count: Int) throws -> [Float] {
        guard count >= 0, count <= capacities.publicationRequests else {
            throw VivoRuntimeError.invalidConfiguration("publication output count is outside the allocated capacity")
        }
        guard count > 0 else { return [] }
        let pointer = publicationOutputBuffer.contents().assumingMemoryBound(to: Float.self)
        return Array(UnsafeBufferPointer(start: pointer, count: count))
    }

    func commit() {
        swap(&currentState, &candidateState)
        swap(&temporalCurrent, &temporalCandidate)
    }

    func replaceTransport(
        _ values: [VivoSpeciesTransportABI],
        commandQueue: MTLCommandQueue
    ) throws {
        let required = Int(pack.runtimeContract.speciesCount)
        guard values.count == required else {
            throw VivoRuntimeError.invalidConfiguration("transport table must contain exactly \(required) records")
        }
        try uploadArray(values, to: transport, commandQueue: commandQueue, label: "NumiVivo.Upload.Transport")
    }

    func replaceVelocity(
        _ values: [SIMD4<Float>],
        commandQueue: MTLCommandQueue
    ) throws {
        guard values.count == Int(configuration.laneCount) else {
            throw VivoRuntimeError.invalidConfiguration("velocity field must contain one float4 per lane")
        }
        try uploadArray(values, to: velocity, commandQueue: commandQueue, label: "NumiVivo.Upload.Velocity")
    }

    func replaceVolumeFraction(
        _ values: [Float],
        commandQueue: MTLCommandQueue
    ) throws {
        guard values.count == Int(configuration.laneCount),
              values.allSatisfy({ $0.isFinite && $0 >= 0 && $0 <= 1 }) else {
            throw VivoRuntimeError.invalidConfiguration("volume fractions must contain one finite [0,1] value per lane")
        }
        try uploadArray(values, to: volumeFraction, commandQueue: commandQueue, label: "NumiVivo.Upload.VolumeFraction")
    }

    private func uploadArray<T>(
        _ values: [T],
        to destination: MTLBuffer,
        commandQueue: MTLCommandQueue,
        label: String
    ) throws {
        let byteCount = values.count * MemoryLayout<T>.stride
        guard byteCount <= destination.length else {
            throw VivoRuntimeError.allocationFailed("upload exceeds destination buffer \(destination.label ?? "unnamed")")
        }
        let staging: MTLBuffer? = values.withUnsafeBytes { bytes in
            guard let base = bytes.baseAddress else { return nil }
            return device.makeBuffer(bytes: base, length: byteCount, options: Self.sharedOptions)
        }
        guard let staging else {
            throw VivoRuntimeError.allocationFailed("Metal could not allocate upload staging for \(label)")
        }
        guard let command = commandQueue.makeCommandBuffer(), let blit = command.makeBlitCommandEncoder() else {
            throw VivoRuntimeError.commandQueueUnavailable
        }
        command.label = label
        blit.copy(from: staging, sourceOffset: 0, to: destination, destinationOffset: 0, size: byteCount)
        blit.endEncoding()
        command.commit()
        command.waitUntilCompleted()
        if let error = command.error {
            throw VivoRuntimeError.commandFailed("\(label): \(error)")
        }
    }

    private static func uploadInitialContents(
        commandQueue: MTLCommandQueue,
        pack: VivoProgramPack,
        configuration: VivoRuntimeConfiguration,
        arena: VivoMetalArena
    ) throws {
        let initialState = try pack.initialState(laneCount: Int(configuration.laneCount))
        let parameterValues = try pack.parameterValues(environmentCount: Int(configuration.environmentCount))
        let transportValues = [VivoSpeciesTransportABI](
            repeating: .init(),
            count: max(Int(pack.runtimeContract.speciesCount), 1)
        )
        let velocityValues = [SIMD4<Float>](
            repeating: .zero,
            count: max(Int(configuration.laneCount), 1)
        )
        let volumeFractions = [Float](
            repeating: 1,
            count: max(Int(configuration.laneCount), 1)
        )

        let programStaging = try makeStagingBuffer(device: arena.device, data: pack.data, label: "NumiVivo.Stage.ProgramPack")
        let stateStaging = try makeStagingBuffer(device: arena.device, values: initialState, label: "NumiVivo.Stage.InitialState")
        let parameterStaging = try makeStagingBuffer(device: arena.device, values: parameterValues, label: "NumiVivo.Stage.Parameters")
        let transportStaging = try makeStagingBuffer(device: arena.device, values: transportValues, label: "NumiVivo.Stage.Transport")
        let velocityStaging = try makeStagingBuffer(device: arena.device, values: velocityValues, label: "NumiVivo.Stage.Velocity")
        let volumeStaging = try makeStagingBuffer(device: arena.device, values: volumeFractions, label: "NumiVivo.Stage.VolumeFraction")

        guard let command = commandQueue.makeCommandBuffer(), let blit = command.makeBlitCommandEncoder() else {
            throw VivoRuntimeError.commandQueueUnavailable
        }
        command.label = "NumiVivo.InitializeArena"
        blit.copy(from: programStaging, sourceOffset: 0, to: arena.programBuffer, destinationOffset: 0, size: pack.data.count)
        blit.copy(from: parameterStaging, sourceOffset: 0, to: arena.parameterBuffer, destinationOffset: 0, size: parameterValues.count * MemoryLayout<Float>.stride)
        for stateBuffer in [arena.currentState, arena.baseState, arena.stageState, arena.candidateState, arena.derivativeK1] {
            blit.copy(from: stateStaging, sourceOffset: 0, to: stateBuffer, destinationOffset: 0, size: initialState.count * MemoryLayout<Float>.stride)
        }
        blit.fill(buffer: arena.temporalCurrent, range: 0..<arena.temporalCurrent.length, value: 0)
        blit.fill(buffer: arena.temporalCandidate, range: 0..<arena.temporalCandidate.length, value: 0)
        blit.fill(buffer: arena.reactionEvents, range: 0..<arena.reactionEvents.length, value: 0)
        blit.copy(from: transportStaging, sourceOffset: 0, to: arena.transport, destinationOffset: 0, size: transportValues.count * MemoryLayout<VivoSpeciesTransportABI>.stride)
        blit.copy(from: velocityStaging, sourceOffset: 0, to: arena.velocity, destinationOffset: 0, size: velocityValues.count * MemoryLayout<SIMD4<Float>>.stride)
        blit.copy(from: volumeStaging, sourceOffset: 0, to: arena.volumeFraction, destinationOffset: 0, size: volumeFractions.count * MemoryLayout<Float>.stride)
        blit.endEncoding()
        command.commit()
        command.waitUntilCompleted()
        if let error = command.error {
            throw VivoRuntimeError.commandFailed("arena initialization: \(error)")
        }
    }

    private static func makeStagingBuffer(
        device: MTLDevice,
        data: Data,
        label: String
    ) throws -> MTLBuffer {
        guard !data.isEmpty else {
            throw VivoRuntimeError.packError("cannot stage an empty ProgramPack")
        }
        let buffer: MTLBuffer? = data.withUnsafeBytes { bytes in
            guard let base = bytes.baseAddress else { return nil }
            return device.makeBuffer(bytes: base, length: bytes.count, options: sharedOptions)
        }
        guard let buffer else {
            throw VivoRuntimeError.allocationFailed("could not allocate \(label)")
        }
        buffer.label = label
        return buffer
    }

    private static func makeStagingBuffer<T>(
        device: MTLDevice,
        values: [T],
        label: String
    ) throws -> MTLBuffer {
        let byteCount = max(values.count * MemoryLayout<T>.stride, 1)
        let buffer: MTLBuffer? = values.withUnsafeBytes { bytes in
            if let base = bytes.baseAddress, bytes.count > 0 {
                return device.makeBuffer(bytes: base, length: bytes.count, options: sharedOptions)
            }
            return device.makeBuffer(length: byteCount, options: sharedOptions)
        }
        guard let buffer else {
            throw VivoRuntimeError.allocationFailed("could not allocate \(label)")
        }
        buffer.label = label
        return buffer
    }

    private static func makePrivateBuffer(heap: MTLHeap, length: Int, label: String) throws -> MTLBuffer {
        guard let buffer = heap.makeBuffer(length: max(length, 1), options: privateOptions) else {
            throw VivoRuntimeError.allocationFailed("private heap could not allocate \(label) (\(length) bytes)")
        }
        buffer.label = label
        return buffer
    }

    private static func makeSharedBuffer(device: MTLDevice, length: Int, label: String) throws -> MTLBuffer {
        guard let buffer = device.makeBuffer(length: max(length, 1), options: sharedOptions) else {
            throw VivoRuntimeError.allocationFailed("device could not allocate \(label) (\(length) bytes)")
        }
        buffer.label = label
        memset(buffer.contents(), 0, buffer.length)
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
                throw VivoRuntimeError.allocationFailed("private heap size overflow while adding \(label)")
            }
            cursor = result.partialValue
        }
        return try align(cursor, to: 4_096, label: "heap final alignment")
    }

    private static func applyHeadroom(_ value: Int, factor: Double) throws -> Int {
        let expanded = Double(value) * factor
        guard expanded.isFinite, expanded <= Double(Int.max) else {
            throw VivoRuntimeError.allocationFailed("private heap headroom calculation overflow")
        }
        return try align(Int(expanded.rounded(.up)), to: 4_096, label: "heap headroom")
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

    private static func checkedProduct(_ lhs: Int, _ rhs: Int, label: String) throws -> Int {
        let result = lhs.multipliedReportingOverflow(by: rhs)
        guard !result.overflow else {
            throw VivoRuntimeError.allocationFailed("\(label) element count overflow")
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
