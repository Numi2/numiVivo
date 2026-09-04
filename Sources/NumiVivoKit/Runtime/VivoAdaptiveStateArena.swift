import Foundation
import Metal
import NumiVivoShaders

public enum VivoAdaptiveScratchRegion: String, Codable, Sendable {
    case randomCounters
    case delayedQueue
    case topology
    case auxiliary
}

private struct VivoAdaptiveChunkKey: Hashable {
    let cohortID: UInt32
    let chunkIndex: UInt32
}

private final class VivoAdaptiveArenaStorage: @unchecked Sendable {
    let heap: MTLHeap
    let chunks: [VivoAdaptiveChunkKey: VivoAdaptiveArenaChunkResources]
    let manifest: VivoAdaptiveArenaManifest

    init(
        heap: MTLHeap,
        chunks: [VivoAdaptiveChunkKey: VivoAdaptiveArenaChunkResources],
        manifest: VivoAdaptiveArenaManifest
    ) {
        self.heap = heap
        self.chunks = chunks
        self.manifest = manifest
    }
}

private struct VivoMigrationCommandABI {
    var elementCount: UInt32
    var reserved0: UInt32 = 0
    var integerTolerance: Float
    var reserved1: Float = 0
}

private struct VivoMigrationStatusABI {
    var flags: UInt32
    var invalidElementCount: UInt32
    var firstInvalidElement: UInt32
    var maximumErrorBits: UInt32
}

private struct VivoMigrationPipelines {
    let clear: MTLComputePipelineState
    let floatToCount: MTLComputePipelineState
    let countToFloat: MTLComputePipelineState
    let validateFloat: MTLComputePipelineState
}

/// Materializes an adaptive-fidelity plan into private Metal heap resources.
/// Reconfiguration builds and validates a complete replacement arena, then swaps
/// ownership only after every conversion and preservation copy has completed.
public actor VivoAdaptiveStateArena {
    public nonisolated let deviceName: String
    public nonisolated let registryID: UInt64

    private let device: MTLDevice
    private let queue: MTLCommandQueue
    private let pipelines: VivoMigrationPipelines
    private let migrationStatus: MTLBuffer
    private var planStorage: VivoAdaptiveFidelityPlan
    private var layoutStorage: [UInt32: VivoAdaptiveCohortRuntimeLayout]
    private var arenaStorage: VivoAdaptiveArenaStorage

    public static func make(
        plan: VivoAdaptiveFidelityPlan,
        layouts: [VivoAdaptiveCohortRuntimeLayout],
        device: MTLDevice
    ) async throws -> VivoAdaptiveStateArena {
        guard let queue = device.makeCommandQueue() else {
            throw VivoAdaptiveArenaError.commandQueueUnavailable
        }
        queue.label = "NumiVivo.AdaptiveArena.Queue"
        let layoutMap = try validate(plan: plan, layouts: layouts)
        let pipelines = try loadPipelines(device: device)
        guard let migrationStatus = device.makeBuffer(
            length: MemoryLayout<VivoMigrationStatusABI>.stride,
            options: .storageModeShared
        ) else {
            throw VivoAdaptiveArenaError.allocationFailed(0, 0, "migration status")
        }
        migrationStatus.label = "NumiVivo.AdaptiveArena.MigrationStatus"
        let storage = try allocate(
            plan: plan,
            layouts: layoutMap,
            device: device,
            queue: queue
        )
        return VivoAdaptiveStateArena(
            plan: plan,
            layouts: layoutMap,
            device: device,
            queue: queue,
            pipelines: pipelines,
            migrationStatus: migrationStatus,
            storage: storage
        )
    }

    private init(
        plan: VivoAdaptiveFidelityPlan,
        layouts: [UInt32: VivoAdaptiveCohortRuntimeLayout],
        device: MTLDevice,
        queue: MTLCommandQueue,
        pipelines: VivoMigrationPipelines,
        migrationStatus: MTLBuffer,
        storage: VivoAdaptiveArenaStorage
    ) {
        self.planStorage = plan
        self.layoutStorage = layouts
        self.device = device
        self.queue = queue
        self.pipelines = pipelines
        self.migrationStatus = migrationStatus
        self.arenaStorage = storage
        self.deviceName = device.name
        self.registryID = device.registryID
    }

    public func plan() -> VivoAdaptiveFidelityPlan { planStorage }
    public func layouts() -> [VivoAdaptiveCohortRuntimeLayout] {
        layoutStorage.values.sorted { $0.cohortID < $1.cohortID }
    }
    public func manifest() -> VivoAdaptiveArenaManifest { arenaStorage.manifest }

    public func resources(
        cohortID: UInt32,
        chunkIndex: UInt32
    ) throws -> VivoAdaptiveArenaChunkResources {
        guard let resources = arenaStorage.chunks[.init(cohortID: cohortID, chunkIndex: chunkIndex)] else {
            throw VivoAdaptiveArenaError.missingChunk(cohortID, chunkIndex)
        }
        return resources
    }

    public func uploadState(
        _ data: Data,
        cohortID: UInt32,
        chunkIndex: UInt32,
        copyIndex: UInt32 = 0,
        clearRemainder: Bool = true
    ) throws {
        let resources = try resources(cohortID: cohortID, chunkIndex: chunkIndex)
        let offset = try resources.descriptor.stateCopyOffset(copyIndex)
        let capacity = Int(resources.descriptor.stateCopyBytes)
        guard data.count <= capacity else {
            throw VivoAdaptiveArenaError.invalidUpload(
                cohortID,
                chunkIndex,
                "state payload has \(data.count) bytes; copy capacity is \(capacity)"
            )
        }
        try upload(
            data,
            destination: resources.state,
            destinationOffset: offset,
            clearRange: clearRemainder ? (offset + data.count)..<(offset + capacity) : nil,
            label: "StateUpload.\(cohortID).\(chunkIndex).\(copyIndex)"
        )
    }

    public func uploadScratch(
        _ data: Data,
        cohortID: UInt32,
        chunkIndex: UInt32,
        region: VivoAdaptiveScratchRegion,
        clearRemainder: Bool = true
    ) throws {
        let resources = try resources(cohortID: cohortID, chunkIndex: chunkIndex)
        guard let layout = layoutStorage[cohortID] else {
            throw VivoAdaptiveArenaError.missingLayout(cohortID)
        }
        let regionInfo = try scratchRegion(
            region,
            layout: layout,
            offsets: resources.descriptor.scratchOffsets
        )
        guard data.count <= regionInfo.length else {
            throw VivoAdaptiveArenaError.invalidUpload(
                cohortID,
                chunkIndex,
                "scratch payload exceeds \(region.rawValue) capacity"
            )
        }
        try upload(
            data,
            destination: resources.scratch,
            destinationOffset: regionInfo.offset,
            clearRange: clearRemainder
                ? (regionInfo.offset + data.count)..<(regionInfo.offset + regionInfo.length)
                : nil,
            label: "ScratchUpload.\(cohortID).\(chunkIndex).\(region.rawValue)"
        )
    }

    public func checkpointSections(
        includeScratch: Bool = true
    ) async throws -> [VivoCheckpointSectionPayload] {
        let ordered = arenaStorage.chunks.values.sorted {
            if $0.descriptor.cohortID != $1.descriptor.cohortID {
                return $0.descriptor.cohortID < $1.descriptor.cohortID
            }
            return $0.descriptor.chunkIndex < $1.descriptor.chunkIndex
        }
        guard let command = queue.makeCommandBuffer(),
              let blit = command.makeBlitCommandEncoder() else {
            throw VivoAdaptiveArenaError.commandQueueUnavailable
        }
        command.label = "NumiVivo.AdaptiveArena.Checkpoint"

        struct Pending {
            let id: String
            let encoding: VivoCheckpointSectionEncoding
            let elementCount: UInt64
            let stride: UInt32
            let buffer: MTLBuffer
            let length: Int
        }
        var pending: [Pending] = []
        for chunk in ordered {
            let stateLength = chunk.state.length
            guard let stateReadback = device.makeBuffer(length: stateLength, options: .storageModeShared) else {
                throw VivoAdaptiveArenaError.allocationFailed(
                    chunk.descriptor.cohortID,
                    chunk.descriptor.chunkIndex,
                    "state checkpoint readback"
                )
            }
            blit.copy(from: chunk.state, sourceOffset: 0, to: stateReadback, destinationOffset: 0, size: stateLength)
            pending.append(.init(
                id: "cohort.\(chunk.descriptor.cohortID).chunk.\(chunk.descriptor.chunkIndex).state",
                encoding: .rawBytes,
                elementCount: UInt64(stateLength),
                stride: 1,
                buffer: stateReadback,
                length: stateLength
            ))

            if includeScratch {
                let scratchLength = chunk.scratch.length
                guard let scratchReadback = device.makeBuffer(length: scratchLength, options: .storageModeShared) else {
                    throw VivoAdaptiveArenaError.allocationFailed(
                        chunk.descriptor.cohortID,
                        chunk.descriptor.chunkIndex,
                        "scratch checkpoint readback"
                    )
                }
                blit.copy(from: chunk.scratch, sourceOffset: 0, to: scratchReadback, destinationOffset: 0, size: scratchLength)
                pending.append(.init(
                    id: "cohort.\(chunk.descriptor.cohortID).chunk.\(chunk.descriptor.chunkIndex).scratch",
                    encoding: .rawBytes,
                    elementCount: UInt64(scratchLength),
                    stride: 1,
                    buffer: scratchReadback,
                    length: scratchLength
                ))
            }
        }
        blit.endEncoding()
        try await complete(command)

        return pending.map { item in
            let bytes = Data(bytes: item.buffer.contents(), count: item.length)
            return .init(
                id: item.id,
                encoding: item.encoding,
                elementCount: item.elementCount,
                elementStride: item.stride,
                data: bytes
            )
        }
    }

    public func reconfigure(
        to newPlan: VivoAdaptiveFidelityPlan,
        layouts newLayouts: [VivoAdaptiveCohortRuntimeLayout]
    ) async throws -> VivoAdaptiveArenaMigrationCertificate {
        let newLayoutMap = try Self.validate(plan: newPlan, layouts: newLayouts)
        let replacement = try Self.allocate(
            plan: newPlan,
            layouts: newLayoutMap,
            device: device,
            queue: queue
        )
        guard let command = queue.makeCommandBuffer() else {
            throw VivoAdaptiveArenaError.commandQueueUnavailable
        }
        command.label = "NumiVivo.AdaptiveArena.Reconfigure"

        try encode(
            pipeline: pipelines.clear,
            commandBuffer: command,
            count: 1,
            label: "AdaptiveArena.ClearMigrationStatus"
        ) { encoder in
            encoder.setBuffer(migrationStatus, offset: 0, index: 0)
        }

        var records: [VivoAdaptiveArenaMigrationRecord] = []
        let orderedDestination = replacement.chunks.values.sorted {
            if $0.descriptor.cohortID != $1.descriptor.cohortID {
                return $0.descriptor.cohortID < $1.descriptor.cohortID
            }
            return $0.descriptor.chunkIndex < $1.descriptor.chunkIndex
        }

        for destination in orderedDestination {
            let cohortID = destination.descriptor.cohortID
            let chunkIndex = destination.descriptor.chunkIndex
            let key = VivoAdaptiveChunkKey(cohortID: cohortID, chunkIndex: chunkIndex)
            guard let newLayout = newLayoutMap[cohortID],
                  let newCohortPlan = newPlan.cohorts.first(where: { $0.cohortID == cohortID }) else {
                throw VivoAdaptiveArenaError.missingLayout(cohortID)
            }
            guard let source = arenaStorage.chunks[key],
                  let oldLayout = layoutStorage[cohortID] else {
                records.append(.init(
                    cohortID: cohortID,
                    chunkIndex: chunkIndex,
                    kind: .resizeAndCopy,
                    sourceRepresentation: newLayout.representation,
                    destinationRepresentation: newLayout.representation,
                    copiedStateBytes: 0,
                    copiedRandomCounterBytes: 0,
                    copiedDelayedQueueBytes: 0,
                    rebuiltTopology: true
                ))
                continue
            }

            let migration = newCohortPlan.migration
            let stateBytes = try migrateState(
                from: source,
                oldLayout: oldLayout,
                to: destination,
                newLayout: newLayout,
                migration: migration,
                commandBuffer: command
            )
            let preserved = try migrateScratch(
                from: source,
                oldLayout: oldLayout,
                to: destination,
                newLayout: newLayout,
                migration: migration,
                commandBuffer: command
            )
            records.append(.init(
                cohortID: cohortID,
                chunkIndex: chunkIndex,
                kind: migration.kind,
                sourceRepresentation: oldLayout.representation,
                destinationRepresentation: newLayout.representation,
                copiedStateBytes: stateBytes,
                copiedRandomCounterBytes: preserved.random,
                copiedDelayedQueueBytes: preserved.delayed,
                rebuiltTopology: migration.kind == .rebuildSpatialTopology || migration.kind == .rebuildTissueCoupling
            ))
        }

        try await complete(command)
        let status = migrationStatus.contents().load(as: VivoMigrationStatusABI.self)
        let migrationStatusValue = VivoAdaptiveMigrationStatus(
            flags: .init(rawValue: status.flags),
            invalidElementCount: status.invalidElementCount,
            firstInvalidElement: status.firstInvalidElement == UInt32.max ? nil : status.firstInvalidElement,
            maximumError: Float(bitPattern: status.maximumErrorBits)
        )
        guard migrationStatusValue.flags.isEmpty else {
            throw VivoAdaptiveArenaError.migrationRejected(migrationStatusValue)
        }

        let sourceFingerprint = arenaStorage.manifest.fingerprint
        planStorage = newPlan
        layoutStorage = newLayoutMap
        arenaStorage = replacement
        return .init(
            sourceManifestFingerprint: sourceFingerprint,
            destinationManifestFingerprint: replacement.manifest.fingerprint,
            committed: true,
            records: records,
            status: migrationStatusValue
        )
    }

    private func migrateState(
        from source: VivoAdaptiveArenaChunkResources,
        oldLayout: VivoAdaptiveCohortRuntimeLayout,
        to destination: VivoAdaptiveArenaChunkResources,
        newLayout: VivoAdaptiveCohortRuntimeLayout,
        migration: VivoBufferMigrationPlan,
        commandBuffer: MTLCommandBuffer
    ) throws -> UInt64 {
        let copyCount = min(oldLayout.stateCopyCount, newLayout.stateCopyCount)
        let oldCopyBytes = try oldLayout.stateCopyBytes()
        let newCopyBytes = try newLayout.stateCopyBytes()
        let oldElements = oldLayout.authoritativeElementsPerChunk
        let newElements = newLayout.authoritativeElementsPerChunk
        let commonElements = min(oldElements, newElements)
        guard commonElements <= UInt64(UInt32.max) else {
            throw VivoAdaptiveArenaError.arithmeticOverflow(newLayout.cohortID)
        }

        if oldLayout.representation == newLayout.representation ||
            migration.kind == .none || migration.kind == .resizeAndCopy ||
            migration.kind == .reinterpretInPlace ||
            migration.kind == .rebuildSpatialTopology ||
            migration.kind == .rebuildTissueCoupling {
            guard let blit = commandBuffer.makeBlitCommandEncoder() else {
                throw VivoAdaptiveArenaError.commandQueueUnavailable
            }
            blit.label = "AdaptiveArena.StateCopy.\(newLayout.cohortID)"
            let bytesPerCopy = min(oldCopyBytes, newCopyBytes)
            for copy in 0..<copyCount {
                let sourceOffset = try source.descriptor.stateCopyOffset(copy)
                let destinationOffset = try destination.descriptor.stateCopyOffset(copy)
                blit.copy(
                    from: source.state,
                    sourceOffset: sourceOffset,
                    to: destination.state,
                    destinationOffset: destinationOffset,
                    size: Int(bytesPerCopy)
                )
            }
            blit.endEncoding()
            return bytesPerCopy * UInt64(copyCount)
        }

        let pipeline: MTLComputePipelineState
        if oldLayout.representation == .continuousFP32,
           newLayout.representation == .discreteUInt32,
           migration.kind == .continuousToDiscrete {
            pipeline = pipelines.floatToCount
        } else if oldLayout.representation == .discreteUInt32,
                  newLayout.representation == .continuousFP32,
                  migration.kind == .discreteToContinuous {
            pipeline = pipelines.countToFloat
        } else {
            throw VivoAdaptiveArenaError.migrationUnsupported(
                newLayout.cohortID,
                oldLayout.representation,
                newLayout.representation
            )
        }

        var migrationCommand = VivoMigrationCommandABI(
            elementCount: UInt32(commonElements),
            integerTolerance: newLayout.integerMigrationTolerance
        )
        for copy in 0..<copyCount {
            let sourceOffset = try source.descriptor.stateCopyOffset(copy)
            let destinationOffset = try destination.descriptor.stateCopyOffset(copy)
            try encode(
                pipeline: pipeline,
                commandBuffer: commandBuffer,
                count: Int(commonElements),
                label: "AdaptiveArena.StateConvert.\(newLayout.cohortID).\(copy)"
            ) { encoder in
                encoder.setBuffer(source.state, offset: sourceOffset, index: 0)
                encoder.setBuffer(destination.state, offset: destinationOffset, index: 1)
                encoder.setBytes(&migrationCommand, length: MemoryLayout<VivoMigrationCommandABI>.stride, index: 2)
                encoder.setBuffer(migrationStatus, offset: 0, index: 3)
            }
        }
        return commonElements * 4 * UInt64(copyCount)
    }

    private func migrateScratch(
        from source: VivoAdaptiveArenaChunkResources,
        oldLayout: VivoAdaptiveCohortRuntimeLayout,
        to destination: VivoAdaptiveArenaChunkResources,
        newLayout: VivoAdaptiveCohortRuntimeLayout,
        migration: VivoBufferMigrationPlan,
        commandBuffer: MTLCommandBuffer
    ) throws -> (random: UInt64, delayed: UInt64) {
        guard let blit = commandBuffer.makeBlitCommandEncoder() else {
            throw VivoAdaptiveArenaError.commandQueueUnavailable
        }
        blit.label = "AdaptiveArena.ScratchMigration.\(newLayout.cohortID)"
        var randomCopied: UInt64 = 0
        var delayedCopied: UInt64 = 0

        if migration.preserveRandomCounter {
            randomCopied = min(oldLayout.randomCounterBytes, newLayout.randomCounterBytes)
            if randomCopied > 0 {
                blit.copy(
                    from: source.scratch,
                    sourceOffset: Int(source.descriptor.scratchOffsets.randomCounterOffset),
                    to: destination.scratch,
                    destinationOffset: Int(destination.descriptor.scratchOffsets.randomCounterOffset),
                    size: Int(randomCopied)
                )
            }
        }
        if migration.preserveDelayedQueue {
            delayedCopied = min(oldLayout.delayedQueueBytes, newLayout.delayedQueueBytes)
            if delayedCopied > 0 {
                blit.copy(
                    from: source.scratch,
                    sourceOffset: Int(source.descriptor.scratchOffsets.delayedQueueOffset),
                    to: destination.scratch,
                    destinationOffset: Int(destination.descriptor.scratchOffsets.delayedQueueOffset),
                    size: Int(delayedCopied)
                )
            }
        }

        let rebuildTopology = migration.kind == .rebuildSpatialTopology ||
                              migration.kind == .rebuildTissueCoupling
        if !rebuildTopology {
            let topologyBytes = min(oldLayout.topologyBytes, newLayout.topologyBytes)
            if topologyBytes > 0 {
                blit.copy(
                    from: source.scratch,
                    sourceOffset: Int(source.descriptor.scratchOffsets.topologyOffset),
                    to: destination.scratch,
                    destinationOffset: Int(destination.descriptor.scratchOffsets.topologyOffset),
                    size: Int(topologyBytes)
                )
            }
        }
        let auxiliaryBytes = min(oldLayout.auxiliaryBytes, newLayout.auxiliaryBytes)
        if auxiliaryBytes > 0, !rebuildTopology {
            blit.copy(
                from: source.scratch,
                sourceOffset: Int(source.descriptor.scratchOffsets.auxiliaryOffset),
                to: destination.scratch,
                destinationOffset: Int(destination.descriptor.scratchOffsets.auxiliaryOffset),
                size: Int(auxiliaryBytes)
            )
        }
        blit.endEncoding()
        return (randomCopied, delayedCopied)
    }

    private func upload(
        _ data: Data,
        destination: MTLBuffer,
        destinationOffset: Int,
        clearRange: Range<Int>?,
        label: String
    ) throws {
        guard destinationOffset >= 0,
              destinationOffset <= destination.length,
              data.count <= destination.length - destinationOffset else {
            throw VivoAdaptiveArenaError.invalidPlan("upload range exceeds destination buffer")
        }
        guard let staging = device.makeBuffer(length: max(1, data.count), options: .storageModeShared),
              let command = queue.makeCommandBuffer(),
              let blit = command.makeBlitCommandEncoder() else {
            throw VivoAdaptiveArenaError.commandQueueUnavailable
        }
        command.label = "NumiVivo.AdaptiveArena.\(label)"
        if !data.isEmpty {
            data.copyBytes(to: staging.contents().assumingMemoryBound(to: UInt8.self), count: data.count)
            blit.copy(
                from: staging,
                sourceOffset: 0,
                to: destination,
                destinationOffset: destinationOffset,
                size: data.count
            )
        }
        if let clearRange, !clearRange.isEmpty {
            blit.fill(buffer: destination, range: clearRange, value: 0)
        }
        blit.endEncoding()
        command.commit()
        command.waitUntilCompleted()
        guard command.status == .completed else {
            throw VivoAdaptiveArenaError.commandFailed(
                command.error?.localizedDescription ?? "upload command failed"
            )
        }
    }

    private func scratchRegion(
        _ region: VivoAdaptiveScratchRegion,
        layout: VivoAdaptiveCohortRuntimeLayout,
        offsets: VivoAdaptiveScratchOffsets
    ) throws -> (offset: Int, length: Int) {
        let raw: (UInt64, UInt64)
        switch region {
        case .randomCounters: raw = (offsets.randomCounterOffset, layout.randomCounterBytes)
        case .delayedQueue: raw = (offsets.delayedQueueOffset, layout.delayedQueueBytes)
        case .topology: raw = (offsets.topologyOffset, layout.topologyBytes)
        case .auxiliary: raw = (offsets.auxiliaryOffset, layout.auxiliaryBytes)
        }
        guard raw.0 <= UInt64(Int.max), raw.1 <= UInt64(Int.max) else {
            throw VivoAdaptiveArenaError.arithmeticOverflow(layout.cohortID)
        }
        return (Int(raw.0), Int(raw.1))
    }

    private func encode(
        pipeline: MTLComputePipelineState,
        commandBuffer: MTLCommandBuffer,
        count: Int,
        label: String,
        bindings: (MTLComputeCommandEncoder) -> Void
    ) throws {
        guard count > 0, let encoder = commandBuffer.makeComputeCommandEncoder() else {
            throw VivoAdaptiveArenaError.commandQueueUnavailable
        }
        encoder.label = label
        encoder.setComputePipelineState(pipeline)
        bindings(encoder)
        let width = Self.threadgroupWidth(pipeline)
        encoder.dispatchThreads(
            MTLSize(width: count, height: 1, depth: 1),
            threadsPerThreadgroup: MTLSize(width: width, height: 1, depth: 1)
        )
        encoder.endEncoding()
    }

    private func complete(_ command: MTLCommandBuffer) async throws {
        try await withCheckedThrowingContinuation { continuation in
            command.addCompletedHandler { completed in
                if completed.status == .completed {
                    continuation.resume()
                } else {
                    continuation.resume(throwing: VivoAdaptiveArenaError.commandFailed(
                        completed.error?.localizedDescription ?? "status \(completed.status.rawValue)"
                    ))
                }
            }
            command.commit()
        }
    }

    private static func validate(
        plan: VivoAdaptiveFidelityPlan,
        layouts: [VivoAdaptiveCohortRuntimeLayout]
    ) throws -> [UInt32: VivoAdaptiveCohortRuntimeLayout] {
        guard plan.schemaVersion == 1,
              !plan.fingerprint.isEmpty,
              plan.cohorts.map(\.cohortID).sorted() == Array(Set(plan.cohorts.map(\.cohortID))).sorted() else {
            throw VivoAdaptiveArenaError.invalidPlan("schema, fingerprint or cohort identity")
        }
        var result: [UInt32: VivoAdaptiveCohortRuntimeLayout] = [:]
        for layout in layouts {
            guard result[layout.cohortID] == nil else {
                throw VivoAdaptiveArenaError.duplicateLayout(layout.cohortID)
            }
            result[layout.cohortID] = layout
        }
        for cohort in plan.cohorts {
            guard cohort.chunkCount > 0,
                  cohort.stateBytesPerChunk > 0,
                  cohort.scratchBytesPerChunk > 0 else {
                throw VivoAdaptiveArenaError.invalidPlan("cohort \(cohort.cohortID) has empty allocation")
            }
            guard let layout = result[cohort.cohortID] else {
                throw VivoAdaptiveArenaError.missingLayout(cohort.cohortID)
            }
            try layout.validate(against: cohort)
        }
        guard result.count == plan.cohorts.count else {
            let known = Set(plan.cohorts.map(\.cohortID))
            if let extra = result.keys.first(where: { !known.contains($0) }) {
                throw VivoAdaptiveArenaError.invalidLayout(extra, "layout has no matching cohort plan")
            }
            throw VivoAdaptiveArenaError.invalidPlan("layout count mismatch")
        }
        guard let actual = plan.actualChunkedArenaBytes else {
            throw VivoAdaptiveArenaError.invalidPlan("chunked byte accounting overflow")
        }
        guard actual <= plan.workingSetBudgetBytes else {
            throw VivoAdaptiveArenaError.workingSetExceeded(
                required: actual,
                available: plan.workingSetBudgetBytes
            )
        }
        return result
    }

    private static func allocate(
        plan: VivoAdaptiveFidelityPlan,
        layouts: [UInt32: VivoAdaptiveCohortRuntimeLayout],
        device: MTLDevice,
        queue: MTLCommandQueue
    ) throws -> VivoAdaptiveArenaStorage {
        struct Request {
            let key: VivoAdaptiveChunkKey
            let plan: VivoFidelityCohortPlan
            let layout: VivoAdaptiveCohortRuntimeLayout
            let stateSizeAndAlign: MTLSizeAndAlign
            let scratchSizeAndAlign: MTLSizeAndAlign
        }
        var requests: [Request] = []
        var heapSize: UInt64 = 0
        var allocatedState: UInt64 = 0
        var allocatedScratch: UInt64 = 0

        for cohort in plan.cohorts.sorted(by: { $0.cohortID < $1.cohortID }) {
            guard let layout = layouts[cohort.cohortID] else {
                throw VivoAdaptiveArenaError.missingLayout(cohort.cohortID)
            }
            guard cohort.stateBytesPerChunk <= UInt64(Int.max),
                  cohort.scratchBytesPerChunk <= UInt64(Int.max) else {
                throw VivoAdaptiveArenaError.arithmeticOverflow(cohort.cohortID)
            }
            let stateRequirements = device.heapBufferSizeAndAlign(
                length: Int(cohort.stateBytesPerChunk),
                options: .storageModePrivate
            )
            let scratchRequirements = device.heapBufferSizeAndAlign(
                length: Int(cohort.scratchBytesPerChunk),
                options: .storageModePrivate
            )
            for chunkIndex in 0..<cohort.chunkCount {
                heapSize = try addHeapResource(
                    current: heapSize,
                    requirements: stateRequirements,
                    cohortID: cohort.cohortID
                )
                heapSize = try addHeapResource(
                    current: heapSize,
                    requirements: scratchRequirements,
                    cohortID: cohort.cohortID
                )
                let stateTotal = allocatedState.addingReportingOverflow(cohort.stateBytesPerChunk)
                let scratchTotal = allocatedScratch.addingReportingOverflow(cohort.scratchBytesPerChunk)
                guard !stateTotal.overflow, !scratchTotal.overflow else {
                    throw VivoAdaptiveArenaError.arithmeticOverflow(cohort.cohortID)
                }
                allocatedState = stateTotal.partialValue
                allocatedScratch = scratchTotal.partialValue
                requests.append(.init(
                    key: .init(cohortID: cohort.cohortID, chunkIndex: chunkIndex),
                    plan: cohort,
                    layout: layout,
                    stateSizeAndAlign: stateRequirements,
                    scratchSizeAndAlign: scratchRequirements
                ))
            }
        }

        guard heapSize <= plan.workingSetBudgetBytes else {
            throw VivoAdaptiveArenaError.workingSetExceeded(
                required: heapSize,
                available: plan.workingSetBudgetBytes
            )
        }
        guard heapSize <= UInt64(Int.max) else {
            throw VivoAdaptiveArenaError.invalidPlan("heap exceeds host address space")
        }
        let descriptor = MTLHeapDescriptor()
        descriptor.storageMode = .private
        descriptor.hazardTrackingMode = .tracked
        descriptor.size = Int(heapSize)
        guard let heap = device.makeHeap(descriptor: descriptor) else {
            throw VivoAdaptiveArenaError.heapCreationFailed
        }
        heap.label = "NumiVivo.AdaptiveArena.\(plan.fingerprint.prefix(12))"

        var chunks: [VivoAdaptiveChunkKey: VivoAdaptiveArenaChunkResources] = [:]
        var chunkDescriptors: [VivoAdaptiveChunkDescriptor] = []
        for request in requests {
            guard let state = heap.makeBuffer(
                length: Int(request.plan.stateBytesPerChunk),
                options: .storageModePrivate
            ) else {
                throw VivoAdaptiveArenaError.allocationFailed(
                    request.key.cohortID,
                    request.key.chunkIndex,
                    "state"
                )
            }
            guard let scratch = heap.makeBuffer(
                length: Int(request.plan.scratchBytesPerChunk),
                options: .storageModePrivate
            ) else {
                throw VivoAdaptiveArenaError.allocationFailed(
                    request.key.cohortID,
                    request.key.chunkIndex,
                    "scratch"
                )
            }
            state.label = "NumiVivo.C\(request.key.cohortID).K\(request.key.chunkIndex).State"
            scratch.label = "NumiVivo.C\(request.key.cohortID).K\(request.key.chunkIndex).Scratch"
            let chunkDescriptor = VivoAdaptiveChunkDescriptor(
                cohortID: request.key.cohortID,
                chunkIndex: request.key.chunkIndex,
                stateBytes: request.plan.stateBytesPerChunk,
                scratchBytes: request.plan.scratchBytesPerChunk,
                stateCopyBytes: try request.layout.stateCopyBytes(),
                stateCopyCount: request.layout.stateCopyCount,
                representation: request.layout.representation,
                scratchOffsets: try VivoAdaptiveScratchOffsets(layout: request.layout)
            )
            chunkDescriptors.append(chunkDescriptor)
            chunks[request.key] = .init(
                descriptor: chunkDescriptor,
                state: state,
                scratch: scratch
            )
        }

        // Fresh and rebuilt regions start from a deterministic zero state.
        guard let command = queue.makeCommandBuffer(),
              let blit = command.makeBlitCommandEncoder() else {
            throw VivoAdaptiveArenaError.commandQueueUnavailable
        }
        command.label = "NumiVivo.AdaptiveArena.Initialize"
        for resources in chunks.values {
            blit.fill(buffer: resources.state, range: 0..<resources.state.length, value: 0)
            blit.fill(buffer: resources.scratch, range: 0..<resources.scratch.length, value: 0)
        }
        blit.endEncoding()
        command.commit()
        command.waitUntilCompleted()
        guard command.status == .completed else {
            throw VivoAdaptiveArenaError.commandFailed(
                command.error?.localizedDescription ?? "arena initialization failed"
            )
        }

        let manifest = try VivoAdaptiveArenaManifest(
            planFingerprint: plan.fingerprint,
            deviceName: device.name,
            registryID: device.registryID,
            heapBytes: heapSize,
            allocatedStateBytes: allocatedState,
            allocatedScratchBytes: allocatedScratch,
            chunkCount: UInt32(requests.count),
            chunks: chunkDescriptors.sorted {
                if $0.cohortID != $1.cohortID { return $0.cohortID < $1.cohortID }
                return $0.chunkIndex < $1.chunkIndex
            }
        )
        return .init(heap: heap, chunks: chunks, manifest: manifest)
    }

    private static func addHeapResource(
        current: UInt64,
        requirements: MTLSizeAndAlign,
        cohortID: UInt32
    ) throws -> UInt64 {
        let alignment = UInt64(requirements.align)
        guard alignment > 0, alignment.nonzeroBitCount == 1 else {
            throw VivoAdaptiveArenaError.invalidAlignment
        }
        let mask = alignment - 1
        let alignedBase = current.addingReportingOverflow(mask)
        guard !alignedBase.overflow else { throw VivoAdaptiveArenaError.arithmeticOverflow(cohortID) }
        let aligned = alignedBase.partialValue & ~mask
        let result = aligned.addingReportingOverflow(UInt64(requirements.size))
        guard !result.overflow else { throw VivoAdaptiveArenaError.arithmeticOverflow(cohortID) }
        return result.partialValue
    }

    private static func loadPipelines(device: MTLDevice) throws -> VivoMigrationPipelines {
        guard let sourceURL = NumiVivoShaderResources.bundle.url(
            forResource: "NumiVivoMigrationKernels",
            withExtension: "metal"
        ) else {
            throw VivoAdaptiveArenaError.shaderResourceUnavailable
        }
        let source = String(decoding: try Data(contentsOf: sourceURL), as: UTF8.self)
        let library: MTLLibrary
        do { library = try device.makeLibrary(source: source, options: nil) }
        catch { throw VivoAdaptiveArenaError.shaderCompilation(error.localizedDescription) }

        func pipeline(_ name: String) throws -> MTLComputePipelineState {
            guard let function = library.makeFunction(name: name) else {
                throw VivoAdaptiveArenaError.pipelineUnavailable(name)
            }
            do { return try device.makeComputePipelineState(function: function) }
            catch { throw VivoAdaptiveArenaError.shaderCompilation("\(name): \(error.localizedDescription)") }
        }
        return .init(
            clear: try pipeline("nvivo_migration_clear_status"),
            floatToCount: try pipeline("nvivo_migration_f32_to_u32"),
            countToFloat: try pipeline("nvivo_migration_u32_to_f32"),
            validateFloat: try pipeline("nvivo_migration_validate_f32")
        )
    }

    private static func threadgroupWidth(_ pipeline: MTLComputePipelineState) -> Int {
        let width = max(1, pipeline.threadExecutionWidth)
        let groups = max(1, min(8, pipeline.maxTotalThreadsPerThreadgroup / width))
        return width * groups
    }
}
