import Foundation
import Metal

public enum VivoCapacityMemoryClass: String, Codable, Sendable {
    case privateGPU
    case sharedCPUAndGPU
    case immutableProgram
}

public struct VivoCapacityComponent: Codable, Sendable, Equatable {
    public let name: String
    public let memoryClass: VivoCapacityMemoryClass
    public let bytes: UInt64
    public let explanation: String
}

public struct VivoCapacityRequest: Codable, Sendable, Equatable {
    public var molecularLaneCount: UInt32
    public var parameterEnvironmentCount: UInt32
    public var molecularEventCapacity: UInt32
    public var molecularCouplingCapacity: UInt32
    public var molecularPublicationCapacity: UInt32
    public var physiologyEnvironmentCount: UInt32
    public var physiologyTransformCapacity: UInt32
    public var physiologyPublicationCapacity: UInt32
    public var workingSetFraction: Double
    public var allocationHeadroom: Double

    public init(
        molecularLaneCount: UInt32,
        parameterEnvironmentCount: UInt32 = 1,
        molecularEventCapacity: UInt32 = 16_384,
        molecularCouplingCapacity: UInt32 = 65_536,
        molecularPublicationCapacity: UInt32 = 65_536,
        physiologyEnvironmentCount: UInt32 = 0,
        physiologyTransformCapacity: UInt32 = 262_144,
        physiologyPublicationCapacity: UInt32 = 65_536,
        workingSetFraction: Double = 0.80,
        allocationHeadroom: Double = 1.15
    ) {
        self.molecularLaneCount = molecularLaneCount
        self.parameterEnvironmentCount = parameterEnvironmentCount
        self.molecularEventCapacity = molecularEventCapacity
        self.molecularCouplingCapacity = molecularCouplingCapacity
        self.molecularPublicationCapacity = molecularPublicationCapacity
        self.physiologyEnvironmentCount = physiologyEnvironmentCount
        self.physiologyTransformCapacity = physiologyTransformCapacity
        self.physiologyPublicationCapacity = physiologyPublicationCapacity
        self.workingSetFraction = workingSetFraction
        self.allocationHeadroom = allocationHeadroom
    }

    public func validate(physiology: PreparedVivoPhysiologyModel?) throws {
        guard molecularLaneCount > 0,
              parameterEnvironmentCount > 0,
              molecularEventCapacity > 0,
              molecularCouplingCapacity > 0,
              molecularPublicationCapacity > 0,
              workingSetFraction.isFinite,
              workingSetFraction > 0,
              workingSetFraction <= 1,
              allocationHeadroom.isFinite,
              allocationHeadroom >= 1 else {
            throw VivoArtifactValidationError.invalid(
                "capacity request contains invalid molecular counts or memory policy"
            )
        }
        if let physiology {
            guard physiologyEnvironmentCount == physiology.environmentCount,
                  physiologyTransformCapacity > 0,
                  physiologyPublicationCapacity > 0 else {
                throw VivoArtifactValidationError.incompatible(
                    "capacity request physiology counts do not match the prepared model"
                )
            }
        } else if physiologyEnvironmentCount != 0 {
            throw VivoArtifactValidationError.invalid(
                "capacity request declares physiology environments without a physiology model"
            )
        }
    }
}

public struct VivoCapacityLimits: Codable, Sendable, Equatable {
    public let maximumMolecularLaneCountKeepingPhysiologyFixed: UInt32
    public let maximumPhysiologyEnvironmentCountKeepingMolecularFixed: UInt32
    public let maximumLockstepScale: UInt32
    public let maximumLockstepMolecularLaneCount: UInt32
    public let maximumLockstepPhysiologyEnvironmentCount: UInt32
}

public struct VivoCapacityPlan: Codable, Sendable, Equatable {
    public let deviceName: String
    public let registryID: UInt64
    public let unifiedMemory: Bool
    public let recommendedWorkingSetBytes: UInt64
    public let currentAllocatedBytes: UInt64
    public let planningBudgetBytes: UInt64
    public let requestedBytesBeforeHeadroom: UInt64
    public let requestedBytesWithHeadroom: UInt64
    public let largestSingleBufferBytes: UInt64
    public let utilizationOfPlanningBudget: Double
    public let fitsPlanningBudget: Bool
    public let fitsMaximumBufferLength: Bool
    public let components: [VivoCapacityComponent]
    public let limits: VivoCapacityLimits
    public let warnings: [String]

    public func fingerprint() throws -> VivoFingerprint {
        try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(self))
    }
}

public struct VivoAppleSiliconCapacityPlanner: Sendable {
    private struct PrivateBufferSpec {
        let name: String
        let length: UInt64
    }

    public init() {}

    public func plan(
        programPack: VivoProgramPack,
        request: VivoCapacityRequest,
        physiology: PreparedVivoPhysiologyModel? = nil,
        device: MTLDevice
    ) throws -> VivoCapacityPlan {
        try request.validate(physiology: physiology)
        let capabilities = VivoMetalCapabilities(device: device)
        let estimate = try estimate(
            programPack: programPack,
            request: request,
            physiology: physiology,
            device: device
        )

        let recommended = device.recommendedMaxWorkingSetSize
        let selectedBudget = try scaled(
            recommended,
            by: request.workingSetFraction,
            label: "working-set fraction"
        )
        let currentAllocated = device.currentAllocatedSize
        let planningBudget = selectedBudget > currentAllocated
            ? selectedBudget - currentAllocated
            : 0
        let fitsBudget = estimate.withHeadroom <= planningBudget
        let fitsBuffer = estimate.largestBuffer <= UInt64(device.maxBufferLength)
        let utilization = planningBudget == 0
            ? Double.greatestFiniteMagnitude
            : Double(estimate.withHeadroom) / Double(planningBudget)

        let limits = try solveLimits(
            programPack: programPack,
            request: request,
            physiology: physiology,
            device: device,
            budget: planningBudget
        )

        var warnings: [String] = []
        if !capabilities.hasUnifiedMemory {
            warnings.append(
                "Selected Metal device does not expose unified memory; production NumiVivo targets Apple silicon unified memory."
            )
        }
        if planningBudget == 0 {
            warnings.append(
                "No planning budget remains after current Metal allocations and the selected working-set fraction."
            )
        }
        if !fitsBudget {
            warnings.append(
                "Requested private heaps and shared transaction buffers exceed the selected planning budget."
            )
        }
        if !fitsBuffer {
            warnings.append(
                "At least one required buffer exceeds device.maxBufferLength."
            )
        }
        if utilization.isFinite, utilization > 0.90 {
            warnings.append(
                "The plan consumes more than 90% of its selected budget; reserve additional space for command buffers, shader libraries, coupled simulators, and operating-system pressure."
            )
        }
        if request.parameterEnvironmentCount != 1,
           request.parameterEnvironmentCount != request.molecularLaneCount {
            warnings.append(
                "Parameter environments are neither global nor one-per-molecular-lane; confirm the runtime parameter mapping."
            )
        }

        return VivoCapacityPlan(
            deviceName: device.name,
            registryID: device.registryID,
            unifiedMemory: device.hasUnifiedMemory,
            recommendedWorkingSetBytes: recommended,
            currentAllocatedBytes: currentAllocated,
            planningBudgetBytes: planningBudget,
            requestedBytesBeforeHeadroom: estimate.raw,
            requestedBytesWithHeadroom: estimate.withHeadroom,
            largestSingleBufferBytes: estimate.largestBuffer,
            utilizationOfPlanningBudget: utilization,
            fitsPlanningBudget: fitsBudget,
            fitsMaximumBufferLength: fitsBuffer,
            components: estimate.components,
            limits: limits,
            warnings: warnings
        )
    }

    private func estimate(
        programPack: VivoProgramPack,
        request: VivoCapacityRequest,
        physiology: PreparedVivoPhysiologyModel?,
        device: MTLDevice
    ) throws -> (
        raw: UInt64,
        withHeadroom: UInt64,
        largestBuffer: UInt64,
        components: [VivoCapacityComponent]
    ) {
        let contract = programPack.runtimeContract
        let lanes = UInt64(request.molecularLaneCount)
        let species = UInt64(contract.speciesCount)
        let reactions = UInt64(contract.reactionCount)
        let temporal = UInt64(contract.temporalStateCount)
        let parameters = UInt64(contract.parameterCount)
        let parameterEnvironments = UInt64(request.parameterEnvironmentCount)

        let stateElements = try product(species, lanes, label: "molecular state elements")
        let stateBytes = try product(stateElements, 4, label: "molecular state bytes")
        let temporalElements = max(
            try product(temporal, lanes, label: "molecular temporal elements"),
            1
        )
        let temporalBytes = try product(temporalElements, 8, label: "molecular temporal bytes")
        let reactionElements = max(
            try product(reactions, lanes, label: "molecular reaction elements"),
            1
        )
        let reactionBytes = try product(reactionElements, 4, label: "molecular reaction bytes")
        let parameterElements = max(
            try product(parameters, parameterEnvironments, label: "molecular parameter elements"),
            1
        )
        let parameterBytes = try product(parameterElements, 4, label: "molecular parameter bytes")
        let transportBytes = try product(max(species, 1), 16, label: "transport bytes")
        let velocityBytes = try product(max(lanes, 1), 16, label: "velocity bytes")
        let volumeBytes = try product(max(lanes, 1), 4, label: "volume-fraction bytes")
        let programBytes = UInt64(max(programPack.data.count, 1))

        var molecularPrivate = [
            PrivateBufferSpec(name: "ProgramPack", length: programBytes),
            PrivateBufferSpec(name: "parameters", length: parameterBytes),
            PrivateBufferSpec(name: "current-state", length: stateBytes),
            PrivateBufferSpec(name: "base-state", length: stateBytes),
            PrivateBufferSpec(name: "stage-state", length: stateBytes),
            PrivateBufferSpec(name: "candidate-state", length: stateBytes),
            PrivateBufferSpec(name: "derivative-k1", length: stateBytes),
            PrivateBufferSpec(name: "temporal-current", length: temporalBytes),
            PrivateBufferSpec(name: "temporal-candidate", length: temporalBytes),
            PrivateBufferSpec(name: "reaction-events", length: reactionBytes),
            PrivateBufferSpec(name: "species-transport", length: transportBytes),
            PrivateBufferSpec(name: "velocity", length: velocityBytes),
            PrivateBufferSpec(name: "volume-fraction", length: volumeBytes)
        ]
        let molecularPrivateHeap = try heapBytes(
            device: device,
            buffers: molecularPrivate,
            headroom: request.allocationHeadroom
        )

        let couplingCapacity = max(
            UInt64(request.molecularCouplingCapacity),
            lanes,
            1
        )
        let molecularShared = try sum([
            UInt64(MemoryLayout<VivoRuntimeCommandABI>.stride),
            UInt64(MemoryLayout<VivoRuntimeStatusABI>.stride),
            try product(UInt64(request.molecularEventCapacity), UInt64(MemoryLayout<VivoEventABI>.stride), label: "event buffer"),
            try product(couplingCapacity, UInt64(MemoryLayout<VivoCouplingUpdateABI>.stride), label: "coupling buffer"),
            try product(UInt64(request.molecularPublicationCapacity), UInt64(MemoryLayout<VivoPublicationRequestABI>.stride), label: "publication request buffer"),
            try product(UInt64(request.molecularPublicationCapacity), 4, label: "publication output buffer"),
            stateBytes
        ], label: "molecular shared buffers")

        var components: [VivoCapacityComponent] = [
            VivoCapacityComponent(
                name: "Molecular private heap",
                memoryClass: .privateGPU,
                bytes: molecularPrivateHeap,
                explanation: "Heap-aligned ProgramPack, parameters, authoritative/candidate state, temporal state, reaction events, and spatial fields."
            ),
            VivoCapacityComponent(
                name: "Molecular shared boundary",
                memoryClass: .sharedCPUAndGPU,
                bytes: molecularShared,
                explanation: "Commands, status, events, coupling updates, publications, and explicit state readback."
            )
        ]
        var raw = try sum([molecularPrivateHeap, molecularShared], label: "molecular runtime")
        var largest = molecularPrivate.map(\.length).max() ?? 0
        largest = max(
            largest,
            try product(UInt64(request.molecularEventCapacity), UInt64(MemoryLayout<VivoEventABI>.stride), label: "event buffer")
        )
        largest = max(largest, stateBytes)

        if let physiology {
            let environments = UInt64(request.physiologyEnvironmentCount)
            let pairCount = UInt64(physiology.pairCount)
            let stateElements = try product(pairCount, environments, label: "physiology state elements")
            let physiologyStateBytes = try product(stateElements, 4, label: "physiology state bytes")
            let offsetsBytes = try product(UInt64(physiology.incidenceOffsets.count), 4, label: "physiology incidence offsets")
            let incidenceBytes = try product(
                UInt64(max(physiology.incidence.count, 1)),
                UInt64(MemoryLayout<VivoPhysiologyIncidenceABI>.stride),
                label: "physiology incidence"
            )
            let clearanceBytes = try product(
                UInt64(physiology.clearances.count),
                UInt64(MemoryLayout<VivoPhysiologyClearanceABI>.stride),
                label: "physiology clearance"
            )
            let boundsBytes = try product(pairCount, UInt64(MemoryLayout<SIMD2<Float>>.stride), label: "physiology bounds")
            let physiologyPrivate = [
                PrivateBufferSpec(name: "current-state", length: physiologyStateBytes),
                PrivateBufferSpec(name: "base-state", length: physiologyStateBytes),
                PrivateBufferSpec(name: "stage-state", length: physiologyStateBytes),
                PrivateBufferSpec(name: "candidate-state", length: physiologyStateBytes),
                PrivateBufferSpec(name: "derivative-k1", length: physiologyStateBytes),
                PrivateBufferSpec(name: "incidence-offsets", length: offsetsBytes),
                PrivateBufferSpec(name: "incidence", length: incidenceBytes),
                PrivateBufferSpec(name: "clearances", length: clearanceBytes),
                PrivateBufferSpec(name: "bounds", length: boundsBytes)
            ]
            let physiologyHeap = try heapBytes(
                device: device,
                buffers: physiologyPrivate,
                headroom: request.allocationHeadroom
            )
            let transformBytes = try product(
                UInt64(request.physiologyTransformCapacity),
                UInt64(MemoryLayout<VivoPhysiologyStateTransformABI>.stride),
                label: "physiology transforms"
            )
            let physiologyShared = try sum([
                UInt64(MemoryLayout<VivoPhysiologyRuntimeCommandABI>.stride),
                UInt64(MemoryLayout<VivoPhysiologyRuntimeStatusABI>.stride),
                transformBytes,
                transformBytes,
                try product(UInt64(request.physiologyPublicationCapacity), UInt64(MemoryLayout<VivoPhysiologyPublicationRequestABI>.stride), label: "physiology publication requests"),
                try product(UInt64(request.physiologyPublicationCapacity), 4, label: "physiology publication output"),
                physiologyStateBytes
            ], label: "physiology shared buffers")
            raw = try sum([raw, physiologyHeap, physiologyShared], label: "coupled runtime")
            largest = max(largest, physiologyPrivate.map(\.length).max() ?? 0, transformBytes)
            components.append(
                VivoCapacityComponent(
                    name: "Physiology private heap",
                    memoryClass: .privateGPU,
                    bytes: physiologyHeap,
                    explanation: "Heap-aligned compartment state, RK2 staging, sparse incidence, clearance, and bound tables."
                )
            )
            components.append(
                VivoCapacityComponent(
                    name: "Physiology shared boundary",
                    memoryClass: .sharedCPUAndGPU,
                    bytes: physiologyShared,
                    explanation: "Pre/post transforms, publications, command/status records, and explicit state readback."
                )
            )
        }

        return (
            raw: raw,
            withHeadroom: raw,
            largestBuffer: largest,
            components: components
        )
    }

    private func solveLimits(
        programPack: VivoProgramPack,
        request: VivoCapacityRequest,
        physiology: PreparedVivoPhysiologyModel?,
        device: MTLDevice,
        budget: UInt64
    ) throws -> VivoCapacityLimits {
        let maximumMolecular = try maximumCount { count in
            var candidate = request
            candidate.molecularLaneCount = count
            if request.parameterEnvironmentCount == request.molecularLaneCount {
                candidate.parameterEnvironmentCount = count
            }
            return try fits(
                programPack: programPack,
                request: candidate,
                physiology: physiology,
                device: device,
                budget: budget
            )
        }

        let maximumPhysiology: UInt32
        if physiology != nil {
            maximumPhysiology = try maximumCount { count in
                var candidate = request
                candidate.physiologyEnvironmentCount = count
                return try fits(
                    programPack: programPack,
                    request: candidate,
                    physiology: physiology,
                    device: device,
                    budget: budget,
                    skipModelEnvironmentValidation: true
                )
            }
        } else {
            maximumPhysiology = 0
        }

        let lockstep = try maximumCount { scale in
            var candidate = request
            candidate.molecularLaneCount = saturatedMultiply(
                request.molecularLaneCount,
                scale
            )
            if request.parameterEnvironmentCount == request.molecularLaneCount {
                candidate.parameterEnvironmentCount = candidate.molecularLaneCount
            }
            if physiology != nil {
                candidate.physiologyEnvironmentCount = saturatedMultiply(
                    request.physiologyEnvironmentCount,
                    scale
                )
            }
            return try fits(
                programPack: programPack,
                request: candidate,
                physiology: physiology,
                device: device,
                budget: budget,
                skipModelEnvironmentValidation: physiology != nil && scale != 1
            )
        }

        return VivoCapacityLimits(
            maximumMolecularLaneCountKeepingPhysiologyFixed: maximumMolecular,
            maximumPhysiologyEnvironmentCountKeepingMolecularFixed: maximumPhysiology,
            maximumLockstepScale: lockstep,
            maximumLockstepMolecularLaneCount: saturatedMultiply(
                request.molecularLaneCount,
                lockstep
            ),
            maximumLockstepPhysiologyEnvironmentCount: physiology == nil
                ? 0
                : saturatedMultiply(request.physiologyEnvironmentCount, lockstep)
        )
    }

    private func fits(
        programPack: VivoProgramPack,
        request: VivoCapacityRequest,
        physiology: PreparedVivoPhysiologyModel?,
        device: MTLDevice,
        budget: UInt64,
        skipModelEnvironmentValidation: Bool = false
    ) throws -> Bool {
        if !skipModelEnvironmentValidation {
            try request.validate(physiology: physiology)
        }
        guard request.molecularLaneCount > 0,
              request.parameterEnvironmentCount > 0 else { return false }
        let result = try estimate(
            programPack: programPack,
            request: request,
            physiology: physiology,
            device: device
        )
        return result.raw <= budget &&
            result.largestBuffer <= UInt64(device.maxBufferLength)
    }

    private func maximumCount(
        fits: (UInt32) throws -> Bool
    ) throws -> UInt32 {
        var low: UInt32 = 0
        var high: UInt32 = 1
        while high < UInt32.max, try fits(high) {
            low = high
            let doubled = high.multipliedReportingOverflow(by: 2)
            high = doubled.overflow ? UInt32.max : doubled.partialValue
            if high == low { break }
        }
        if high == UInt32.max, try fits(high) { return high }
        while UInt64(high) - UInt64(low) > 1 {
            let middle = UInt32(UInt64(low) + (UInt64(high) - UInt64(low)) / 2)
            if try fits(middle) {
                low = middle
            } else {
                high = middle
            }
        }
        return low
    }

    private func heapBytes(
        device: MTLDevice,
        buffers: [PrivateBufferSpec],
        headroom: Double
    ) throws -> UInt64 {
        var cursor: UInt64 = 0
        for buffer in buffers {
            guard buffer.length > 0,
                  buffer.length <= UInt64(Int.max) else {
                throw VivoArtifactValidationError.invalid(
                    "private buffer \(buffer.name) has an invalid length"
                )
            }
            let alignment = device.heapBufferSizeAndAlign(
                length: Int(buffer.length),
                options: .storageModePrivate
            )
            let align = UInt64(max(alignment.align, 1))
            cursor = try alignUp(cursor, alignment: align, label: buffer.name)
            let addition = cursor.addingReportingOverflow(UInt64(alignment.size))
            guard !addition.overflow else {
                throw VivoArtifactValidationError.invalid(
                    "private heap sizing overflow for \(buffer.name)"
                )
            }
            cursor = addition.partialValue
        }
        return try scaled(cursor, by: headroom, label: "private heap headroom")
    }

    private func alignUp(
        _ value: UInt64,
        alignment: UInt64,
        label: String
    ) throws -> UInt64 {
        guard alignment > 0,
              alignment & (alignment - 1) == 0 else {
            throw VivoArtifactValidationError.invalid(
                "Metal returned a non-power-of-two heap alignment for \(label)"
            )
        }
        let mask = alignment - 1
        let addition = value.addingReportingOverflow(mask)
        guard !addition.overflow else {
            throw VivoArtifactValidationError.invalid(
                "heap alignment overflow for \(label)"
            )
        }
        return addition.partialValue & ~mask
    }

    private func product(
        _ left: UInt64,
        _ right: UInt64,
        label: String
    ) throws -> UInt64 {
        let result = left.multipliedReportingOverflow(by: right)
        guard !result.overflow else {
            throw VivoArtifactValidationError.invalid("\(label) overflow")
        }
        return result.partialValue
    }

    private func sum(_ values: [UInt64], label: String) throws -> UInt64 {
        var total: UInt64 = 0
        for value in values {
            let result = total.addingReportingOverflow(value)
            guard !result.overflow else {
                throw VivoArtifactValidationError.invalid("\(label) overflow")
            }
            total = result.partialValue
        }
        return total
    }

    private func scaled(
        _ value: UInt64,
        by factor: Double,
        label: String
    ) throws -> UInt64 {
        let result = Double(value) * factor
        guard result.isFinite,
              result >= 0,
              result <= Double(UInt64.max) else {
            throw VivoArtifactValidationError.invalid("\(label) overflow")
        }
        return UInt64(result.rounded(.up))
    }

    private func saturatedMultiply(_ left: UInt32, _ right: UInt32) -> UInt32 {
        let result = left.multipliedReportingOverflow(by: right)
        return result.overflow ? UInt32.max : result.partialValue
    }
}
