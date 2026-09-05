import Foundation
@preconcurrency import Metal

public struct VivoMetalCapabilities: Sendable, Codable, Equatable {
    public enum Family: String, Sendable, Codable {
        case apple9OrNewer
        case apple8
        case apple7
        case apple6
        case apple5
        case apple4
        case apple3
        case apple2
        case apple1
        case mac2
        case mac1
        case unknown
    }

    public let deviceName: String
    public let registryID: UInt64
    public let family: Family
    public let hasUnifiedMemory: Bool
    public let isLowPower: Bool
    public let isRemovable: Bool
    public let supportsMetal3: Bool
    public let supportsNonUniformThreadgroups: Bool
    public let supportsSIMDGroupReduction: Bool
    public let supportsSIMDGroupMatrix: Bool
    public let supportsIndirectCommandBuffers: Bool
    public let recommendedMaximumWorkingSet: UInt64
    public let maximumBufferLength: UInt64
    public let defaultExecutionWidth: Int
    public let recommendedThreadsPerThreadgroup: Int

    public init(device: MTLDevice) {
        self.deviceName = device.name
        self.registryID = device.registryID
        self.family = Self.detectFamily(device)
        self.hasUnifiedMemory = device.hasUnifiedMemory
        self.isLowPower = device.isLowPower
        self.isRemovable = device.isRemovable
        self.supportsMetal3 = device.supportsFamily(.metal3)
        self.supportsNonUniformThreadgroups = Self.detectFamily(device) != .unknown
        self.supportsSIMDGroupReduction = Self.familyRank(Self.detectFamily(device)) >= Self.familyRank(.apple6)
        self.supportsSIMDGroupMatrix = Self.familyRank(Self.detectFamily(device)) >= Self.familyRank(.apple7)
        self.supportsIndirectCommandBuffers = Self.familyRank(Self.detectFamily(device)) >= Self.familyRank(.apple4)
        self.recommendedMaximumWorkingSet = device.recommendedMaxWorkingSetSize
        self.maximumBufferLength = UInt64(device.maxBufferLength)
        self.defaultExecutionWidth = 32
        self.recommendedThreadsPerThreadgroup = Self.recommendedThreads(for: Self.detectFamily(device))
    }

    public var isAppleSiliconOptimizedTarget: Bool {
        hasUnifiedMemory && family != .mac1 && family != .mac2 && family != .unknown
    }

    public func validate(programBytes: UInt64, privateHeapBytes: UInt64) throws {
        guard hasUnifiedMemory else {
            throw VivoRuntimeError.incompatibleDevice("NumiVivo requires an Apple unified-memory Metal device for its production runtime")
        }
        guard programBytes <= maximumBufferLength else {
            throw VivoRuntimeError.incompatibleDevice("ProgramPack exceeds Metal maximumBufferLength")
        }
        guard privateHeapBytes <= maximumBufferLength else {
            throw VivoRuntimeError.incompatibleDevice("a required runtime buffer exceeds Metal maximumBufferLength")
        }
        if recommendedMaximumWorkingSet > 0 {
            let total = programBytes.addingReportingOverflow(privateHeapBytes)
            guard !total.overflow, total.partialValue <= recommendedMaximumWorkingSet else {
                throw VivoRuntimeError.allocationFailed(
                    "ProgramPack plus private runtime state exceeds recommendedMaxWorkingSetSize"
                )
            }
        }
    }

    public func threadgroupSize(
        elementCount: Int,
        executionWidth: Int,
        pipelineMaximum: Int,
        preferred: Int? = nil
    ) -> MTLSize {
        guard elementCount > 0 else { return MTLSize(width: 1, height: 1, depth: 1) }
        let width = max(executionWidth, executionWidth)
        let target = min(preferred ?? recommendedThreadsPerThreadgroup, pipelineMaximum)
        let aligned = max(width, target / width * width)
        return MTLSize(width: min(aligned, pipelineMaximum), height: 1, depth: 1)
    }

    private static func detectFamily(_ device: MTLDevice) -> Family {
        if #available(macOS 15.0, iOS 18.0, *), device.supportsFamily(.apple9) { return .apple9OrNewer }
        if device.supportsFamily(.apple8) { return .apple8 }
        if device.supportsFamily(.apple7) { return .apple7 }
        if device.supportsFamily(.apple6) { return .apple6 }
        if device.supportsFamily(.apple5) { return .apple5 }
        if device.supportsFamily(.apple4) { return .apple4 }
        if device.supportsFamily(.apple3) { return .apple3 }
        if device.supportsFamily(.apple2) { return .apple2 }
        if device.supportsFamily(.apple1) { return .apple1 }
        if device.supportsFamily(.mac2) { return .mac2 }
        if device.supportsFamily(.mac1) { return .mac1 }
        return .unknown
    }

    private static func familyRank(_ family: Family) -> Int {
        switch family {
        case .unknown: 0
        case .mac1: 1
        case .mac2: 2
        case .apple1: 10
        case .apple2: 20
        case .apple3: 30
        case .apple4: 40
        case .apple5: 50
        case .apple6: 60
        case .apple7: 70
        case .apple8: 80
        case .apple9OrNewer: 90
        }
    }

    private static func recommendedThreads(for family: Family) -> Int {
        switch family {
        case .apple9OrNewer, .apple8, .apple7: 256
        case .apple6, .apple5, .apple4: 128
        default: 64
        }
    }
}

public enum VivoMetalDeviceSelector {
    public static func productionDevice() throws -> MTLDevice {
        #if os(macOS)
        let devices = MTLCopyAllDevices()
        if let unified = devices
            .filter({ $0.hasUnifiedMemory && !$0.isRemovable })
            .sorted(by: preferredDeviceOrder)
            .first {
            return unified
        }
        #endif
        guard let device = MTLCreateSystemDefaultDevice() else {
            throw VivoRuntimeError.metalUnavailable
        }
        guard device.hasUnifiedMemory else {
            throw VivoRuntimeError.incompatibleDevice("default Metal device does not expose unified memory")
        }
        return device
    }

    private static func preferredDeviceOrder(_ lhs: MTLDevice, _ rhs: MTLDevice) -> Bool {
        if lhs.isLowPower != rhs.isLowPower { return !lhs.isLowPower }
        if lhs.recommendedMaxWorkingSetSize != rhs.recommendedMaxWorkingSetSize {
            return lhs.recommendedMaxWorkingSetSize > rhs.recommendedMaxWorkingSetSize
        }
        return lhs.registryID < rhs.registryID
    }
}
