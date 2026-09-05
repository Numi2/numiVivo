import Foundation
import Metal

public enum NumiVivoShaderError: Error, Sendable, CustomStringConvertible {
    case sourceResourceMissing
    case functionMissing(String)
    case compilationFailed(String)
    case pipelineFailed(String, String)
    public var description: String {
        switch self {
        case .sourceResourceMissing: return "The selected NumiVivo shader source resource is missing."
        case .functionMissing(let name): return "Metal function '\(name)' is missing from its selected module."
        case .compilationFailed(let message): return "NumiVivo Metal compilation failed: \(message)"
        case .pipelineFailed(let name, let message): return "NumiVivo pipeline '\(name)' could not be created: \(message)"
        }
    }
}

public enum NumiVivoKernel: String, CaseIterable, Sendable {
    case clearStatus = "nvivo_clear_status"
    case prepareTransaction = "nvivo_prepare_transaction"
    case applyCouplingUpdates = "nvivo_apply_coupling_updates"
    case f1HeunPredict = "nvivo_f1_heun_predict"
    case f1HeunCorrect = "nvivo_f1_heun_correct"
    case f2SampleReactions = "nvivo_f2_sample_reactions"
    case f2ApplyReactions = "nvivo_f2_apply_reactions"
    case f3Transport = "nvivo_f3_transport"
    case executeRules = "nvivo_execute_rules"
    case evaluateMonitors = "nvivo_evaluate_monitors"
    case validateShadow = "nvivo_validate_shadow"
    case publish = "nvivo_publish"

    case physiologyClearStatus = "nvivo_phys_clear_status"
    case physiologyPrepareTransaction = "nvivo_phys_prepare_transaction"
    case physiologyApplyTransforms = "nvivo_phys_apply_transforms"
    case physiologyHeunPredict = "nvivo_phys_heun_predict"
    case physiologyHeunCorrect = "nvivo_phys_heun_correct"
    case physiologyValidateCandidate = "nvivo_phys_validate_candidate"
    case physiologyPublish = "nvivo_phys_publish"

    case mdClearForce = "nvivo_md_clear_force"
    case mdClearStatus = "nvivo_md_clear_status"
    case mdUpdateVirtualPosition = "nvivo_md_update_virtual_position"
    case mdUpdateVirtualVelocity = "nvivo_md_update_virtual_velocity"
    case mdRedistributeVirtualForce = "nvivo_md_redistribute_virtual_force"
    case mdBonded = "nvivo_md_bonded"
    case mdBuildNeighborList = "nvivo_md_build_neighbor_list"
    case mdValidateNeighborDisplacement = "nvivo_md_validate_neighbor_displacement"
    case mdNonbondedNeighbor = "nvivo_md_nonbonded_neighbor"
    case mdNonbondedDirect = "nvivo_md_nonbonded_direct"
    case mdGridClear = "nvivo_md_grid_clear"
    case mdGridBin = "nvivo_md_grid_bin"
    case mdGridBuildNeighbors = "nvivo_md_grid_build_neighbors"
    case mdHalfKick = "nvivo_md_half_kick"
    case mdDrift = "nvivo_md_drift"
    case mdLangevin = "nvivo_md_langevin"
    case mdConstraintPosition = "nvivo_md_constraint_position"
    case mdConstraintVelocity = "nvivo_md_constraint_velocity"
    case mdValidateConstraints = "nvivo_md_validate_constraints"
    case mdKinetic = "nvivo_md_kinetic"
    case mdValidate = "nvivo_md_validate"

    fileprivate var sourceModule: String {
        switch self {
        case .physiologyClearStatus, .physiologyPrepareTransaction, .physiologyApplyTransforms,
             .physiologyHeunPredict, .physiologyHeunCorrect, .physiologyValidateCandidate,
             .physiologyPublish:
            return "NumiVivoPhysiologyKernels"
        case .mdUpdateVirtualPosition, .mdUpdateVirtualVelocity, .mdRedistributeVirtualForce:
            return "NumiVivoMDVirtualSites"
        case .mdGridClear, .mdGridBin, .mdGridBuildNeighbors:
            return "NumiVivoMDNeighborGrid"
        case .mdClearForce, .mdClearStatus, .mdBonded, .mdBuildNeighborList,
             .mdValidateNeighborDisplacement, .mdNonbondedNeighbor, .mdNonbondedDirect,
             .mdHalfKick, .mdDrift, .mdLangevin, .mdConstraintPosition,
             .mdConstraintVelocity, .mdValidateConstraints, .mdKinetic, .mdValidate:
            return "NumiVivoMDKernels"
        default:
            return "NumiVivoProgramPackRuntime"
        }
    }
}

public struct NumiVivoPipeline: @unchecked Sendable {
    public let state: MTLComputePipelineState
    public let executionWidth: Int
    public let maximumThreadsPerThreadgroup: Int
    public init(state: MTLComputePipelineState) {
        self.state = state
        executionWidth = state.threadExecutionWidth
        maximumThreadsPerThreadgroup = state.maxTotalThreadsPerThreadgroup
    }
    public func threadgroupSize(for elementCount: Int, preferred: Int? = nil) -> MTLSize {
        let width = max(executionWidth, 1)
        let requested = max(width, preferred ?? width * 4)
        let capped = min(requested, maximumThreadsPerThreadgroup)
        let aligned = max(1, capped >= width ? (capped / width) * width : capped)
        return .init(width: aligned, height: 1, depth: 1)
    }
    public func gridSize(for elementCount: Int) -> MTLSize {
        .init(width: max(elementCount, 1), height: 1, depth: 1)
    }
}

private final class NumiVivoRuntimeLibraryCache: @unchecked Sendable {
    static let shared = NumiVivoRuntimeLibraryCache()
    private let lock = NSLock()
    private var libraries: [String: MTLLibrary] = [:]

    func library(device: MTLDevice, module: String) throws -> MTLLibrary {
        lock.lock(); defer { lock.unlock() }
        let key = "\(device.registryID)/\(module)/v1"
        if let value = libraries[key] { return value }
        let url = Bundle.module.url(forResource: module, withExtension: "metal")
            ?? Bundle.module.url(forResource: module, withExtension: "metal", subdirectory: "Resources")
        guard let url else { throw NumiVivoShaderError.sourceResourceMissing }
        do {
            let source = try String(contentsOf: url, encoding: .utf8)
            let options = MTLCompileOptions()
            options.fastMathEnabled = false
            let library = try device.makeLibrary(source: source, options: options)
            libraries[key] = library
            return library
        } catch {
            throw NumiVivoShaderError.compilationFailed(String(describing: error))
        }
    }
}

public actor NumiVivoPipelineCatalog {
    private let device: MTLDevice
    private var cache: [NumiVivoKernel: NumiVivoPipeline] = [:]
    public init(device: MTLDevice) throws { self.device = device }

    public func pipeline(_ kernel: NumiVivoKernel) throws -> NumiVivoPipeline {
        if let value = cache[kernel] { return value }
        let library = try NumiVivoRuntimeLibraryCache.shared.library(device: device,
                                                                     module: kernel.sourceModule)
        guard let function = library.makeFunction(name: kernel.rawValue) else {
            throw NumiVivoShaderError.functionMissing(kernel.rawValue)
        }
        do {
            let descriptor = MTLComputePipelineDescriptor()
            descriptor.label = "NumiVivo.\(kernel.rawValue)"
            descriptor.computeFunction = function
            descriptor.threadGroupSizeIsMultipleOfThreadExecutionWidth = false
            let state = try device.makeComputePipelineState(descriptor: descriptor,
                                                            options: [], reflection: nil)
            let value = NumiVivoPipeline(state: state)
            cache[kernel] = value
            return value
        } catch {
            throw NumiVivoShaderError.pipelineFailed(kernel.rawValue, String(describing: error))
        }
    }

    public func preloadAll() throws {
        for kernel in NumiVivoKernel.allCases { _ = try pipeline(kernel) }
    }
}
