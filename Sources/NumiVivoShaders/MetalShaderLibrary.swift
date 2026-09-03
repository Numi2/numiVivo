import Foundation
import Metal

public enum NumiVivoShaderError: Error, Sendable, CustomStringConvertible {
    case sourceResourceMissing
    case functionMissing(String)
    case compilationFailed(String)
    case pipelineFailed(String, String)

    public var description: String {
        switch self {
        case .sourceResourceMissing:
            return "NumiVivoKernels.metal is missing from the NumiVivoShaders resource bundle."
        case .functionMissing(let name):
            return "Metal function '\(name)' is missing from the NumiVivo shader library."
        case .compilationFailed(let message):
            return "NumiVivo Metal compilation failed: \(message)"
        case .pipelineFailed(let name, let message):
            return "NumiVivo pipeline '\(name)' could not be created: \(message)"
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
}

public struct NumiVivoPipeline: @unchecked Sendable {
    public let state: MTLComputePipelineState
    public let executionWidth: Int
    public let maximumThreadsPerThreadgroup: Int

    public init(state: MTLComputePipelineState) {
        self.state = state
        self.executionWidth = state.threadExecutionWidth
        self.maximumThreadsPerThreadgroup = state.maxTotalThreadsPerThreadgroup
    }

    public func threadgroupSize(for elementCount: Int, preferred: Int? = nil) -> MTLSize {
        guard elementCount > 0 else { return MTLSize(width: 1, height: 1, depth: 1) }
        let requested = preferred ?? max(executionWidth, min(maximumThreadsPerThreadgroup, executionWidth * 4))
        let aligned = max(executionWidth, (requested / executionWidth) * executionWidth)
        return MTLSize(width: min(aligned, maximumThreadsPerThreadgroup), height: 1, depth: 1)
    }

    public func gridSize(for elementCount: Int) -> MTLSize {
        MTLSize(width: max(elementCount, 1), height: 1, depth: 1)
    }
}

public actor NumiVivoPipelineCatalog {
    private let device: MTLDevice
    private let library: MTLLibrary
    private var cache: [NumiVivoKernel: NumiVivoPipeline] = [:]

    public init(device: MTLDevice) throws {
        self.device = device
        self.library = try Self.loadLibrary(device: device)
    }

    public func pipeline(_ kernel: NumiVivoKernel) throws -> NumiVivoPipeline {
        if let cached = cache[kernel] { return cached }
        guard let function = library.makeFunction(name: kernel.rawValue) else {
            throw NumiVivoShaderError.functionMissing(kernel.rawValue)
        }
        do {
            let descriptor = MTLComputePipelineDescriptor()
            descriptor.label = "NumiVivo.\(kernel.rawValue)"
            descriptor.computeFunction = function
            descriptor.threadGroupSizeIsMultipleOfThreadExecutionWidth = true
            let state = try device.makeComputePipelineState(descriptor: descriptor, options: [], reflection: nil)
            let pipeline = NumiVivoPipeline(state: state)
            cache[kernel] = pipeline
            return pipeline
        } catch {
            throw NumiVivoShaderError.pipelineFailed(kernel.rawValue, String(describing: error))
        }
    }

    public func preloadAll() throws {
        for kernel in NumiVivoKernel.allCases {
            _ = try pipeline(kernel)
        }
    }

    private static func loadLibrary(device: MTLDevice) throws -> MTLLibrary {
        if let library = try? device.makeDefaultLibrary(bundle: .module) {
            return library
        }
        guard let url = Bundle.module.url(forResource: "NumiVivoKernels", withExtension: "metal") else {
            throw NumiVivoShaderError.sourceResourceMissing
        }
        do {
            let source = try String(contentsOf: url, encoding: .utf8)
            let options = MTLCompileOptions()
            options.fastMathEnabled = false
            return try device.makeLibrary(source: source, options: options)
        } catch let error as NumiVivoShaderError {
            throw error
        } catch {
            throw NumiVivoShaderError.compilationFailed(String(describing: error))
        }
    }
}
