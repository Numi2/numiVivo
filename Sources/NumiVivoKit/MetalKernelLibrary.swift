@preconcurrency import Metal
import Foundation
import NumiVivoShaders

final class VivoMetalKernelLibrary: @unchecked Sendable {
    enum Kernel: String, CaseIterable, Sendable {
        case initializeProgram = "nvivo_initialize_program"
        case prepareStep = "nvivo_prepare_step"
        case stageInputUpdates = "nvivo_stage_input_updates"
        case evaluateReactionCohort = "nvivo_evaluate_reaction_cohort"
        case applyIncidence = "nvivo_apply_incidence"
        case applyRules = "nvivo_apply_rules"
        case validateState = "nvivo_validate_state"
        case evaluateMonitors = "nvivo_evaluate_monitors"
        case commitIfValid = "nvivo_commit_if_valid"
        case diffuseCSR = "nvivo_diffuse_csr"
        case reduceSumStage1 = "nvivo_reduce_sum_stage1"
    }

    struct Pipeline: @unchecked Sendable {
        let state: MTLComputePipelineState
        let executionWidth: Int
        let maximumThreads: Int

        func threadsPerThreadgroup(preferred: Int? = nil) -> MTLSize {
            let requested = preferred ?? max(executionWidth, 64)
            return MTLSize(width: max(1, min(requested, maximumThreads)), height: 1, depth: 1)
        }
    }

    let device: MTLDevice
    let library: MTLLibrary
    let argumentEncoder: MTLArgumentEncoder
    private let pipelines: [Kernel: Pipeline]

    init(device: MTLDevice) throws {
        self.device = device
        let loadedLibrary: MTLLibrary
        if let precompiledURL = NumiVivoShaderResources.libraryURL(),
           FileManager.default.fileExists(atPath: precompiledURL.path) {
            do {
                loadedLibrary = try device.makeLibrary(URL: precompiledURL)
            } catch {
                throw VivoRuntimeError.metalLibrary(
                    "Unable to load precompiled NumiVivo Metal library: \(error)"
                )
            }
        } else {
            let source: String
            do {
                source = try NumiVivoShaderResources.completeMetalSource()
            } catch {
                throw VivoRuntimeError.missingShaderResource(String(describing: error))
            }
            let options = MTLCompileOptions()
            options.fastMathEnabled = false
            options.languageVersion = .version3_1
            do {
                loadedLibrary = try device.makeLibrary(source: source, options: options)
            } catch {
                throw VivoRuntimeError.metalLibrary(
                    "Unable to compile NumiVivo Metal kernels: \(error)"
                )
            }
        }
        self.library = loadedLibrary

        func function(for kernel: Kernel) -> MTLFunction? {
            loadedLibrary.makeFunction(name: kernel.rawValue) ??
            loadedLibrary.makeFunction(name: "numivivo::\(kernel.rawValue)")
        }

        var built: [Kernel: Pipeline] = [:]
        built.reserveCapacity(Kernel.allCases.count)
        var encoderFunction: MTLFunction?
        for kernel in Kernel.allCases {
            guard let function = function(for: kernel) else {
                throw VivoRuntimeError.missingKernel(kernel.rawValue)
            }
            if kernel == .initializeProgram { encoderFunction = function }
            do {
                let state = try device.makeComputePipelineState(function: function)
                built[kernel] = Pipeline(
                    state: state,
                    executionWidth: state.threadExecutionWidth,
                    maximumThreads: state.maxTotalThreadsPerThreadgroup
                )
            } catch {
                throw VivoRuntimeError.pipeline(
                    "Unable to create pipeline for \(kernel.rawValue): \(error)"
                )
            }
        }
        guard let encoderFunction else {
            throw VivoRuntimeError.missingKernel(Kernel.initializeProgram.rawValue)
        }
        self.argumentEncoder = encoderFunction.makeArgumentEncoder(bufferIndex: 0)
        self.pipelines = built
    }

    subscript(kernel: Kernel) -> Pipeline {
        pipelines[kernel]!
    }
}

public struct VivoMetalDeviceInfo: Sendable, Hashable, Codable {
    public var name: String
    public var registryID: UInt64
    public var hasUnifiedMemory: Bool
    public var maximumBufferLength: UInt64
    public var recommendedMaximumWorkingSetSize: UInt64
    public var isLowPower: Bool
    public var isRemovable: Bool
}

extension MTLDevice {
    var nvivoInfo: VivoMetalDeviceInfo {
        #if os(macOS)
        let recommended = UInt64(recommendedMaxWorkingSetSize)
        let lowPower = isLowPower
        let removable = isRemovable
        #else
        let recommended: UInt64 = 0
        let lowPower = true
        let removable = false
        #endif
        return VivoMetalDeviceInfo(
            name: name,
            registryID: registryID,
            hasUnifiedMemory: hasUnifiedMemory,
            maximumBufferLength: UInt64(maxBufferLength),
            recommendedMaximumWorkingSetSize: recommended,
            isLowPower: lowPower,
            isRemovable: removable
        )
    }
}
