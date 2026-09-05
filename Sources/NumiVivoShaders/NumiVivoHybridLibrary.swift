import Foundation
@preconcurrency import Metal

/// Independent source loading prevents legacy ProgramPack argument-buffer
/// declarations from being concatenated into the executable hybrid library.
public final class NumiVivoHybridPipelineSet: @unchecked Sendable {
    private let pipelines: [String: MTLComputePipelineState]
    fileprivate init(device: MTLDevice) throws {
        let names = ["clear", "reset_continuation", "rk_rates", "rk_predict", "rk_correct",
                     "tau_sample", "tau_apply", "exact_advance", "validate", "publish"]
        let url = Bundle.module.url(forResource: "NumiVivoHybridExecution", withExtension: "metal")
            ?? Bundle.module.url(forResource: "NumiVivoHybridExecution", withExtension: "metal", subdirectory: "Resources")
        let library: MTLLibrary
        if let compiled = try? device.makeDefaultLibrary(bundle: .module),
           names.allSatisfy({ compiled.makeFunction(name: "nvivo_hybrid_" + $0) != nil }) {
            library = compiled
        } else {
            guard let url else { throw NumiVivoShaderError.sourceResourceMissing }
            let source = try String(contentsOf: url, encoding: .utf8)
            let options = MTLCompileOptions()
            options.fastMathEnabled = false
            library = try device.makeLibrary(source: source, options: options)
        }
        var result: [String: MTLComputePipelineState] = [:]
        for name in names {
            let symbol = "nvivo_hybrid_" + name
            guard let function = library.makeFunction(name: symbol) else {
                throw NumiVivoShaderError.functionMissing(symbol)
            }
            let descriptor = MTLComputePipelineDescriptor()
            descriptor.computeFunction = function
            descriptor.label = symbol
            // dispatchThreads may produce a partial final SIMD group.
            descriptor.threadGroupSizeIsMultipleOfThreadExecutionWidth = false
            result[name] = try device.makeComputePipelineState(descriptor: descriptor, options: [], reflection: nil)
        }
        pipelines = result
    }
    public func pipeline(_ name: String) throws -> MTLComputePipelineState {
        guard let pipeline = pipelines[name] else { throw NumiVivoShaderError.functionMissing(name) }
        return pipeline
    }
}

public actor NumiVivoHybridPipelineCache {
    public static let shared = NumiVivoHybridPipelineCache()
    private var values: [UInt64: NumiVivoHybridPipelineSet] = [:]
    public func pipelines(device: MTLDevice) throws -> NumiVivoHybridPipelineSet {
        if let existing = values[device.registryID] { return existing }
        let value = try NumiVivoHybridPipelineSet(device: device)
        values[device.registryID] = value
        return value
    }
}
