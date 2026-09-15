import Foundation
import Dispatch

/// Wall-clock accounting for the opt-in Metal NB objective within one native
/// NB cohort evaluation. These fields deliberately describe execution only:
/// they do not qualify GPU performance, statistical agreement, or biological
/// validity.
public struct VivoMetalNBExecutionProfile: Codable, Sendable, Equatable {
    public let schemaVersion: String
    public let cohortEvaluationWallClockNanoseconds: UInt64
    public let fixedDispersionFitAttempts: Int
    public let completedFixedDispersionFits: Int
    public let fixedDispersionFitWallClockNanoseconds: UInt64
    public let metalObjectiveAttempts: Int
    public let completedMetalObjectiveEvaluations: Int
    public let metalObjectiveWallClockNanoseconds: UInt64
    public let metalObjectiveBatches: Int
    public let metalObjectiveObservations: UInt64
    public let fp64ToFP32ConversionNanoseconds: UInt64
    public let sharedBufferWriteNanoseconds: UInt64
    public let commandEncodingAndSubmissionNanoseconds: UInt64
    public let commandCompletionNanoseconds: UInt64
    public let cpuReductionNanoseconds: UInt64
    public let cpuExactLikelihoodEvaluations: Int
    public let cpuExactLikelihoodObservations: UInt64
    public let cpuExactLikelihoodNanoseconds: UInt64
    public let qualification: String
}

/// The collector is dynamically scoped to one Metal-backed cohort evaluation.
/// Its lock keeps an accidental concurrent caller from corrupting diagnostics;
/// it does not serialize the actual statistical work.
final class VivoMetalNBExecutionProfileCollector: @unchecked Sendable {
    private let lock = NSLock()
    private var cohortEvaluationWallClockNanoseconds: UInt64 = 0
    private var fixedDispersionFitAttempts = 0
    private var completedFixedDispersionFits = 0
    private var fixedDispersionFitWallClockNanoseconds: UInt64 = 0
    private var metalObjectiveAttempts = 0
    private var completedMetalObjectiveEvaluations = 0
    private var metalObjectiveWallClockNanoseconds: UInt64 = 0
    private var metalObjectiveBatches = 0
    private var metalObjectiveObservations: UInt64 = 0
    private var fp64ToFP32ConversionNanoseconds: UInt64 = 0
    private var sharedBufferWriteNanoseconds: UInt64 = 0
    private var commandEncodingAndSubmissionNanoseconds: UInt64 = 0
    private var commandCompletionNanoseconds: UInt64 = 0
    private var cpuReductionNanoseconds: UInt64 = 0
    private var cpuExactLikelihoodEvaluations = 0
    private var cpuExactLikelihoodObservations: UInt64 = 0
    private var cpuExactLikelihoodNanoseconds: UInt64 = 0

    func recordFitStarted() {
        lock.lock(); defer { lock.unlock() }
        fixedDispersionFitAttempts += 1
    }

    func recordFitFinished(completed: Bool, elapsedNanoseconds: UInt64) {
        lock.lock(); defer { lock.unlock() }
        if completed { completedFixedDispersionFits += 1 }
        fixedDispersionFitWallClockNanoseconds += elapsedNanoseconds
    }

    func recordObjectiveStarted() {
        lock.lock(); defer { lock.unlock() }
        metalObjectiveAttempts += 1
    }

    func recordObjectiveFinished(completed: Bool, elapsedNanoseconds: UInt64) {
        lock.lock(); defer { lock.unlock() }
        if completed { completedMetalObjectiveEvaluations += 1 }
        metalObjectiveWallClockNanoseconds += elapsedNanoseconds
    }

    func recordBatch(observations: Int, bufferWriteNanoseconds: UInt64,
                     commandEncodingAndSubmissionNanoseconds: UInt64,
                     commandCompletionNanoseconds: UInt64,
                     cpuReductionNanoseconds: UInt64) {
        lock.lock(); defer { lock.unlock() }
        metalObjectiveBatches += 1
        metalObjectiveObservations += UInt64(observations)
        sharedBufferWriteNanoseconds += bufferWriteNanoseconds
        self.commandEncodingAndSubmissionNanoseconds += commandEncodingAndSubmissionNanoseconds
        self.commandCompletionNanoseconds += commandCompletionNanoseconds
        self.cpuReductionNanoseconds += cpuReductionNanoseconds
    }

    func recordFP64ToFP32Conversion(elapsedNanoseconds: UInt64) {
        lock.lock(); defer { lock.unlock() }
        fp64ToFP32ConversionNanoseconds += elapsedNanoseconds
    }

    func recordCPUExactLikelihood(observations: Int, elapsedNanoseconds: UInt64) {
        lock.lock(); defer { lock.unlock() }
        cpuExactLikelihoodEvaluations += 1
        cpuExactLikelihoodObservations += UInt64(observations)
        cpuExactLikelihoodNanoseconds += elapsedNanoseconds
    }

    func recordCohortEvaluation(elapsedNanoseconds: UInt64) {
        lock.lock(); defer { lock.unlock() }
        cohortEvaluationWallClockNanoseconds = elapsedNanoseconds
    }

    func snapshot() -> VivoMetalNBExecutionProfile {
        lock.lock(); defer { lock.unlock() }
        return .init(schemaVersion: "numi.vivo.metal-nb-execution-profile.v1",
            cohortEvaluationWallClockNanoseconds: cohortEvaluationWallClockNanoseconds,
            fixedDispersionFitAttempts: fixedDispersionFitAttempts,
            completedFixedDispersionFits: completedFixedDispersionFits,
            fixedDispersionFitWallClockNanoseconds: fixedDispersionFitWallClockNanoseconds,
            metalObjectiveAttempts: metalObjectiveAttempts,
            completedMetalObjectiveEvaluations: completedMetalObjectiveEvaluations,
            metalObjectiveWallClockNanoseconds: metalObjectiveWallClockNanoseconds,
            metalObjectiveBatches: metalObjectiveBatches,
            metalObjectiveObservations: metalObjectiveObservations,
            fp64ToFP32ConversionNanoseconds: fp64ToFP32ConversionNanoseconds,
            sharedBufferWriteNanoseconds: sharedBufferWriteNanoseconds,
            commandEncodingAndSubmissionNanoseconds: commandEncodingAndSubmissionNanoseconds,
            commandCompletionNanoseconds: commandCompletionNanoseconds,
            cpuReductionNanoseconds: cpuReductionNanoseconds,
            cpuExactLikelihoodEvaluations: cpuExactLikelihoodEvaluations,
            cpuExactLikelihoodObservations: cpuExactLikelihoodObservations,
            cpuExactLikelihoodNanoseconds: cpuExactLikelihoodNanoseconds,
            qualification: "Experimental wall-clock accounting for the opt-in Metal FP32 fixed-dispersion NB path only. This is not a GPU-profiler trace or a performance, numerical-agreement, statistical, or biological qualification. CPU gene-wise dispersion profiling, import/export, aggregation and other analysis work are outside the component totals.")
    }
}

enum VivoMetalNBExecutionProfileContext {
    @TaskLocal static var collector: VivoMetalNBExecutionProfileCollector?
}

@inline(__always) func vivoMetalNBClock() -> UInt64 {
    DispatchTime.now().uptimeNanoseconds
}
#if canImport(Metal)
@preconcurrency import Metal

/// The NB objective is evaluated many times during a dispersion profile. Keep
/// the compiled pipeline per process so the opt-in cohort backend does not
/// recompile the same shader for every line-search evaluation. The cache is
/// protected because cohort callers may run independent contrasts concurrently.
private struct VivoMetalNBResources {
    let device: MTLDevice
    let pipeline: MTLComputePipelineState
    let queue: MTLCommandQueue
    let countBuffer: MTLBuffer
    let meanBuffer: MTLBuffer
    let termBuffer: MTLBuffer
}

private final class VivoMetalNBPipelineCache: @unchecked Sendable {
    private let lock = NSLock()
    private let source: String
    private var cached: VivoMetalNBResources?

    init(source: String) { self.source = source }

    func withResources<T>(_ body: (VivoMetalNBResources) throws -> T) throws -> T {
        lock.lock(); defer { lock.unlock() }
        if let cached { return try body(cached) }
        guard let device = MTLCreateSystemDefaultDevice(), device.hasUnifiedMemory,
              device.supportsFamily(.apple7), !device.name.lowercased().contains("paravirtual") else {
            throw VivoOmicsStatisticsError.invalid("Metal NB objective requires a physical Apple GPU")
        }
        let options = MTLCompileOptions()
        options.fastMathEnabled = false
        let library = try device.makeLibrary(source: source, options: options)
        guard let function = library.makeFunction(name: "vivo_nb_variable_objective") else {
            throw VivoOmicsStatisticsError.invalid("missing Metal NB objective kernel")
        }
        let pipeline = try device.makeComputePipelineState(function: function)
        guard let queue = device.makeCommandQueue(),
              let countBuffer = device.makeBuffer(length: VivoMetalNegativeBinomialLikelihood.batchObservations * MemoryLayout<UInt32>.stride,
                                                  options: .storageModeShared),
              let meanBuffer = device.makeBuffer(length: VivoMetalNegativeBinomialLikelihood.batchObservations * MemoryLayout<Float>.stride,
                                                 options: .storageModeShared),
              let termBuffer = device.makeBuffer(length: VivoMetalNegativeBinomialLikelihood.batchObservations * MemoryLayout<Float>.stride,
                                                 options: .storageModeShared) else {
            throw VivoOmicsStatisticsError.invalid("Metal NB objective buffers")
        }
        let result = VivoMetalNBResources(device: device, pipeline: pipeline, queue: queue,
                                          countBuffer: countBuffer, meanBuffer: meanBuffer, termBuffer: termBuffer)
        cached = result
        return try body(result)
    }
}
#endif

/// Bounded FP32 Metal evaluation of the mean-dependent part of an NB2 log
/// likelihood. Count-only log-gamma terms stay on the exact CPU owner, so this
/// operator is an explicit objective accelerator rather than a second
/// statistical authority.
final class VivoMetalNegativeBinomialLikelihood {
    static let method = "metal-FP32-batched-NB-objective-v1"
    static let maximumCount = 16_777_216
    static let maximumObservations = 4_000_000
    static let batchObservations = 262_144
    private let counts: [UInt32]
    static let source = """
    #include <metal_stdlib>
    using namespace metal;
    kernel void vivo_nb_variable_objective(
        device const uint *counts [[buffer(0)]],
        device const float *means [[buffer(1)]],
        device float *terms [[buffer(2)]],
        constant float &dispersion [[buffer(3)]],
        constant uint &count [[buffer(4)]],
        uint i [[thread_position_in_grid]]) {
        if (i >= count) return;
        float y = float(counts[i]);
        float mu = means[i];
        float x = dispersion * mu;
        // Separate the 1/a term before evaluating it. Directly multiplying
        // (y + 1/a) by log1p(a*mu) loses the Poisson-limit objective in FP32
        // when a is small. The short series keeps log1p(a*mu)/a resolved.
        float scaledLog;
        if (fabs(x) < 1.0e-3f) {
            float x2 = x * x;
            float series = 1.0f - x * 0.5f + x2 / 3.0f - x2 * x * 0.25f;
            scaledLog = mu * series;
        } else {
            scaledLog = log(1.0f + x) / dispersion;
        }
        float logOnePlus = fabs(x) < 1.0e-3f
            ? x - x * x * 0.5f + x * x * x / 3.0f - x * x * x * x * 0.25f
            : log(1.0f + x);
        terms[i] = y * (log(mu) - logOnePlus) - scaledLog;
    }
    """

    #if canImport(Metal)
    private static let pipelineCache = VivoMetalNBPipelineCache(source: source)
    private let dispersion: Float
    #endif

    init(counts: [UInt64], dispersion: Double) throws {
        guard !counts.isEmpty, counts.count <= Self.maximumObservations,
              counts.allSatisfy({ $0 <= UInt64(Self.maximumCount) }),
              dispersion.isFinite, (1e-8...100).contains(dispersion),
              Float(dispersion).isFinite else {
            throw VivoOmicsStatisticsError.invalid("Metal NB objective count or dispersion domain")
        }
        self.counts = counts.map(UInt32.init)
        #if canImport(Metal)
        try Task.checkCancellation()
        try Self.pipelineCache.withResources { _ in }
        self.dispersion = Float(dispersion)
        #else
        throw VivoOmicsStatisticsError.invalid("Metal NB objective is unavailable")
        #endif
    }

    /// Evaluates only y*log(mu) - (y + 1/a)*log(1+a*mu), which is the part
    /// that changes during coefficient line search. Results are widened and
    /// compensated on the CPU after each bounded dispatch.
    func evaluate(means: [Double]) throws -> Double {
        guard means.count == counts.count else {
            throw VivoOmicsStatisticsError.invalid("Metal NB objective mean dimensions")
        }
        #if canImport(Metal)
        let execution = VivoMetalNBExecutionProfileContext.collector
        let evaluationStarted = execution == nil ? nil : vivoMetalNBClock()
        execution?.recordObjectiveStarted()
        var completed = false
        defer {
            if let execution, let evaluationStarted {
                execution.recordObjectiveFinished(completed: completed,
                    elapsedNanoseconds: vivoMetalNBClock() - evaluationStarted)
            }
        }
        try Task.checkCancellation()
        let conversionStarted = execution == nil ? nil : vivoMetalNBClock()
        let fp32Means: [Float]
        do {
            fp32Means = try means.map { value -> Float in
                guard value.isFinite, value > 0, value >= Double(Float.leastNormalMagnitude),
                      value <= Double(Float.greatestFiniteMagnitude) else {
                    throw VivoOmicsStatisticsError.invalid("Metal NB objective mean outside FP32 domain")
                }
                let converted = Float(value)
                guard converted.isFinite, converted > 0,
                      (converted * dispersion).isFinite else {
                    throw VivoOmicsStatisticsError.invalid("Metal NB objective mean or product outside FP32 domain")
                }
                return converted
            }
        } catch {
            if let execution, let conversionStarted {
                execution.recordFP64ToFP32Conversion(elapsedNanoseconds: vivoMetalNBClock() - conversionStarted)
            }
            throw error
        }
        if let execution, let conversionStarted {
            execution.recordFP64ToFP32Conversion(elapsedNanoseconds: vivoMetalNBClock() - conversionStarted)
        }
        let result = try Self.pipelineCache.withResources { resources in
            let countPointer = resources.countBuffer.contents().bindMemory(to: UInt32.self, capacity: Self.batchObservations)
            let meanPointer = resources.meanBuffer.contents().bindMemory(to: Float.self, capacity: Self.batchObservations)
            let termPointer = resources.termBuffer.contents().bindMemory(to: Float.self, capacity: Self.batchObservations)
            var total = 0.0
            var compensation = 0.0
            for start in stride(from: 0, to: counts.count, by: Self.batchObservations) {
                try Task.checkCancellation()
                let count = min(Self.batchObservations, counts.count - start)
                let bufferWriteStarted = execution == nil ? nil : vivoMetalNBClock()
                counts[start..<start + count].withUnsafeBufferPointer { values in
                    countPointer.update(from: values.baseAddress!, count: count)
                }
                fp32Means[start..<start + count].withUnsafeBufferPointer { values in
                    meanPointer.update(from: values.baseAddress!, count: count)
                }
                let bufferWriteElapsed = bufferWriteStarted.map { vivoMetalNBClock() - $0 } ?? 0
                let encodingStarted = execution == nil ? nil : vivoMetalNBClock()
                guard let command = resources.queue.makeCommandBuffer(), let encoder = command.makeComputeCommandEncoder() else {
                    throw VivoOmicsStatisticsError.invalid("Metal NB objective command")
                }
                encoder.setComputePipelineState(resources.pipeline)
                encoder.setBuffer(resources.countBuffer, offset: 0, index: 0)
                encoder.setBuffer(resources.meanBuffer, offset: 0, index: 1)
                encoder.setBuffer(resources.termBuffer, offset: 0, index: 2)
                var a = dispersion, n = UInt32(count)
                encoder.setBytes(&a, length: MemoryLayout<Float>.size, index: 3)
                encoder.setBytes(&n, length: MemoryLayout<UInt32>.size, index: 4)
                let width = min(resources.pipeline.threadExecutionWidth, resources.pipeline.maxTotalThreadsPerThreadgroup)
                encoder.dispatchThreads(MTLSize(width: count, height: 1, depth: 1),
                                        threadsPerThreadgroup: MTLSize(width: width, height: 1, depth: 1))
                encoder.endEncoding()
                command.commit()
                let encodingElapsed = encodingStarted.map { vivoMetalNBClock() - $0 } ?? 0
                let completionStarted = execution == nil ? nil : vivoMetalNBClock()
                command.waitUntilCompleted()
                let completionElapsed = completionStarted.map { vivoMetalNBClock() - $0 } ?? 0
                guard command.status == .completed else {
                    throw VivoOmicsStatisticsError.invalid("Metal NB objective failed: " +
                        (command.error?.localizedDescription ?? "unknown"))
                }
                let reductionStarted = execution == nil ? nil : vivoMetalNBClock()
                for i in 0..<count {
                    let term = Double(termPointer[i])
                    guard term.isFinite else { throw VivoOmicsStatisticsError.invalid("nonfinite Metal NB objective term") }
                    let corrected = term - compensation
                    let next = total + corrected
                    compensation = (next - total) - corrected
                    total = next
                }
                let reductionElapsed = reductionStarted.map { vivoMetalNBClock() - $0 } ?? 0
                execution?.recordBatch(observations: count, bufferWriteNanoseconds: bufferWriteElapsed,
                    commandEncodingAndSubmissionNanoseconds: encodingElapsed,
                    commandCompletionNanoseconds: completionElapsed, cpuReductionNanoseconds: reductionElapsed)
            }
            guard total.isFinite else { throw VivoOmicsStatisticsError.invalid("nonfinite Metal NB objective") }
            return total
        }
        completed = true
        return result
        #else
        throw VivoOmicsStatisticsError.invalid("Metal NB objective is unavailable")
        #endif
    }
}
