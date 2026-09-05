import Foundation
@preconcurrency import Metal

/// Swift 6 imports several Metal creation APIs as async-only even though the
/// underlying completion-handler operations remain appropriate for one-time
/// initialization. These blocking wrappers are deliberately confined to
/// construction paths; command submission and simulation execution remain async.
private final class VivoMetalCompletionBox<Value>: @unchecked Sendable {
    private let lock = NSLock()
    private var storage: Result<Value, Error>?

    func store(_ result: Result<Value, Error>) {
        lock.lock(); defer { lock.unlock() }
        storage = result
    }

    func load() -> Result<Value, Error>? {
        lock.lock(); defer { lock.unlock() }
        return storage
    }
}

extension MTLDevice {
    /// Synchronous construction compatibility for Swift 6 SDK imports.
    /// Do not use this from a simulation hot path.
    func makeLibrary(source: String, options: MTLCompileOptions?) throws -> MTLLibrary {
        let semaphore = DispatchSemaphore(value: 0)
        let box = VivoMetalCompletionBox<MTLLibrary>()
        self.makeLibrary(source: source, options: options) { library, error in
            if let library {
                box.store(.success(library))
            } else {
                box.store(.failure(error ?? VivoRuntimeError.incompatibleDevice("Metal library creation returned neither library nor error")))
            }
            semaphore.signal()
        }
        semaphore.wait()
        guard let result = box.load() else {
            throw VivoRuntimeError.incompatibleDevice("Metal library completion was not published")
        }
        return try result.get()
    }

    /// Synchronous construction compatibility for Swift 6 SDK imports.
    /// Pipeline compilation remains outside steady-state execution.
    func makeComputePipelineState(function: MTLFunction) throws -> MTLComputePipelineState {
        let semaphore = DispatchSemaphore(value: 0)
        let box = VivoMetalCompletionBox<MTLComputePipelineState>()
        self.makeComputePipelineState(function: function) { state, error in
            if let state {
                box.store(.success(state))
            } else {
                box.store(.failure(error ?? VivoRuntimeError.incompatibleDevice("Metal pipeline creation returned neither state nor error")))
            }
            semaphore.signal()
        }
        semaphore.wait()
        guard let result = box.load() else {
            throw VivoRuntimeError.incompatibleDevice("Metal pipeline completion was not published")
        }
        return try result.get()
    }
}

/// `MTLHeapDescriptor` has no label property in the macOS 15 SDK used by CI.
/// Existing assignments are diagnostic-only. Keep source compatibility without
/// pretending the descriptor owns a label; runtime heaps may still be labeled.
extension MTLHeapDescriptor {
    var label: String? {
        get { nil }
        set { _ = newValue }
    }
}
