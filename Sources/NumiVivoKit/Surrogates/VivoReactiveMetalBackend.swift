import Foundation
#if canImport(Metal)
@preconcurrency import Metal

public extension VivoReactiveEnergyForceBackend {
    /// Batched FP32 proposal inference. Training, checkpoint accounting and exact
    /// endpoint corrections remain FP64. No speedup is asserted without timing.
    static func metal(model: VivoReactiveSurrogateModel, maximumBatchSize: Int = 1024,
                      maximumBytes: Int = 256*1024*1024, device: MTLDevice? = nil) throws -> Self {
        let backend = try VivoReactiveMetalInference(model: model,maximumBatchSize: maximumBatchSize,
            maximumBytes: maximumBytes,device: device)
        return .init(modelFingerprint: model.fingerprint,numericalProfile: "metal-fp32-rbf-analytic-derivative/v1",
            predict: { try await backend.predict($0) })
    }
}
private actor VivoReactiveMetalInference {
    let model: VivoReactiveSurrogateModel
    let maximumBatchSize: Int
    let maximumBytes: Int
    let device: MTLDevice
    let queue: MTLCommandQueue
    let pipeline: MTLComputePipelineState
    let centers: MTLBuffer
    let coefficients: MTLBuffer
    let scales: MTLBuffer
    var inFlight = false
    struct Uniforms { var frames: UInt32; var features: UInt32; var centers: UInt32; var members: UInt32 }
    init(model: VivoReactiveSurrogateModel, maximumBatchSize: Int, maximumBytes: Int, device requested: MTLDevice?) throws {
        try model.validate()
        guard (1...65536).contains(maximumBatchSize), maximumBytes > 0,
              let device = requested ?? MTLCreateSystemDefaultDevice(), let queue = device.makeCommandQueue() else {
            throw VivoChemistryError.unsupported("reactive Metal device/queue or batch capacity")
        }
        let payload = model.payload
        let n = payload.configuration.features.count, c = payload.centersNM.count, k = payload.coefficientsKJPerMol.count
        let reserved = Double(maximumBatchSize)*Double(n+k*(n+2))*4+Double(c*n+k*(c+1)+n)*4
        guard reserved <= Double(maximumBytes), reserved <= Double(device.maxBufferLength) else {
            throw VivoChemistryError.resourceLimit("reactive Metal batch allocation")
        }
        func buffer(_ values: [Double]) throws -> MTLBuffer {
            let converted = values.map(Float.init)
            guard converted.allSatisfy(\.isFinite), let buffer = converted.withUnsafeBytes({ raw in
                device.makeBuffer(bytes: raw.baseAddress!,length: raw.count,options: .storageModeShared)
            }) else { throw VivoChemistryError.resourceLimit("reactive Metal model coefficients or buffer") }
            return buffer
        }
        self.model = model; self.maximumBatchSize = maximumBatchSize; self.maximumBytes = maximumBytes
        self.device = device; self.queue = queue
        centers = try buffer(payload.centersNM.flatMap { $0 })
        coefficients = try buffer(payload.coefficientsKJPerMol.flatMap { $0 })
        scales = try buffer(payload.configuration.features.map { 1/($0.scaleNM*$0.scaleNM) })
        let options = MTLCompileOptions(); options.fastMathEnabled = false
        let library = try device.makeLibrary(source: Self.source,options: options)
        guard let function = library.makeFunction(name: "vivo_reactive_delta") else { throw VivoChemistryError.invalid("missing reactive Metal kernel") }
        pipeline = try device.makeComputePipelineState(function: function)
    }
    func predict(_ positions: [[VivoVector3D]]) async throws -> [VivoReactiveSurrogatePrediction] {
        guard !inFlight else { throw VivoChemistryError.invalid("reactive Metal backend is already in flight") }
        inFlight = true; defer { inFlight = false }
        let p = model.payload, d = p.configuration.features.count, k = p.coefficientsKJPerMol.count, n = p.authorityDefinition.atomIndices.count
        guard !positions.isEmpty, positions.count <= maximumBatchSize,
              positions.allSatisfy({ $0.count == n && $0.allSatisfy(\.isFinite) }) else {
            throw VivoChemistryError.invalid("reactive Metal input geometry")
        }
        let geometries = try positions.map { try VivoReactiveDeltaSurrogate.geometry($0,features: p.configuration.features) }
        let input = geometries.flatMap { $0.distances.map(Float.init) }
        guard input.allSatisfy(\.isFinite) else { throw VivoChemistryError.invalid("reactive descriptor FP32 overflow") }
        let outputCount = positions.count*k*(d+2)
        guard let inputBuffer = input.withUnsafeBytes({ device.makeBuffer(bytes: $0.baseAddress!,length: $0.count,options: .storageModeShared) }),
              let outputBuffer = device.makeBuffer(length: outputCount*MemoryLayout<Float>.stride,options: .storageModeShared),
              let command = queue.makeCommandBuffer(), let encoder = command.makeComputeCommandEncoder() else {
            throw VivoChemistryError.resourceLimit("reactive Metal execution buffers")
        }
        var uniforms = Uniforms(frames: UInt32(positions.count),features: UInt32(d),centers: UInt32(p.centersNM.count),members: UInt32(k))
        encoder.setComputePipelineState(pipeline)
        for (i,buffer) in [inputBuffer,centers,coefficients,scales,outputBuffer].enumerated() { encoder.setBuffer(buffer,offset: 0,index: i) }
        encoder.setBytes(&uniforms,length: MemoryLayout<Uniforms>.stride,index: 5)
        let width = min(pipeline.maxTotalThreadsPerThreadgroup,max(1,pipeline.threadExecutionWidth*4))
        encoder.dispatchThreads(.init(width: positions.count*k,height: 1,depth: 1),threadsPerThreadgroup: .init(width: width,height: 1,depth: 1))
        encoder.endEncoding()
        await withCheckedContinuation { (continuation: CheckedContinuation<Void,Never>) in
            command.addCompletedHandler { _ in continuation.resume() }; command.commit()
        }
        try Task.checkCancellation()
        guard command.status == .completed else { throw VivoChemistryError.convergence("reactive Metal command failed: \(command.error?.localizedDescription ?? "unknown")") }
        let raw = outputBuffer.contents().bindMemory(to: Float.self,capacity: outputCount)
        var result: [VivoReactiveSurrogatePrediction] = []
        for frame in positions.indices {
            var energies: [Double] = [], forces: [[VivoVector3D]] = [], nearest = Double.infinity
            for member in 0..<k {
                let offset = (frame*k+member)*(d+2)
                energies.append(Double(raw[offset])); nearest = min(nearest,Double(raw[offset+1]))
                var force = [VivoVector3D](repeating: .zero,count: n)
                for feature in 0..<d {
                    let f = p.configuration.features[feature], derivative = Double(raw[offset+2+feature])
                    let contribution = geometries[frame].directions[feature]*(-derivative)
                    force[f.atomA] = force[f.atomA]+contribution; force[f.atomB] = force[f.atomB]-contribution
                }
                forces.append(force)
            }
            result.append(try VivoReactiveDeltaSurrogate.summarize(model: model,positions: positions[frame],
                energies: energies,forces: forces,nearest: nearest))
        }
        return result
    }
    private static let source = """
    #include <metal_stdlib>
    using namespace metal;
    struct Uniforms { uint frames; uint features; uint centers; uint members; };
    kernel void vivo_reactive_delta(device const float *x [[buffer(0)]],
        device const float *centers [[buffer(1)]], device const float *coefficients [[buffer(2)]],
        device const float *inverseScale2 [[buffer(3)]], device float *output [[buffer(4)]],
        constant Uniforms &u [[buffer(5)]], uint tid [[thread_position_in_grid]]) {
        if (tid >= u.frames*u.members) return;
        uint frame = tid/u.members, member = tid%u.members, out = tid*(u.features+2);
        float energy = coefficients[member*(u.centers+1)], nearest2 = INFINITY;
        for (uint j=0; j<u.features; ++j) output[out+2+j] = 0;
        for (uint c=0; c<u.centers; ++c) {
            float distance2 = 0;
            for (uint j=0; j<u.features; ++j) {
                float delta = x[frame*u.features+j]-centers[c*u.features+j];
                distance2 += delta*delta*inverseScale2[j];
            }
            nearest2 = min(nearest2,distance2);
            float weighted = exp(-0.5f*distance2)*coefficients[member*(u.centers+1)+c+1];
            energy += weighted;
            for (uint j=0; j<u.features; ++j) {
                float delta = x[frame*u.features+j]-centers[c*u.features+j];
                output[out+2+j] -= weighted*delta*inverseScale2[j];
            }
        }
        output[out] = energy; output[out+1] = sqrt(nearest2);
    }
    """
}
#endif
