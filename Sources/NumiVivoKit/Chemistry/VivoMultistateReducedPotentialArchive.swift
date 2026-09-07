import Foundation

public struct VivoMultistateReducedPotentialArchiveRequest: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/multistate-reduced-potential-archive/v1"
    public var schema: String
    public var identifier: String
    public var temperatureK: Double
    public var stateIdentifiers: [String]
    public var physicalManifoldFingerprint: VivoFingerprint
    public var hamiltonianFingerprints: [VivoFingerprint]
    public var initialPhysicalStateFingerprints: [VivoFingerprint]
    public var equilibrationSteps: UInt64
    public var productionSteps: UInt64
    public var sampleEvery: UInt64
    public var checkpointBlockSteps: UInt64
    public var maximumMatrixElements: Int
    public var mbar: VivoMBARConfiguration

    public init(identifier: String, temperatureK: Double,
                stateIdentifiers: [String],
                physicalManifoldFingerprint: VivoFingerprint,
                hamiltonianFingerprints: [VivoFingerprint],
                initialPhysicalStateFingerprints: [VivoFingerprint],
                equilibrationSteps: UInt64,
                productionSteps: UInt64,
                sampleEvery: UInt64,
                checkpointBlockSteps: UInt64 = 5_000,
                maximumMatrixElements: Int = 100_000_000,
                mbar: VivoMBARConfiguration = .init()) {
        self.schema = Self.schema
        self.identifier = identifier
        self.temperatureK = temperatureK
        self.stateIdentifiers = stateIdentifiers
        self.physicalManifoldFingerprint = physicalManifoldFingerprint
        self.hamiltonianFingerprints = hamiltonianFingerprints
        self.initialPhysicalStateFingerprints = initialPhysicalStateFingerprints
        self.equilibrationSteps = equilibrationSteps
        self.productionSteps = productionSteps
        self.sampleEvery = sampleEvery
        self.checkpointBlockSteps = checkpointBlockSteps
        self.maximumMatrixElements = maximumMatrixElements
        self.mbar = mbar
    }

    public func validate() throws {
        guard schema == Self.schema,
              !identifier.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty,
              identifier.utf8.count <= 512,
              temperatureK.isFinite, temperatureK > 0,
              stateIdentifiers.count >= 2, stateIdentifiers.count <= 4096,
              Set(stateIdentifiers).count == stateIdentifiers.count,
              hamiltonianFingerprints.count == stateIdentifiers.count,
              initialPhysicalStateFingerprints.count == stateIdentifiers.count,
              productionSteps > 0, sampleEvery > 0,
              productionSteps % sampleEvery == 0,
              checkpointBlockSteps > 0,
              checkpointBlockSteps % sampleEvery == 0,
              maximumMatrixElements > 0,
              maximumMatrixElements <= 500_000_000 else {
            throw VivoChemistryError.invalid("multistate reduced-potential archive request")
        }
        try mbar.validate()
        let samplesPerOrigin = productionSteps / sampleEvery
        let totalSamples = samplesPerOrigin.multipliedReportingOverflow(by: UInt64(stateIdentifiers.count))
        guard !totalSamples.overflow, totalSamples.partialValue <= UInt64(Int.max) else {
            throw VivoChemistryError.resourceLimit("multistate reduced-potential sample count overflow")
        }
        let matrix = Int(totalSamples.partialValue).multipliedReportingOverflow(by: stateIdentifiers.count)
        guard !matrix.overflow, matrix.partialValue <= maximumMatrixElements,
              matrix.partialValue <= mbar.maximumWorkElements else {
            throw VivoChemistryError.resourceLimit("multistate reduced-potential matrix exceeds declared capacity")
        }
    }

    public func fingerprint() throws -> VivoFingerprint {
        try validate()
        return try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(self))
    }
}

public struct VivoMultistateReducedPotentialChunk: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/multistate-reduced-potential-chunk/v1"
    public let schema: String
    public let requestFingerprint: VivoFingerprint
    public let originStateIndex: Int
    public let firstSampleIndex: Int
    public let samples: [VivoMBARReducedPotentialSample]

    public init(requestFingerprint: VivoFingerprint, originStateIndex: Int,
                firstSampleIndex: Int, samples: [VivoMBARReducedPotentialSample]) {
        self.schema = Self.schema
        self.requestFingerprint = requestFingerprint
        self.originStateIndex = originStateIndex
        self.firstSampleIndex = firstSampleIndex
        self.samples = samples
    }
}

private struct VivoMultistateReducedPotentialCursor: Codable, Sendable, Equatable {
    static let schema = "numivivo.org/multistate-reduced-potential-checkpoint/v1"
    var schema: String
    var requestFingerprint: VivoFingerprint
    var completedOrigins: Int
    var currentOriginProductionSteps: UInt64
    var currentOriginSampleCount: Int
    var currentPhysicalState: VivoFingerprint?
    var chunks: [VivoFingerprint]
}

public enum VivoMultistateReducedPotentialArchiveStatus: String, Codable, Sendable {
    case converged
    case analysisNotConverged
    case rejected
    case cancelled
}

public struct VivoMultistateReducedPotentialArchiveResult: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/multistate-reduced-potential-archive-result/v1"
    public let schema: String
    public let requestFingerprint: VivoFingerprint
    public let status: VivoMultistateReducedPotentialArchiveStatus
    public let checkpointFingerprint: VivoFingerprint
    public let completedOrigins: Int
    public let currentOriginProductionSteps: UInt64
    public let chunkFingerprints: [VivoFingerprint]
    public let mbar: VivoMBARResult?
    public let diagnostic: String?
}

/// Durable sampling of an explicit multistate Hamiltonian set. Every accepted
/// frame from origin state i is evaluated under every state Hamiltonian at the
/// exact same physical conformation, producing the reduced-potential matrix MBAR
/// requires. Only complete checkpoint blocks are committed, so cancellation can
/// resume without interpreting a partial block as durable sampling evidence.
public enum VivoMultistateReducedPotentialArchiveRunner {
    private static let gasConstantKJ = 0.00831446261815324
    public static let interpretation = "Artifact-backed multistate equilibrium sampling over an explicit common physical manifold. Each decorrelation interval is propagated under its origin Hamiltonian and the resulting conformation is evaluated under every declared Hamiltonian. The archive supplies physical reduced potentials only; proton-reservoir terms, state enumeration and claims of decorrelation/convergence remain separate qualifications."

    public static func checkpointReferenceName(requestFingerprint: VivoFingerprint) -> String {
        "multistate-reduced-potential/\(requestFingerprint.hex)/checkpoint"
    }

    private static func put<T: Encodable>(_ value: T, kind: String,
                                           store: VivoArtifactStore) async throws -> VivoFingerprint {
        try await store.put(data: VivoCanonicalJSON.encode(value), kind: kind,
                            mediaType: "application/json").fingerprint
    }

    private static func read<T: Decodable>(_ type: T.Type, id: VivoFingerprint,
                                            kind: String, store: VivoArtifactStore) async throws -> T {
        let descriptor = try await store.descriptor(for: id)
        guard descriptor.kind == kind else {
            throw VivoChemistryError.invalid("multistate reduced-potential artifact kind mismatch")
        }
        return try VivoCanonicalJSON.decode(type, from: await store.data(for: id, verify: true))
    }

    public static func run(_ request: VivoMultistateReducedPotentialArchiveRequest,
                           initialPhysicalStates: [VivoConstantPHPhysicalState],
                           executableStates: [VivoConstantPHExecutableState],
                           store: VivoArtifactStore,
                           resumeFrom: VivoFingerprint? = nil) async throws -> VivoMultistateReducedPotentialArchiveResult {
        try request.validate()
        guard initialPhysicalStates.count == request.stateIdentifiers.count,
              executableStates.count == request.stateIdentifiers.count else {
            throw VivoChemistryError.invalid("multistate archive state/resource count")
        }
        for state in initialPhysicalStates { try state.validate() }
        let executableIDs = executableStates.map(\.identifier)
        guard executableIDs == request.stateIdentifiers,
              executableStates.allSatisfy({ $0.physicalManifoldFingerprint == request.physicalManifoldFingerprint }),
              executableStates.map(\.hamiltonianFingerprint) == request.hamiltonianFingerprints else {
            throw VivoChemistryError.invalid("multistate archive executable identity or physical manifold")
        }
        let initialIDs = try initialPhysicalStates.map { try $0.fingerprint() }
        guard initialIDs == request.initialPhysicalStateFingerprints else {
            throw VivoChemistryError.invalid("multistate archive initial physical-state identity")
        }
        let requestID = try await put(request, kind: "multistate-reduced-potential-request", store: store)
        let samplesPerOrigin = Int(request.productionSteps / request.sampleEvery)
        let rt = gasConstantKJ * request.temperatureK

        var cursor: VivoMultistateReducedPotentialCursor
        if let resumeFrom {
            cursor = try await read(VivoMultistateReducedPotentialCursor.self, id: resumeFrom,
                                    kind: "multistate-reduced-potential-checkpoint", store: store)
            try validate(cursor: cursor, request: request, requestID: requestID,
                         samplesPerOrigin: samplesPerOrigin)
            if let physicalID = cursor.currentPhysicalState {
                let physical = try await read(VivoConstantPHPhysicalState.self, id: physicalID,
                                              kind: "multistate-physical-state", store: store)
                try physical.validate()
            }
            for chunkID in cursor.chunks {
                let chunk = try await read(VivoMultistateReducedPotentialChunk.self, id: chunkID,
                                           kind: "multistate-reduced-potential-chunk", store: store)
                try validate(chunk: chunk, request: request, requestID: requestID,
                             samplesPerOrigin: samplesPerOrigin)
            }
        } else {
            cursor = .init(schema: VivoMultistateReducedPotentialCursor.schema,
                           requestFingerprint: requestID, completedOrigins: 0,
                           currentOriginProductionSteps: 0, currentOriginSampleCount: 0,
                           currentPhysicalState: nil, chunks: [])
        }
        var durable = try await put(cursor, kind: "multistate-reduced-potential-checkpoint", store: store)

        func receipt(_ status: VivoMultistateReducedPotentialArchiveStatus,
                     mbar: VivoMBARResult? = nil,
                     diagnostic: String? = nil) -> VivoMultistateReducedPotentialArchiveResult {
            .init(schema: VivoMultistateReducedPotentialArchiveResult.schema,
                  requestFingerprint: requestID, status: status,
                  checkpointFingerprint: durable,
                  completedOrigins: cursor.completedOrigins,
                  currentOriginProductionSteps: cursor.currentOriginProductionSteps,
                  chunkFingerprints: cursor.chunks, mbar: mbar, diagnostic: diagnostic)
        }

        do {
            while cursor.completedOrigins < request.stateIdentifiers.count {
                try Task.checkCancellation()
                let origin = cursor.completedOrigins
                let runtime = executableStates[origin]
                var physical: VivoConstantPHPhysicalState
                if let physicalID = cursor.currentPhysicalState {
                    physical = try await read(VivoConstantPHPhysicalState.self, id: physicalID,
                                              kind: "multistate-physical-state", store: store)
                } else {
                    physical = initialPhysicalStates[origin]
                    if request.equilibrationSteps > 0 {
                        physical = try await runtime.propagate(physical, request.equilibrationSteps)
                        try physical.validate()
                    }
                    cursor.currentPhysicalState = try await put(physical, kind: "multistate-physical-state", store: store)
                }

                let remaining = request.productionSteps - cursor.currentOriginProductionSteps
                let block = min(remaining, request.checkpointBlockSteps)
                guard block > 0, block % request.sampleEvery == 0 else {
                    throw VivoChemistryError.invalid("multistate archive checkpoint block alignment")
                }
                let blockSamples = Int(block / request.sampleEvery)
                var samples: [VivoMBARReducedPotentialSample] = []
                samples.reserveCapacity(blockSamples)
                let firstSample = cursor.currentOriginSampleCount
                for local in 0..<blockSamples {
                    try Task.checkCancellation()
                    physical = try await runtime.propagate(physical, request.sampleEvery)
                    try physical.validate()
                    let energies = try await evaluateAll(executableStates, at: physical)
                    guard energies.allSatisfy(\.isFinite) else {
                        throw VivoChemistryError.convergence("nonfinite multistate cross-Hamiltonian energy")
                    }
                    let reduced = energies.map { $0 / rt }
                    let globalSample = firstSample + local
                    samples.append(.init(identifier: "\(request.identifier)/\(origin)/\(globalSample)",
                                         originStateIndex: origin,
                                         reducedPotentials: reduced))
                }
                let chunk = VivoMultistateReducedPotentialChunk(requestFingerprint: requestID,
                    originStateIndex: origin, firstSampleIndex: firstSample, samples: samples)
                let chunkID = try await put(chunk, kind: "multistate-reduced-potential-chunk", store: store)
                cursor.chunks.append(chunkID)
                cursor.currentOriginProductionSteps += block
                cursor.currentOriginSampleCount += blockSamples

                if cursor.currentOriginProductionSteps == request.productionSteps {
                    guard cursor.currentOriginSampleCount == samplesPerOrigin else {
                        throw VivoChemistryError.invalid("multistate archive completed-origin sample accounting")
                    }
                    cursor.completedOrigins += 1
                    cursor.currentOriginProductionSteps = 0
                    cursor.currentOriginSampleCount = 0
                    cursor.currentPhysicalState = nil
                } else {
                    cursor.currentPhysicalState = try await put(physical, kind: "multistate-physical-state", store: store)
                }
                try validate(cursor: cursor, request: request, requestID: requestID,
                             samplesPerOrigin: samplesPerOrigin)
                durable = try await put(cursor, kind: "multistate-reduced-potential-checkpoint", store: store)
            }

            var allSamples: [VivoMBARReducedPotentialSample] = []
            allSamples.reserveCapacity(samplesPerOrigin * request.stateIdentifiers.count)
            for chunkID in cursor.chunks {
                let chunk = try await read(VivoMultistateReducedPotentialChunk.self, id: chunkID,
                                           kind: "multistate-reduced-potential-chunk", store: store)
                try validate(chunk: chunk, request: request, requestID: requestID,
                             samplesPerOrigin: samplesPerOrigin)
                allSamples.append(contentsOf: chunk.samples)
            }
            guard allSamples.count == samplesPerOrigin * request.stateIdentifiers.count,
                  Set(allSamples.map(\.identifier)).count == allSamples.count else {
                throw VivoChemistryError.invalid("multistate archive final sample accounting")
            }
            let mbar = try VivoMultistateMBAR.solve(samples: allSamples,
                                                    stateCount: request.stateIdentifiers.count,
                                                    configuration: request.mbar)
            _ = try await put(mbar, kind: "multistate-mbar-result", store: store)
            return receipt(mbar.converged ? .converged : .analysisNotConverged, mbar: mbar)
        } catch is CancellationError {
            return receipt(.cancelled, diagnostic: "cancelled; resume from the last complete checkpoint block")
        } catch {
            return receipt(.rejected, diagnostic: String(describing: error))
        }
    }

    private static func evaluateAll(_ states: [VivoConstantPHExecutableState],
                                    at physical: VivoConstantPHPhysicalState) async throws -> [Double] {
        var energies = [Double](repeating: 0, count: states.count)
        try await withThrowingTaskGroup(of: (Int, Double).self) { group in
            for (index, state) in states.enumerated() {
                group.addTask {
                    (index, try await state.potentialEnergyKJPerMol(physical))
                }
            }
            for try await (index, energy) in group { energies[index] = energy }
        }
        return energies
    }

    private static func validate(cursor: VivoMultistateReducedPotentialCursor,
                                 request: VivoMultistateReducedPotentialArchiveRequest,
                                 requestID: VivoFingerprint,
                                 samplesPerOrigin: Int) throws {
        guard cursor.schema == VivoMultistateReducedPotentialCursor.schema,
              cursor.requestFingerprint == requestID,
              cursor.completedOrigins >= 0,
              cursor.completedOrigins <= request.stateIdentifiers.count,
              cursor.currentOriginProductionSteps <= request.productionSteps,
              cursor.currentOriginProductionSteps % request.sampleEvery == 0,
              cursor.currentOriginSampleCount == Int(cursor.currentOriginProductionSteps / request.sampleEvery),
              cursor.currentOriginSampleCount <= samplesPerOrigin,
              (cursor.completedOrigins < request.stateIdentifiers.count || cursor.currentOriginProductionSteps == 0),
              (cursor.currentOriginProductionSteps == 0) == (cursor.currentPhysicalState == nil) else {
            throw VivoChemistryError.invalid("multistate reduced-potential checkpoint identity or accounting")
        }
    }

    private static func validate(chunk: VivoMultistateReducedPotentialChunk,
                                 request: VivoMultistateReducedPotentialArchiveRequest,
                                 requestID: VivoFingerprint,
                                 samplesPerOrigin: Int) throws {
        guard chunk.schema == VivoMultistateReducedPotentialChunk.schema,
              chunk.requestFingerprint == requestID,
              chunk.originStateIndex >= 0,
              chunk.originStateIndex < request.stateIdentifiers.count,
              chunk.firstSampleIndex >= 0,
              chunk.firstSampleIndex + chunk.samples.count <= samplesPerOrigin,
              !chunk.samples.isEmpty,
              chunk.samples.allSatisfy({ sample in
                  sample.originStateIndex == chunk.originStateIndex &&
                  sample.reducedPotentials.count == request.stateIdentifiers.count &&
                  sample.reducedPotentials.allSatisfy(\.isFinite)
              }) else {
            throw VivoChemistryError.invalid("multistate reduced-potential chunk identity or dimensions")
        }
    }
}
