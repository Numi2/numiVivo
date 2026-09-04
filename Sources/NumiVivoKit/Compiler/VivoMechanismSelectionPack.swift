import CryptoKit
import Foundation

public enum VivoMechanismRole: String, Codable, CaseIterable, Sendable {
    case sensor
    case transducer
    case logic
    case temporal
    case memory
    case effector
    case communication
    case containment
    case shutdown
    case monitor
}

public enum VivoMechanismEvidenceTier: String, Codable, CaseIterable, Sendable {
    case observedTargetContext = "observed-target-context"
    case observedRelatedContext = "observed-related-context"
    case calibrated
    case inferred
    case assumed
    case hypothetical
}

public enum VivoMechanismReversibility: String, Codable, CaseIterable, Sendable {
    case reversible
    case resettable
    case conditionallyIrreversible = "conditionally-irreversible"
    case irreversible
}

public struct VivoMechanismPerformanceEnvelope: Codable, Equatable, Sendable {
    public let payloadBytes: Double
    public let cellularBurden: Double
    public let latencySeconds: Double
    public let leakProbability: Double
    public let dynamicRange: Double
    public let specificity: Double
    public let robustness: Double
    public let relativeUncertainty: Double

    public func validate() throws {
        let values = [
            payloadBytes, cellularBurden, latencySeconds, leakProbability,
            dynamicRange, specificity, robustness, relativeUncertainty
        ]
        guard values.allSatisfy(\.isFinite),
              payloadBytes >= 0,
              cellularBurden >= 0,
              latencySeconds >= 0,
              leakProbability >= 0,
              leakProbability <= 1,
              dynamicRange >= 1,
              specificity >= 0,
              specificity <= 1,
              robustness >= 0,
              robustness <= 1,
              relativeUncertainty >= 0 else {
            throw VivoMechanismSelectionError.invalidPart("performance envelope")
        }
    }
}

public struct VivoAbstractMechanismPart: Codable, Equatable, Sendable {
    public let id: String
    public let role: VivoMechanismRole
    public let acceptedInputs: [String]?
    public let producedOutputs: [String]?
    public let supportedHosts: [String]?
    public let supportedDeliveryModes: [String]?
    public let dependencies: [String]?
    public let incompatibilities: [String]?
    public let orthogonalityGroup: String?
    public let evidence: VivoMechanismEvidenceTier?
    public let reversibility: VivoMechanismReversibility?
    public let independentlyControlled: Bool?
    public let contextInsulated: Bool?
    public let resourceBuffered: Bool?
    public let performance: VivoMechanismPerformanceEnvelope

    public func validate() throws {
        guard Self.validIdentifier(id),
              id.utf8.count <= 128,
              (orthogonalityGroup?.utf8.count ?? 0) <= 1024 else {
            throw VivoMechanismSelectionError.invalidPart(id)
        }
        try performance.validate()
        for collection in [
            acceptedInputs, producedOutputs, supportedHosts,
            supportedDeliveryModes, dependencies, incompatibilities
        ] {
            guard let collection else { continue }
            guard collection.count <= 65_536,
                  Set(collection).count == collection.count,
                  collection.allSatisfy({ !$0.isEmpty && $0.utf8.count <= 1024 }) else {
                throw VivoMechanismSelectionError.invalidPart(id)
            }
        }
        if dependencies?.contains(id) == true || incompatibilities?.contains(id) == true {
            throw VivoMechanismSelectionError.invalidPart(id)
        }
    }

    private static func validIdentifier(_ value: String) -> Bool {
        guard let first = value.unicodeScalars.first,
              CharacterSet.letters.union(CharacterSet(charactersIn: "_")).contains(first) else {
            return false
        }
        let allowed = CharacterSet.alphanumerics.union(CharacterSet(charactersIn: "_.-"))
        return value.unicodeScalars.dropFirst().allSatisfy(allowed.contains)
    }
}

public struct VivoMechanismLibrarySource: Codable, Equatable, Sendable {
    public let schemaVersion: UInt32
    public let parts: [VivoAbstractMechanismPart]

    public func validate() throws {
        guard schemaVersion == 1,
              !parts.isEmpty,
              parts.count <= 1_000_000,
              Set(parts.map(\.id)).count == parts.count else {
            throw VivoMechanismSelectionError.invalidLibrary
        }
        for part in parts { try part.validate() }
        let known = Set(parts.map(\.id))
        for part in parts {
            for dependency in part.dependencies ?? [] where !known.contains(dependency) {
                throw VivoMechanismSelectionError.unknownDependency(part.id, dependency)
            }
            for incompatibility in part.incompatibilities ?? [] where !known.contains(incompatibility) {
                throw VivoMechanismSelectionError.unknownIncompatibility(part.id, incompatibility)
            }
        }
    }
}

public struct VivoMechanismProblemSummary: Codable, Equatable, Sendable {
    public let id: String
    public let host: String
    public let deliveryMode: String
    public let requireIndependentShutdown: Bool
    public let requireMonitor: Bool
    public let requireContextInsulation: Bool
    public let requireResourceBuffering: Bool
    public let requireDistinctOrthogonalityGroups: Bool
}

public struct VivoMechanismCandidateMetrics: Codable, Equatable, Sendable {
    public let payloadBytes: Double
    public let cellularBurden: Double
    public let latencySeconds: Double
    public let leakProbability: Double
    public let dynamicRange: Double
    public let specificity: Double
    public let robustness: Double
    public let relativeUncertainty: Double
}

public struct VivoMechanismSynthesisCandidate: Codable, Equatable, Sendable {
    public let fingerprint: String
    public let objective: Double
    public let weakestEvidence: VivoMechanismEvidenceTier
    public let worstReversibility: VivoMechanismReversibility
    public let metrics: VivoMechanismCandidateMetrics
    public let selectedPartIDs: [String]
    public let slotAssignments: [String: String?]
}

public struct VivoMechanismRejectionCount: Codable, Equatable, Sendable {
    public let reason: String
    public let count: UInt64
}

public struct VivoMechanismSynthesisResult: Codable, Equatable, Sendable {
    public let schemaVersion: UInt32
    public let visitedNodes: UInt64
    public let completedAssignments: UInt64
    public let searchBudgetExhausted: Bool
    public let candidates: [VivoMechanismSynthesisCandidate]
    public let rejectionCounts: [VivoMechanismRejectionCount]
}

public struct VivoMechanismSelectionPack: Codable, Equatable, Sendable {
    public let schemaVersion: UInt32
    public let problemSourceFingerprint: String
    public let librarySourceFingerprint: String
    public let diagnosticsFingerprint: String
    public let problem: VivoMechanismProblemSummary
    public let candidate: VivoMechanismSynthesisCandidate
    public let selectedParts: [VivoAbstractMechanismPart]
    public let visitedNodes: UInt64
    public let completedAssignments: UInt64
    public let searchBudgetExhausted: Bool
    public let rejectionCounts: [VivoMechanismRejectionCount]
    public let sequenceFree: Bool
    public let fingerprint: String

    public func validate() throws {
        guard schemaVersion == 1,
              Self.isFingerprint(problemSourceFingerprint),
              Self.isFingerprint(librarySourceFingerprint),
              Self.isFingerprint(diagnosticsFingerprint),
              Self.isFingerprint(candidate.fingerprint),
              Self.isFingerprint(fingerprint),
              sequenceFree,
              !selectedParts.isEmpty,
              Set(selectedParts.map(\.id)).count == selectedParts.count,
              selectedParts.map(\.id).sorted() == candidate.selectedPartIDs.sorted() else {
            throw VivoMechanismSelectionError.invalidPack("identity or part membership")
        }
        for part in selectedParts { try part.validate() }
        let selected = Set(selectedParts.map(\.id))
        for (slot, part) in candidate.slotAssignments {
            guard !slot.isEmpty else {
                throw VivoMechanismSelectionError.invalidPack("empty slot identifier")
            }
            if let part, !selected.contains(part) {
                throw VivoMechanismSelectionError.invalidPack(
                    "slot \(slot) references unselected part \(part)"
                )
            }
        }
        if problem.requireIndependentShutdown {
            guard selectedParts.contains(where: {
                $0.role == .shutdown &&
                $0.independentlyControlled == true &&
                ($0.reversibility ?? .reversible) != .irreversible &&
                ($0.reversibility ?? .reversible) != .conditionallyIrreversible
            }) else {
                throw VivoMechanismSelectionError.invalidPack("independent shutdown is absent")
            }
        }
        if problem.requireMonitor,
           !selectedParts.contains(where: { $0.role == .monitor }) {
            throw VivoMechanismSelectionError.invalidPack("monitor is absent")
        }
        let unsigned = VivoMechanismSelectionPack(
            schemaVersion: schemaVersion,
            problemSourceFingerprint: problemSourceFingerprint,
            librarySourceFingerprint: librarySourceFingerprint,
            diagnosticsFingerprint: diagnosticsFingerprint,
            problem: problem,
            candidate: candidate,
            selectedParts: selectedParts,
            visitedNodes: visitedNodes,
            completedAssignments: completedAssignments,
            searchBudgetExhausted: searchBudgetExhausted,
            rejectionCounts: rejectionCounts,
            sequenceFree: sequenceFree,
            fingerprint: ""
        )
        let actual = try Self.canonicalFingerprint(unsigned)
        guard actual == fingerprint else {
            throw VivoMechanismSelectionError.fingerprintMismatch(
                expected: fingerprint,
                actual: actual
            )
        }
    }

    fileprivate static func isFingerprint(_ value: String) -> Bool {
        value.utf8.count == 64 && value.allSatisfy {
            $0.isHexDigit && !$0.isUppercase
        }
    }

    fileprivate static func canonicalFingerprint<T: Encodable>(_ value: T) throws -> String {
        let encoder = JSONEncoder()
        encoder.outputFormatting = [.sortedKeys, .withoutEscapingSlashes]
        do {
            return SHA256.hash(data: try encoder.encode(value)).map {
                String(format: "%02x", $0)
            }.joined()
        } catch {
            throw VivoMechanismSelectionError.encoding(error.localizedDescription)
        }
    }
}

public struct VivoMechanismSelectionCompilation: Sendable {
    public let pack: VivoMechanismSelectionPack
    public let synthesis: VivoMechanismSynthesisResult
    public let diagnostics: Data
}

public enum VivoMechanismSelectionError: Error, LocalizedError, Sendable {
    case invalidLibrary
    case invalidPart(String)
    case unknownDependency(String, String)
    case unknownIncompatibility(String, String)
    case invalidProblem(String)
    case synthesisFailed(VivoNativeStatus, String)
    case candidateUnavailable(Int)
    case invalidCandidate(String)
    case invalidPack(String)
    case fingerprintMismatch(expected: String, actual: String)
    case decoding(String)
    case encoding(String)

    public var errorDescription: String? {
        switch self {
        case .invalidLibrary: return "abstract mechanism library is invalid"
        case .invalidPart(let part): return "abstract mechanism part is invalid: \(part)"
        case .unknownDependency(let part, let dependency): return "mechanism part \(part) depends on unknown part \(dependency)"
        case .unknownIncompatibility(let part, let other): return "mechanism part \(part) conflicts with unknown part \(other)"
        case .invalidProblem(let reason): return "mechanism synthesis problem is invalid: \(reason)"
        case .synthesisFailed(let status, let diagnostics): return "mechanism synthesis failed with \(status): \(diagnostics)"
        case .candidateUnavailable(let index): return "mechanism candidate index \(index) is unavailable"
        case .invalidCandidate(let reason): return "mechanism candidate is invalid: \(reason)"
        case .invalidPack(let reason): return "mechanism selection pack is invalid: \(reason)"
        case .fingerprintMismatch(let expected, let actual): return "mechanism pack fingerprint mismatch: expected \(expected), found \(actual)"
        case .decoding(let reason): return "mechanism document decoding failed: \(reason)"
        case .encoding(let reason): return "mechanism document encoding failed: \(reason)"
        }
    }
}

public struct VivoMechanismSelectionCompiler: Sendable {
    public init() {}

    public func compile(
        problemSource: Data,
        librarySource: Data,
        candidateIndex: Int = 0,
        maximumSolutions: UInt32? = nil,
        maximumVisitedNodes: UInt64? = nil
    ) throws -> VivoMechanismSelectionCompilation {
        let library: VivoMechanismLibrarySource
        do {
            library = try JSONDecoder().decode(
                VivoMechanismLibrarySource.self,
                from: librarySource
            )
        } catch {
            throw VivoMechanismSelectionError.decoding(error.localizedDescription)
        }
        try library.validate()
        let problem = try decodeProblemSummary(problemSource)

        let invocation = VivoNativeCompilerBridge.synthesizeMechanisms(
            problem: problemSource,
            library: librarySource,
            maximumSolutions: maximumSolutions,
            maximumVisitedNodes: maximumVisitedNodes
        )
        guard invocation.succeeded else {
            throw VivoMechanismSelectionError.synthesisFailed(
                invocation.status,
                invocation.diagnosticsString()
            )
        }
        let synthesis: VivoMechanismSynthesisResult
        do {
            synthesis = try JSONDecoder().decode(
                VivoMechanismSynthesisResult.self,
                from: invocation.primary
            )
        } catch {
            throw VivoMechanismSelectionError.decoding(error.localizedDescription)
        }
        guard synthesis.schemaVersion == 1,
              candidateIndex >= 0,
              candidateIndex < synthesis.candidates.count else {
            throw VivoMechanismSelectionError.candidateUnavailable(candidateIndex)
        }
        let candidate = synthesis.candidates[candidateIndex]
        guard VivoMechanismSelectionPack.isFingerprint(candidate.fingerprint),
              candidate.objective.isFinite,
              Set(candidate.selectedPartIDs).count == candidate.selectedPartIDs.count else {
            throw VivoMechanismSelectionError.invalidCandidate("identity, objective or part list")
        }
        let partByID = Dictionary(uniqueKeysWithValues: library.parts.map { ($0.id, $0) })
        let selectedParts = try candidate.selectedPartIDs.map { id -> VivoAbstractMechanismPart in
            guard let part = partByID[id] else {
                throw VivoMechanismSelectionError.invalidCandidate("unknown selected part \(id)")
            }
            return part
        }.sorted { $0.id < $1.id }

        let unsigned = VivoMechanismSelectionPack(
            schemaVersion: 1,
            problemSourceFingerprint: fingerprint(problemSource),
            librarySourceFingerprint: fingerprint(librarySource),
            diagnosticsFingerprint: fingerprint(invocation.diagnostics),
            problem: problem,
            candidate: candidate,
            selectedParts: selectedParts,
            visitedNodes: synthesis.visitedNodes,
            completedAssignments: synthesis.completedAssignments,
            searchBudgetExhausted: synthesis.searchBudgetExhausted,
            rejectionCounts: synthesis.rejectionCounts,
            sequenceFree: true,
            fingerprint: ""
        )
        let pack = VivoMechanismSelectionPack(
            schemaVersion: unsigned.schemaVersion,
            problemSourceFingerprint: unsigned.problemSourceFingerprint,
            librarySourceFingerprint: unsigned.librarySourceFingerprint,
            diagnosticsFingerprint: unsigned.diagnosticsFingerprint,
            problem: unsigned.problem,
            candidate: unsigned.candidate,
            selectedParts: unsigned.selectedParts,
            visitedNodes: unsigned.visitedNodes,
            completedAssignments: unsigned.completedAssignments,
            searchBudgetExhausted: unsigned.searchBudgetExhausted,
            rejectionCounts: unsigned.rejectionCounts,
            sequenceFree: unsigned.sequenceFree,
            fingerprint: try VivoMechanismSelectionPack.canonicalFingerprint(unsigned)
        )
        try pack.validate()
        return .init(pack: pack, synthesis: synthesis, diagnostics: invocation.diagnostics)
    }

    private func decodeProblemSummary(_ data: Data) throws -> VivoMechanismProblemSummary {
        struct Source: Decodable {
            let id: String
            let constraints: Constraints
            struct Constraints: Decodable {
                let host: String
                let deliveryMode: String
                let requireIndependentShutdown: Bool?
                let requireMonitor: Bool?
                let requireContextInsulation: Bool?
                let requireResourceBuffering: Bool?
                let requireDistinctOrthogonalityGroups: Bool?
            }
        }
        let source: Source
        do { source = try JSONDecoder().decode(Source.self, from: data) }
        catch { throw VivoMechanismSelectionError.decoding(error.localizedDescription) }
        guard !source.id.isEmpty,
              !source.constraints.host.isEmpty,
              !source.constraints.deliveryMode.isEmpty else {
            throw VivoMechanismSelectionError.invalidProblem("identity, host or delivery mode")
        }
        return .init(
            id: source.id,
            host: source.constraints.host,
            deliveryMode: source.constraints.deliveryMode,
            requireIndependentShutdown: source.constraints.requireIndependentShutdown ?? true,
            requireMonitor: source.constraints.requireMonitor ?? true,
            requireContextInsulation: source.constraints.requireContextInsulation ?? false,
            requireResourceBuffering: source.constraints.requireResourceBuffering ?? false,
            requireDistinctOrthogonalityGroups: source.constraints.requireDistinctOrthogonalityGroups ?? false
        )
    }

    private func fingerprint(_ data: Data) -> String {
        SHA256.hash(data: data).map { String(format: "%02x", $0) }.joined()
    }
}
