import Foundation
import NumiVivoKit

/// CLI document only; the library's immutable admission limits stay unchanged.
struct VivoMolecularSamplingCLILimits: Decodable {
    static let schema = "numivivo.org/molecular-sampling-read-limits/v1"
    let value: VivoMolecularSamplingReadLimits

    private struct Key: CodingKey {
        let stringValue: String
        var intValue: Int? { nil }
        init(stringValue: String) { self.stringValue = stringValue }
        init?(intValue: Int) { return nil }
    }

    init(from decoder: Decoder) throws {
        let input = try decoder.container(keyedBy: Key.self)
        let fields: Set<String> = ["schema", "maximumRequestBytes", "maximumCursorBytes", "maximumCheckpointBytes",
            "maximumDiagnosticBytes", "maximumReadBytes", "maximumReplicas", "maximumParticles", "maximumScalarElements",
            "maximumAcceptedSteps", "maximumDiagnosticEvaluations", "maximumAutocorrelationWork", "maximumTrajectoryChunks"]
        guard Set(input.allKeys.map(\.stringValue)).isSubset(of: fields),
              try input.decode(String.self, forKey: Key(stringValue: "schema")) == Self.schema else {
            throw VivoChemistryError.invalid("sampling read limits schema or unknown field")
        }
        let defaults = VivoMolecularSamplingReadLimits()
        func positive(_ name: String, _ fallback: Int) throws -> Int {
            let key = Key(stringValue: name)
            let result: Int
            if input.contains(key) { result = try input.decode(Int.self, forKey: key) } else { result = fallback }
            guard result > 0 else { throw VivoChemistryError.invalid("sampling read limit \(name) must be positive") }
            return result
        }
        func positive64(_ name: String, _ fallback: UInt64) throws -> UInt64 {
            let key = Key(stringValue: name)
            let result: UInt64
            if input.contains(key) { result = try input.decode(UInt64.self, forKey: key) } else { result = fallback }
            guard result > 0 else { throw VivoChemistryError.invalid("sampling read limit \(name) must be positive") }
            return result
        }
        value = try .init(
            maximumRequestBytes: positive("maximumRequestBytes", defaults.maximumRequestBytes),
            maximumCursorBytes: positive("maximumCursorBytes", defaults.maximumCursorBytes),
            maximumCheckpointBytes: positive("maximumCheckpointBytes", defaults.maximumCheckpointBytes),
            maximumDiagnosticBytes: positive("maximumDiagnosticBytes", defaults.maximumDiagnosticBytes),
            maximumReadBytes: positive64("maximumReadBytes", defaults.maximumReadBytes),
            maximumReplicas: positive("maximumReplicas", defaults.maximumReplicas),
            maximumParticles: positive("maximumParticles", defaults.maximumParticles),
            maximumScalarElements: positive("maximumScalarElements", defaults.maximumScalarElements),
            maximumAcceptedSteps: positive64("maximumAcceptedSteps", defaults.maximumAcceptedSteps),
            maximumDiagnosticEvaluations: positive("maximumDiagnosticEvaluations", defaults.maximumDiagnosticEvaluations),
            maximumAutocorrelationWork: positive64("maximumAutocorrelationWork", defaults.maximumAutocorrelationWork),
            maximumTrajectoryChunks: input.decodeIfPresent(UInt64.self, forKey: Key(stringValue: "maximumTrajectoryChunks")))
    }

    static func load(_ path: String?) throws -> VivoMolecularSamplingReadLimits {
        guard let path else { return .init() }
        return try VivoKineticsDocumentIO.read(Self.self, from: URL(fileURLWithPath: path), maximumBytes: 64 * 1024).value
    }

    static let help = """
    --read-limits FILE overrides reader admission limits using a bounded JSON file:
      {"schema":"numivivo.org/molecular-sampling-read-limits/v1","maximumReadBytes":8589934592}
    Optional fields (omitted fields use the library defaults): maximumRequestBytes,
    maximumCursorBytes, maximumCheckpointBytes, maximumDiagnosticBytes (each 134217728),
    maximumReadBytes (2147483648), maximumReplicas (64), maximumParticles (1000000),
    maximumScalarElements (2000000), maximumAcceptedSteps (10000000000),
    maximumDiagnosticEvaluations (101), maximumAutocorrelationWork (1000000000),
    maximumTrajectoryChunks (null, no separate chunk ceiling; zero rejects nonempty
    archives). Unknown fields, negative limits and zero positive-only limits reject.
    The read budget includes conservative trajectory object ceilings, not measured
    RSS; store descriptor limits remain separate. Limits do not change MD physics.
    """
}
