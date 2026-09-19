import Foundation

/// Original UTF-8 PDB/mmCIF bytes, not caller-constructed molecular coordinates.
/// Parsing reuses the existing molecular owners. A checksum binds bytes, not
/// model correctness, upstream authenticity, or biological target identity.
public enum VivoBinderStructureSources {
    public enum Format: String, Codable, Sendable { case pdb, mmcif }

    public struct Source: Codable, Sendable {
        public let candidateID: String
        public let target: String
        public let sourceLabel: String
        public let format: Format
        public let contents: String
        public let sha256: String
        /// Exact selected chain -> expected canonical sequence. This is a
        /// declared target reference, not proof that the reference is biological truth.
        public let targetChainSequences: [String: String]
        public let interfacePlan: VivoMolecularInterface.Plan

        public init(candidateID: String, target: String, sourceLabel: String,
                    format: Format, contents: String, sha256: String,
                    targetChainSequences: [String: String], interfacePlan: VivoMolecularInterface.Plan) {
            self.candidateID = candidateID; self.target = target; self.sourceLabel = sourceLabel
            self.format = format; self.contents = contents; self.sha256 = sha256
            self.targetChainSequences = targetChainSequences; self.interfacePlan = interfacePlan
        }
    }

    public struct Input: Codable, Sendable {
        public let schemaVersion: Int
        public let sources: [Source]
        public init(sources: [Source]) { schemaVersion = 1; self.sources = sources }
    }

    /// Materializes the existing structural input without dropping alternate
    /// locations, unit cells, missing atoms, unsupported chemistry or conformers.
    /// Those remain subject to the existing interface/protein admission rules.
    public static func reconstruct(_ input: Input) throws -> VivoBinderStructuralFeatures.Input {
        guard input.schemaVersion == 1, !input.sources.isEmpty, input.sources.count <= 10_000 else {
            throw invalid("unsupported raw structure source schema/count")
        }
        var totalBytes = 0, totalAtoms = 0
        var candidates = Set<String>()
        var observations: [VivoBinderStructuralFeatures.Observation] = []
        for source in input.sources {
            guard !source.candidateID.isEmpty, candidates.insert(source.candidateID).inserted,
                  !source.target.isEmpty, !source.sourceLabel.isEmpty,
                  source.sourceLabel.utf8.count <= 1024 else {
                throw invalid("empty or duplicate candidate/target/source identity")
            }
            let bytes = Data(source.contents.utf8)
            totalBytes += bytes.count
            guard !bytes.isEmpty, bytes.count <= 16 * 1024 * 1024,
                  totalBytes <= 64 * 1024 * 1024, !bytes.contains(0) else {
                throw invalid("raw structure text exceeds bounded UTF-8 source capacity")
            }
            guard source.sha256.utf8.count == 64,
                  source.sha256.utf8.allSatisfy({ (48...57).contains($0) || (97...102).contains($0) }),
                  try VivoCanonicalJSON.fingerprint(bytes).hex == source.sha256 else {
                throw invalid("raw structure SHA-256 mismatch: \(source.candidateID)")
            }
            let selectedTargets = Set(source.interfacePlan.targetChains)
            guard !selectedTargets.isEmpty,
                  selectedTargets.count == source.interfacePlan.targetChains.count,
                  Set(source.targetChainSequences.keys) == selectedTargets,
                  source.targetChainSequences.values.allSatisfy({ !$0.isEmpty && $0.utf8.count <= 100_000 }) else {
                throw invalid("exact expected sequence required for every selected target chain")
            }
            let structure: VivoMolecularStructure
            switch source.format {
            case .pdb:
                structure = try VivoPDB.read(source.contents,
                    options: .init(identifier: source.candidateID, keepAlternateLocations: true))
            case .mmcif:
                structure = try VivoMMCIF.read(source.contents, identifier: source.candidateID)
            }
            _ = try VivoStructureValidator.validate(structure)
            totalAtoms += structure.atoms.count
            guard totalAtoms <= 1_000_000 else { throw invalid("raw structure archive atom capacity exceeded") }
            for chain in source.interfacePlan.targetChains {
                let sequence = try VivoBinderStructuralFeatures.proteinSequence(structure, chain: chain)
                guard sequence == source.targetChainSequences[chain] else {
                    throw invalid("target structure sequence differs from declared reference: \(source.candidateID)/\(chain)")
                }
            }
            observations.append(.init(candidateID: source.candidateID, target: source.target,
                sourceLabel: source.sourceLabel, structure: structure, interfacePlan: source.interfacePlan))
        }
        return .init(observations: observations)
    }

    private static func invalid(_ reason: String) -> VivoBinderBenchmark.Failure { .invalid(reason) }
}
