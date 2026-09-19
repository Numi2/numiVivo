import Foundation

/// Joins computed geometry to experimental candidate identity. No arbitrary numeric
/// feature map is admitted, and no outcome is used in geometry computation.
public enum VivoBinderStructuralFeatures {
    public struct Observation: Codable, Sendable {
        public let candidateID: String
        public let target: String
        /// Caller-declared origin; not an authenticated model execution receipt.
        public let sourceLabel: String
        public let structure: VivoMolecularStructure
        public let interfacePlan: VivoMolecularInterface.Plan
        public init(candidateID: String, target: String, sourceLabel: String,
                    structure: VivoMolecularStructure, interfacePlan: VivoMolecularInterface.Plan) {
            self.candidateID = candidateID; self.target = target; self.sourceLabel = sourceLabel
            self.structure = structure; self.interfacePlan = interfacePlan
        }
    }
    public struct Input: Codable, Sendable {
        public let schemaVersion: Int
        public let observations: [Observation]
        public init(observations: [Observation]) { schemaVersion = 1; self.observations = observations }
    }
    public struct Augmentation: Codable, Sendable {
        public let schemaVersion: Int
        public let dataset: VivoBinderBenchmark.Dataset
        public let interfaces: [String: VivoMolecularInterface.Report]
        public let sourceLabels: [String: String]
        public let unavailableCandidateIDs: [String]
        public let limitations: [String]
    }
    public static let featureNames: Set<String> = [
        "numi.interface.contactAtomPairs", "numi.interface.contactResiduePairs",
        "numi.interface.binderContactFraction", "numi.interface.targetContactFraction",
        "numi.interface.shortDistanceFraction", "numi.interface.minimumDistanceNM",
        "numi.interface.centroidDistanceNM", "numi.interface.binderGeometricRgNM"
    ]

    public static func augment(_ dataset: VivoBinderBenchmark.Dataset, input: Input) throws -> Augmentation {
        try VivoBinderBenchmark.validate(dataset)
        guard input.schemaVersion == 1, !input.observations.isEmpty, input.observations.count <= 10_000 else {
            throw invalid("unsupported structural input schema/count")
        }
        let records = Dictionary(uniqueKeysWithValues: dataset.records.map { ($0.id, $0) })
        var reports: [String: VivoMolecularInterface.Report] = [:], labels: [String: String] = [:]
        var totalAtoms = 0, totalWork = 0
        let profile = input.observations[0].interfacePlan
        var targetSequences: [String: [String]] = [:]
        for observation in input.observations {
            guard let row = records[observation.candidateID], row.target == observation.target,
                  reports[row.id] == nil, !observation.sourceLabel.isEmpty,
                  observation.sourceLabel.utf8.count <= 1024 else { throw invalid("unknown/duplicate candidate or target/source mismatch") }
            guard observation.interfacePlan.contactDistanceNM == profile.contactDistanceNM,
                  observation.interfacePlan.shortDistanceNM == profile.shortDistanceNM else {
                throw invalid("all candidates require the same interface-distance profile")
            }
            totalAtoms += observation.structure.atoms.count
            guard totalAtoms <= 1_000_000 else { throw invalid("structural archive atom capacity exceeded") }
            guard observation.interfacePlan.binderChains.count == 1 else {
                throw invalid("published candidate sequence requires exactly one binder chain")
            }
            let report = try VivoMolecularInterface.analyze(observation.structure, plan: observation.interfacePlan)
            totalWork += report.evaluatedAtomPairs
            guard totalWork <= 100_000_000 else { throw invalid("structural archive pair-evaluation budget exceeded") }
            let sequence = try proteinSequence(observation.structure, chain: observation.interfacePlan.binderChains[0])
            guard let expected = row.sourceFields["sequence"], sequence == expected else {
                throw invalid("binder structure sequence does not match published candidate: \(row.id)")
            }
            // Check the selected target's representation too; its identity remains caller-declared.
            let sequences = try observation.interfacePlan.targetChains.map {
                try proteinSequence(observation.structure, chain: $0)
            }.sorted()
            if let previous = targetSequences[row.target], previous != sequences {
                throw invalid("inconsistent target sequences within one target cohort")
            }
            targetSequences[row.target] = sequences
            reports[row.id] = report; labels[row.id] = observation.sourceLabel
        }
        let augmented = try dataset.records.map { row -> VivoBinderBenchmark.Record in
            guard Set(row.features.keys).isDisjoint(with: featureNames) else { throw invalid("structural feature names already present") }
            var features = row.features
            if let report = reports[row.id] {
                features["numi.interface.contactAtomPairs"] = Double(report.contactAtomPairs)
                features["numi.interface.contactResiduePairs"] = Double(report.residueContacts.count)
                features["numi.interface.binderContactFraction"] = Double(report.binderContactAtomIndices.count) / Double(report.binderHeavyAtomIndices.count)
                features["numi.interface.targetContactFraction"] = Double(report.targetContactAtomIndices.count) / Double(report.targetHeavyAtomIndices.count)
                features["numi.interface.shortDistanceFraction"] = Double(report.shortDistanceAtomPairs) / Double(report.evaluatedAtomPairs)
                features["numi.interface.minimumDistanceNM"] = report.minimumDistanceNM
                features["numi.interface.centroidDistanceNM"] = report.geometricCentroidDistanceNM
                features["numi.interface.binderGeometricRgNM"] = report.binderGeometricRadiusOfGyrationNM
            }
            return .init(id: row.id, target: row.target, leakageGroup: row.leakageGroup,
                outcome: row.outcome, rawOutcome: row.rawOutcome, features: features, sourceFields: row.sourceFields)
        }
        return Augmentation(schemaVersion: 1,
            dataset: .init(sourceSHA256: dataset.sourceSHA256, assay: dataset.assay,
                           groupingMethod: dataset.groupingMethod, records: augmented),
            interfaces: reports, sourceLabels: labels,
            unavailableCandidateIDs: dataset.records.filter { reports[$0.id] == nil }.map(\.id).sorted(),
            limitations: ["Static geometry, not molecular preparation, dynamics or binding evidence.",
                "Binder sequence is checked; target sequence and model origin remain caller-declared.",
                "Canonical protein heavy-atom names are required; no atoms, charges or protonation states are repaired.",
                "Missing structures remain missing features, never zero observations."])
    }

    /// Same candidate population, same training-only scaling, same ridge penalty.
    /// Keeping excluded rows with empty features preserves the original test-group purge.
    public static func matchedScoreOnly(_ augmented: VivoBinderBenchmark.Dataset,
                                         plan: VivoBinderBenchmark.Plan) throws -> VivoBinderBenchmark.Report {
        let allowed = VivoBinderAnthropicImport.allowedFeatures
        let scoreFeatures = plan.modelFeatures.filter { allowed.contains($0) }
        guard !scoreFeatures.isEmpty, plan.modelFeatures.contains(where: { featureNames.contains($0) }),
              Set(plan.modelFeatures).isSubset(of: allowed.union(featureNames)),
              allowed.contains(plan.baselineFeature) else { throw invalid("structural plan must combine declared in-silico and geometry features") }
        let required = Set(plan.modelFeatures + [plan.baselineFeature])
        let rows = augmented.records.map { row in
            VivoBinderBenchmark.Record(id: row.id, target: row.target, leakageGroup: row.leakageGroup,
                outcome: row.outcome, rawOutcome: row.rawOutcome,
                features: required.isSubset(of: Set(row.features.keys)) ? row.features : [:], sourceFields: row.sourceFields)
        }
        return try VivoBinderBenchmark.evaluate(.init(sourceSHA256: augmented.sourceSHA256, assay: augmented.assay,
            groupingMethod: augmented.groupingMethod, records: rows),
            plan: .init(sourceSHA256: plan.sourceSHA256, trainingTargets: plan.trainingTargets,
                testTargets: plan.testTargets, baselineFeature: plan.baselineFeature,
                modelFeatures: scoreFeatures, topK: plan.topK, ridgePenalty: plan.ridgePenalty))
    }

    static func proteinSequence(_ structure: VivoMolecularStructure, chain: String) throws -> String {
        guard let chain = structure.chains.first(where: { $0.identifier == chain }), !chain.residueIndices.isEmpty else {
            throw invalid("empty/missing protein chain")
        }
        var sequence = ""
        for index in chain.residueIndices {
            let residue = structure.residues[Int(index)]
            guard let entry = standardResidues[residue.name] else { throw invalid("noncanonical protein residue: \(residue.name)") }
            let atoms = residue.atomIndices.map { structure.atoms[Int($0)] }.filter { $0.element.atomicNumber != 1 }
            let names = Set(atoms.map(\.name)), expected = Set(entry.1.split(separator: " ").map(String.init))
            guard expected.isSubset(of: names), names.isSubset(of: expected.union(["OXT"])) else {
                throw invalid("missing or unsupported heavy atoms in \(residue.name) residue \(index)")
            }
            for atom in atoms {
                let element: UInt16 = atom.name.hasPrefix("N") ? 7 : atom.name.hasPrefix("O") ? 8 : atom.name.hasPrefix("S") ? 16 : 6
                guard atom.element.atomicNumber == element else { throw invalid("protein atom name/element mismatch") }
            }
            sequence += entry.0
        }
        return sequence
    }
    private static let standardResidues: [String: (String, String)] = [
        "ALA": ("A", "N CA C O CB"), "ARG": ("R", "N CA C O CB CG CD NE CZ NH1 NH2"),
        "ASN": ("N", "N CA C O CB CG OD1 ND2"), "ASP": ("D", "N CA C O CB CG OD1 OD2"),
        "CYS": ("C", "N CA C O CB SG"), "GLN": ("Q", "N CA C O CB CG CD OE1 NE2"),
        "GLU": ("E", "N CA C O CB CG CD OE1 OE2"), "GLY": ("G", "N CA C O"),
        "HIS": ("H", "N CA C O CB CG ND1 CD2 CE1 NE2"), "ILE": ("I", "N CA C O CB CG1 CG2 CD1"),
        "LEU": ("L", "N CA C O CB CG CD1 CD2"), "LYS": ("K", "N CA C O CB CG CD CE NZ"),
        "MET": ("M", "N CA C O CB CG SD CE"), "PHE": ("F", "N CA C O CB CG CD1 CD2 CE1 CE2 CZ"),
        "PRO": ("P", "N CA C O CB CG CD"), "SER": ("S", "N CA C O CB OG"),
        "THR": ("T", "N CA C O CB OG1 CG2"), "TRP": ("W", "N CA C O CB CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2"),
        "TYR": ("Y", "N CA C O CB CG CD1 CD2 CE1 CE2 CZ OH"), "VAL": ("V", "N CA C O CB CG1 CG2")
    ]
    private static func invalid(_ reason: String) -> VivoBinderBenchmark.Failure { .invalid(reason) }
}
