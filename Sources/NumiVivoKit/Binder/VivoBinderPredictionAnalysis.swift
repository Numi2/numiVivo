import Foundation

/// One outcome-blind prediction at a time. Source parsing, chain/sequence checks,
/// interface geometry and SASA use the existing molecular owners.
public enum VivoBinderPredictionAnalysis {
    public struct Identity: Codable, Sendable, Hashable {
        public let candidateID: String
        public let predictor: String
        /// Retain a source sample/seed identifier; do not invent an integer seed.
        public let sampleID: String
        public let targetConstructID: String
        public init(candidateID: String, predictor: String, sampleID: String, targetConstructID: String) {
            self.candidateID = candidateID; self.predictor = predictor
            self.sampleID = sampleID; self.targetConstructID = targetConstructID
        }
    }
    public struct Input: Codable, Sendable {
        public let schemaVersion: Int
        public let identity: Identity
        public let binderSequence: String
        public let source: VivoBinderStructureSources.Source
        public let surfacePlan: VivoMolecularInterface.SurfacePlan
        public init(identity: Identity, binderSequence: String, source: VivoBinderStructureSources.Source,
                    surfacePlan: VivoMolecularInterface.SurfacePlan) {
            schemaVersion = 1; self.identity = identity; self.binderSequence = binderSequence
            self.source = source; self.surfacePlan = surfacePlan
        }
    }
    public struct Contact: Codable, Sendable, Hashable {
        /// Zero-based positions in the checked complete single-chain sequences.
        public let binderOffset: Int
        public let targetOffset: Int
    }
    public struct Report: Codable, Sendable {
        public let schemaVersion: Int
        public let identity: Identity
        public let inputSHA256: String
        public let sourceSHA256: String
        public let binderSequenceSHA256: String
        public let targetSequenceSHA256: String
        public let target: String
        public let interface: VivoMolecularInterface.Report
        public let surface: VivoMolecularInterface.SurfaceReport
        public let contacts: [Contact]
        public let features: [String: Double]
        public let limitations: [String]
    }

    public static func decodeInput(_ data: Data) throws -> Input {
        guard data.count <= 20 * 1024 * 1024 else { throw invalid("prediction input byte capacity") }
        func object(_ value: Any?, keys: Set<String>) throws -> [String: Any] {
            guard let d = value as? [String: Any], Set(d.keys) == keys else {
                throw invalid("prediction input contains unknown or outcome-bearing fields")
            }
            return d
        }
        let root = try object(JSONSerialization.jsonObject(with: data),
            keys: ["schemaVersion", "identity", "binderSequence", "source", "surfacePlan"])
        _ = try object(root["identity"], keys: ["candidateID", "predictor", "sampleID", "targetConstructID"])
        let source = try object(root["source"], keys: ["candidateID", "target", "sourceLabel", "format", "contents",
            "sha256", "targetChainSequences", "interfacePlan"])
        _ = try object(source["interfacePlan"], keys: ["schemaVersion", "binderChains", "targetChains", "conformerID",
            "contactDistanceNM", "shortDistanceNM", "maximumPairEvaluations"])
        _ = try object(root["surfacePlan"], keys: ["schemaVersion", "radiusProfile", "radiiNM", "probeRadiusNM",
            "pointsPerAtom", "overlapThresholdNM", "maximumNeighborPairs", "maximumPointTests", "maximumNeighborSearchTests"])
        return try VivoCanonicalJSON.decode(Input.self, from: data)
    }

    public static func analyze(_ input: Input) throws -> Report {
        let id = input.identity, source = input.source
        guard input.schemaVersion == 1, id.candidateID == source.candidateID,
              [id.candidateID, id.predictor, id.sampleID, id.targetConstructID].allSatisfy({ !$0.isEmpty && $0.utf8.count <= 512 }),
              !input.binderSequence.isEmpty, input.binderSequence.utf8.count <= 5000,
              input.binderSequence.allSatisfy({ "ACDEFGHIKLMNPQRSTVWY".contains($0) }),
              source.interfacePlan.binderChains.count == 1, source.interfacePlan.targetChains.count == 1 else {
            throw invalid("prediction requires bounded identity, canonical sequence and explicit 1:1 chain selection")
        }
        let parsed = try VivoBinderStructureSources.reconstruct(.init(sources: [source])).observations[0]
        let structure = parsed.structure
        let binderChain = source.interfacePlan.binderChains[0], targetChain = source.interfacePlan.targetChains[0]
        guard try VivoBinderStructuralFeatures.proteinSequence(structure, chain: binderChain) == input.binderSequence else {
            throw invalid("prediction binder sequence differs from declared candidate")
        }
        let targetSequence = try VivoBinderStructuralFeatures.proteinSequence(structure, chain: targetChain)
        let interface = try VivoMolecularInterface.analyze(structure, plan: source.interfacePlan)
        let surface = try VivoMolecularInterface.analyzeSurface(structure, interfacePlan: source.interfacePlan, plan: input.surfacePlan)
        func offsets(_ chain: String) -> [UInt32: Int] {
            Dictionary(uniqueKeysWithValues: structure.chains.first { $0.identifier == chain }!
                .residueIndices.enumerated().map { ($0.element, $0.offset) })
        }
        let b = offsets(binderChain), t = offsets(targetChain)
        var contacts: [Contact] = interface.residueContacts.map { row in
            Contact(binderOffset: b[row.binderResidueIndex]!, targetOffset: t[row.targetResidueIndex]!)
        }
        contacts.sort { a, b in
            if a.binderOffset == b.binderOffset { return a.targetOffset < b.targetOffset }
            return a.binderOffset < b.binderOffset
        }
        var features: [String: Double] = [
            "numi.surface.buriedAreaSumNM2": surface.buriedAreaSumNM2,
            "numi.surface.halfBuriedAreaSumNM2": surface.halfBuriedAreaSumNM2,
            "numi.surface.overlapPairs": Double(surface.crossPartnerOverlapPairs),
            "numi.surface.maximumOverlapNM": surface.maximumCrossPartnerOverlapNM,
            "numi.interface.contactResiduePairs": Double(contacts.count)]
        if surface.binderIsolatedAreaNM2 > 0 {
            features["numi.surface.binderBuriedFraction"] = surface.binderBuriedAreaNM2 / surface.binderIsolatedAreaNM2
        }
        if surface.buriedAreaSumNM2 > 0 {
            features["numi.surface.contactResiduePairsPerNM2"] = Double(contacts.count) / surface.buriedAreaSumNM2
        }
        return Report(schemaVersion: 1, identity: id, inputSHA256: try hash(input), sourceSHA256: source.sha256,
            binderSequenceSHA256: try VivoCanonicalJSON.fingerprint(Data(input.binderSequence.utf8)).hex,
            targetSequenceSHA256: try VivoCanonicalJSON.fingerprint(Data(targetSequence.utf8)).hex,
            target: source.target, interface: interface, surface: surface, contacts: contacts, features: features,
            limitations: ["Outcome-blind structural observations, not a fitted ranking or affinity estimate.",
                "Declared model/sample/construct provenance is retained, not authenticated model execution.",
                "Complete canonical 1:1 protein chains only; no cropped-sequence alignment or automatic repair.",
                "Undefined ratios are omitted, not replaced with zero; absolute areas and counts remain available."])
    }

    public struct Distribution: Codable, Sendable {
        public let availableCount: Int
        public let minimum: Double
        public let maximum: Double
        public let mean: Double
        public let populationSD: Double
    }
    public struct Series: Codable, Sendable {
        public let candidateID: String
        public let targetConstructID: String
        public let predictor: String
        public let sampleIDs: [String]
        public let sourceSHA256: [String]
        public let reportSHA256: [String]
        public let features: [String: Distribution]
        public let meanContactJaccard: Double?
        public let contactComparisons: Int
        public let bothEmptyContactComparisons: Int
        public let limitations: [String]
    }
    /// Retains compact observations only; each append computes from a single original
    /// source. Failed append leaves state unchanged. No assay labels can enter.
    public struct Accumulator: Sendable {
        private struct Sample: Sendable {
            let reportSHA: String; let sourceSHA: String; let identity: Identity
            let binderSHA: String; let targetSHA: String; let profileSHA: String
            let target: String; let features: [String: Double]; let contacts: Set<Contact>
        }
        private var samples: [Identity: Sample] = [:]
        private var totalContacts = 0
        private var candidateSignatures: [String: String] = [:]
        private var constructSignatures: [String: String] = [:]
        public init() {}
        @discardableResult
        public mutating func append(_ input: Input) throws -> Report {
            guard samples.count < 10_000, samples[input.identity] == nil else { throw invalid("duplicate prediction identity or series capacity") }
            let report = try analyze(input)
            let candidateSignature = try hash([report.binderSequenceSHA256, report.target])
            let constructSignature = try hash([report.targetSequenceSHA256, report.target])
            if let previous = candidateSignatures[input.identity.candidateID], previous != candidateSignature {
                throw invalid("candidate sequence or target identity changed across predictions")
            }
            if let previous = constructSignatures[input.identity.targetConstructID], previous != constructSignature {
                throw invalid("target construct sequence or identity changed across predictions")
            }
            guard totalContacts + report.contacts.count <= 1_000_000 else { throw invalid("series contact capacity") }
            // Chain names may differ across predictors. Sequence offsets are canonical;
            // only one target and binder chain are admitted in this version.
            let profileSHA = try hash(Profile(surface: input.surfacePlan,
                contact: input.source.interfacePlan.contactDistanceNM, short: input.source.interfacePlan.shortDistanceNM))
            samples[input.identity] = Sample(reportSHA: try hash(report), sourceSHA: report.sourceSHA256,
                identity: input.identity, binderSHA: report.binderSequenceSHA256, targetSHA: report.targetSequenceSHA256,
                profileSHA: profileSHA, target: report.target, features: report.features, contacts: Set(report.contacts))
            totalContacts += report.contacts.count
            candidateSignatures[input.identity.candidateID] = candidateSignature
            constructSignatures[input.identity.targetConstructID] = constructSignature
            return report
        }
        public func finish() throws -> [Series] {
            let groups = Dictionary(grouping: samples.values, by: { Group(candidate: $0.identity.candidateID,
                construct: $0.identity.targetConstructID, predictor: $0.identity.predictor) })
            var output: [Series] = []
            var contactWork = 0
            for key in groups.keys.sorted(by: { ($0.candidate, $0.construct, $0.predictor) < ($1.candidate, $1.construct, $1.predictor) }) {
                let rows = groups[key]!.sorted { $0.identity.sampleID < $1.identity.sampleID }, first = rows[0]
                guard rows.count <= 64, rows.allSatisfy({ $0.binderSHA == first.binderSHA && $0.targetSHA == first.targetSHA
                    && $0.profileSHA == first.profileSHA && $0.target == first.target }) else { throw invalid("inconsistent sequence/profile or too many samples in series") }
                var distributions: [String: Distribution] = [:]
                for feature in Set(rows.flatMap { $0.features.keys }).sorted() {
                    let values = rows.compactMap { $0.features[feature] }, anchor = values[0]
                    let mean = anchor + values.reduce(0) { $0 + ($1 - anchor) } / Double(values.count)
                    let sd = sqrt(values.reduce(0) { $0 + ($1 - mean) * ($1 - mean) } / Double(values.count))
                    distributions[feature] = .init(availableCount: values.count, minimum: values.min()!, maximum: values.max()!, mean: mean, populationSD: sd)
                }
                var agreement = 0.0, comparisons = 0, bothEmpty = 0
                for i in rows.indices { for j in rows.indices where j > i {
                    contactWork += rows[i].contacts.count + rows[j].contacts.count
                    guard contactWork <= 10_000_000 else { throw invalid("series contact-comparison budget exceeded") }
                    let union = rows[i].contacts.union(rows[j].contacts).count
                    if union == 0 { bothEmpty += 1; continue }
                    agreement += Double(rows[i].contacts.intersection(rows[j].contacts).count) / Double(union); comparisons += 1
                } }
                output.append(.init(candidateID: key.candidate, targetConstructID: key.construct, predictor: key.predictor,
                    sampleIDs: rows.map { $0.identity.sampleID }, sourceSHA256: rows.map(\.sourceSHA), reportSHA256: rows.map(\.reportSHA),
                    features: distributions, meanContactJaccard: comparisons == 0 ? nil : agreement / Double(comparisons),
                    contactComparisons: comparisons, bothEmptyContactComparisons: bothEmpty,
                    limitations: ["Sample dispersion is descriptive, not confidence intervals or independent experimental evidence.",
                        "Both-empty contact sets are unavailable agreement, not perfect interface agreement.",
                        "Predictors and target constructs are summarized separately; no consensus filter is applied."]))
            }
            return output
        }
    }
    private struct Group: Hashable { let candidate: String; let construct: String; let predictor: String }
    private struct Profile: Codable { let surface: VivoMolecularInterface.SurfacePlan; let contact: Double; let short: Double }
    private static func hash<T: Encodable>(_ value: T) throws -> String { try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(value)).hex }
    private static func invalid(_ reason: String) -> VivoBinderBenchmark.Failure { .invalid(reason) }
}
