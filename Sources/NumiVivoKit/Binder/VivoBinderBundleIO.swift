import Foundation

/// Immutable source-reconstructable bundles using the suite's canonical JSON/SHA-256 owner.
public enum VivoBinderBundleIO {
    public struct Receipt: Codable, Sendable {
        public let schemaVersion: Int
        public let kind: String
        public let implementationSHA256: String
        public let files: [String: String]
    }
    private static let inputNames: Set<String> = ["source.csv", "import.json", "imported.json"]
    private static let evaluationNames: Set<String> = inputNames.union(["plan.json", "report.json"])
    private static let structuralNames: Set<String> = evaluationNames.union(["structures.json", "geometry.json", "score-only-report.json"])

    private static let sourceStructuralNames: Set<String> = structuralNames.union(["structure-sources.json"])

    private static let rankingQueryNames: Set<String> = ["query.json"]
    private static let rankingNames: Set<String> = ["query.json", "ranking-plan.json", "ranking.json"]
    private static let rankingAssessmentNames: Set<String> = inputNames.union(rankingNames).union(["assessment.json"])

    /// Publishes no raw CSV, assay, outcomes or source fields to the ranking process.
    public static func rankingQuery(bundle: URL, to output: URL, implementationSHA256: String) throws -> Receipt {
        let receipt = try verify(bundle, implementationSHA256: implementationSHA256)
        guard receipt.kind == "import" else { throw invalid("ranking query requires an import bundle") }
        let imported = try VivoCanonicalJSON.decode(VivoBinderAnthropicImport.Result.self,
            from: read(bundle.appendingPathComponent("imported.json")))
        let query = try VivoBinderRanking.query(from: imported.dataset)
        return try publish(["query.json": VivoCanonicalJSON.encode(query)], kind: "rankingQuery",
                           to: output, implementationSHA256: implementationSHA256)
    }

    public static func rank(bundle: URL, plan: URL, to output: URL, implementationSHA256: String) throws -> Receipt {
        let receipt = try verify(bundle, implementationSHA256: implementationSHA256)
        guard receipt.kind == "rankingQuery" else { throw invalid("rank requires an outcome-blind query bundle") }
        let queryData = try read(bundle.appendingPathComponent("query.json")), planData = try read(plan)
        let query = try VivoBinderRanking.decodeQuery(queryData), request = try VivoBinderRanking.decodePlan(planData)
        let result = try VivoBinderRanking.rank(query, plan: request)
        return try publish(["query.json": queryData, "ranking-plan.json": planData,
                            "ranking.json": VivoCanonicalJSON.encode(result)], kind: "ranking",
                           to: output, implementationSHA256: implementationSHA256)
    }

    public static func assessRanking(bundle: URL, ranking: URL, to output: URL,
                                     implementationSHA256: String) throws -> Receipt {
        let input = try verify(bundle, implementationSHA256: implementationSHA256)
        let ranked = try verify(ranking, implementationSHA256: implementationSHA256)
        guard input.kind == "import", ranked.kind == "ranking" else {
            throw invalid("assessment requires an import bundle and a completed ranking bundle")
        }
        var files: [String: Data] = [:]
        for name in inputNames { files[name] = try read(bundle.appendingPathComponent(name)) }
        for name in rankingNames { files[name] = try read(ranking.appendingPathComponent(name)) }
        let imported = try VivoCanonicalJSON.decode(VivoBinderAnthropicImport.Result.self, from: files["imported.json"]!)
        files["assessment.json"] = try rankingAssessment(imported.dataset, files: files)
        return try publish(files, kind: "rankingAssessment", to: output, implementationSHA256: implementationSHA256)
    }

    public static func importCSV(source: URL, configuration: URL, to output: URL,
                                 implementationSHA256: String) throws -> Receipt {
        try checkDigest(implementationSHA256)
        let raw = try read(source), configData = try read(configuration)
        let config = try VivoCanonicalJSON.decode(VivoBinderAnthropicImport.Configuration.self, from: configData)
        let imported = try VivoBinderAnthropicImport.parse(raw, sourceSHA256: digest(raw), configuration: config)
        return try publish(["source.csv": raw, "import.json": configData,
                            "imported.json": VivoCanonicalJSON.encode(imported)],
                           kind: "import", to: output, implementationSHA256: implementationSHA256)
    }

    public static func evaluate(bundle: URL, plan: URL, to output: URL,
                                implementationSHA256: String) throws -> Receipt {
        let receipt = try verify(bundle, implementationSHA256: implementationSHA256)
        guard receipt.kind == "import" else { throw invalid("evaluation requires an import bundle") }
        var files: [String: Data] = [:]
        for name in inputNames { files[name] = try read(bundle.appendingPathComponent(name)) }
        let imported = try VivoCanonicalJSON.decode(VivoBinderAnthropicImport.Result.self, from: files["imported.json"]!)
        let planData = try read(plan), request = try VivoCanonicalJSON.decode(VivoBinderBenchmark.Plan.self, from: planData)
        let report = try VivoBinderBenchmark.evaluate(imported.dataset, plan: request)
        files["plan.json"] = planData; files["report.json"] = try VivoCanonicalJSON.encode(report)
        return try publish(files, kind: "evaluation", to: output, implementationSHA256: implementationSHA256)
    }

    public static func evaluateStructures(bundle: URL, plan: URL, structures: URL, to output: URL,
                                          implementationSHA256: String) throws -> Receipt {
        let receipt = try verify(bundle, implementationSHA256: implementationSHA256)
        guard receipt.kind == "import" else { throw invalid("structural evaluation requires an import bundle") }
        var files: [String: Data] = [:]
        for name in inputNames { files[name] = try read(bundle.appendingPathComponent(name)) }
        let imported = try VivoCanonicalJSON.decode(VivoBinderAnthropicImport.Result.self, from: files["imported.json"]!)
        files["plan.json"] = try read(plan); files["structures.json"] = try read(structures)
        let computed = try structuralReports(imported.dataset, plan: files["plan.json"]!, structures: files["structures.json"]!)
        files.merge(computed) { _, new in new }
        return try publish(files, kind: "structuralEvaluation", to: output, implementationSHA256: implementationSHA256)
    }

    /// Parses original retained PDB/mmCIF bytes through the shared molecular owners.
    /// Raw source text, chain references, parsed structures and both matched fits replay.
    public static func evaluateStructureSources(bundle: URL, plan: URL, sources: URL, to output: URL,
                                                implementationSHA256: String) throws -> Receipt {
        let receipt = try verify(bundle, implementationSHA256: implementationSHA256)
        guard receipt.kind == "import" else { throw invalid("source structural evaluation requires an import bundle") }
        var files: [String: Data] = [:]
        for name in inputNames { files[name] = try read(bundle.appendingPathComponent(name)) }
        let imported = try VivoCanonicalJSON.decode(VivoBinderAnthropicImport.Result.self, from: files["imported.json"]!)
        files["plan.json"] = try read(plan)
        files["structure-sources.json"] = try read(sources)
        let sourceInput = try VivoCanonicalJSON.decode(VivoBinderStructureSources.Input.self,
                                                       from: files["structure-sources.json"]!)
        files["structures.json"] = try VivoCanonicalJSON.encode(VivoBinderStructureSources.reconstruct(sourceInput))
        files.merge(try structuralReports(imported.dataset, plan: files["plan.json"]!,
                                          structures: files["structures.json"]!)) { _, new in new }
        return try publish(files, kind: "sourceStructuralEvaluation", to: output,
                           implementationSHA256: implementationSHA256)
    }

    /// Reconstructs imported records from original bytes; evaluation bundles also recompute
    /// fitting and every held-out score. Verification requires the original executable identity.
    @discardableResult
    public static func verify(_ bundle: URL, implementationSHA256: String) throws -> Receipt {
        try checkDigest(implementationSHA256)
        let r = try VivoCanonicalJSON.decode(Receipt.self, from: read(bundle.appendingPathComponent("receipt.json")))
        guard r.schemaVersion == 1, ["import", "evaluation", "structuralEvaluation", "sourceStructuralEvaluation", "rankingQuery", "ranking", "rankingAssessment"].contains(r.kind),
              r.implementationSHA256 == implementationSHA256 else { throw invalid("receipt schema/kind/executable mismatch") }
        let fileSets = ["import": inputNames, "evaluation": evaluationNames, "structuralEvaluation": structuralNames,
                        "sourceStructuralEvaluation": sourceStructuralNames, "rankingQuery": rankingQueryNames, "ranking": rankingNames, "rankingAssessment": rankingAssessmentNames]
        let expected = fileSets[r.kind]!
        let names = try FileManager.default.contentsOfDirectory(atPath: bundle.path)
        guard Set(r.files.keys) == expected, Set(names) == expected.union(["receipt.json"]) else {
            throw invalid("unexpected/missing bundle files")
        }
        var files: [String: Data] = [:]
        for name in expected.sorted() {
            let data = try read(bundle.appendingPathComponent(name))
            guard try digest(data) == r.files[name] else { throw invalid("source/artifact hash mismatch: \(name)") }
            files[name] = data
        }
        if r.kind == "rankingQuery" || r.kind == "ranking" {
            let query = try VivoBinderRanking.decodeQuery(files["query.json"]!)
            if r.kind == "ranking" {
                let plan = try VivoBinderRanking.decodePlan(files["ranking-plan.json"]!)
                let result = try VivoBinderRanking.rank(query, plan: plan)
                guard try VivoCanonicalJSON.encode(result) == files["ranking.json"]! else {
                    throw invalid("outcome-blind ranking does not replay")
                }
            }
            return r
        }
        let config = try VivoCanonicalJSON.decode(VivoBinderAnthropicImport.Configuration.self, from: files["import.json"]!)
        let reconstructed = try VivoBinderAnthropicImport.parse(files["source.csv"]!,
            sourceSHA256: digest(files["source.csv"]!), configuration: config)
        guard try VivoCanonicalJSON.encode(reconstructed) == files["imported.json"]! else {
            throw invalid("imported records do not reconstruct from source")
        }
        if r.kind == "rankingAssessment" {
            guard try rankingAssessment(reconstructed.dataset, files: files) == files["assessment.json"]! else {
                throw invalid("ranking assessment does not reconstruct")
            }
        } else if r.kind == "evaluation" {
            let plan = try VivoCanonicalJSON.decode(VivoBinderBenchmark.Plan.self, from: files["plan.json"]!)
            let recomputed = try VivoBinderBenchmark.evaluate(reconstructed.dataset, plan: plan)
            guard try VivoCanonicalJSON.encode(recomputed) == files["report.json"]! else {
                throw invalid("evaluation report does not replay")
            }
        } else if r.kind == "structuralEvaluation" || r.kind == "sourceStructuralEvaluation" {
            if r.kind == "sourceStructuralEvaluation" {
                let sourceInput = try VivoCanonicalJSON.decode(VivoBinderStructureSources.Input.self,
                                                               from: files["structure-sources.json"]!)
                let parsed = try VivoBinderStructureSources.reconstruct(sourceInput)
                guard try VivoCanonicalJSON.encode(parsed) == files["structures.json"]! else {
                    throw invalid("structures do not reconstruct from original PDB/mmCIF sources")
                }
            }
            let recomputed = try structuralReports(reconstructed.dataset, plan: files["plan.json"]!, structures: files["structures.json"]!)
            for (name, data) in recomputed where data != files[name] {
                throw invalid("structural evaluation does not reconstruct: \(name)")
            }
        }
        return r
    }

    private static func rankingAssessment(_ dataset: VivoBinderBenchmark.Dataset, files: [String: Data]) throws -> Data {
        let query = try VivoBinderRanking.decodeQuery(files["query.json"]!)
        let plan = try VivoBinderRanking.decodePlan(files["ranking-plan.json"]!)
        let result = try VivoCanonicalJSON.decode(VivoBinderRanking.Report.self, from: files["ranking.json"]!)
        return try VivoCanonicalJSON.encode(VivoBinderRanking.assess(dataset, query: query, plan: plan, ranking: result))
    }

    private static func structuralReports(_ dataset: VivoBinderBenchmark.Dataset, plan: Data,
                                           structures: Data) throws -> [String: Data] {
        let request = try VivoCanonicalJSON.decode(VivoBinderBenchmark.Plan.self, from: plan)
        let input = try VivoCanonicalJSON.decode(VivoBinderStructuralFeatures.Input.self, from: structures)
        let geometry = try VivoBinderStructuralFeatures.augment(dataset, input: input)
        let control = try VivoBinderStructuralFeatures.matchedScoreOnly(geometry.dataset, plan: request)
        let enhanced = try VivoBinderBenchmark.evaluate(geometry.dataset, plan: request)
        guard control.trainingIDs == enhanced.trainingIDs,
              control.targets.map(\.target) == enhanced.targets.map(\.target),
              control.targets.map(\.candidateIDs) == enhanced.targets.map(\.candidateIDs) else {
            throw invalid("score-only and structural models have different populations")
        }
        return ["geometry.json": try VivoCanonicalJSON.encode(geometry),
                "score-only-report.json": try VivoCanonicalJSON.encode(control),
                "report.json": try VivoCanonicalJSON.encode(enhanced)]
    }

    private static func publish(_ files: [String: Data], kind: String, to output: URL,
                                implementationSHA256: String) throws -> Receipt {
        let fm = FileManager.default, destination = output.standardizedFileURL
        guard !fm.fileExists(atPath: destination.path) else { throw invalid("output already exists") }
        let parent = destination.deletingLastPathComponent()
        var isDirectory: ObjCBool = false
        guard fm.fileExists(atPath: parent.path, isDirectory: &isDirectory), isDirectory.boolValue else {
            throw invalid("output parent directory must exist")
        }
        let temp = parent.appendingPathComponent(".numivivo-binder-\(UUID().uuidString)", isDirectory: true)
        try fm.createDirectory(at: temp, withIntermediateDirectories: false)
        defer { try? fm.removeItem(at: temp) }
        var hashes: [String: String] = [:]
        for name in files.keys.sorted() {
            let data = files[name]!
            hashes[name] = try digest(data)
            try data.write(to: temp.appendingPathComponent(name), options: .withoutOverwriting)
        }
        let receipt = Receipt(schemaVersion: 1, kind: kind, implementationSHA256: implementationSHA256, files: hashes)
        try VivoCanonicalJSON.encode(receipt).write(to: temp.appendingPathComponent("receipt.json"), options: .withoutOverwriting)
        // Validate before the atomic same-parent rename. A failed run never publishes a bundle.
        try verify(temp, implementationSHA256: implementationSHA256)
        try fm.moveItem(at: temp, to: destination)
        return receipt
    }
    private static func read(_ url: URL) throws -> Data {
        let info = try url.resourceValues(forKeys: [.isRegularFileKey, .isSymbolicLinkKey, .fileSizeKey])
        guard info.isRegularFile == true, info.isSymbolicLink != true,
              let size = info.fileSize, size <= 128 * 1024 * 1024 else { throw invalid("unsupported file or size: \(url.lastPathComponent)") }
        let data = try Data(contentsOf: url)
        guard data.count == size else { throw invalid("file changed while reading") }
        return data
    }
    private static func digest(_ data: Data) throws -> String { try VivoCanonicalJSON.fingerprint(data).hex }
    private static func checkDigest(_ text: String) throws {
        guard text.utf8.count == 64, text.utf8.allSatisfy({ (48...57).contains($0) || (97...102).contains($0) }) else {
            throw invalid("implementation SHA-256 required")
        }
    }
    private static func invalid(_ reason: String) -> VivoBinderBenchmark.Failure { .invalid(reason) }
}
