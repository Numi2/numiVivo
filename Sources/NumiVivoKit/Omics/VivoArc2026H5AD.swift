import Foundation

/// Frozen, control-only inputs and raw-count submission artifacts. This owner
/// reuses NumiVivo's snapshot, HDF5 frame/count readers and canonical JSON; no
/// inference engine, expression transformation or competing count store is added.
public enum VivoArc2026H5AD {
    public static func readPlan<T: Decodable>(_ type: T.Type, at url: URL) throws -> T {
        try VivoCanonicalJSON.decode(type, from: VivoSingleCellCampaignIO.readDocument(url, maximumBytes: VivoArc2026.maximumJSONBytes))
    }
    private static func read(_ root: URL, _ path: String) throws -> Data {
        try VivoSingleCellCampaignIO.readDocument(root.appendingPathComponent(path), maximumBytes: VivoArc2026.maximumJSONBytes)
    }
    private static func hash(_ data: Data) throws -> String { try VivoCanonicalJSON.fingerprint(data).hex }
    private static func snapshot(_ source: URL, _ target: URL?, maximum: Int = VivoArc2026.maximumFileBytes) throws -> String {
        try VivoOmicsFileSnapshot.fingerprint(source, copyTo: target, maximumBytes: maximum).hex
    }
    private static func mkdir(_ url: URL) throws {
        try FileManager.default.createDirectory(at: url, withIntermediateDirectories: true, attributes: [.posixPermissions: 0o700])
    }
    private static func stage(_ parent: URL) throws -> URL {
        let url = parent.appendingPathComponent(".numivivo-arc2026-" + UUID().uuidString)
        try FileManager.default.createDirectory(at: url, withIntermediateDirectories: false, attributes: [.posixPermissions: 0o700])
        return url
    }
    private static func fresh(_ destination: URL) throws {
        guard destination.isFileURL, (try? FileManager.default.attributesOfItem(atPath: destination.path)) == nil else {
            throw VivoOmicsError.invalid("Arc output exists or is not local")
        }
    }
    private static func write<T: Encodable>(_ value: T, _ root: URL, _ path: String) throws {
        try VivoCanonicalJSON.encode(value).write(to: root.appendingPathComponent(path), options: .withoutOverwriting)
    }
    private static func inspect(_ source: URL, column: String, expectedGenes: [String]?,
                                targets: [String], controlsOnly: Bool) throws -> ([String], VivoArc2026CountInspection) {
        var limits = VivoOmicsLimits()
        limits.maximumCells = VivoArc2026.maximumCells; limits.maximumFeatures = VivoArc2026.featureCount
        limits.maximumNonzeros = VivoArc2026.maximumEntries; limits.maximumInputBytes = VivoArc2026.maximumFileBytes
        return try VivoHDF5.lock.withLock {
            let h = try VivoHDF5(), file = try h.file(source.path)
            defer { h.close(file, "H5Fclose") }
            guard try h.text(file, "encoding-type") == "anndata", try h.text(file, "encoding-version") == "0.1.0" else {
                throw VivoOmicsError.invalid("Arc input requires supported AnnData encoding")
            }
            let reader = VivoH5ADFrameReader(h: h, file: file)
            let genes = try reader.index("var", maximum: limits.maximumFeatures)
            try VivoArc2026.labels(genes, count: VivoArc2026.featureCount)
            if let expectedGenes, genes != expectedGenes { throw VivoOmicsError.invalid("Arc gene axis/order mismatch; no intersection or reordering is permitted") }
            let cells = try reader.index("obs", maximum: limits.maximumCells)
            guard Set(cells).count == cells.count else { throw VivoOmicsError.invalid("Arc duplicate observation identity") }
            let labels = try reader.required(reader.column("obs", column, maximum: limits.maximumCells))
            guard labels.count == cells.count else { throw VivoOmicsError.invalid("Arc label/cell shape mismatch") }
            let expected = controlsOnly ? Set([VivoArc2026.control]) : Set(targets + [VivoArc2026.control])
            var accumulator = try VivoArc2026CountAccumulator(labels: labels, featureCount: genes.count, expectedLabels: expected)
            try VivoH5ADCountReader.scan(h, file: file, path: "X", rows: cells.count, features: genes.count, limits: limits) {
                try accumulator.append(row: $0, feature: $1, count: $2)
            }
            return (genes, try accumulator.finish())
        }
    }

    public static func prepare(_ plan: VivoArc2026PreparePlan, to destination: URL) throws -> VivoArc2026Query {
        try plan.validate(); try fresh(destination)
        let temp = try stage(destination.deletingLastPathComponent())
        defer { try? FileManager.default.removeItem(at: temp) }
        var contexts: [VivoArc2026QueryContext] = []
        for c in plan.contexts {
            let dir = temp.appendingPathComponent("contexts/" + c.id); try mkdir(dir)
            let targetHash = try snapshot(URL(fileURLWithPath: c.targetsJSON), dir.appendingPathComponent("targets.json"), maximum: VivoArc2026.maximumJSONBytes)
            let targets = try readPlan([String].self, at: dir.appendingPathComponent("targets.json"))
            try VivoArc2026.labels(targets, count: VivoArc2026.targetCount)
            guard !targets.contains(VivoArc2026.control) else { throw VivoOmicsError.invalid("Arc targets include the control") }
            let controlFile = dir.appendingPathComponent("controls.h5ad")
            let controlsHash = try snapshot(URL(fileURLWithPath: c.controlsH5AD), controlFile)
            let (genes, result) = try inspect(controlFile, column: c.perturbationColumn,
                expectedGenes: contexts.first?.genes, targets: targets, controlsOnly: true)
            contexts.append(.init(id: c.id, perturbationColumn: c.perturbationColumn, genes: genes,
                targets: targets, controlsSHA256: controlsHash, targetsSHA256: targetHash, inspection: result))
        }
        let query = VivoArc2026Query(schemaVersion: 1, phase: plan.phase, evaluatorCommit: VivoArc2026.evaluatorCommit, contexts: contexts)
        try query.validate(); try write(query, temp, "query.json")
        try Task.checkCancellation(); try FileManager.default.moveItem(at: temp, to: destination)
        return query
    }

    // Verification always reads immutable copies, not a manifest followed by a
    // later independent read of changing source data. Copy only named payloads.
    private static func copyQuery(_ source: URL, into destination: URL) throws -> (VivoArc2026Query, String) {
        try mkdir(destination)
        let bytes = try read(source, "query.json")
        let query = try VivoCanonicalJSON.decode(VivoArc2026Query.self, from: bytes); try query.validate()
        try bytes.write(to: destination.appendingPathComponent("query.json"), options: .withoutOverwriting)
        for c in query.contexts {
            let path = "contexts/" + c.id, dir = destination.appendingPathComponent(path); try mkdir(dir)
            let th = try snapshot(source.appendingPathComponent(path + "/targets.json"), dir.appendingPathComponent("targets.json"), maximum: VivoArc2026.maximumJSONBytes)
            let ch = try snapshot(source.appendingPathComponent(path + "/controls.h5ad"), dir.appendingPathComponent("controls.h5ad"))
            guard th == c.targetsSHA256, ch == c.controlsSHA256,
                  try readPlan([String].self, at: dir.appendingPathComponent("targets.json")) == c.targets else {
                throw VivoOmicsError.invalid("Arc query payload fingerprint/target mismatch")
            }
            let (_, actual) = try inspect(dir.appendingPathComponent("controls.h5ad"), column: c.perturbationColumn,
                expectedGenes: c.genes, targets: c.targets, controlsOnly: true)
            guard actual == c.inspection else { throw VivoOmicsError.invalid("Arc query count reconstruction mismatch") }
        }
        return (query, try hash(bytes))
    }
    public static func verifyQuery(_ source: URL) throws -> VivoArc2026Query {
        let temp = try stage(FileManager.default.temporaryDirectory); defer { try? FileManager.default.removeItem(at: temp) }
        return try copyQuery(source, into: temp.appendingPathComponent("query")).0
    }

    public static func pack(query source: URL, plan: VivoArc2026PackPlan, to destination: URL) throws -> VivoArc2026Submission {
        guard plan.schemaVersion == 1, vivoOmicsID(plan.modelID), !plan.modelArtifact.isEmpty,
              !plan.trainingDataManifest.isEmpty else { throw VivoOmicsError.invalid("Arc prediction plan/model identity") }
        try VivoArc2026.contextIDs(plan.contexts.map(\.id)); try fresh(destination)
        let temp = try stage(destination.deletingLastPathComponent()); defer { try? FileManager.default.removeItem(at: temp) }
        let (query, queryHash) = try copyQuery(source, into: temp.appendingPathComponent("query"))
        let ids = Set(query.contexts.map(\.id))
        guard Set(plan.contexts.map(\.id)) == ids, Set(plan.excludedPerturbedContexts) == ids,
              plan.excludedPerturbedContexts.count == ids.count else { throw VivoOmicsError.invalid("Arc contexts or declared zero-shot exclusion differ") }
        // A declared model identity is not proof that this model emitted a file.
        // Retain that distinction; the competition's code review is independent.
        let modelHash = try snapshot(URL(fileURLWithPath: plan.modelArtifact), nil)
        let trainingHash = try snapshot(URL(fileURLWithPath: plan.trainingDataManifest), temp.appendingPathComponent("training-data.json"), maximum: VivoArc2026.maximumJSONBytes)
        guard try JSONSerialization.jsonObject(with: read(temp, "training-data.json")) is [String: Any] else {
            throw VivoOmicsError.invalid("Arc training-data provenance must be a JSON object")
        }
        var receipts: [VivoArc2026PredictionReceipt] = []
        for c in query.contexts {
            let input = plan.contexts.first { $0.id == c.id }!
            guard !input.predictionH5AD.isEmpty else { throw VivoOmicsError.invalid("Arc prediction path missing") }
            let dir = temp.appendingPathComponent("contexts/" + c.id); try mkdir(dir)
            let predicted = dir.appendingPathComponent("prediction.h5ad")
            let digest = try snapshot(URL(fileURLWithPath: input.predictionH5AD), predicted)
            let (_, inspection) = try inspect(predicted, column: c.perturbationColumn,
                expectedGenes: c.genes, targets: c.targets, controlsOnly: false)
            receipts.append(.init(id: c.id, predictionSHA256: digest, inspection: inspection))
        }
        let result = VivoArc2026Submission(schemaVersion: 1, evidenceClass: "input-conformance-only", phase: query.phase,
            evaluatorCommit: VivoArc2026.evaluatorCommit, querySHA256: queryHash,
            model: .init(id: plan.modelID, artifactSHA256: modelHash, trainingManifestSHA256: trainingHash,
                         excludedPerturbedContexts: plan.excludedPerturbedContexts), contexts: receipts)
        try write(result, temp, "submission.json")
        try Task.checkCancellation(); try FileManager.default.moveItem(at: temp, to: destination)
        return result
    }

    public static func verifySubmission(_ source: URL) throws -> VivoArc2026Submission {
        let temp = try stage(FileManager.default.temporaryDirectory); defer { try? FileManager.default.removeItem(at: temp) }
        let result = try readPlan(VivoArc2026Submission.self, at: source.appendingPathComponent("submission.json"))
        let (query, queryHash) = try copyQuery(source.appendingPathComponent("query"), into: temp.appendingPathComponent("query"))
        try VivoArc2026.contextIDs(result.contexts.map(\.id))
        let ids = Set(query.contexts.map(\.id))
        guard result.schemaVersion == 1, result.evidenceClass == "input-conformance-only", result.phase == query.phase,
              result.evaluatorCommit == VivoArc2026.evaluatorCommit, result.querySHA256 == queryHash,
              Set(result.contexts.map(\.id)) == ids, Set(result.model.excludedPerturbedContexts) == ids,
              result.model.excludedPerturbedContexts.count == ids.count, vivoOmicsID(result.model.id),
              VivoArc2026.sha256(result.model.artifactSHA256) else { throw VivoOmicsError.invalid("Arc submission identity mismatch") }
        let trainingHash = try snapshot(source.appendingPathComponent("training-data.json"), temp.appendingPathComponent("training-data.json"), maximum: VivoArc2026.maximumJSONBytes)
        guard trainingHash == result.model.trainingManifestSHA256,
              try JSONSerialization.jsonObject(with: read(temp, "training-data.json")) is [String: Any] else {
            throw VivoOmicsError.invalid("Arc training provenance mismatch")
        }
        for c in query.contexts {
            let receipt = result.contexts.first { $0.id == c.id }!
            let path = "contexts/" + c.id, dir = temp.appendingPathComponent(path); try mkdir(dir)
            let predicted = dir.appendingPathComponent("prediction.h5ad")
            guard try snapshot(source.appendingPathComponent(path + "/prediction.h5ad"), predicted) == receipt.predictionSHA256 else {
                throw VivoOmicsError.invalid("Arc prediction fingerprint mismatch")
            }
            let (_, actual) = try inspect(predicted, column: c.perturbationColumn, expectedGenes: c.genes, targets: c.targets, controlsOnly: false)
            guard actual == receipt.inspection else { throw VivoOmicsError.invalid("Arc prediction count reconstruction mismatch") }
        }
        return result
    }
}
