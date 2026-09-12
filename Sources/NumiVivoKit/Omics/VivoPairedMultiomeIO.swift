import Foundation

public struct VivoPairedMultiomePlan: Codable, Sendable {
    public let schemaVersion: Int
    public let mapping: VivoTenXMultiAssayPlan
    public let rnaAssayID: String
    public let atacAssayID: String
    public let rnaPCA: VivoSingleCellReductionOptions
    public let atacLSI: VivoAccessibilityLSIPlan
    public let neighborsIncludingSelf: Int
    private enum CodingKeys: String, CodingKey { case schemaVersion, mapping, rnaAssayID, atacAssayID, rnaPCA, atacLSI, neighborsIncludingSelf }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["schemaVersion","mapping","rnaAssayID","atacAssayID","rnaPCA","atacLSI","neighborsIncludingSelf"])
        let c = try decoder.container(keyedBy: CodingKeys.self)
        schemaVersion = try c.decode(Int.self, forKey: .schemaVersion)
        mapping = try c.decode(VivoTenXMultiAssayPlan.self, forKey: .mapping)
        rnaAssayID = try c.decode(String.self, forKey: .rnaAssayID)
        atacAssayID = try c.decode(String.self, forKey: .atacAssayID)
        rnaPCA = try c.decode(VivoSingleCellReductionOptions.self, forKey: .rnaPCA)
        atacLSI = try c.decode(VivoAccessibilityLSIPlan.self, forKey: .atacLSI)
        neighborsIncludingSelf = try c.decode(Int.self, forKey: .neighborsIncludingSelf)
    }
    public func validate() throws {
        try rnaPCA.validate(); try atacLSI.validate()
        guard schemaVersion == 1, rnaAssayID != atacAssayID,
              mapping.assays.contains(where: { $0.id == rnaAssayID && $0.kind == .rna }),
              mapping.assays.contains(where: { $0.id == atacAssayID && $0.kind == .chromatinAccessibility }),
              rnaPCA.components + atacLSI.components <= 64 else { throw VivoOmicsError.invalid("paired multiome plan") }
        try VivoMultiAssayTenX.validate(mapping); try rnaPCA.validate(); try atacLSI.validate()
        var options = VivoSingleCellNeighborOptions(); options.neighbors = neighborsIncludingSelf; try options.validate()
    }
}

public struct VivoPairedMultiomeReceipt: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let source: VivoFingerprint
    public let plan: VivoFingerprint
    public let result: VivoFingerprint
    public let implementation: VivoFingerprint
}

public enum VivoPairedMultiomeIO {
    private struct Result: Codable {
        let method: String
        let rnaMetadata: VivoSingleCellCountMetadata
        let atacMetadata: VivoAccessibilityTFIDFIO.Axes
        let rna: VivoSingleCellReductionResult
        let atac: VivoAccessibilityLSIResult
        let atacStatistics: VivoAccessibilityTFIDFReport
        let weighted: VivoWeightedNeighborResult
    }
    public struct Verification: Codable, Sendable {
        public let receipt: VivoPairedMultiomeReceipt
        public let verifier: VivoFingerprint
    }
    private static func compute(_ source: URL, _ plan: VivoPairedMultiomePlan) throws -> Data {
        try plan.validate()
        let data = try VivoMultiAssayTenX.readSnapshot(source, plan: plan.mapping)
        guard let rna = data.assays.first(where: { $0.id == plan.rnaAssayID }),
              let atac = data.assays.first(where: { $0.id == plan.atacAssayID }),
              rna.observationIndices == atac.observationIndices else { throw VivoOmicsError.invalid("paired assay observations differ") }
        let countUnit: VivoOmicsCountUnit
        switch rna.countUnit {
        case .umiCount: countUnit = .umiCount
        case .readCount: countUnit = .readCount
        default: throw VivoOmicsError.invalid("RNA count units")
        }
        let cells = rna.observationIndices.map { data.observations[$0].identity }
        let identities = cells.map { VivoOmicsCellIdentity(sampleID: $0.sampleID, barcode: $0.barcode) }
        let counts = VivoSingleCellDataset(id: data.id, evidence: data.evidence, sourceDescription: data.sourceDescription,
            countUnit: countUnit, samples: data.samples, features: rna.features.map { .init(id: $0.id, name: $0.name) },
            cells: cells, matrix: rna.matrix)
        var limits = VivoOmicsLimits(); limits.maximumNonzeros = 32_000_000
        let processed = try VivoSingleCellProcessing.run(counts, limits: limits)
        guard processed.sourceCellIndices == Array(cells.indices) else { throw VivoOmicsError.invalid("paired workflow would filter RNA cells") }
        let pca = try VivoSingleCellReduction.run(processed, options: plan.rnaPCA)
        let matrix = atac.matrix
        var values = [Double](); values.reserveCapacity(matrix.counts.count)
        let statistics = try VivoAccessibilityTFIDF.normalize(cells: cells.count, features: atac.features.count,
            maximumRecords: matrix.counts.count, scan: { accept in
                for row in cells.indices {
                    try Task.checkCancellation()
                    for k in matrix.rowOffsets[row]..<matrix.rowOffsets[row+1] {
                        try accept(row, matrix.featureIndices[k], matrix.counts[k])
                    }
                }
            }, emit: { _, _, value in values.append(value) })
        let lsi = try VivoAccessibilityLSI.fit(cells: identities, featureIDs: atac.features.map(\.id),
            components: plan.atacLSI.components, maximumBasis: plan.atacLSI.maximumBasis,
            relativeResidualTolerance: plan.atacLSI.relativeResidualTolerance, seed: plan.atacLSI.seed,
            scan: { accept in
                for row in cells.indices {
                    try Task.checkCancellation()
                    for k in matrix.rowOffsets[row]..<matrix.rowOffsets[row+1] {
                        try accept(row, matrix.featureIndices[k], values[k])
                    }
                }
            })
        var rnaScores = pca.scores
        for k in 0..<plan.rnaPCA.components {
            var mean = 0.0, m2 = 0.0
            for i in cells.indices { let delta = rnaScores[i][k]-mean; mean += delta/Double(i+1); m2 += delta*(rnaScores[i][k]-mean) }
            let sd = sqrt(m2/Double(cells.count-1))
            guard sd.isFinite, sd > 0 else { throw VivoOmicsError.invalid("constant RNA component") }
            for i in cells.indices { rnaScores[i][k] = (rnaScores[i][k]-mean)/sd }
        }
        func normalize(_ scores: [[Double]]) throws -> [[Double]] {
            try scores.map { row in
                let norm = sqrt(row.reduce(0) { $0+$1*$1 })
                guard norm.isFinite, norm > 0 else { throw VivoOmicsError.invalid("zero paired modality row") }
                return row.map { $0/norm }
            }
        }
        let weighted = try VivoWeightedNeighbors.fit(first: normalize(rnaScores), second: normalize(lsi.standardizedEmbeddings),
            cells: identities, neighborsIncludingSelf: plan.neighborsIncludingSelf)
        let axes = VivoAccessibilityTFIDFIO.Axes(assayID: atac.id, kind: atac.kind, featureNamespace: atac.featureNamespace,
            countUnit: atac.countUnit, genomeAssembly: atac.genomeAssembly, features: atac.features,
            sourceObservationIndices: atac.observationIndices, observations: atac.observationIndices.map { data.observations[$0] }, samples: data.samples)
        let result = try VivoCanonicalJSON.encode(Result(method: "native-paired-rna-pca-atac-lsi-weighted-graph-v1",
            rnaMetadata: counts.metadata, atacMetadata: axes, rna: pca, atac: lsi, atacStatistics: statistics, weighted: weighted))
        guard result.count <= 536_870_912 else { throw VivoOmicsError.limit("paired result bytes") }
        return result
    }
    private static func staging(_ parent: URL) throws -> URL {
        let temp = parent.appendingPathComponent(".numivivo-paired-" + UUID().uuidString)
        try FileManager.default.createDirectory(at: temp, withIntermediateDirectories: false, attributes: [.posixPermissions: 0o700])
        return temp
    }
    public static func publish(source: URL, plan: VivoPairedMultiomePlan, implementation: VivoFingerprint,
                               to destination: URL) throws -> VivoPairedMultiomeReceipt {
        try plan.validate(); try VivoH5ADCountStore.requireNew(destination)
        let temp = try staging(destination.deletingLastPathComponent()); defer { try? FileManager.default.removeItem(at: temp) }
        let snapshot = temp.appendingPathComponent("original.h5")
        let sourceID = try VivoH5ADPseudobulk.fingerprint(source, copyTo: snapshot)
        let planData = try VivoCanonicalJSON.encode(plan), result = try compute(snapshot, plan)
        let receipt = try VivoPairedMultiomeReceipt(schemaVersion: 1, source: sourceID,
            plan: VivoCanonicalJSON.fingerprint(planData), result: VivoCanonicalJSON.fingerprint(result), implementation: implementation)
        try planData.write(to: temp.appendingPathComponent("plan.json"), options: .withoutOverwriting)
        try result.write(to: temp.appendingPathComponent("result.json"), options: .withoutOverwriting)
        try VivoCanonicalJSON.encode(receipt).write(to: temp.appendingPathComponent("receipt.json"), options: .withoutOverwriting)
        try Task.checkCancellation(); try FileManager.default.moveItem(at: temp, to: destination)
        return receipt
    }
    public static func verify(_ directory: URL, implementation: VivoFingerprint) throws -> Verification {
        func read(_ name: String, _ limit: Int) throws -> Data {
            try VivoSingleCellCampaignIO.readDocument(directory.appendingPathComponent(name), maximumBytes: limit)
        }
        let receipt = try VivoCanonicalJSON.decode(VivoPairedMultiomeReceipt.self, from: read("receipt.json", 65_536))
        guard receipt.schemaVersion == 1 else { throw VivoOmicsError.invalid("paired receipt schema") }
        let planData = try read("plan.json", 2_097_152), result = try read("result.json", 536_870_912)
        guard try VivoCanonicalJSON.fingerprint(planData) == receipt.plan,
              try VivoCanonicalJSON.fingerprint(result) == receipt.result else { throw VivoOmicsError.invalid("paired artifact hash") }
        let plan = try VivoCanonicalJSON.decode(VivoPairedMultiomePlan.self, from: planData)
        let temp = try staging(FileManager.default.temporaryDirectory); defer { try? FileManager.default.removeItem(at: temp) }
        let snapshot = temp.appendingPathComponent("original.h5")
        guard try VivoH5ADPseudobulk.fingerprint(directory.appendingPathComponent("original.h5"), copyTo: snapshot) == receipt.source else {
            throw VivoOmicsError.invalid("paired source hash")
        }
        guard try compute(snapshot, plan) == result else { throw VivoOmicsError.invalid("paired result differs from source reconstruction") }
        return .init(receipt: receipt, verifier: implementation)
    }
}
