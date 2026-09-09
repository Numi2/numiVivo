import Foundation

/// MNN publication shares the verified original-PCA/metadata contract with ridge
/// integration, while retaining its own anchor witness instead of dummy memberships.
enum VivoPCAMNNIntegration {
    static func publish(input: URL, plan: VivoPCAIntegrationPlan, implementation: VivoFingerprint, to destination: URL) throws -> VivoPCAIntegrationReceipt {
        try plan.validate(); guard let options = plan.mnn else { throw VivoOmicsError.invalid("MNN plan") }
        try VivoH5ADCountStore.requireNew(destination)
        let temp = try VivoPCAIntegration.staging(destination.deletingLastPathComponent())
        defer { try? FileManager.default.removeItem(at: temp) }
        let source = temp.appendingPathComponent("input")
        try VivoPCANeighborBundle.snapshot(input, kind: plan.inputKind, to: source)
        let parent: VivoFingerprint, dimensions: Int
        switch plan.inputKind {
        case .fitted:
            parent = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(VivoH5ADPCA.verify(source, implementation: implementation)))
            dimensions = try VivoPCAIntegration.read(VivoH5ADPCAModel.self, source, "model.json", maximum: 67_108_864).options.components
        case .query:
            parent = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(VivoH5ADPCAQuery.verify(source, implementation: implementation)))
            dimensions = try VivoPCAIntegration.read(VivoH5ADPCAQueryReport.self, source, "report.json", maximum: 1_048_576).components
        case .integrated: throw VivoOmicsError.invalid("recursive MNN input")
        }
        let metadata = try VivoPCAIntegration.read(VivoSingleCellCountMetadata.self, source, "metadata.json", maximum: 536_870_912)
        let n = metadata.cells.count
        _ = try VivoMNNIntegration.admittedBytes(rows: n, dimensions: dimensions, options: options)
        let reader = try VivoPCAScoreReader(source.appendingPathComponent("scores.bin"), rows: n, columns: dimensions)
        var values = [Double](repeating: 0, count: n*dimensions)
        for start in stride(from: 0, to: n, by: 2_048) {
            let end = min(n,start+2_048), tile = try reader.readRows(start..<end)
            values.replaceSubrange((start*dimensions)..<(end*dimensions), with: tile)
        }
        let cells = metadata.cells.map { VivoOmicsCellIdentity(sampleID: $0.sampleID, barcode: $0.barcode) }
        let result = try VivoMNNIntegration.run(cells: cells, scores: values, dimensions: dimensions, samples: metadata.samples, options: options)
        let scoreWriter = try VivoCountRecordWriter(temp.appendingPathComponent("scores.bin"))
        for i in 0..<n { for j in 0..<dimensions { try scoreWriter.append(row: i, feature: j, bits: result.scores[i*dimensions+j].bitPattern) } }
        _ = try scoreWriter.finish()
        let anchorWriter = try VivoCountRecordWriter(temp.appendingPathComponent("anchors.bin"))
        for (i,anchor) in result.anchors.enumerated() {
            try anchorWriter.append(row: anchor.first, feature: anchor.second, bits: result.anchorDistances[i].bitPattern)
        }
        _ = try anchorWriter.finish()
        let receipt = try VivoPCAIntegrationReceipt(schemaVersion: 1, input: parent,
            plan: VivoPCAIntegration.write(plan, temp, "plan.json", maximum: 65_536),
            metadata: VivoOmicsFileSnapshot.fingerprint(source.appendingPathComponent("metadata.json"), copyTo: temp.appendingPathComponent("metadata.json"), maximumBytes: 536_870_912),
            scores: VivoOmicsFileSnapshot.fingerprint(temp.appendingPathComponent("scores.bin"), maximumBytes: 1_024_000_000),
            memberships: nil, assignmentScores: nil,
            report: VivoPCAIntegration.write(result.report, temp, "report.json", maximum: 16_777_216), implementation: implementation,
            anchors: VivoOmicsFileSnapshot.fingerprint(temp.appendingPathComponent("anchors.bin"), maximumBytes: 1_600_000_000))
        _ = try VivoPCAIntegration.write(receipt, temp, "receipt.json", maximum: 65_536)
        try Task.checkCancellation(); try FileManager.default.moveItem(at: temp, to: destination)
        return receipt
    }
    static func verify(_ root: URL, plan: VivoPCAIntegrationPlan, receipt: VivoPCAIntegrationReceipt, implementation: VivoFingerprint) throws -> VivoPCAIntegrationReceipt {
        guard let anchors = receipt.anchors, receipt.memberships == nil, receipt.assignmentScores == nil else { throw VivoOmicsError.invalid("MNN integration witnesses") }
        let temp = try VivoPCAIntegration.staging(FileManager.default.temporaryDirectory)
        defer { try? FileManager.default.removeItem(at: temp) }
        let rebuilt = try publish(input: root.appendingPathComponent("input"), plan: plan, implementation: implementation, to: temp.appendingPathComponent("rebuilt"))
        guard rebuilt == receipt else { throw VivoOmicsError.invalid("MNN integration reconstruction differs") }
        for (name, hash, limit) in [("plan.json",receipt.plan,65_536), ("metadata.json",receipt.metadata,536_870_912),
            ("report.json",receipt.report,16_777_216), ("scores.bin",receipt.scores,1_024_000_000), ("anchors.bin",anchors,1_600_000_000)] {
            guard try VivoOmicsFileSnapshot.fingerprint(root.appendingPathComponent(name), maximumBytes: limit) == hash else { throw VivoOmicsError.invalid("MNN artifact fingerprint differs") }
        }
        return receipt
    }
}
