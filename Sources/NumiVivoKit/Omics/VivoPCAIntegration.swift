import Foundation

public struct VivoPCAIntegrationPlan: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let inputKind: VivoPCANeighborPlan.InputKind
    public let integration: VivoSingleCellIntegrationOptions?
    public let mnn: VivoMNNIntegrationOptions?
    public init(inputKind: VivoPCANeighborPlan.InputKind = .fitted, integration: VivoSingleCellIntegrationOptions = .init()) {
        schemaVersion = 1; self.inputKind = inputKind; self.integration = integration; mnn = nil
    }
    public init(inputKind: VivoPCANeighborPlan.InputKind = .fitted, mnn: VivoMNNIntegrationOptions) {
        schemaVersion = 1; self.inputKind = inputKind; integration = nil; self.mnn = mnn
    }
    private enum CodingKeys: String, CodingKey { case schemaVersion, inputKind, integration, mnn }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["schemaVersion", "inputKind", "integration", "mnn"])
        let c = try decoder.container(keyedBy: CodingKeys.self)
        schemaVersion = try c.decode(Int.self, forKey: .schemaVersion)
        inputKind = try c.decodeIfPresent(VivoPCANeighborPlan.InputKind.self, forKey: .inputKind) ?? .fitted
        mnn = try c.decodeIfPresent(VivoMNNIntegrationOptions.self, forKey: .mnn)
        integration = try c.decodeIfPresent(VivoSingleCellIntegrationOptions.self, forKey: .integration) ?? (mnn == nil ? .init() : nil)
    }
    public func validate() throws {
        guard (integration == nil) != (mnn == nil) else { throw VivoOmicsError.invalid("select exactly one integration method") }
        try integration?.validate(); try mnn?.validate()
        guard schemaVersion == 1, inputKind != .integrated else { throw VivoOmicsError.invalid("integration needs original fitted or query PCA") }
    }
}
public struct VivoPCAIntegrationReport: Codable, Sendable, Equatable {
    public let method: String
    public let cells: Int
    public let components: Int
    public let clusters: Int
    public let levels: [String]
    public let cellLevels: [Int]
    public let assignmentCenters: [[Double]]
    public let objectives: [Double]
    public let relativeImprovements: [Double]
    public let stoppingReason: String
    public let maximumRidgeResidual: Double
    public let scratchBytes: Int
    public let maximumMappedBytesPerMatrix: Int
    public let maximumSimultaneousMatrices: Int
    public let qualification: String
    public var ridgePenalties: [[Double]]? = nil
}
public struct VivoPCAIntegrationReceipt: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let input: VivoFingerprint
    public let plan: VivoFingerprint
    public let metadata: VivoFingerprint
    public let scores: VivoFingerprint
    public let memberships: VivoFingerprint?
    public let assignmentScores: VivoFingerprint?
    public let report: VivoFingerprint
    public let implementation: VivoFingerprint
    public var anchors: VivoFingerprint? = nil
}

public enum VivoPCAIntegration {
    static func components(_ root: URL) throws -> Int {
        struct Dimensions: Decodable { let components: Int }
        return try read(Dimensions.self, root, "report.json", maximum: 16_777_216).components
    }
    static func read<T: Decodable>(_ type: T.Type, _ root: URL, _ name: String, maximum: Int) throws -> T {
        try VivoCanonicalJSON.decode(type, from: VivoH5ADCountStore.read(root, name, maximum: maximum))
    }
    static func write<T: Encodable>(_ value: T, _ root: URL, _ name: String, maximum: Int) throws -> VivoFingerprint {
        let bytes = try VivoCanonicalJSON.encode(value)
        guard bytes.count <= maximum else { throw VivoOmicsError.limit("integration artifact bytes") }
        try bytes.write(to: root.appendingPathComponent(name), options: .withoutOverwriting)
        return try VivoCanonicalJSON.fingerprint(bytes)
    }
    static func staging(_ parent: URL) throws -> URL {
        let root = parent.appendingPathComponent(".numivivo-integration-" + UUID().uuidString)
        try FileManager.default.createDirectory(at: root, withIntermediateDirectories: false, attributes: [.posixPermissions: 0o700]); return root
    }
    static func snapshot(_ input: URL, to output: URL) throws {
        let plan = try read(VivoPCAIntegrationPlan.self, input, "plan.json", maximum: 65_536)
        try plan.validate()
        try FileManager.default.createDirectory(at: output, withIntermediateDirectories: false, attributes: [.posixPermissions: 0o700])
        try VivoPCANeighborBundle.snapshot(input.appendingPathComponent("input"), kind: plan.inputKind, to: output.appendingPathComponent("input"))
        let witnesses = plan.mnn == nil ? [("assignment-scores.bin", 1_024_000_000), ("memberships.bin", 1_600_000_000)] : [("anchors.bin", 1_600_000_000)]
        for (name, limit) in [("plan.json", 65_536), ("receipt.json", 65_536), ("report.json", 16_777_216), ("metadata.json", 536_870_912),
            ("scores.bin", 1_024_000_000)] + witnesses {
            _ = try VivoOmicsFileSnapshot.fingerprint(input.appendingPathComponent(name), copyTo: output.appendingPathComponent(name), maximumBytes: limit)
        }
    }
    public static func publish(input: URL, plan: VivoPCAIntegrationPlan, implementation: VivoFingerprint, to destination: URL) throws -> VivoPCAIntegrationReceipt {
        try plan.validate()
        if plan.mnn != nil { return try VivoPCAMNNIntegration.publish(input: input, plan: plan, implementation: implementation, to: destination) }
        try VivoH5ADCountStore.requireNew(destination)
        let temp = try staging(destination.deletingLastPathComponent()); defer { try? FileManager.default.removeItem(at: temp) }
        let source = temp.appendingPathComponent("input")
        try VivoPCANeighborBundle.snapshot(input, kind: plan.inputKind, to: source)
        let parent: VivoFingerprint, dimensions: Int
        switch plan.inputKind {
        case .fitted:
            parent = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(VivoH5ADPCA.verify(source, implementation: implementation)))
            dimensions = try read(VivoH5ADPCAModel.self, source, "model.json", maximum: 67_108_864).options.components
        case .query:
            parent = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(VivoH5ADPCAQuery.verify(source, implementation: implementation)))
            dimensions = try read(VivoH5ADPCAQueryReport.self, source, "report.json", maximum: 1_048_576).components
        case .integrated: throw VivoOmicsError.invalid("recursive integration input")
        }
        let metadata = try read(VivoSingleCellCountMetadata.self, source, "metadata.json", maximum: 536_870_912)
        let n = metadata.cells.count
        try VivoSingleCellIntegration.validateAxes(rows: n, columns: dimensions, options: plan.integration!)
        var matrices: [VivoIntegrationMatrix] = []
        defer { for matrix in matrices { try? matrix.remove() } }
        func matrix(_ rows: Int, _ columns: Int) throws -> VivoIntegrationMatrix {
            let value = try VivoIntegrationMatrix(rows: rows, columns: columns, scratch: temp)
            matrices.append(value); return value
        }
        let x = try matrix(n, dimensions)
        let reader = try VivoPCAScoreReader(source.appendingPathComponent("scores.bin"), rows: n, columns: dimensions)
        for start in stride(from: 0, to: n, by: 2_048) {
            let end = min(n, start + 2_048), tile = try reader.readRows(start..<end)
            for i in start..<end {
                let offset = (i - start) * dimensions
                try x.setRow(i, Array(tile[offset..<(offset + dimensions)]))
            }
        }
        let cells = metadata.cells.map { VivoOmicsCellIdentity(sampleID: $0.sampleID, barcode: $0.barcode) }
        let result = try VivoSingleCellIntegration.run(cells: cells, x: x, samples: metadata.samples, options: plan.integration!, matrix: matrix)
        let report = VivoPCAIntegrationReport(method: VivoIntegrationSolution.method(options: plan.integration!), cells: n, components: dimensions, clusters: plan.integration!.clusters,
            levels: result.levels, cellLevels: result.cellLevels, assignmentCenters: result.assignmentCenters, objectives: result.objectives,
            relativeImprovements: result.relativeImprovements, stoppingReason: result.stoppingReason, maximumRidgeResidual: result.maximumRidgeResidual,
            scratchBytes: matrices.reduce(0) { $0 + $1.fileBytes }, maximumMappedBytesPerMatrix: VivoIntegrationMatrix.maximumWindowBytes,
            maximumSimultaneousMatrices: matrices.count,
            qualification: VivoIntegrationSolution.qualification + " Latent matrices use private row-major f64 scratch with one 64 MiB mapping per matrix; cell identities, covariate indices, permutations and cluster statistics remain resident. Count and PCA parent reconstruction is included. Final matrices use complete row-major u32-row/u32-column/f64 little-endian records. Work budget is an admission index, not an operation counter. No million-cell, Metal or biological qualification follows from file storage.", ridgePenalties: result.ridgePenalties)
        let receipt = try VivoPCAIntegrationReceipt(schemaVersion: 1, input: parent,
            plan: write(plan, temp, "plan.json", maximum: 65_536),
            metadata: VivoOmicsFileSnapshot.fingerprint(source.appendingPathComponent("metadata.json"), copyTo: temp.appendingPathComponent("metadata.json"), maximumBytes: 536_870_912),
            scores: result.scores.writeRecords(to: temp.appendingPathComponent("scores.bin")),
            memberships: result.memberships.writeRecords(to: temp.appendingPathComponent("memberships.bin")),
            assignmentScores: result.assignmentScores.writeRecords(to: temp.appendingPathComponent("assignment-scores.bin")),
            report: write(report, temp, "report.json", maximum: 16_777_216), implementation: implementation)
        for matrix in matrices { try matrix.remove() }
        _ = try write(receipt, temp, "receipt.json", maximum: 65_536)
        try Task.checkCancellation(); try FileManager.default.moveItem(at: temp, to: destination); return receipt
    }
    public static func verify(_ root: URL, implementation: VivoFingerprint) throws -> VivoPCAIntegrationReceipt {
        let bytes = try VivoH5ADCountStore.read(root, "receipt.json", maximum: 65_536)
        let receipt = try VivoCanonicalJSON.decode(VivoPCAIntegrationReceipt.self, from: bytes)
        guard receipt.schemaVersion == 1, receipt.implementation == implementation, try VivoCanonicalJSON.encode(receipt) == bytes else { throw VivoOmicsError.invalid("integration receipt") }
        let plan = try read(VivoPCAIntegrationPlan.self, root, "plan.json", maximum: 65_536)
        try plan.validate()
        if plan.mnn != nil { return try VivoPCAMNNIntegration.verify(root, plan: plan, receipt: receipt, implementation: implementation) }
        guard receipt.anchors == nil, receipt.memberships != nil, receipt.assignmentScores != nil else { throw VivoOmicsError.invalid("ridge integration witnesses") }
        let temp = try staging(FileManager.default.temporaryDirectory); defer { try? FileManager.default.removeItem(at: temp) }
        let rebuilt = try publish(input: root.appendingPathComponent("input"), plan: plan, implementation: implementation, to: temp.appendingPathComponent("rebuilt"))
        guard rebuilt == receipt else { throw VivoOmicsError.invalid("integration reconstruction differs") }
        for (name, hash, limit) in [("plan.json", receipt.plan, 65_536), ("metadata.json", receipt.metadata, 536_870_912),
            ("report.json", receipt.report, 16_777_216), ("scores.bin", receipt.scores, 1_024_000_000),
            ("assignment-scores.bin", receipt.assignmentScores!, 1_024_000_000), ("memberships.bin", receipt.memberships!, 1_600_000_000)] {
            guard try VivoOmicsFileSnapshot.fingerprint(root.appendingPathComponent(name), maximumBytes: limit) == hash else { throw VivoOmicsError.invalid("integration artifact fingerprint differs") }
        }
        return receipt
    }
}
