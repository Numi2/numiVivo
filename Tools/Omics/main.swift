import Foundation

@main struct OmicsChecks {
    static func main() {
        do { try run() }
        catch { FileHandle.standardError.write(Data("\(error.localizedDescription)\n".utf8)); exit(1) }
    }
    static func run() throws {
        var passed = 0
        func check(_ condition: @autoclosure () throws -> Bool, _ name: String) throws {
            guard try condition() else { throw VivoOmicsError.invalid("CHECK FAILED: \(name)") }
            passed += 1
        }
        func rejects(_ name: String, _ body: () throws -> Void) throws {
            do { try body() } catch { passed += 1; return }
            throw VivoOmicsError.invalid("CHECK DID NOT REJECT: \(name)")
        }
        let featureText = "g1\tMT-A\tGene Expression\ng2\tA\tGene Expression\ng3\tA\tGene Expression\n"
        let barcodeText = "B1\nB2\nEMPTY\n"
        let matrixText = "%%MatrixMarket matrix coordinate integer general\n% unordered fixture\n3 3 4\n3 2 7\n2 1 3\n1 1 2\n2 2 1\n"
        let sample = VivoOmicsSample(id: "s1", biologicalReplicateID: "r1", donorID: "d1", condition: "control", batchID: "b1", organism: "NCBITaxon:9606")
        func metadata(_ sample: VivoOmicsSample = sample, mitochondrial: [String] = ["g1"], groups: [String: String] = [:]) -> VivoSingleCellImport {
            .init(datasetID: sample.id, evidence: .synthetic, sourceDescription: "Synthetic integer fixture, not observations",
                  countUnit: .umiCount, sample: sample, mitochondrialFeatureIDs: mitochondrial, cellGroups: groups)
        }
        func decode(_ matrix: String = matrixText, _ features: String = featureText, _ barcodes: String = barcodeText,
                    _ config: VivoSingleCellImport? = nil, _ limits: VivoOmicsLimits = .init()) throws -> VivoSingleCellDataset {
            try VivoMatrixMarketCounts.decode(matrix: Data(matrix.utf8), features: Data(features.utf8),
                barcodes: Data(barcodes.utf8), metadata: config ?? metadata(), limits: limits)
        }
        let dataset = try decode()
        try check(dataset.matrix.rowOffsets == [0, 2, 4, 4], "feature-by-cell transpose and empty cell")
        try check(dataset.matrix.featureIndices == [0, 1, 1, 2] && dataset.matrix.counts == [2, 3, 1, 7], "sorted exact CSR")
        try check(dataset.features[1].name == dataset.features[2].name, "duplicate symbols retain distinct IDs")
        let qc = try VivoSingleCellAnalysis.quality(dataset)
        try check(qc.map(\.totalCounts) == [5, 8, 0], "quality exact totals")
        try check(qc[0].mitochondrialFraction == 0.4 && qc[2].mitochondrialFraction == nil, "mitochondrial fraction and empty-cell missingness")
        let unannotated = try decode(matrixText, featureText, barcodeText, metadata(mitochondrial: []))
        try check(VivoSingleCellAnalysis.quality(unannotated)[0].mitochondrialFraction == nil, "absent mitochondrial annotation is missing")
        let normalized = try VivoSingleCellAnalysis.logNormalize(dataset)
        try check(abs(normalized.values[0] - log1p(4_000.0)) < 1e-12, "normalization reference")
        try check(dataset.matrix.counts == [2, 3, 1, 7], "normalization preserves authoritative counts")
        try rejects("nan target") { _ = try VivoSingleCellAnalysis.logNormalize(dataset, targetSum: .nan) }
        try rejects("zero target") { _ = try VivoSingleCellAnalysis.logNormalize(dataset, targetSum: 0) }
        try rejects("real matrix header") { _ = try decode(matrixText.replacingOccurrences(of: "integer", with: "real")) }
        try rejects("negative count") { _ = try decode(matrixText.replacingOccurrences(of: "3 2 7", with: "3 2 -7")) }
        try rejects("fractional count") { _ = try decode(matrixText.replacingOccurrences(of: "3 2 7", with: "3 2 7.1")) }
        try rejects("zero stored count") { _ = try decode(matrixText.replacingOccurrences(of: "3 2 7", with: "3 2 0")) }
        try rejects("out-of-range row") { _ = try decode(matrixText.replacingOccurrences(of: "3 2 7", with: "4 2 7")) }
        try rejects("out-of-range cell") { _ = try decode(matrixText.replacingOccurrences(of: "3 2 7", with: "3 4 7")) }
        try rejects("truncated nonzeros") { _ = try decode(matrixText.replacingOccurrences(of: "3 3 4", with: "3 3 5")) }
        try rejects("excess records") { _ = try decode(matrixText + "1 3 1\n") }
        try rejects("duplicate coordinate") { _ = try decode(matrixText.replacingOccurrences(of: "3 2 7", with: "2 2 7")) }
        try rejects("duplicate feature ID") { _ = try decode(matrixText, featureText.replacingOccurrences(of: "g3", with: "g2")) }
        try rejects("duplicate barcode") { _ = try decode(matrixText, featureText, "B1\nB1\nEMPTY\n") }
        try rejects("wrong feature modality") { _ = try decode(matrixText, featureText.replacingOccurrences(of: "Gene Expression", with: "Antibody Capture")) }
        try rejects("unknown mitochondrial feature") { _ = try decode(matrixText, featureText, barcodeText, metadata(mitochondrial: ["typo"])) }
        try rejects("unknown group barcode") { _ = try decode(matrixText, featureText, barcodeText, metadata(groups: ["missing": "T"])) }
        var bytesLimit = VivoOmicsLimits(); bytesLimit.maximumInputBytes = 10
        try rejects("source byte budget") { _ = try decode(matrixText, featureText, barcodeText, nil, bytesLimit) }
        var nnzLimit = VivoOmicsLimits(); nnzLimit.maximumNonzeros = 3
        try rejects("nonzero budget") { _ = try decode(matrixText, featureText, barcodeText, nil, nnzLimit) }
        try rejects("Int.max dimensions") { _ = try decode(matrixText.replacingOccurrences(of: "3 3 4", with: "3 \(Int.max) 4")) }
        try rejects("UInt64 overflow") { _ = try decode(matrixText.replacingOccurrences(of: "1 1 2", with: "1 1 \(UInt64.max)")) }
        let largeExact = try decode(matrixText.replacingOccurrences(of: "1 1 2", with: "1 1 9007199254740993"))
        try check(largeExact.matrix.counts[0] == 9_007_199_254_740_993, "counts above Double integer precision remain exact")
        let encoded = try JSONEncoder().encode(largeExact)
        try check(JSONDecoder().decode(VivoSingleCellDataset.self, from: encoded) == largeExact, "exact JSON count round trip")
        let badVersion = String(decoding: encoded, as: UTF8.self).replacingOccurrences(of: "\"schemaVersion\":1", with: "\"schemaVersion\":99")
        try rejects("decoded unsupported schema") { try JSONDecoder().decode(VivoSingleCellDataset.self, from: Data(badVersion.utf8)).validate() }
        try rejects("malformed CSR offsets") { try VivoSparseCounts(cellCount: 1, featureCount: 2, rowOffsets: [0, -1], featureIndices: [], counts: []).validate() }
        try rejects("unsorted CSR") { try VivoSparseCounts(cellCount: 1, featureCount: 2, rowOffsets: [0, 2], featureIndices: [1, 0], counts: [1, 2]).validate() }
        let s2 = VivoOmicsSample(id: "s2", biologicalReplicateID: "r2", donorID: "d1", condition: "control", batchID: "b2", organism: sample.organism)
        let second = try decode(matrixText, featureText, barcodeText, metadata(s2))
        let combined = try VivoSingleCellAnalysis.concatenate([dataset, second], id: "combined", sourceDescription: "Synthetic libraries")
        try check(combined.cells.count == 6, "duplicate barcode across samples is valid")
        let pb = try VivoSingleCellAnalysis.pseudobulk(combined)
        try check(pb.groups.count == 2 && pb.matrix.counts == [2, 4, 7, 2, 4, 7], "replicates not pooled because they share donor ID")
        try check(pb.groups[0].sourceCellIndices == [0, 1, 2] && pb.groups[1].sourceCellIndices == [3, 4, 5], "aggregation lineage")
        let technical = VivoOmicsSample(id: "technical", biologicalReplicateID: "r1", donorID: "d1", condition: "control", batchID: "b2", organism: sample.organism)
        let third = try decode(matrixText, featureText, barcodeText, metadata(technical))
        let technicalPool = try VivoSingleCellAnalysis.concatenate([dataset, third], id: "technical", sourceDescription: "Synthetic technical replicates")
        let pooled = try VivoSingleCellAnalysis.pseudobulk(technicalPool)
        try check(pooled.groups.count == 1 && pooled.matrix.counts == [4, 8, 14], "technical libraries pool by explicit replicate")
        try check(pooled.groups[0].sampleIDs == ["s1", "technical"] && pooled.groups[0].batchIDs == ["b1", "b2"], "technical library and batch provenance")
        try rejects("duplicate sample identity during concatenate") { _ = try VivoSingleCellAnalysis.concatenate([dataset, dataset], id: "bad", sourceDescription: "bad") }
        try rejects("feature annotation mismatch during concatenate") { _ = try VivoSingleCellAnalysis.concatenate([dataset, unannotated], id: "bad", sourceDescription: "bad") }
        let grouped = try decode(matrixText, featureText, barcodeText, metadata(groups: ["B1": "T", "B2": "B"]))
        let groupedPB = try VivoSingleCellAnalysis.pseudobulk(grouped)
        try check(groupedPB.groups.count == 3 && groupedPB.groups.first?.cellGroup == nil, "unassigned cells kept distinct from cell groups")
        try check(groupedPB.groups.first?.sourceCellIndices == [2], "zero-count group retained")
        let badDonor = VivoOmicsSample(id: "bad", biologicalReplicateID: "r1", donorID: "different", condition: "control", batchID: "b1", organism: sample.organism)
        let badLibrary = try decode(matrixText, featureText, barcodeText, metadata(badDonor))
        try rejects("inconsistent replicate donor") { _ = try VivoSingleCellAnalysis.concatenate([dataset, badLibrary], id: "bad", sourceDescription: "bad") }
        try check(decode(matrixText.replacingOccurrences(of: "\n", with: "\r\n")) == dataset, "CRLF input")
        // Deterministic matrix panel checked against a separate dense reference.
        for seed in 0..<25 {
            let featureCount = 7, cellCount = 11
            var dense = [[UInt64]](repeating: [UInt64](repeating: 0, count: featureCount), count: cellCount)
            var records: [String] = []
            for cell in 0..<cellCount { for feature in 0..<featureCount {
                if (cell * 17 + feature * 3 + seed) % 5 < 2 {
                    let value = UInt64(1 + (seed + cell + feature) % 13)
                    dense[cell][feature] = value; records.append("\(feature + 1) \(cell + 1) \(value)")
                }
            } }
            let text = "%%MatrixMarket matrix coordinate integer general\n\(featureCount) \(cellCount) \(records.count)\n" + records.reversed().joined(separator: "\n") + "\n"
            let ft = (0..<featureCount).map { "g\($0)\tgene\($0)\tGene Expression" }.joined(separator: "\n")
            let bc = (0..<cellCount).map { "c\($0)" }.joined(separator: "\n")
            let imported = try decode(text, ft, bc, metadata(mitochondrial: ["g0"]))
            let metrics = try VivoSingleCellAnalysis.quality(imported)
            try check(metrics.map(\.totalCounts) == dense.map { $0.reduce(0, +) }, "dense independent cell totals \(seed)")
            let bulk = try VivoSingleCellAnalysis.pseudobulk(imported)
            var expected = [UInt64](repeating: 0, count: featureCount)
            for row in dense { for feature in 0..<featureCount { expected[feature] += row[feature] } }
            let actual = Dictionary(uniqueKeysWithValues: zip(bulk.matrix.featureIndices, bulk.matrix.counts))
            try check((0..<featureCount).allSatisfy { actual[$0, default: 0] == expected[$0] }, "dense independent pseudobulk \(seed)")
        }
        print("Omics portable native checks passed: \(passed)")
    }
}
