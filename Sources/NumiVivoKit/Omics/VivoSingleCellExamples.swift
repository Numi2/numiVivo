import Foundation

public enum VivoSingleCellExamples {
    /// Six declared paired donors, two nonempty cells and one empty barcode in
    /// each condition/library. All values are synthetic and explicitly labelled.
    public static func pairedCounts() throws -> VivoSingleCellDataset {
        var samples: [VivoOmicsSample] = [], cells: [VivoOmicsCell] = []
        var offsets = [0], columns: [Int] = [], counts: [UInt64] = []
        let features = (0..<32).map { VivoOmicsFeature(id: "g\($0)", name: "synthetic-gene-\($0)", mitochondrial: $0 == 31) }
        for donor in 0..<6 { for condition in 0..<2 {
            let name = condition == 0 ? "control" : "treated", sample = "sample-\(donor)-\(condition)"
            samples.append(.init(id: sample, biologicalReplicateID: "replicate-\(donor)", donorID: "donor-\(donor)",
                condition: name, batchID: "shared-batch", organism: "synthetic-organism"))
            for cell in 0..<3 {
                cells.append(.init(barcode: "cell-\(cell)", sampleID: sample, group: "declared-group"))
                if cell < 2 { for gene in 0..<32 {
                    let raw = 50 + gene * 3 + donor * 7 + cell * 2 + (donor * 11 + gene * 5 + condition * gene * 3 + cell * 7) % 17
                    columns.append(gene); counts.append(UInt64(raw * (condition == 1 && gene < 3 ? 4 : 1)))
                } }
                offsets.append(counts.count)
            }
        } }
        let dataset = VivoSingleCellDataset(id: "synthetic-paired-donor-example", evidence: .synthetic,
            sourceDescription: "Synthetic six-donor count fixture, not experimental observations or statistical calibration",
            countUnit: .umiCount, samples: samples, features: features, cells: cells,
            matrix: .init(cellCount: cells.count, featureCount: features.count, rowOffsets: offsets, featureIndices: columns, counts: counts))
        try dataset.validate(); return dataset
    }
    public static func pairedPlan() -> VivoSingleCellAnalysisPlan {
        var contrast = VivoOmicsExpressionContrast(id: "treated-minus-control", controlCondition: "control", treatmentCondition: "treated",
            cellGroup: "declared-group", design: .pairedDonors)
        contrast.minimumCellsPerPseudobulk = 2
        return .init(id: "synthetic-paired-donor-analysis", contrasts: [contrast])
    }
}
