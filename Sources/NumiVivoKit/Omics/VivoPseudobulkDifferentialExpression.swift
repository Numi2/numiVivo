import Foundation

public enum VivoOmicsExpressionModel: String, Codable, Sendable { case logLinear, negativeBinomial }
public enum VivoOmicsReplicationDesign: String, Codable, Sendable { case independentReplicates, pairedDonors }
public enum VivoOmicsSizeFactorMethod: String, Codable, Sendable { case librarySize, medianRatio }
public enum VivoOmicsVarianceMethod: String, Codable, Sendable { case ordinary, empiricalBayes }
public struct VivoOmicsExpressionContrast: Codable, Sendable, Equatable {
    public var id: String
    public var model: VivoOmicsExpressionModel = .logLinear
    public var negativeBinomialOptions: VivoOmicsNBCohortOptions?
    public var controlCondition: String
    public var treatmentCondition: String
    public var cellGroup: String?
    public var design: VivoOmicsReplicationDesign
    public var sizeFactors: VivoOmicsSizeFactorMethod = .medianRatio
    public var variance: VivoOmicsVarianceMethod = .empiricalBayes
    public var adjustForBatch: Bool = true
    public var minimumCellsPerPseudobulk: Int = 10
    public var minimumReplicatesPerCondition: Int = 3
    public var minimumFeatureCounts: UInt64 = 10
    public var minimumExpressingPseudobulks: Int = 3
    public var minimumReferenceFeatures: Int = 10
    public var priorCount: Double = 0.5
    public var intervalCoverage: Double = 0.95
    public var includedDonorIDs: [String]?
    public init(id: String, controlCondition: String, treatmentCondition: String,
                cellGroup: String? = nil, design: VivoOmicsReplicationDesign) {
        self.id = id; self.controlCondition = controlCondition; self.treatmentCondition = treatmentCondition
        self.cellGroup = cellGroup; self.design = design
    }
    private enum CodingKeys: String, CodingKey {
        case id, model, negativeBinomialOptions, controlCondition, treatmentCondition, cellGroup, design, sizeFactors, variance, adjustForBatch
        case minimumCellsPerPseudobulk, minimumReplicatesPerCondition, minimumFeatureCounts
        case minimumExpressingPseudobulks, minimumReferenceFeatures, priorCount, intervalCoverage, includedDonorIDs
    }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["id", "model", "negativeBinomialOptions", "controlCondition", "treatmentCondition", "cellGroup", "design",
            "sizeFactors", "variance", "adjustForBatch", "minimumCellsPerPseudobulk", "minimumReplicatesPerCondition",
            "minimumFeatureCounts", "minimumExpressingPseudobulks", "minimumReferenceFeatures", "priorCount", "intervalCoverage", "includedDonorIDs"])
        let c = try decoder.container(keyedBy: CodingKeys.self)
        id = try c.decode(String.self, forKey: .id)
        model = try c.decodeIfPresent(VivoOmicsExpressionModel.self, forKey: .model) ?? .logLinear
        negativeBinomialOptions = try c.decodeIfPresent(VivoOmicsNBCohortOptions.self, forKey: .negativeBinomialOptions)
        controlCondition = try c.decode(String.self, forKey: .controlCondition)
        treatmentCondition = try c.decode(String.self, forKey: .treatmentCondition)
        cellGroup = try c.decodeIfPresent(String.self, forKey: .cellGroup)
        design = try c.decode(VivoOmicsReplicationDesign.self, forKey: .design)
        sizeFactors = try c.decodeIfPresent(VivoOmicsSizeFactorMethod.self, forKey: .sizeFactors) ?? .medianRatio
        variance = try c.decodeIfPresent(VivoOmicsVarianceMethod.self, forKey: .variance) ?? .empiricalBayes
        adjustForBatch = try c.decodeIfPresent(Bool.self, forKey: .adjustForBatch) ?? true
        minimumCellsPerPseudobulk = try c.decodeIfPresent(Int.self, forKey: .minimumCellsPerPseudobulk) ?? 10
        minimumReplicatesPerCondition = try c.decodeIfPresent(Int.self, forKey: .minimumReplicatesPerCondition) ?? 3
        minimumFeatureCounts = try c.decodeIfPresent(UInt64.self, forKey: .minimumFeatureCounts) ?? 10
        minimumExpressingPseudobulks = try c.decodeIfPresent(Int.self, forKey: .minimumExpressingPseudobulks) ?? 3
        minimumReferenceFeatures = try c.decodeIfPresent(Int.self, forKey: .minimumReferenceFeatures) ?? 10
        priorCount = try c.decodeIfPresent(Double.self, forKey: .priorCount) ?? 0.5
        intervalCoverage = try c.decodeIfPresent(Double.self, forKey: .intervalCoverage) ?? 0.95
        includedDonorIDs = try c.decodeIfPresent([String].self, forKey: .includedDonorIDs)
    }
    public func validate() throws {
        guard [id, controlCondition, treatmentCondition].allSatisfy(vivoOmicsID), controlCondition != treatmentCondition,
              cellGroup.map(vivoOmicsID) ?? true, minimumCellsPerPseudobulk > 0,
              (3...256).contains(minimumReplicatesPerCondition), minimumExpressingPseudobulks > 0,
              minimumReferenceFeatures >= 3, minimumReferenceFeatures <= 100_000,
              priorCount.isFinite, priorCount > 0, priorCount <= 1_000,
              intervalCoverage.isFinite, intervalCoverage >= 0.5, intervalCoverage <= 0.9999 else {
            throw VivoOmicsError.invalid("expression contrast identity, replication or numerical settings")
        }
        if model == .negativeBinomial {
            guard variance == .empiricalBayes, priorCount == 0.5 else {
                throw VivoOmicsError.invalid("NB does not accept overridden log-linear variance/priorCount settings")
            }
            try (negativeBinomialOptions ?? .init()).validate()
        } else if negativeBinomialOptions != nil {
            throw VivoOmicsError.invalid("NB options require negativeBinomial model")
        }
        if let ids = includedDonorIDs {
            guard !ids.isEmpty, ids.count <= 512, Set(ids).count == ids.count, ids.allSatisfy(vivoOmicsID) else {
                throw VivoOmicsError.invalid("contrast donor subset")
            }
        }
    }
}
public struct VivoOmicsDesignMatrix: Codable, Sendable, Equatable {
    public let columnNames: [String]
    public let rows: [[Double]]
    public let contrast: [Double]
    public let sourcePseudobulkIndices: [Int]
    public let observations: [VivoPseudobulkGroup]
    public let excludedSmallPseudobulkIndices: [Int]
    public let controlReplicates: Int
    public let treatmentReplicates: Int
    public let residualDegreesOfFreedom: Int
    public let sizeFactorValues: [Double]
    public let libraryCounts: [UInt64]
    public let referenceFeatureIndices: [Int]
}
public enum VivoOmicsExpressionStatus: String, Codable, Sendable { case tested, filteredLowExpression, zeroResidualVariance, rankDeficientSupport, numericalFailure, dispersionBoundary, influentialObservation }
public struct VivoOmicsExpressionFeature: Codable, Sendable, Equatable {
    public let featureIndex: Int
    public let featureID: String
    public let status: VivoOmicsExpressionStatus
    public let totalCounts: UInt64
    public let expressingPseudobulks: Int
    public let meanNormalizedCount: Double
    public let log2FoldChange: Double?
    public let residualVariance: Double?
    public let posteriorVariance: Double?
    public let standardError: Double?
    public let tStatistic: Double?
    public let degreesOfFreedom: Double?
    public let intervalLower: Double?
    public let intervalUpper: Double?
    public let pValue: Double?
    public var adjustedPValue: Double?
    public var zStatistic: Double? = nil
}
public struct VivoOmicsExpressionResult: Codable, Sendable, Equatable {
    public let method: String
    public let request: VivoOmicsExpressionContrast
    public let evidence: VivoOmicsEvidence
    public let design: VivoOmicsDesignMatrix
    public let variancePrior: VivoOmicsVariancePrior?
    public let features: [VivoOmicsExpressionFeature]
    public let testedFeatures: Int
    public let multiplicityScope: String
    public var negativeBinomial: VivoOmicsNBCohortDiagnostics? = nil
}

public enum VivoPseudobulkDifferentialExpression {
    public static func run(_ dataset: VivoSingleCellDataset, contrast: VivoOmicsExpressionContrast,
                           limits: VivoOmicsLimits = .init()) throws -> VivoOmicsExpressionResult {
        try dataset.validate(limits: limits)
        return try evaluate(dataset, bulk: VivoSingleCellAnalysis.pseudobulk(dataset, limits: limits), contrast: contrast)
    }
    /// The shared processing route calls this only with a freshly computed bulk.
    static func evaluate(_ dataset: VivoSingleCellDataset, bulk: VivoPseudobulkCounts,
                         contrast request: VivoOmicsExpressionContrast) throws -> VivoOmicsExpressionResult {
        try request.validate()
        let donorSubset = request.includedDonorIDs.map { Set($0) }
        let knownDonors = Set(bulk.groups.compactMap(\.donorID))
        guard donorSubset.map({ $0.isSubset(of: knownDonors) }) ?? true else { throw VivoOmicsError.invalid("unknown selected donor") }
        var indices: [Int] = [], excluded: [Int] = []
        for (i, group) in bulk.groups.enumerated() where group.cellGroup == request.cellGroup &&
            [request.controlCondition, request.treatmentCondition].contains(group.condition) {
            if let donors = donorSubset, !(group.donorID.map(donors.contains) ?? false) { continue }
            if group.sourceCellIndices.count < request.minimumCellsPerPseudobulk { excluded.append(i) }
            else { indices.append(i) }
        }
        guard indices.count <= 512 else { throw VivoOmicsError.limit("at most 512 pseudobulk observations per contrast") }
        let observations = indices.map { bulk.groups[$0] }, n = observations.count
        guard Set(observations.map(\.organism)).count == 1 else { throw VivoOmicsError.invalid("contrast must contain exactly one organism") }
        let controlCount = observations.filter { $0.condition == request.controlCondition }.count, treatmentCount = n - controlCount
        guard controlCount >= request.minimumReplicatesPerCondition, treatmentCount >= request.minimumReplicatesPerCondition,
              request.minimumExpressingPseudobulks <= n else {
            throw VivoOmicsError.invalid("insufficient biological replication after selection; cells are not independent replicates")
        }
        var names = ["intercept", "treatment-minus-control"]
        var design = observations.map { [1.0, $0.condition == request.treatmentCondition ? 1.0 : 0.0] }
        if request.design == .pairedDonors {
            guard observations.allSatisfy({ $0.donorID != nil }) else { throw VivoOmicsError.invalid("paired design requires donor IDs") }
            let donors = Set(observations.compactMap(\.donorID)).sorted()
            for donor in donors {
                let rows = observations.filter { $0.donorID == donor }
                guard rows.count == 2, Set(rows.map(\.condition)).count == 2 else {
                    throw VivoOmicsError.invalid("each paired donor must have exactly one pseudobulk in both conditions; incomplete pairs are not silently removed")
                }
            }
            for donor in donors.dropFirst() {
                names.append("donor:" + donor)
                for row in 0..<n { design[row].append(observations[row].donorID == donor ? 1 : 0) }
            }
        } else {
            let donorIDs = observations.compactMap(\.donorID)
            guard Set(donorIDs).count == donorIDs.count,
                  Set(observations.map(\.biologicalReplicateID)).count == n else {
                throw VivoOmicsError.invalid("repeated donors or biological replicates require an explicit paired design")
            }
        }
        if request.adjustForBatch {
            guard observations.allSatisfy({ $0.batchIDs.count == 1 }) else {
                throw VivoOmicsError.invalid("a pooled row spans multiple batches; categorical batch adjustment is not identifiable from this aggregation")
            }
            for batch in Set(observations.flatMap(\.batchIDs)).sorted().dropFirst() {
                names.append("batch:" + batch)
                for row in 0..<n { design[row].append(observations[row].batchIDs[0] == batch ? 1 : 0) }
            }
        }
        let qr = try VivoOmicsQR(design: design)
        var c = [Double](repeating: 0, count: names.count); c[1] = 1
        // Sparse feature-major access, built once. No genes-by-cells dense matrix.
        var entries = [[(row: Int, count: UInt64)]](repeating: [], count: dataset.features.count)
        var libraries = [UInt64](repeating: 0, count: n)
        for (row, sourceRow) in indices.enumerated() {
            for k in bulk.matrix.rowOffsets[sourceRow]..<bulk.matrix.rowOffsets[sourceRow + 1] {
                let count = bulk.matrix.counts[k]
                guard count <= 9_007_199_254_740_992 else { throw VivoOmicsError.invalid("inference rejects counts outside exact FP64 integer range; raw import/export remains UInt64") }
                libraries[row] = try vivoOmicsSum(libraries[row], count)
                entries[bulk.matrix.featureIndices[k]].append((row, count))
            }
            guard libraries[row] > 0, libraries[row] <= 9_007_199_254_740_992 else {
                throw VivoOmicsError.invalid("zero or unsupported-size pseudobulk library")
            }
        }
        let referenceFeatures: [Int]
        var logFactors: [Double]
        if request.sizeFactors == .medianRatio {
            referenceFeatures = entries.indices.filter { entries[$0].count == n }
            guard referenceFeatures.count >= request.minimumReferenceFeatures else {
                throw VivoOmicsError.invalid("too few all-positive reference genes for median-ratio factors; choose librarySize explicitly rather than silently change normalization")
            }
            var ratios = [[Double]](repeating: [], count: n)
            for feature in referenceFeatures {
                let logged = entries[feature].map { log(Double($0.count)) }, center = logged.reduce(0, +) / Double(n)
                for entry in entries[feature] { ratios[entry.row].append(exp(log(Double(entry.count)) - center)) }
            }
            logFactors = try ratios.map { log(try VivoOmicsLinearStatistics.median($0)) }
        } else { referenceFeatures = []; logFactors = libraries.map { log(Double($0)) } }
        let factorCenter = logFactors.reduce(0, +) / Double(n)
        let factors = logFactors.map { exp($0 - factorCenter) }
        guard factors.allSatisfy({ $0.isFinite && $0 > 0 }) else { throw VivoOmicsError.invalid("nonfinite library normalization") }
        let df = n - names.count
        let matrix = VivoOmicsDesignMatrix(columnNames: names, rows: design, contrast: c, sourcePseudobulkIndices: indices,
            observations: observations, excludedSmallPseudobulkIndices: excluded, controlReplicates: controlCount,
            treatmentReplicates: treatmentCount, residualDegreesOfFreedom: df, sizeFactorValues: factors,
            libraryCounts: libraries, referenceFeatureIndices: referenceFeatures)
        if request.model == .negativeBinomial {
            return try VivoOmicsNBCohort.evaluate(dataset: dataset, entries: entries, design: matrix, request: request)
        }
        var fits = [VivoOmicsLinearFit?](repeating: nil, count: entries.count)
        var totals = [UInt64](repeating: 0, count: entries.count), means = [Double](repeating: 0, count: entries.count)
        var modelWork = 0
        for feature in entries.indices {
            if feature % 128 == 0 { try Task.checkCancellation() }
            for entry in entries[feature] {
                totals[feature] = try vivoOmicsSum(totals[feature], entry.count)
                means[feature] += Double(entry.count) / factors[entry.row] / Double(n)
            }
            if totals[feature] < request.minimumFeatureCounts || entries[feature].count < request.minimumExpressingPseudobulks { continue }
            modelWork += n * names.count
            guard modelWork <= 250_000_000 else { throw VivoOmicsError.limit("contrast exceeds 250 million design-response products") }
            var response = [Double](repeating: log2(request.priorCount), count: n)
            for entry in entries[feature] { response[entry.row] = log2(Double(entry.count) / factors[entry.row] + request.priorCount) }
            fits[feature] = try qr.fit(response, contrast: c)
        }
        let eligible = fits.compactMap { $0 }
        guard !eligible.isEmpty else { throw VivoOmicsError.invalid("no features satisfy the declared expression filter") }
        let prior: VivoOmicsVariancePrior? = request.variance == .empiricalBayes ?
            try VivoOmicsLinearStatistics.variancePrior(eligible.map { $0.residualVariance <= 1e-24 ? 0 : $0.residualVariance }, residualDF: df) : nil
        let totalDF = min(Double(df) + (prior?.degreesOfFreedom ?? 0), Double(df) * Double(eligible.count))
        let critical = try VivoOmicsLinearStatistics.studentCriticalValue(degreesOfFreedom: totalDF, coverage: request.intervalCoverage)
        var features: [VivoOmicsExpressionFeature] = [], tested: [Int] = [], probabilities: [Double] = []
        for feature in entries.indices {
            guard let fit = fits[feature] else {
                features.append(.init(featureIndex: feature, featureID: dataset.features[feature].id, status: .filteredLowExpression,
                    totalCounts: totals[feature], expressingPseudobulks: entries[feature].count, meanNormalizedCount: means[feature],
                    log2FoldChange: nil, residualVariance: nil, posteriorVariance: nil, standardError: nil, tStatistic: nil,
                    degreesOfFreedom: nil, intervalLower: nil, intervalUpper: nil, pValue: nil, adjustedPValue: nil)); continue
            }
            let posterior = prior.map { (Double(df) * fit.residualVariance + $0.degreesOfFreedom * $0.variance) / (Double(df) + $0.degreesOfFreedom) } ?? fit.residualVariance
            if posterior <= 1e-24 {
                features.append(.init(featureIndex: feature, featureID: dataset.features[feature].id, status: .zeroResidualVariance,
                    totalCounts: totals[feature], expressingPseudobulks: entries[feature].count, meanNormalizedCount: means[feature],
                    log2FoldChange: fit.effect, residualVariance: fit.residualVariance, posteriorVariance: posterior, standardError: nil,
                    tStatistic: nil, degreesOfFreedom: totalDF, intervalLower: nil, intervalUpper: nil, pValue: nil, adjustedPValue: nil)); continue
            }
            let error = sqrt(posterior * fit.contrastVarianceScale), t = fit.effect / error
            let probability = try VivoOmicsLinearStatistics.studentTwoSidedP(t: t, degreesOfFreedom: totalDF)
            tested.append(feature); probabilities.append(probability)
            features.append(.init(featureIndex: feature, featureID: dataset.features[feature].id, status: .tested,
                totalCounts: totals[feature], expressingPseudobulks: entries[feature].count, meanNormalizedCount: means[feature],
                log2FoldChange: fit.effect, residualVariance: fit.residualVariance, posteriorVariance: posterior, standardError: error,
                tStatistic: t, degreesOfFreedom: totalDF, intervalLower: fit.effect - critical * error, intervalUpper: fit.effect + critical * error,
                pValue: probability, adjustedPValue: nil))
        }
        let adjusted = try VivoOmicsLinearStatistics.benjaminiHochberg(probabilities)
        for (index, feature) in tested.enumerated() { features[feature].adjustedPValue = adjusted[index] }
        return .init(method: "donor-aware-log2-normalized-pseudobulk-linear-model-v1", request: request,
            evidence: dataset.evidence, design: matrix, variancePrior: prior, features: features, testedFeatures: tested.count,
            multiplicityScope: "Benjamini-Hochberg across tested features in this contrast only; no across-contrast or selective-inference guarantee")
    }
}

