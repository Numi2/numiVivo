import Foundation

struct Training: Codable {
    let bulk: VivoPseudobulkCounts
    let request: VivoOmicsExpressionContrast
}
struct Fit: Codable {
    let coefficients: [Double]
    let standardDeviation: Double
    let scaledScore: Double
    let converged: Bool
    init(_ fit: VivoOmicsNBFit) {
        coefficients=fit.coefficients;standardDeviation=fit.standardError!
        scaledScore=fit.maximumScaledScore;converged=fit.converged
    }
    init(_ fit: VivoOmicsNBContrastMAPFit) {
        coefficients=fit.coefficients;standardDeviation=fit.posteriorStandardDeviation
        scaledScore=fit.maximumScaledScore;converged=fit.converged
    }
}
struct Gene: Codable {
    let featureIndex: Int
    let status: String
    let mean: Double
    let trendDispersion: Double?
    let dispersion: Double?
    let mle: Fit?
    let empirical: Fit?
    let fixed: Fit?
    let error: String?
}
struct Model: Codable {
    let request: VivoOmicsExpressionContrast
    let design: VivoOmicsDesignMatrix
    let trend: VivoOmicsNBTrend
    let prior: VivoOmicsNBEffectPriorEstimate?
    let priorError: String?
    let trainingLogLibraryCenter: Double
    let genes: [Gene]
}
struct TestCounts: Codable {
    let sampleIDs: [String]
    let conditions: [String]
    let counts: [[UInt64]]
}
struct ObservationScore: Codable {
    let sampleID: String
    let condition: String
    let libraryCounts: UInt64
    let sizeFactor: Double
    let method: String
    let genes: Int
    let meanLogScore: Double
    let meanSquaredLog1pError: Double
}
struct Scores: Codable {
    let primaryAvailable: Bool
    let testedTrainingGenes: Int
    let unavailableMAPFeatureIndices: [Int]
    let observations: [ObservationScore]
}

@main struct Main {
    static func read<T: Decodable>(_ type: T.Type,_ path: String) throws -> T {
        let raw=try VivoOmicsSourceDecoder.decode(Data(contentsOf: URL(fileURLWithPath: path)),maximumExpandedBytes: 128*1024*1024)
        return try JSONDecoder().decode(type,from: raw)
    }
    static func emit<T: Encodable>(_ value: T) throws {
        let encoder=JSONEncoder();encoder.outputFormatting=[.sortedKeys]
        FileHandle.standardOutput.write(try encoder.encode(value))
    }
    static func main() throws {
        let args=CommandLine.arguments
        guard args.count==4 else { throw VivoOmicsError.invalid("usage: nb-heldout fit metadata training | score model test") }
        if args[1]=="fit" {
            let metadata=try read(VivoSingleCellCountMetadata.self,args[2]), input=try read(Training.self,args[3])
            try metadata.validate();try input.bulk.matrix.validate()
            guard input.bulk.featureIDs==metadata.features.map(\.id),input.bulk.groups.count==input.bulk.matrix.cellCount,
                  input.request.model == .negativeBinomial,input.request.design == .independentReplicates,
                  input.request.negativeBinomialOptions?.effectPriorEstimation == .weightedUpperQuantile else {
                throw VivoOmicsError.invalid("training source axes or empirical NB design")
            }
            let result=try VivoPseudobulkDifferentialExpression.evaluate(metadata: metadata,bulk: input.bulk,contrast: input.request)
            let diagnostics=result.negativeBinomial!, design=result.design
            var counts=Array(repeating: [UInt64](repeating: 0,count: design.rows.count),count: metadata.features.count)
            for (i,row) in design.sourcePseudobulkIndices.enumerated() {
                for k in input.bulk.matrix.rowOffsets[row]..<input.bulk.matrix.rowOffsets[row+1] {
                    counts[input.bulk.matrix.featureIndices[k]][i]=input.bulk.matrix.counts[k]
                }
            }
            var genes: [Gene]=[]
            for feature in result.features {
                let diagnostic=diagnostics.features[feature.featureIndex]
                var fixed: Fit?, failure=diagnostic.effectShrinkageError
                if feature.status == .tested {
                    do {
                        fixed=Fit(try VivoOmicsNegativeBinomial.fitContrastMAP(counts: counts[feature.featureIndex],
                            design: design.rows,offsets: design.sizeFactorValues.map(log),contrast: design.contrast,
                            dispersion: diagnostic.finalDispersion!,priorStandardDeviation: log(2)))
                    } catch { failure=error.localizedDescription }
                }
                genes.append(.init(featureIndex: feature.featureIndex,status: feature.status.rawValue,mean: feature.meanNormalizedCount,
                    trendDispersion: diagnostic.trendDispersion,dispersion: diagnostic.finalDispersion,
                    mle: feature.status == .tested ? diagnostic.finalFit.map(Fit.init) : nil,
                    empirical: diagnostic.effectShrinkageFit.map(Fit.init),fixed: fixed,error: failure))
            }
            try emit(Model(request: input.request,design: design,trend: diagnostics.trend,prior: diagnostics.effectPriorEstimate,
                priorError: diagnostics.effectPriorEstimationError,
                trainingLogLibraryCenter: design.libraryCounts.map { log(Double($0)) }.reduce(0,+)/Double(design.rows.count),genes: genes))
        } else if args[1]=="score" {
            let model=try read(Model.self,args[2]), test=try read(TestCounts.self,args[3])
            guard model.request.sizeFactors == .librarySize,test.sampleIDs.count==2,test.conditions.count==2,test.counts.count==2,
                  Set(test.sampleIDs).count==2,Set(test.conditions)==Set([model.request.controlCondition,model.request.treatmentCondition]),
                  test.counts.allSatisfy({ $0.count==model.genes.count }) else { throw VivoOmicsError.invalid("held-out score axes or exposure policy") }
            let training=Set(model.design.observations.flatMap(\.sampleIDs))
            guard Set(test.sampleIDs).isDisjoint(with: training) else { throw VivoOmicsError.invalid("held-out sample overlaps training") }
            let tested=model.genes.filter { $0.status=="tested" }
            let unavailable=tested.filter { $0.empirical?.converged != true || $0.fixed?.converged != true || $0.mle?.converged != true }
            var observations: [ObservationScore]=[]
            if unavailable.isEmpty && !tested.isEmpty {
                for i in test.sampleIDs.indices {
                    let library=try test.counts[i].reduce(UInt64(0),vivoOmicsSum)
                    guard library>0 else { throw VivoOmicsError.invalid("empty held-out library") }
                    let factor=exp(log(Double(library))-model.trainingLogLibraryCenter)
                    let arm=test.conditions[i]==model.request.treatmentCondition ? 1.0 : 0.0
                    for name in ["MLE","fixed","empirical"] {
                        var score=0.0,squared=0.0
                        for gene in tested {
                            let fit=name=="MLE" ? gene.mle! : name=="fixed" ? gene.fixed! : gene.empirical!
                            guard fit.coefficients.count==2 else { throw VivoOmicsError.invalid("independent-animal design coefficients") }
                            let prediction=exp(fit.coefficients[0]+arm*fit.coefficients[1]), y=test.counts[i][gene.featureIndex]
                            score += try VivoOmicsNegativeBinomial.logMass(count: y,mean: factor*prediction,dispersion: gene.dispersion!)
                            squared += pow(log1p(Double(y)/factor)-log1p(prediction),2)
                        }
                        observations.append(.init(sampleID: test.sampleIDs[i],condition: test.conditions[i],libraryCounts: library,
                            sizeFactor: factor,method: name,genes: tested.count,meanLogScore: score/Double(tested.count),
                            meanSquaredLog1pError: squared/Double(tested.count)))
                    }
                }
            }
            try emit(Scores(primaryAvailable: unavailable.isEmpty && !tested.isEmpty,testedTrainingGenes: tested.count,
                unavailableMAPFeatureIndices: unavailable.map(\.featureIndex),observations: observations))
        } else { throw VivoOmicsError.invalid("unknown measurement command") }
    }
}
