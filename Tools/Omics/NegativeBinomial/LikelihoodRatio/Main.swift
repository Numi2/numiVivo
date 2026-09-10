import Foundation

struct Baseline: Decodable {
    let metadata: VivoSingleCellCountMetadata
    let pseudobulk: VivoPseudobulkCounts
    let contrasts: [VivoOmicsExpressionResult]
}
struct Gene: Encodable {
    let feature: VivoOmicsExpressionFeature
    let likelihoodRatio: VivoOmicsNBLikelihoodRatioFit?
    let error: String?
    let fullFitExactlyEqual: Bool
}
struct Measurement: Encodable {
    let request: VivoOmicsExpressionContrast
    let method: String
    let originalFeaturesExactlyEqual: Bool
    let originalDesignExactlyEqual: Bool
    let originalDiagnosticsExactlyEqual: Bool
    let lrtDesignExactlyEqual: Bool
    let lrtTrendExactlyEqual: Bool
    let qualification: String
    let genes: [Gene]
}
@main struct Main {
    static func main() throws {
        guard CommandLine.arguments.count == 2 else { throw VivoOmicsError.invalid("usage: nb-lrt BASELINE_REPORT_GZ") }
        let bytes = try VivoOmicsSourceDecoder.decode(Data(contentsOf: URL(fileURLWithPath: CommandLine.arguments[1])),maximumExpandedBytes: 256*1024*1024)
        let source = try JSONDecoder().decode(Baseline.self,from: bytes)
        guard source.contrasts.count == 1 else { throw VivoOmicsError.invalid("one source contrast required") }
        let baseline = source.contrasts[0]
        let old = try VivoPseudobulkDifferentialExpression.evaluate(metadata: source.metadata,bulk: source.pseudobulk,contrast: baseline.request)
        var request = baseline.request
        var options = request.negativeBinomialOptions ?? .init()
        options.testMethod = .likelihoodRatio; request.negativeBinomialOptions = options
        let lrt = try VivoPseudobulkDifferentialExpression.evaluate(metadata: source.metadata,bulk: source.pseudobulk,contrast: request)
        let genes = lrt.features.map { f in
            let d = lrt.negativeBinomial!.features[f.featureIndex]
            return Gene(feature: f,likelihoodRatio: d.likelihoodRatioFit,error: d.error,
                fullFitExactlyEqual: d.finalFit == old.negativeBinomial!.features[f.featureIndex].finalFit)
        }
        let output = Measurement(request: request,method: lrt.method,
            originalFeaturesExactlyEqual: old.features == baseline.features,
            originalDesignExactlyEqual: old.design == baseline.design,
            originalDiagnosticsExactlyEqual: old.negativeBinomial == baseline.negativeBinomial,
            lrtDesignExactlyEqual: lrt.design == old.design,lrtTrendExactlyEqual: lrt.negativeBinomial!.trend == old.negativeBinomial!.trend,
            qualification: lrt.negativeBinomial!.qualification,genes: genes)
        let encoder = JSONEncoder(); encoder.outputFormatting = [.sortedKeys]
        FileHandle.standardOutput.write(try encoder.encode(output))
    }
}
