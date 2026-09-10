import Foundation
// Exercise the same internal pseudobulk owner used by the validated public
// dataset and H5AD routes, using the already-qualified source aggregation.
@testable import NumiVivoKit
struct Baseline: Decodable {
    let metadata: VivoSingleCellCountMetadata
    let pseudobulk: VivoPseudobulkCounts
    let contrasts: [VivoOmicsExpressionResult]
}
struct ProductOutput: Encodable {
    let sourceDesignExactlyEqual: Bool
    let sourceTrendExactlyEqual: Bool
    let sourceEvidenceExactlyEqual: Bool
    let expression: VivoOmicsExpressionResult
    let seconds: Double
}
@main struct Product {
    static func main() throws {
        guard CommandLine.arguments.count == 2 else { throw VivoOmicsError.invalid("usage: ql-product BASELINE_REPORT_GZ") }
        let data = try VivoOmicsSourceDecoder.decode(Data(contentsOf: URL(fileURLWithPath: CommandLine.arguments[1])),maximumExpandedBytes: 256*1024*1024)
        let source = try JSONDecoder().decode(Baseline.self,from: data)
        guard source.contrasts.count == 1 else { throw VivoOmicsError.invalid("one source contrast required") }
        let baseline = source.contrasts[0]
        var request = baseline.request, options = request.negativeBinomialOptions ?? .init()
        // This harness explicitly requests unshrunk adjusted QL inference.
        // The product rejects unsupported priors instead of silently ignoring them.
        options.testMethod = .quasiLikelihoodAdjusted
        options.effectPriorEstimation = nil; options.effectPriorStandardDeviationLog2 = nil
        request.negativeBinomialOptions = options
        let start = Date()
        let expression = try VivoPseudobulkDifferentialExpression.evaluate(metadata: source.metadata,bulk: source.pseudobulk,contrast: request)
        let output = ProductOutput(sourceDesignExactlyEqual: expression.design == baseline.design,
            sourceTrendExactlyEqual: expression.negativeBinomial?.trend == baseline.negativeBinomial?.trend,
            sourceEvidenceExactlyEqual: expression.evidence == baseline.evidence,expression: expression,
            seconds: Date().timeIntervalSince(start))
        let encoder = JSONEncoder(); encoder.outputFormatting = [.sortedKeys]
        FileHandle.standardOutput.write(try encoder.encode(output))
    }
}
