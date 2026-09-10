import Foundation
import CryptoKit

struct Source: Decodable {
    let metadata: VivoSingleCellCountMetadata
    let pseudobulk: VivoPseudobulkCounts
}
struct Output: Encodable {
    let sourceReportSHA256: String
    let request: VivoOmicsExpressionContrast
    let result: VivoOmicsExpressionResult?
    let error: String?
    let seconds: Double
}
@main struct Main {
    static func main() throws {
        guard CommandLine.arguments.count==3 else { throw VivoOmicsError.invalid("independent-null REPORT_GZ REQUEST_JSON") }
        let raw=try Data(contentsOf: URL(fileURLWithPath: CommandLine.arguments[1]))
        let identity=SHA256.hash(data: raw).map { String(format: "%02x",$0) }.joined()
        let source=try JSONDecoder().decode(Source.self,from: VivoOmicsSourceDecoder.decode(raw,maximumExpandedBytes: 512*1024*1024))
        let request=try JSONDecoder().decode(VivoOmicsExpressionContrast.self,from: Data(contentsOf: URL(fileURLWithPath: CommandLine.arguments[2])))
        try request.validate()
        let start=Date();let result: VivoOmicsExpressionResult?,error: String?
        do {
            result=try VivoPseudobulkDifferentialExpression.evaluate(metadata: source.metadata,bulk: source.pseudobulk,contrast: request)
            error=nil
        } catch let failure { result=nil;error=String(describing: failure) }
        let out=Output(sourceReportSHA256: identity,request: request,result: result,error: error,seconds: Date().timeIntervalSince(start))
        let encoder=JSONEncoder();encoder.outputFormatting=[.sortedKeys]
        FileHandle.standardOutput.write(try encoder.encode(out))
    }
}
