import Foundation
import NumiVivoKit

struct VivoMDBenchmarkCLI {
    static func handles(_ name:String?) -> Bool { ["md-benchmark","md-benchmark-help"].contains(name ?? "") }
    func run(arguments:[String]) async -> Int32 {
        do {
            if arguments == ["md-benchmark-help"] {
                print("numivivo md-benchmark <request.json> --output <new-report.json>\nCompares native forces and energies with supplied reference observations.\nShort dynamics are execution evidence; ensemble qualification remains inconclusive.")
                return 0
            }
            guard arguments.count==4,arguments[0]=="md-benchmark",arguments[2]=="--output" else {
                throw VivoChemistryError.invalid("usage: md-benchmark <request.json> --output <new-report.json>")
            }
            let input=URL(fileURLWithPath:arguments[1]),output=URL(fileURLWithPath:arguments[3])
            guard try !VivoWorkflowCLIDocumentPaths.aliases(input,output),!FileManager.default.fileExists(atPath:output.path) else {
                throw VivoChemistryError.invalid("benchmark output exists or aliases input")
            }
            let destination=try VivoKineticsDocumentIO.prepareNoClobberOutput(to:output)
            let request=try VivoKineticsDocumentIO.read(VivoMDBenchmarkRequest.self,from:input,maximumBytes:536_870_912)
            let report=try await VivoMDBenchmark.run(request)
            try Task.checkCancellation();try destination.write(VivoCanonicalJSON.encode(report))
            return report.outcome == .passed ? 0:75
        } catch is CancellationError { return 130 }
        catch { FileHandle.standardError.write(Data("numivivo md-benchmark: \(error)\n".utf8));return 65 }
    }
}
