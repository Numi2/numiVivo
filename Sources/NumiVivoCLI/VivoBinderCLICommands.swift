import Foundation
import NumiVivoKit

struct VivoBinderCLICommands {
    static func handles(_ command: String?) -> Bool { command?.hasPrefix("binder-") ?? false }
    func run(arguments: [String]) -> Int32 {
        do {
            let command = arguments.first ?? "binder-help"
            if command == "binder-help" {
                print("""
                Retrospective binder scoring (not biological qualification):
                  binder-import SOURCE.csv IMPORT.json NEW_BUNDLE
                  binder-ranking-query IMPORT_BUNDLE NEW_QUERY
                  binder-rank QUERY_BUNDLE RANKING_PLAN.json NEW_RANKING
                  binder-assess-ranking IMPORT_BUNDLE RANKING_BUNDLE NEW_ASSESSMENT
                  binder-evaluate IMPORT_BUNDLE PLAN.json NEW_RESULT
                  binder-evaluate-supported IMPORT_BUNDLE PLAN.json POLICY.json NEW_RESULT
                  binder-evaluate-structures IMPORT_BUNDLE PLAN.json STRUCTURES.json NEW_RESULT
                  binder-evaluate-structure-sources IMPORT_BUNDLE PLAN.json SOURCES.json NEW_RESULT
                  binder-analyze-prediction PREDICTION.json NEW_BUNDLE
                  binder-verify BUNDLE
                Structural evaluation derives heavy-atom geometry and compares with a score-only
                ensemble on identical candidates. No structure repair, MD, affinity or inference.
                Raw-source evaluation reparses retained PDB/mmCIF text and checks declared target
                chain sequences. Source and configuration are retained unchanged; bundles never overwrite.
                Verification reimports source and, for results, repeats the complete fit/evaluation.
                Use the same executable to create/verify a bundle. No network, GPU or paid API.
                Examples and source retrieval: Tools/BinderBenchmark/README.md
                """)
                return 0
            }
            guard let executable = Bundle.main.executableURL else {
                throw VivoBinderBenchmark.Failure.invalid("cannot identify current executable")
            }
            let identity = try VivoCanonicalJSON.fingerprint(Data(contentsOf: executable)).hex
            func path(_ i: Int) -> URL { URL(fileURLWithPath: arguments[i]) }
            var status: Int32 = 0
            let receipt: VivoBinderBundleIO.Receipt
            switch command {
            case "binder-analyze-prediction":
                guard arguments.count == 3 else { return usage() }
                receipt = try VivoBinderBundleIO.analyzePrediction(input: path(1), to: path(2), implementationSHA256: identity)

            case "binder-ranking-query":
                guard arguments.count == 3 else { return usage() }
                receipt = try VivoBinderBundleIO.rankingQuery(bundle: path(1), to: path(2), implementationSHA256: identity)
            case "binder-rank":
                guard arguments.count == 4 else { return usage() }
                receipt = try VivoBinderBundleIO.rank(bundle: path(1), plan: path(2), to: path(3), implementationSHA256: identity)
            case "binder-assess-ranking":
                guard arguments.count == 4 else { return usage() }
                receipt = try VivoBinderBundleIO.assessRanking(bundle: path(1), ranking: path(2), to: path(3), implementationSHA256: identity)
            case "binder-import":
                guard arguments.count == 4 else { return usage() }
                receipt = try VivoBinderBundleIO.importCSV(source: path(1), configuration: path(2),
                    to: path(3), implementationSHA256: identity)
            case "binder-evaluate":
                guard arguments.count == 4 else { return usage() }
                receipt = try VivoBinderBundleIO.evaluate(bundle: path(1), plan: path(2),
                    to: path(3), implementationSHA256: identity)
            case "binder-evaluate-supported":
                guard arguments.count == 5 else { return usage() }
                receipt = try VivoBinderBundleIO.evaluateSupported(bundle: path(1), plan: path(2), policy: path(3),
                    to: path(4), implementationSHA256: identity)
                let report = try VivoCanonicalJSON.decode(VivoBinderTrainingSupport.Evaluation.self,
                    from: Data(contentsOf: path(4).appendingPathComponent("report.json")))
                if !report.support.eligibleForExperimentalFit { status = 2 }
            case "binder-evaluate-structures":
                guard arguments.count == 5 else { return usage() }
                receipt = try VivoBinderBundleIO.evaluateStructures(bundle: path(1), plan: path(2), structures: path(3),
                    to: path(4), implementationSHA256: identity)
            case "binder-evaluate-structure-sources":
                guard arguments.count == 5 else { return usage() }
                receipt = try VivoBinderBundleIO.evaluateStructureSources(bundle: path(1), plan: path(2), sources: path(3),
                    to: path(4), implementationSHA256: identity)
            case "binder-verify":
                guard arguments.count == 2 else { return usage() }
                receipt = try VivoBinderBundleIO.verify(path(1), implementationSHA256: identity)
            default: return usage()
            }
            FileHandle.standardOutput.write(try VivoCanonicalJSON.encode(receipt))
            FileHandle.standardOutput.write(Data("\n".utf8))
            return status
        } catch {
            FileHandle.standardError.write(Data("Binder workflow rejected: \(error)\n".utf8)); return 65
        }
    }
    private func usage() -> Int32 {
        FileHandle.standardError.write(Data("Use binder-help for the exact syntax.\n".utf8)); return 64
    }
}
