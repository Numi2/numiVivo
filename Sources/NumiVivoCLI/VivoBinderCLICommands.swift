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
                  binder-evaluate IMPORT_BUNDLE PLAN.json NEW_RESULT
                  binder-verify BUNDLE
                Source and configuration are retained unchanged; bundles never overwrite.
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
            let receipt: VivoBinderBundleIO.Receipt
            switch command {
            case "binder-import":
                guard arguments.count == 4 else { return usage() }
                receipt = try VivoBinderBundleIO.importCSV(source: path(1), configuration: path(2),
                    to: path(3), implementationSHA256: identity)
            case "binder-evaluate":
                guard arguments.count == 4 else { return usage() }
                receipt = try VivoBinderBundleIO.evaluate(bundle: path(1), plan: path(2),
                    to: path(3), implementationSHA256: identity)
            case "binder-verify":
                guard arguments.count == 2 else { return usage() }
                receipt = try VivoBinderBundleIO.verify(path(1), implementationSHA256: identity)
            default: return usage()
            }
            FileHandle.standardOutput.write(try VivoCanonicalJSON.encode(receipt))
            FileHandle.standardOutput.write(Data("\n".utf8))
            return 0
        } catch {
            FileHandle.standardError.write(Data("Binder workflow rejected: \(error)\n".utf8)); return 65
        }
    }
    private func usage() -> Int32 {
        FileHandle.standardError.write(Data("Use binder-help for the exact syntax.\n".utf8)); return 64
    }
}
