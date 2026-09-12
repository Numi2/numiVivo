import Foundation
import NumiVivoKit

/// Native input/submission boundary only. Arc's official Python scorer runs
/// separately under Tools/Omics/Arc2026; it never enters production inference.
struct VivoArc2026CLICommands {
    static func handles(_ command: String?) -> Bool { command?.hasPrefix("arc2026-") ?? false }
    func run(arguments: [String]) -> Int32 {
        do {
            let command = arguments.first ?? "arc2026-help"
            func emit<T: Encodable>(_ value: T) throws {
                FileHandle.standardOutput.write(try VivoCanonicalJSON.encode(value))
                FileHandle.standardOutput.write(Data("\n".utf8))
            }
            func local(_ value: String) -> URL { URL(fileURLWithPath: value) }
            switch command {
            case "arc2026-help":
                print("""
                Arc VCC 2026 native interface (input conformance, not biological qualification):
                  arc2026-prepare PLAN.json --output NEW_QUERY_DIR
                  arc2026-pack QUERY_DIR --plan PLAN.json --output NEW_SUBMISSION_DIR
                  arc2026-verify-query QUERY_DIR
                  arc2026-verify SUBMISSION_DIR
                Scoring: Tools/Omics/Arc2026/evaluate.py, pinned cell-eval2 0.16.0/rule 3.
                All three contexts, all 300 targets and the ordered 18,533 genes are required.
                Raw integer counts only. Predicted cell counts need not equal the reference.
                No reference outcomes, score anchors, model fitting or upload in these commands.
                """)
            case "arc2026-prepare":
                guard arguments.count == 4, arguments[2] == "--output" else { return usage() }
                let plan = try VivoArc2026H5AD.readPlan(VivoArc2026PreparePlan.self, at: local(arguments[1]))
                try emit(VivoArc2026H5AD.prepare(plan, to: local(arguments[3])))
            case "arc2026-pack":
                guard arguments.count == 6, arguments[2] == "--plan", arguments[4] == "--output" else { return usage() }
                let plan = try VivoArc2026H5AD.readPlan(VivoArc2026PackPlan.self, at: local(arguments[3]))
                try emit(VivoArc2026H5AD.pack(query: local(arguments[1]), plan: plan, to: local(arguments[5])))
            case "arc2026-verify-query":
                guard arguments.count == 2 else { return usage() }
                try emit(VivoArc2026H5AD.verifyQuery(local(arguments[1])))
            case "arc2026-verify":
                guard arguments.count == 2 else { return usage() }
                try emit(VivoArc2026H5AD.verifySubmission(local(arguments[1])))
            default: return usage()
            }
            return 0
        } catch {
            FileHandle.standardError.write(Data("Arc 2026 rejected: \(error.localizedDescription)\n".utf8))
            return 65
        }
    }
    private func usage() -> Int32 {
        FileHandle.standardError.write(Data("Use arc2026-help for the exact native command syntax.\n".utf8)); return 64
    }
}
