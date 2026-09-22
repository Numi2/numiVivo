import Foundation
import NumiVivoKit
import NumiVivoLearning

/// Commands for the source-bound MLX response learner. These commands remain
/// separate from the pseudobulk perturbation baseline commands in the
/// single-cell router.
struct VivoCellResponseCLICommands {
    static func handles(_ name: String?) -> Bool {
        [
            "cell-response-help", "cell-response-prepare", "cell-response-verify",
            "cell-response-cohort-admit", "cell-response-cohort-verify",
            "cell-response-train", "cell-response-cohort-train",
            "cell-response-model-verify", "cell-response-resume",
            "cell-response-evaluate", "cell-response-evaluation-verify", "cell-response-evaluation-qualify",
            "cell-response-predict", "cell-response-prediction-verify",
            "cell-response-prediction-verify-bound"
        ].contains(name ?? "")
    }

    private static let help = """
    cell-response-prepare <count-store> --plan <corpus-plan.json> --descriptor-source <target-descriptors.json> --output <new-corpus>
    cell-response-verify <corpus> --store <count-store>
    cell-response-cohort-admit <corpus> --store <count-store> [--requirements <cohort-requirements.json>] --output <new-admission.json>
    cell-response-cohort-verify <admission.json> --corpus <corpus> --store <count-store>
    cell-response-train <corpus> --store <count-store> --plan <training-plan.json> --cohort <admission.json> --output <new-model>
    cell-response-cohort-train <corpus> --store <count-store> --plan <training-plan.json> --cohort <admission.json> --output <new-model>
    cell-response-model-verify <model>
    cell-response-resume <model> --corpus <corpus> --store <count-store> --plan <resume-plan.json> --output <new-model>
    cell-response-evaluate <model> --corpus <corpus> --store <count-store> --partition training|validation|test --max-examples <1...4096> --output <new-evaluation>
    cell-response-evaluation-verify <evaluation> --model <model> --corpus <corpus> --store <count-store>
    cell-response-evaluation-qualify <evaluation> --model <model> --corpus <corpus> --store <count-store>
    cell-response-predict <model> --corpus <corpus> --store <count-store> --plan <prediction-plan.json> --output <new-prediction>
    cell-response-prediction-verify <prediction>
    cell-response-prediction-verify-bound <prediction> --model <model> --corpus <corpus> --store <count-store>

    """

    private func canonicalURL(_ path: String) throws -> URL {
        try VivoWorkflowCLIDocumentPaths.canonicalURL(URL(fileURLWithPath: path))
    }

    private func load<T: Decodable>(_ type: T.Type, _ path: String, maximumBytes: Int = 4_194_304) throws -> T {
        let url = try canonicalURL(path)
        return try VivoCanonicalJSON.decode(type, from: VivoSingleCellCampaignIO.readDocument(url, maximumBytes: maximumBytes))
    }

    private func loadCanonical<T: Codable>(_ type: T.Type, _ path: String,
                                           maximumBytes: Int = 4_194_304) throws -> T {
        let bytes = try VivoSingleCellCampaignIO.readDocument(try canonicalURL(path), maximumBytes: maximumBytes)
        let value = try VivoCanonicalJSON.decode(type, from: bytes)
        guard try VivoCanonicalJSON.encode(value) == bytes else {
            throw VivoOmicsError.invalid("cell-response CLI document is not canonical")
        }
        return value
    }

    private func writeNew<T: Encodable>(_ value: T, _ path: String, maximumBytes: Int) throws {
        let output = try canonicalURL(path)
        guard !FileManager.default.fileExists(atPath: output.path) else {
            throw VivoOmicsError.invalid("cell-response CLI output already exists")
        }
        let bytes = try VivoCanonicalJSON.encode(value)
        guard bytes.count <= maximumBytes else {
            throw VivoOmicsError.limit("cell-response CLI output bytes")
        }
        try bytes.write(to: output, options: .withoutOverwriting)
    }

    private func printJSON<T: Encodable>(_ value: T) throws {
        FileHandle.standardOutput.write(try VivoCanonicalJSON.encode(value))
        FileHandle.standardOutput.write(Data("\n".utf8))
    }

    private func corpus(_ path: String, store: String, implementation: VivoFingerprint) throws -> VivoCellResponseCorpusReader {
        try VivoCellResponseCorpus.open(try canonicalURL(path), sourceStore: try canonicalURL(store), implementation: implementation)
    }

    func run(arguments: [String]) async -> Int32 {
        do {
            guard let command = arguments.first else { throw VivoOmicsError.invalid(Self.help) }
            if command == "cell-response-help" {
                guard arguments.count == 1 else { throw VivoOmicsError.invalid(Self.help) }
                FileHandle.standardOutput.write(Data(Self.help.utf8)); return 0
            }
            let implementation = try VivoWorkflowCLIImplementation.fingerprint()
            switch command {
            case "cell-response-prepare":
                guard arguments.count == 8, arguments[2] == "--plan", arguments[4] == "--descriptor-source", arguments[6] == "--output" else {
                    throw VivoOmicsError.invalid(Self.help)
                }
                let plan = try load(VivoCellResponseCorpusPlan.self, arguments[3], maximumBytes: 67_108_864)
                let descriptorSource = try VivoSingleCellCampaignIO.readDocument(try canonicalURL(arguments[5]), maximumBytes: 67_108_864)
                try printJSON(VivoCellResponseCorpus.prepare(sourceStore: try canonicalURL(arguments[1]), plan: plan,
                                                              descriptorSourceBytes: descriptorSource,
                                                              implementation: implementation, to: try canonicalURL(arguments[7])))
                return 0
            case "cell-response-verify":
                guard arguments.count == 4, arguments[2] == "--store" else { throw VivoOmicsError.invalid(Self.help) }
                try printJSON(VivoCellResponseCorpus.verify(try canonicalURL(arguments[1]), sourceStore: try canonicalURL(arguments[3]), implementation: implementation))
                return 0
            case "cell-response-cohort-admit":
                guard (arguments.count == 6 && arguments[2] == "--store" && arguments[4] == "--output") ||
                      (arguments.count == 8 && arguments[2] == "--store" && arguments[4] == "--requirements" &&
                       arguments[6] == "--output") else {
                    throw VivoOmicsError.invalid(Self.help)
                }
                let reader = try corpus(arguments[1], store: arguments[3], implementation: implementation)
                let requirements: VivoCellResponseCohortRequirements
                let output: String
                if arguments.count == 6 {
                    requirements = .init()
                    output = arguments[5]
                } else {
                    requirements = try loadCanonical(VivoCellResponseCohortRequirements.self, arguments[5])
                    output = arguments[7]
                }
                let admission = try VivoCellResponseCohort.admit(readers: [reader], requirements: requirements)
                try writeNew(admission, output, maximumBytes: 268_435_456)
                try printJSON(admission)
                return 0
            case "cell-response-cohort-verify":
                guard arguments.count == 6, arguments[2] == "--corpus", arguments[4] == "--store" else {
                    throw VivoOmicsError.invalid(Self.help)
                }
                let admission = try loadCanonical(VivoCellResponseCohortAdmission.self, arguments[1], maximumBytes: 268_435_456)
                let reader = try corpus(arguments[3], store: arguments[5], implementation: implementation)
                try VivoCellResponseCohort.verify(admission, readers: [reader])
                try printJSON(admission)
                return 0
            case "cell-response-train", "cell-response-cohort-train":
                guard arguments.count == 10, arguments[2] == "--store", arguments[4] == "--plan",
                      arguments[6] == "--cohort", arguments[8] == "--output" else {
                    throw VivoOmicsError.invalid(Self.help)
                }
                let reader = try corpus(arguments[1], store: arguments[3], implementation: implementation)
                let plan = try load(VivoCellResponseTrainingPlan.self, arguments[5])
                let admission = try loadCanonical(VivoCellResponseCohortAdmission.self, arguments[7], maximumBytes: 268_435_456)
                try printJSON(VivoCellResponseLearning.train(corpus: reader, plan: plan, cohort: admission,
                                                              implementation: implementation,
                                                              to: try canonicalURL(arguments[9])))
                return 0
            case "cell-response-model-verify":
                guard arguments.count == 2 else { throw VivoOmicsError.invalid(Self.help) }
                try printJSON(VivoCellResponseLearning.verifyModel(try canonicalURL(arguments[1]), implementation: implementation))
                return 0
            case "cell-response-resume":
                guard arguments.count == 10, arguments[2] == "--corpus", arguments[4] == "--store",
                      arguments[6] == "--plan", arguments[8] == "--output" else {
                    throw VivoOmicsError.invalid(Self.help)
                }
                let reader = try corpus(arguments[3], store: arguments[5], implementation: implementation)
                let plan = try load(VivoCellResponseResumePlan.self, arguments[7])
                try printJSON(VivoCellResponseLearning.resume(model: try canonicalURL(arguments[1]), corpus: reader,
                                                               plan: plan, implementation: implementation,
                                                               to: try canonicalURL(arguments[9])))
                return 0
            case "cell-response-evaluate":
                guard arguments.count == 12, arguments[2] == "--corpus", arguments[4] == "--store",
                      arguments[6] == "--partition", let partition = VivoCellResponsePartition(rawValue: arguments[7]),
                      arguments[8] == "--max-examples", let maximum = Int(arguments[9]), (1...4_096).contains(maximum),
                      arguments[10] == "--output" else {
                    throw VivoOmicsError.invalid(Self.help)
                }
                let reader = try corpus(arguments[3], store: arguments[5], implementation: implementation)
                try printJSON(VivoCellResponseLearning.evaluate(model: try canonicalURL(arguments[1]), corpus: reader,
                                                                 partition: partition, maximumExamples: maximum,
                                                                 implementation: implementation, to: try canonicalURL(arguments[11])))
                return 0
            case "cell-response-evaluation-verify":
                guard arguments.count == 8, arguments[2] == "--model", arguments[4] == "--corpus",
                      arguments[6] == "--store" else { throw VivoOmicsError.invalid(Self.help) }
                let reader = try corpus(arguments[5], store: arguments[7], implementation: implementation)
                try printJSON(VivoCellResponseLearning.verifyEvaluation(try canonicalURL(arguments[1]),
                                                                         model: try canonicalURL(arguments[3]), corpus: reader,
                                                                         implementation: implementation))
                return 0
            case "cell-response-evaluation-qualify":
                guard arguments.count == 8, arguments[2] == "--model", arguments[4] == "--corpus",
                      arguments[6] == "--store" else { throw VivoOmicsError.invalid(Self.help) }
                let reader = try corpus(arguments[5], store: arguments[7], implementation: implementation)
                try printJSON(VivoCellResponseLearning.qualifyEvaluation(try canonicalURL(arguments[1]),
                                                                          model: try canonicalURL(arguments[3]), corpus: reader,
                                                                          implementation: implementation))
                return 0
            case "cell-response-predict":
                guard arguments.count == 10, arguments[2] == "--corpus", arguments[4] == "--store",
                      arguments[6] == "--plan", arguments[8] == "--output" else {
                    throw VivoOmicsError.invalid(Self.help)
                }
                let reader = try corpus(arguments[3], store: arguments[5], implementation: implementation)
                let plan = try load(VivoCellResponsePredictionPlan.self, arguments[7])
                try printJSON(VivoCellResponseLearning.predict(model: try canonicalURL(arguments[1]), corpus: reader,
                                                                plan: plan, implementation: implementation,
                                                                to: try canonicalURL(arguments[9])))
                return 0
            case "cell-response-prediction-verify":
                guard arguments.count == 2 else { throw VivoOmicsError.invalid(Self.help) }
                try printJSON(VivoCellResponseLearning.verifyPrediction(try canonicalURL(arguments[1]), implementation: implementation))
                return 0
            case "cell-response-prediction-verify-bound":
                guard arguments.count == 8, arguments[2] == "--model", arguments[4] == "--corpus",
                      arguments[6] == "--store" else { throw VivoOmicsError.invalid(Self.help) }
                let reader = try corpus(arguments[5], store: arguments[7], implementation: implementation)
                try printJSON(VivoCellResponseLearning.verifyPrediction(try canonicalURL(arguments[1]),
                                                                        model: try canonicalURL(arguments[3]), corpus: reader,
                                                                        implementation: implementation))
                return 0
            default:
                throw VivoOmicsError.invalid(Self.help)
            }
        } catch {
            // Several native HDF5 and fingerprint errors do not conform to
            // LocalizedError. Preserve their concrete description so a failed
            // source-preflight cannot be mistaken for a corpus failure.
            let detail = String(describing: error)
            let message = detail == "" ? error.localizedDescription : detail
            FileHandle.standardError.write(Data((message + "\n").utf8))
            return 65
        }
    }
}
