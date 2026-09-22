import Foundation
import NumiVivoKit
import NumiVivoLearning

/// File-system locations are an invocation input, never training evidence.
/// The model records only the receipt-derived composite identity resolved from
/// this manifest, so moving a store does not change model provenance.
private struct VivoCellResponseCLISource: Codable, Hashable {
    let corpus: String
    let store: String

    private enum CodingKeys: String, CodingKey { case corpus, store }

    init(corpus: String, store: String) {
        self.corpus = corpus
        self.store = store
    }

    init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["corpus", "store"])
        let values = try decoder.container(keyedBy: CodingKeys.self)
        corpus = try values.decode(String.self, forKey: .corpus)
        store = try values.decode(String.self, forKey: .store)
    }
}

/// Bounded input for a source-qualified response-learning run. The reader
/// resolves and receipt-sorts entries before training, so this document's path
/// order cannot change row IDs, sampler order, or checkpoint provenance.
private struct VivoCellResponseCLISources: Codable, Equatable {
    let schemaVersion: Int
    let format: String
    let sources: [VivoCellResponseCLISource]

    private enum CodingKeys: String, CodingKey { case schemaVersion, format, sources }

    init(sources: [VivoCellResponseCLISource]) {
        schemaVersion = 1
        format = "numivivo-cell-response-cli-sources/v1"
        self.sources = sources
    }

    init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["schemaVersion", "format", "sources"])
        let values = try decoder.container(keyedBy: CodingKeys.self)
        schemaVersion = try values.decode(Int.self, forKey: .schemaVersion)
        format = try values.decode(String.self, forKey: .format)
        sources = try values.decode([VivoCellResponseCLISource].self, forKey: .sources)
    }

    func validate() throws {
        guard schemaVersion == 1, format == "numivivo-cell-response-cli-sources/v1",
              (1...64).contains(sources.count),
              sources.allSatisfy({ source in
                  source.corpus.hasPrefix("/") && source.corpus.utf8.count <= 8_192 &&
                  source.corpus.unicodeScalars.allSatisfy({ $0.value >= 0x20 && $0.value != 0x7f }) &&
                  source.store.hasPrefix("/") && source.store.utf8.count <= 8_192 &&
                  source.store.unicodeScalars.allSatisfy({ $0.value >= 0x20 && $0.value != 0x7f })
              }),
              Set(sources).count == sources.count else {
            throw VivoOmicsError.invalid("cell-response source manifest")
        }
    }
}

/// Commands for the source-bound MLX response learner. These commands remain
/// separate from the pseudobulk perturbation baseline commands in the
/// single-cell router.
struct VivoCellResponseCLICommands {
    static func handles(_ name: String?) -> Bool {
        [
            "cell-response-help", "cell-response-prepare", "cell-response-verify",
            "cell-response-sources-verify",
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
    cell-response-sources-verify <corpus-sources.json>
    cell-response-cohort-admit <corpus> --store <count-store> [--requirements <cohort-requirements.json>] --output <new-admission.json>
    cell-response-cohort-admit --sources <corpus-sources.json> [--requirements <cohort-requirements.json>] --output <new-admission.json>
    cell-response-cohort-verify <admission.json> --corpus <corpus> --store <count-store>
    cell-response-cohort-verify <admission.json> --sources <corpus-sources.json>
    cell-response-train <corpus> --store <count-store> --plan <training-plan.json> --cohort <admission.json> --output <new-model>
    cell-response-train --sources <corpus-sources.json> --plan <training-plan.json> --cohort <admission.json> --output <new-model>
    cell-response-cohort-train <corpus> --store <count-store> --plan <training-plan.json> --cohort <admission.json> --output <new-model>
    cell-response-cohort-train --sources <corpus-sources.json> --plan <training-plan.json> --cohort <admission.json> --output <new-model>
    cell-response-model-verify <model>
    cell-response-resume <model> --corpus <corpus> --store <count-store> --plan <resume-plan.json> --output <new-model>
    cell-response-resume <model> --sources <corpus-sources.json> --plan <resume-plan.json> --output <new-model>
    cell-response-evaluate <model> --corpus <corpus> --store <count-store> --partition training|validation|test --max-examples <1...4096> --output <new-evaluation>
    cell-response-evaluate <model> --sources <corpus-sources.json> --partition training|validation|test --max-examples <1...4096> --output <new-evaluation>
    cell-response-evaluation-verify <evaluation> --model <model> --corpus <corpus> --store <count-store>
    cell-response-evaluation-verify <evaluation> --model <model> --sources <corpus-sources.json>
    cell-response-evaluation-qualify <evaluation> --model <model> --corpus <corpus> --store <count-store>
    cell-response-evaluation-qualify <evaluation> --model <model> --sources <corpus-sources.json>
    cell-response-predict <model> --corpus <corpus> --store <count-store> --plan <prediction-plan.json> --output <new-prediction>
    cell-response-predict <model> --sources <corpus-sources.json> --plan <prediction-plan.json> --output <new-prediction>
    cell-response-prediction-verify <prediction>
    cell-response-prediction-verify-bound <prediction> --model <model> --corpus <corpus> --store <count-store>
    cell-response-prediction-verify-bound <prediction> --model <model> --sources <corpus-sources.json>

    Source manifests use absolute corpus/store paths. `cell-response-sources-verify`
    prints the member corpus fingerprints required by a multi-source prediction
    plan's schemaVersion 2 `contextCorpus` field.

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

    private func composite(_ path: String, implementation: VivoFingerprint) throws -> VivoCellResponseCompositeCorpusReader {
        let manifest = try load(VivoCellResponseCLISources.self, path, maximumBytes: 4_194_304)
        try manifest.validate()
        let readers = try manifest.sources.map { source in
            try corpus(source.corpus, store: source.store, implementation: implementation)
        }
        return try VivoCellResponseCompositeCorpusReader(readers: readers)
    }

    private func singleComposite(_ path: String, store: String,
                                 implementation: VivoFingerprint) throws -> VivoCellResponseCompositeCorpusReader {
        try VivoCellResponseCompositeCorpusReader(readers: [try corpus(path, store: store, implementation: implementation)])
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
            case "cell-response-sources-verify":
                guard arguments.count == 2 else { throw VivoOmicsError.invalid(Self.help) }
                let reader = try composite(arguments[1], implementation: implementation)
                try printJSON(reader.identity)
                return 0
            case "cell-response-cohort-admit":
                let reader: VivoCellResponseCompositeCorpusReader
                let requirements: VivoCellResponseCohortRequirements
                let output: String
                if arguments.count == 6, arguments[2] == "--store", arguments[4] == "--output" {
                    reader = try singleComposite(arguments[1], store: arguments[3], implementation: implementation)
                    requirements = .init()
                    output = arguments[5]
                } else if arguments.count == 8, arguments[2] == "--store", arguments[4] == "--requirements",
                          arguments[6] == "--output" {
                    reader = try singleComposite(arguments[1], store: arguments[3], implementation: implementation)
                    requirements = try loadCanonical(VivoCellResponseCohortRequirements.self, arguments[5])
                    output = arguments[7]
                } else if arguments.count == 5, arguments[1] == "--sources", arguments[3] == "--output" {
                    reader = try composite(arguments[2], implementation: implementation)
                    requirements = .init()
                    output = arguments[4]
                } else if arguments.count == 7, arguments[1] == "--sources", arguments[3] == "--requirements",
                          arguments[5] == "--output" {
                    reader = try composite(arguments[2], implementation: implementation)
                    requirements = try loadCanonical(VivoCellResponseCohortRequirements.self, arguments[4])
                    output = arguments[6]
                } else {
                    throw VivoOmicsError.invalid(Self.help)
                }
                let admission = try VivoCellResponseCohort.admit(readers: reader.sourceReaders, requirements: requirements)
                try writeNew(admission, output, maximumBytes: 268_435_456)
                try printJSON(admission)
                return 0
            case "cell-response-cohort-verify":
                guard arguments.count >= 2 else { throw VivoOmicsError.invalid(Self.help) }
                let admission = try loadCanonical(VivoCellResponseCohortAdmission.self, arguments[1], maximumBytes: 268_435_456)
                let reader: VivoCellResponseCompositeCorpusReader
                if arguments.count == 6, arguments[2] == "--corpus", arguments[4] == "--store" {
                    reader = try singleComposite(arguments[3], store: arguments[5], implementation: implementation)
                } else if arguments.count == 4, arguments[2] == "--sources" {
                    reader = try composite(arguments[3], implementation: implementation)
                } else {
                    throw VivoOmicsError.invalid(Self.help)
                }
                try VivoCellResponseCohort.verify(admission, readers: reader.sourceReaders)
                try printJSON(admission)
                return 0
            case "cell-response-train", "cell-response-cohort-train":
                let reader: VivoCellResponseCompositeCorpusReader
                let plan: VivoCellResponseTrainingPlan
                let admission: VivoCellResponseCohortAdmission
                let output: String
                if arguments.count == 10, arguments[2] == "--store", arguments[4] == "--plan",
                   arguments[6] == "--cohort", arguments[8] == "--output" {
                    reader = try singleComposite(arguments[1], store: arguments[3], implementation: implementation)
                    plan = try load(VivoCellResponseTrainingPlan.self, arguments[5])
                    admission = try loadCanonical(VivoCellResponseCohortAdmission.self, arguments[7], maximumBytes: 268_435_456)
                    output = arguments[9]
                } else if arguments.count == 9, arguments[1] == "--sources", arguments[3] == "--plan",
                          arguments[5] == "--cohort", arguments[7] == "--output" {
                    reader = try composite(arguments[2], implementation: implementation)
                    plan = try load(VivoCellResponseTrainingPlan.self, arguments[4])
                    admission = try loadCanonical(VivoCellResponseCohortAdmission.self, arguments[6], maximumBytes: 268_435_456)
                    output = arguments[8]
                } else {
                    throw VivoOmicsError.invalid(Self.help)
                }
                try printJSON(VivoCellResponseLearning.train(corpus: reader, plan: plan, cohort: admission,
                                                              implementation: implementation,
                                                              to: try canonicalURL(output)))
                return 0
            case "cell-response-model-verify":
                guard arguments.count == 2 else { throw VivoOmicsError.invalid(Self.help) }
                try printJSON(VivoCellResponseLearning.verifyModel(try canonicalURL(arguments[1]), implementation: implementation))
                return 0
            case "cell-response-resume":
                let reader: VivoCellResponseCompositeCorpusReader
                let plan: VivoCellResponseResumePlan
                let output: String
                if arguments.count == 10, arguments[2] == "--corpus", arguments[4] == "--store",
                   arguments[6] == "--plan", arguments[8] == "--output" {
                    reader = try singleComposite(arguments[3], store: arguments[5], implementation: implementation)
                    plan = try load(VivoCellResponseResumePlan.self, arguments[7])
                    output = arguments[9]
                } else if arguments.count == 8, arguments[2] == "--sources", arguments[4] == "--plan",
                          arguments[6] == "--output" {
                    reader = try composite(arguments[3], implementation: implementation)
                    plan = try load(VivoCellResponseResumePlan.self, arguments[5])
                    output = arguments[7]
                } else {
                    throw VivoOmicsError.invalid(Self.help)
                }
                try printJSON(VivoCellResponseLearning.resume(model: try canonicalURL(arguments[1]), corpus: reader,
                                                               plan: plan, implementation: implementation,
                                                               to: try canonicalURL(output)))
                return 0
            case "cell-response-evaluate":
                let reader: VivoCellResponseCompositeCorpusReader
                let partition: VivoCellResponsePartition
                let maximum: Int
                let output: String
                if arguments.count == 12, arguments[2] == "--corpus", arguments[4] == "--store",
                   arguments[6] == "--partition", let parsed = VivoCellResponsePartition(rawValue: arguments[7]),
                   arguments[8] == "--max-examples", let count = Int(arguments[9]), (1...4_096).contains(count),
                   arguments[10] == "--output" {
                    reader = try singleComposite(arguments[3], store: arguments[5], implementation: implementation)
                    partition = parsed; maximum = count; output = arguments[11]
                } else if arguments.count == 10, arguments[2] == "--sources", arguments[4] == "--partition",
                          let parsed = VivoCellResponsePartition(rawValue: arguments[5]),
                          arguments[6] == "--max-examples", let count = Int(arguments[7]), (1...4_096).contains(count),
                          arguments[8] == "--output" {
                    reader = try composite(arguments[3], implementation: implementation)
                    partition = parsed; maximum = count; output = arguments[9]
                } else {
                    throw VivoOmicsError.invalid(Self.help)
                }
                try printJSON(VivoCellResponseLearning.evaluate(model: try canonicalURL(arguments[1]), corpus: reader,
                                                                 partition: partition, maximumExamples: maximum,
                                                                 implementation: implementation, to: try canonicalURL(output)))
                return 0
            case "cell-response-evaluation-verify":
                let reader: VivoCellResponseCompositeCorpusReader
                if arguments.count == 8, arguments[2] == "--model", arguments[4] == "--corpus",
                   arguments[6] == "--store" {
                    reader = try singleComposite(arguments[5], store: arguments[7], implementation: implementation)
                } else if arguments.count == 6, arguments[2] == "--model", arguments[4] == "--sources" {
                    reader = try composite(arguments[5], implementation: implementation)
                } else {
                    throw VivoOmicsError.invalid(Self.help)
                }
                try printJSON(VivoCellResponseLearning.verifyEvaluation(try canonicalURL(arguments[1]),
                                                                         model: try canonicalURL(arguments[3]), corpus: reader,
                                                                         implementation: implementation))
                return 0
            case "cell-response-evaluation-qualify":
                let reader: VivoCellResponseCompositeCorpusReader
                if arguments.count == 8, arguments[2] == "--model", arguments[4] == "--corpus",
                   arguments[6] == "--store" {
                    reader = try singleComposite(arguments[5], store: arguments[7], implementation: implementation)
                } else if arguments.count == 6, arguments[2] == "--model", arguments[4] == "--sources" {
                    reader = try composite(arguments[5], implementation: implementation)
                } else {
                    throw VivoOmicsError.invalid(Self.help)
                }
                try printJSON(VivoCellResponseLearning.qualifyEvaluation(try canonicalURL(arguments[1]),
                                                                          model: try canonicalURL(arguments[3]), corpus: reader,
                                                                          implementation: implementation))
                return 0
            case "cell-response-predict":
                let reader: VivoCellResponseCompositeCorpusReader
                let plan: VivoCellResponsePredictionPlan
                let output: String
                if arguments.count == 10, arguments[2] == "--corpus", arguments[4] == "--store",
                   arguments[6] == "--plan", arguments[8] == "--output" {
                    reader = try singleComposite(arguments[3], store: arguments[5], implementation: implementation)
                    plan = try load(VivoCellResponsePredictionPlan.self, arguments[7])
                    output = arguments[9]
                } else if arguments.count == 8, arguments[2] == "--sources", arguments[4] == "--plan",
                          arguments[6] == "--output" {
                    reader = try composite(arguments[3], implementation: implementation)
                    plan = try load(VivoCellResponsePredictionPlan.self, arguments[5])
                    output = arguments[7]
                } else {
                    throw VivoOmicsError.invalid(Self.help)
                }
                try printJSON(VivoCellResponseLearning.predict(model: try canonicalURL(arguments[1]), corpus: reader,
                                                                plan: plan, implementation: implementation,
                                                                to: try canonicalURL(output)))
                return 0
            case "cell-response-prediction-verify":
                guard arguments.count == 2 else { throw VivoOmicsError.invalid(Self.help) }
                try printJSON(VivoCellResponseLearning.verifyPrediction(try canonicalURL(arguments[1]), implementation: implementation))
                return 0
            case "cell-response-prediction-verify-bound":
                let reader: VivoCellResponseCompositeCorpusReader
                if arguments.count == 8, arguments[2] == "--model", arguments[4] == "--corpus",
                   arguments[6] == "--store" {
                    reader = try singleComposite(arguments[5], store: arguments[7], implementation: implementation)
                } else if arguments.count == 6, arguments[2] == "--model", arguments[4] == "--sources" {
                    reader = try composite(arguments[5], implementation: implementation)
                } else {
                    throw VivoOmicsError.invalid(Self.help)
                }
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
