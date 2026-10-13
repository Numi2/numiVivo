import Foundation
import NumiVivoKit

// Dated experiment runner: the native prediction owner receives only training
// labels and typed outcome-free query fields. It never opens the source CSV.
let args = Array(CommandLine.arguments.dropFirst())
if args.first != "heldout-predict" {
    exit(VivoBinderCLICommands().run(arguments: args))
}
do {
    guard args.count == 3 else { throw VivoBinderBenchmark.Failure.invalid("heldout-predict INPUT_DIRECTORY NEW_OUTPUT_DIRECTORY") }
    let input = URL(fileURLWithPath: args[1]), output = URL(fileURLWithPath: args[2])
    guard !FileManager.default.fileExists(atPath: output.path) else { throw VivoBinderBenchmark.Failure.invalid("output exists") }
    func read<T: Decodable>(_ name: String, _ type: T.Type) throws -> T {
        try JSONDecoder().decode(type, from: Data(contentsOf: input.appendingPathComponent(name)))
    }
    let query = try read("query.json", [VivoBinderBenchmark.QueryCandidate].self)
    let plan = try read("plan.json", VivoBinderBenchmark.Plan.self)
    let rankingQuery = VivoBinderRanking.Query(sourceSHA256: plan.sourceSHA256,
        candidates: query.map { .init(id: $0.id, target: $0.target, scores: $0.features) })
    var reports: [String: Data] = [:]
    for assay in ["adaptyv", "twist"] {
        let training = try read("training-\(assay).json", VivoBinderBenchmark.Dataset.self)
        let prediction = try VivoBinderBenchmark.predict(training: training, query: query, plan: plan)
        reports["prediction-\(assay).json"] = try VivoCanonicalJSON.encode(prediction)
    }
    for method in ["boltz2", "confidence3"] {
        let fixedPlan = try read("ranking-\(method).json", VivoBinderRanking.Plan.self)
        reports["ranking-\(method).json"] = try VivoCanonicalJSON.encode(VivoBinderRanking.rank(rankingQuery, plan: fixedPlan))
    }
    let temporary = output.deletingLastPathComponent().appendingPathComponent(".heldout-" + UUID().uuidString)
    try FileManager.default.createDirectory(at: temporary, withIntermediateDirectories: false)
    defer { try? FileManager.default.removeItem(at: temporary) }
    for (name, data) in reports { try data.write(to: temporary.appendingPathComponent(name), options: .withoutOverwriting) }
    try FileManager.default.moveItem(at: temporary, to: output)
    print("Native predictions and fixed rankings completed; no held-out outcomes read.")
} catch {
    FileHandle.standardError.write(Data("Held-out prediction failed: \(error)\n".utf8))
    exit(65)
}
