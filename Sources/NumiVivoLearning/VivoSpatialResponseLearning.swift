import Foundation
import MLX
import MLXNN
import MLXOptimizers

/// Spatial conditioning of the existing cell-response network. This is an
/// experimental endpoint-distribution learner, not a cell trajectory model.
public enum VivoSpatialResponseLearning {
    struct Plan: Codable {
        let featureCount: Int
        let targetCount: Int
        let descriptorCount: Int
        let hiddenWidth: Int
        let seed: UInt64
        let steps: [Int]
        let batchSize: Int
        let learningRate: Float
        let weightDecay: Float
    }
    static func read(_ url: URL) throws -> [String: MLXArray] {
        try MLX.loadArraysAndMetadata(data: Data(contentsOf: url)).0
    }
    static func write(_ arrays: [String: MLXArray], _ url: URL) throws {
        try MLX.saveToData(arrays: arrays).write(to: url, options: .withoutOverwriting)
    }
    static func output(_ model: VivoCellResponseMLXModel, _ a: [String: MLXArray], _ rows: MLXArray) -> (mean: MLXArray, variance: MLXArray) {
        let c = a["context"]![rows]
        let o = model.outputs(context: c.reshaped([c.shape[0], 1, c.shape[1]]),
                              targetIDs: a["target"]![rows], knownTargetMask: a["known"]![rows],
                              descriptors: a["descriptor"]![rows])
        return (o.mean, o.variance)
    }
    static func validate(_ a: [String: MLXArray], _ p: Plan, observations: Bool) throws {
        func check(_ test: Bool, _ message: String) throws {
            if !test { throw NSError(domain: "spatial-response", code: 1, userInfo: [NSLocalizedDescriptionKey: message]) }
        }
        guard let c = a["context"], c.ndim == 2 else { throw NSError(domain: "spatial-response", code: 2) }
        let n = c.shape[0]
        try check(n > 0 && c.shape[1] == p.featureCount, "context shape")
        try check(a["descriptor"]?.shape == [n, p.descriptorCount], "descriptor shape")
        try check(a["target"]?.shape == [n] && a["known"]?.shape == [n, 1], "target shape")
        try check(a["target"]!.asArray(Int32.self).allSatisfy { $0 >= 0 && $0 < p.targetCount }, "target range")
        for key in ["context", "descriptor", "known"] { try check(a[key]!.asArray(Float.self).allSatisfy(\.isFinite), "nonfinite input") }
        if observations {
            for key in ["observed", "observedVariance", "mask"] {
                try check(a[key]?.shape == [n, p.featureCount], "observation shape")
                try check(a[key]!.asArray(Float.self).allSatisfy(\.isFinite), "nonfinite observation")
            }
            try check(a["mask"]!.sum().item(Float.self) > 0, "empty measured feature mask")
        }
    }
    /// Arguments: train|predict plan.json input.safetensors output-directory [weights.safetensors]
    public static func run(_ arguments: [String]) throws {
        guard arguments.count >= 4 else { throw NSError(domain: "spatial-response: train|predict plan inputs output [weights]", code: 1) }
        let mode = arguments[0], planURL = URL(fileURLWithPath: arguments[1])
        let p = try JSONDecoder().decode(Plan.self, from: Data(contentsOf: planURL))
        guard p.featureCount > 0, p.featureCount <= 50000, p.targetCount > 0, p.descriptorCount > 0,
              p.descriptorCount <= 2048, p.hiddenWidth == 64, p.batchSize > 0, p.batchSize <= 64,
              p.steps == [240,720,1440], p.learningRate > 0, p.learningRate.isFinite,
              p.weightDecay >= 0, p.weightDecay.isFinite else { throw NSError(domain: "invalid spatial learning plan", code: 1) }
        let a = try read(URL(fileURLWithPath: arguments[2]))
        try validate(a, p, observations: mode == "train")
        let root = URL(fileURLWithPath: arguments[3], isDirectory: true)
        guard !FileManager.default.fileExists(atPath: root.path) else { throw NSError(domain: "output exists", code: 1) }
        try FileManager.default.createDirectory(at: root, withIntermediateDirectories: false)
        MLXRandom.seed(p.seed)
        let model = VivoCellResponseMLXModel(featureCount: p.featureCount, targetCount: p.targetCount,
             descriptorCount: p.descriptorCount, architecture: .init(hiddenWidth: p.hiddenWidth, contextCells: 1),
             targetResponsePriorValues: a["prior"])
        if mode == "train" {
            guard let prior = a["prior"], prior.shape == [p.targetCount, p.featureCount],
                  prior.asArray(Float.self).allSatisfy(\.isFinite) else { throw NSError(domain: "invalid prior", code: 1) }
            let opt = MultiOptimizer(optimizers: [
                SGD(learningRate: p.learningRate * Float(p.featureCount), momentum: 0, weightDecay: p.weightDecay / Float(p.featureCount)),
                SGD(learningRate: p.learningRate, momentum: 0, weightDecay: p.weightDecay)
            ], filters: [{ key, _ in key == "meanHead.weight" || key == "meanHead.bias" }])
            let grad = valueAndGrad(model: model) { m, arrays in
                let idx = arrays[0], o = output(m, a, idx), mask = a["mask"]![idx]
                let residual = o.mean - a["observed"]![idx]
                let loss = (log(o.variance) + (residual * residual + a["observedVariance"]![idx]) / o.variance) * mask
                return [0.5 * loss.sum() / maximum(mask.sum(), MLXArray(Float(1)))]
            }
            var state = p.seed
            func next() -> UInt64 { state &+= 0x9e3779b97f4a7c15; var z = state; z = (z ^ (z >> 30)) &* 0xbf58476d1ce4e5b9; z = (z ^ (z >> 27)) &* 0x94d049bb133111eb; return z ^ (z >> 31) }
            // Balance target and response role strata, not the number of cells.
            let strata = a["stratum"]!.asArray(Int32.self)
            let groups = Dictionary(grouping: strata.indices, by: { strata[$0] }).sorted { $0.key < $1.key }.map(\.value)
            for step in 1...p.steps.last! {
                let indices = (0..<p.batchSize).map { _ -> Int32 in let g = groups[Int(next() % UInt64(groups.count))]; return Int32(g[Int(next() % UInt64(g.count))]) }
                let (loss, grads) = grad(model, [MLXArray(indices)])
                opt.update(model: model, gradients: grads); eval(model, opt)
                guard loss[0].item(Float.self).isFinite else { throw NSError(domain: "nonfinite spatial loss", code: step) }
                if p.steps.contains(step) {
                    try write(Dictionary(uniqueKeysWithValues: model.parameters().flattened()), root.appendingPathComponent("weights-\(step).safetensors"))
                    print("step=\(step) nll=\(loss[0].item(Float.self))")
                    fflush(stdout)
                }
            }
        } else if mode == "predict", arguments.count == 5 {
            try model.update(parameters: ModuleParameters.unflattened(try read(URL(fileURLWithPath: arguments[4]))), verify: .all)
            model.train(false)
            var means: [MLXArray] = [], variances: [MLXArray] = []
            for start in stride(from: 0, to: a["context"]!.shape[0], by: p.batchSize) {
                let end = min(start + p.batchSize, a["context"]!.shape[0])
                let o = output(model, a, MLXArray((start..<end).map(Int32.init)))
                eval(o.mean, o.variance); means.append(o.mean); variances.append(o.variance)
            }
            try write(["mean": concatenated(means, axis: 0), "variance": concatenated(variances, axis: 0)], root.appendingPathComponent("prediction.safetensors"))
        } else { throw NSError(domain: "unknown spatial command", code: 1) }
    }
}
