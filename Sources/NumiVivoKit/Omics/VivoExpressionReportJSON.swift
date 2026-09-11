import Foundation
import CryptoKit

/// Canonical expression-report output without a whole-report JSONEncoder tree or
/// Data buffer. The model still owns its fit arrays. Each feature uses the same
/// canonical Foundation encoder; envelopes enumerate the versioned report schema.
/// Exact-byte regression checks against synthesized Codable protect every field.
enum VivoExpressionReportJSON {
    static func fingerprint(_ result: VivoOmicsExpressionResult, maximumBytes: Int,
                            sink: ((Data) throws -> Void)? = nil) throws -> VivoFingerprint {
        let output = Writer(maximumBytes: maximumBytes, sink: sink)
        try output.object { o in
            o.field("design", result.design); o.field("evidence", result.evidence)
            o.array("features", result.features); o.field("method", result.method)
            o.field("multiplicityScope", result.multiplicityScope)
            if let nb = result.negativeBinomial { o.nested("negativeBinomial") { try output.negativeBinomial(nb) } }
            o.field("request", result.request); o.field("testedFeatures", result.testedFeatures)
            o.optional("variancePrior", result.variancePrior)
        }
        return try output.finish()
    }
    private final class Object {
        let output: Writer
        var fields: [String: () throws -> Void] = [:]
        init(_ output: Writer) { self.output = output }
        func field<T: Encodable>(_ key: String, _ value: T) {
            let output = output; fields[key] = { try output.value(value) }
        }
        func optional<T: Encodable>(_ key: String, _ value: T?) { if let value { field(key, value) } }
        func array<T: Encodable>(_ key: String, _ values: [T]) {
            let output = output; fields[key] = { try output.array(values) }
        }
        func nested(_ key: String, _ body: @escaping () throws -> Void) { fields[key] = body }
    }
    private final class Writer {
        let maximumBytes: Int
        let sink: ((Data) throws -> Void)?
        var bytes = 0, hash = SHA256(), pending = Data()
        init(maximumBytes: Int, sink: ((Data) throws -> Void)?) { self.maximumBytes = maximumBytes; self.sink = sink }
        func append(_ data: Data) throws {
            try Task.checkCancellation()
            guard maximumBytes >= 0, bytes <= maximumBytes, data.count <= maximumBytes - bytes else {
                throw VivoOmicsError.limit("file expression document size")
            }
            bytes += data.count; hash.update(data: data)
            if let sink {
                if pending.count + data.count > 65_536 { try sink(pending); pending.removeAll(keepingCapacity: true) }
                if data.count >= 65_536 { try sink(data) } else { pending.append(data) }
            }
        }
        func raw(_ text: String) throws { try append(Data(text.utf8)) }
        func value<T: Encodable>(_ value: T) throws { try vivoAxisPool { try append(VivoCanonicalJSON.encode(value)) } }
        func array<T: Encodable>(_ values: [T]) throws {
            try raw("[")
            for (index, value) in values.enumerated() { if index > 0 { try raw(",") }; try self.value(value) }
            try raw("]")
        }
        func object(_ configure: (Object) -> Void) throws {
            let o = Object(self); configure(o); try raw("{")
            // These envelope keys are fixed ASCII schema names. Value/key string
            // escaping and number encoding remain owned by VivoCanonicalJSON.
            for (index, key) in o.fields.keys.sorted().enumerated() {
                if index > 0 { try raw(",") }; try value(key); try raw(":"); try o.fields[key]!()
            }
            try raw("}")
        }
        func finish() throws -> VivoFingerprint {
            try Task.checkCancellation(); if let sink, !pending.isEmpty { try sink(pending) }; pending.removeAll()
            return try VivoFingerprint(bytes: Array(hash.finalize()))
        }
        func negativeBinomial(_ v: VivoOmicsNBCohortDiagnostics) throws {
            try object { o in
                o.optional("effectPriorEstimate", v.effectPriorEstimate); o.optional("effectPriorEstimationError", v.effectPriorEstimationError)
                o.optional("effectShrinkage", v.effectShrinkage); o.array("features", v.features); o.field("qualification", v.qualification)
                if let ql = v.quasiLikelihood { o.nested("quasiLikelihood") { try self.quasiLikelihood(ql) } }
                o.field("trend", v.trend)
            }
        }
        func quasiLikelihood(_ v: VivoOmicsNBQLCohortDiagnostics) throws {
            try object { o in
                o.array("abundanceFits", v.abundanceFits); o.optional("averageQLDispersion", v.averageQLDispersion)
                if let failed = v.failedUpstreamFit { o.nested("failedUpstreamFit") { try self.upstream(failed) } }
                o.array("failures", v.failures); o.field("featureIndices", v.featureIndices)
                if let inference = v.inference { o.nested("inference") { try self.inference(inference) } }
                if let moderation = v.moderation { o.nested("moderation") { try self.moderation(moderation) } }
            }
        }
        func inference(_ v: VivoOmicsNBQLTestFamily) throws {
            try object { o in
                o.field("completed", v.completed); o.field("excludedInfluentialIndices", v.excludedInfluentialIndices)
                o.array("failedNullFits", v.failedNullFits); o.array("failures", v.failures)
                o.field("ordinaryResidualDFCap", v.ordinaryResidualDFCap); o.field("poissonBound", v.poissonBound)
                o.field("testedIndices", v.testedIndices); o.array("tests", v.tests)
            }
        }
        func moderation(_ v: VivoOmicsQLModeratedVariances) throws {
            try object { o in
                o.field("commonPriorDegreesOfFreedom", v.commonPriorDegreesOfFreedom); o.field("excludedLowDFIndices", v.excludedLowDFIndices)
                o.field("featureEvaluations", v.featureEvaluations); o.field("flooredVarianceIndices", v.flooredVarianceIndices)
                o.field("informativeIndices", v.informativeIndices); o.optional("notOutlierProbabilities", v.notOutlierProbabilities)
                o.optional("outlierDegreesOfFreedom", v.outlierDegreesOfFreedom); o.field("posteriorVariances", v.posteriorVariances)
                o.field("priorDegreesOfFreedom", v.priorDegreesOfFreedom); o.field("priorScales", v.priorScales); o.array("profiles", v.profiles)
                o.optional("screeningFDRWeights", v.screeningFDRWeights); o.optional("screeningRightProbabilities", v.screeningRightProbabilities)
            }
        }
        func upstream(_ v: VivoOmicsNBQLNativeAbundanceFit) throws {
            try object { o in
                o.array("abundanceFits", v.abundanceFits); o.field("completed", v.completed); o.array("failures", v.failures)
                if let global = v.globalFit { o.nested("globalFit") { try self.global(global) } }
            }
        }
        func global(_ v: VivoOmicsNBQLGlobalFit) throws {
            try object { o in
                o.array("adjustedResiduals", v.adjustedResiduals); o.optional("averageQuasiDispersion", v.averageQuasiDispersion)
                o.field("completed", v.completed); o.array("failures", v.failures); o.array("initialFits", v.initialFits)
                o.array("refittedFits", v.refittedFits); o.array("updates", v.updates)
            }
        }
    }
}
