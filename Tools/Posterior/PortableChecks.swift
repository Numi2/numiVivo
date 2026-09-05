import Foundation
#if NUMIVIVO_PORTABLE_SEPARATE
@testable import NumiVivoPortable
#endif

struct CheckFailure: Error, LocalizedError { let message: String; var errorDescription: String? { message } }
struct Pause: Error, Sendable {}
struct References: Decodable {
    let normalLogCDF: [[Double]]
    let intervalLogProbability: [[Double]]
    let correlatedLogLikelihood: Double
}
struct SBCOutcome: Codable {
    let truth: Double
    let observation: Double
    let posteriorMean: Double
    let exactMean: Double
    let rank: Int
    let covered90: Bool
}
@main struct PortableChecks {
    static func main() async throws {
        var names: [String] = []
        func check(_ condition: Bool, _ label: String) throws {
            guard condition else { throw CheckFailure(message: label) }
            names.append(label)
        }
        func near(_ x: Double, _ y: Double, _ label: String, tolerance: Double = 1e-10) throws {
            try check(x.isFinite && abs(x-y) <= tolerance, "\(label): \(x) versus \(y)")
        }
        func rejects(_ label: String, _ operation: () throws -> Void) throws {
            do { try operation() } catch { names.append(label); return }
            throw CheckFailure(message: "unexpected acceptance: " + label)
        }
        let digest = try VivoCanonicalJSON.fingerprint(Data("abc".utf8))
        try check(digest.hex == "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad", "portable SHA-256 known answer")
        let reference = try JSONDecoder().decode(References.self, from: Data(contentsOf: URL(fileURLWithPath: "Tools/Posterior/reference-data.json")))
        for r in reference.normalLogCDF {
            let actual = try VivoGaussianAssay.logCDF(r[0])
            try near(actual, r[1], "normal logCDF \(r[0])", tolerance: max(abs(r[1])*3e-14, 1e-30))
        }
        for r in reference.intervalLogProbability {
            let actual = try VivoGaussianAssay.logContribution(mean: 0, sd: 1, value: 123,
                support: .interval(lower: r[0], upper: r[1]))
            try near(actual, r[2], "interval probability \(r[0]),\(r[1])")
        }
        let lower = try VivoGaussianAssay.logContribution(mean: 0, sd: 1, value: -999, support: .below(limit: 0))
        let upper = try VivoGaussianAssay.logContribution(mean: 0, sd: 1, value: 999, support: .above(limit: 0))
        try near(lower, -log(2), "left censoring probability")
        try near(upper, lower, "right censoring and placeholder independence")
        let data: [VivoGaussianAssayDatum] = [.init(timeSeconds: 0, value: 0.4, reportedStandardDeviation: 0.1),
            .init(timeSeconds: 1, value: 0.2, reportedStandardDeviation: 0.2),
            .init(timeSeconds: 1, value: 0.7, reportedStandardDeviation: 0.3)]
        let correlated = try VivoGaussianAssay.logLikelihood(predictions: [0.3,0.1,0.5], observations: data,
            noise: .init(scale: 1.2, additionalStandardDeviation: 0.03, bias: 0.02, correlatedFraction: 0.7, correlationTimeSeconds: 2))
        try near(correlated, reference.correlatedLogLikelihood, "correlated likelihood against independent NumPy covariance solve")
        try rejects("missing precision remains unknown") { _ = try VivoGaussianAssayNoise().standardDeviation(reported: nil) }
        try near(try VivoGaussianAssayNoise(additionalStandardDeviation: 0.2).standardDeviation(reported: nil), 0.2, "explicit residual SD")
        try rejects("singular correlation rejected") { try VivoGaussianAssayNoise(correlatedFraction: 1).validate() }
        try rejects("correlated censoring rejected") {
            _ = try VivoGaussianAssay.logLikelihood(predictions: [0], observations: [.init(timeSeconds: 0,value: 0,reportedStandardDeviation: 1,support: .below(limit: 0))], noise: .init(correlatedFraction: 0.5))
        }
        try rejects("inverted interval rejected") { try VivoAssaySupport.interval(lower: 1, upper: 0).validate() }
        try rejects("unrepresentable mixture bracket rejected") {
            _ = try VivoGaussianAssay.mixtureQuantile(means: [-1e308,1e308],standardDeviations: [1,1],probability: 0.5)
        }
        try near(try VivoGaussianAssay.mixtureQuantile(means: [0,0],standardDeviations: [1,1],probability: 0.975), 1.959963984540054, "predictive normal quantile")
        let parameters = [VivoPosteriorParameter(identifier: "x", unit: "1", lower: 0, upper: 1, prior: .uniformPhysical),
                          VivoPosteriorParameter(identifier: "y", unit: "1", lower: 0, upper: 1, prior: .uniformPhysical)]
        let evaluate: VivoPosteriorBatchEvaluator = { batch in batch.reversed().map { c in
            let x=(c.values[0]-0.35)/0.08, y=(c.values[1]-0.65)/0.05, rho=0.75
            return .init(candidate: c, logLikelihood: -0.5*(x*x-2*rho*x*y+y*y)/(1-rho*rho))
        } }
        var means: [[Double]] = []
        for seed: UInt64 in [1234567, 12, 9831, 7755] {
            let plan = VivoPosteriorPlan(likelihoodFingerprint: digest,parameters: parameters,
                configuration: .init(particleCount: 256,mutationSweeps: 8,seed: seed))
            let engine = try VivoTemperedPosteriorSampler(plan: plan)
            let result = try await engine.run(evaluate: evaluate)
            try result.validate(requireComplete: true)
            let cp=result.checkpoint!, n=Double(cp.particles.count)
            let mx=cp.particles.reduce(0) { $0+$1.coordinates[0]/n }, my=cp.particles.reduce(0) { $0+$1.coordinates[1]/n }
            try check(abs(mx-0.35)<0.04 && abs(my-0.65)<0.03,"correlated Gaussian seed \(seed)")
            means.append([mx,my])
            if seed == 1234567 {
                let stopped = try VivoTemperedPosteriorSampler(plan: plan)
                let partial = try await stopped.run(evaluate: evaluate,progress: { if $0.stages.count == 1 { throw Pause() } })
                try check(!partial.completed && partial.failure != nil, "partial run is not completed")
                let saved=try VivoCanonicalJSON.decode(VivoPosteriorCheckpoint.self,from: VivoCanonicalJSON.encode(partial.checkpoint!))
                let resumed=try VivoTemperedPosteriorSampler(plan: plan,checkpoint: saved)
                let rerun=try await resumed.run(evaluate: evaluate)
                try check(rerun.checkpoint==cp, "exact stage checkpoint and RNG continuation")
                let changed=VivoPosteriorPlan(likelihoodFingerprint: digest,parameters: parameters,configuration: .init(seed: 8))
                try rejects("changed-plan resume rejected") { try saved.validate(for: changed) }
                let failing=try VivoTemperedPosteriorSampler(plan: plan)
                let failed=try await failing.run(evaluate: { $0.map { .init(candidate: $0, failure: "numerical failure fixture") } })
                try check(!failed.completed && failed.checkpoint == nil && failed.failure?.failedEvaluations.count == 256,"failed evaluations retained")
                let missing=try VivoTemperedPosteriorSampler(plan: plan)
                let malformed=try await missing.run(evaluate: { $0.dropFirst().map { .init(candidate: $0,logLikelihood: 0) } })
                try check(malformed.failure != nil && malformed.checkpoint == nil,"omitted candidate rejected")
            }
        }
        var random=VivoSplitMix64(state: 20260904), sbc: [SBCOutcome]=[]
        for trial in 0..<24 {
            let truth=VivoPosteriorNumerics.openUnit(&random), observation=truth+0.15*random.normal()
            let key=try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode([observation]))
            let plan=VivoPosteriorPlan(likelihoodFingerprint: key,parameters: [parameters[0]],
                configuration: .init(particleCount: 256,mutationSweeps: 8,seed: UInt64(9000+trial)))
            let sampler=try VivoTemperedPosteriorSampler(plan: plan)
            let result=try await sampler.run(evaluate: { $0.map { candidate in
                .init(candidate: candidate,logLikelihood: -0.5*pow((observation-candidate.values[0])/0.15,2))
            } })
            try result.validate(requireComplete: true)
            let draws=result.checkpoint!.particles.map { $0.coordinates[0] }.sorted()
            let a = -observation/0.15, b=(1-observation)/0.15
            let phi: (Double)->Double = { exp(-0.5*$0*$0)/sqrt(2*Double.pi) }
            let normalization=0.5*(erfc(-b/sqrt(2))-erfc(-a/sqrt(2)))
            let exact=observation+0.15*(phi(a)-phi(b))/normalization
            let average=draws.reduce(0,+)/Double(draws.count)
            try check(abs(average-exact)<0.055,"prior-predictive trial \(trial) mean against analytic truncated normal")
            sbc.append(.init(truth: truth,observation: observation,posteriorMean: average,exactMean: exact,
                rank: draws.filter { $0<truth }.count,
                covered90: truth>=VivoPosteriorNumerics.quantile(sorted: draws,probability: 0.05) && truth<=VivoPosteriorNumerics.quantile(sorted: draws,probability: 0.95)))
        }
        let problem = try VivoCanonicalJSON.decode(VivoTargetPosteriorProblem.self,
            from: Data(contentsOf: URL(fileURLWithPath: "Examples/target-engagement/synthetic-assay-inference.json")))
        let prepared = try VivoPreparedTargetPosterior(problem)
        let nominal = try prepared.logLikelihood(values: [100000,0.2,0.01])
        try check(nominal.isFinite,"kinetic likelihood with unknown reported SD and explicit inferred error")
        var j=try JSONSerialization.jsonObject(with: VivoCanonicalJSON.encode(problem)) as! [String:Any]
        var study=j["study"] as! [String:Any], cases=study["cases"] as! [[String:Any]]
        var observations=cases[2]["observations"] as! [[String:Any]]
        observations[0]["value"]=123.0;cases[2]["observations"]=observations;study["cases"]=cases;j["study"]=study
        let heldout=try VivoCanonicalJSON.decode(VivoTargetPosteriorProblem.self,from: JSONSerialization.data(withJSONObject:j))
        let other=try VivoPreparedTargetPosterior(heldout)
        try check(prepared.plan==other.plan,"held-out values excluded from assay likelihood identity")
        try near(nominal,try other.logLikelihood(values: [100000,0.2,0.01]),"held-out values excluded from assay likelihood")
        let fit=try await VivoTargetPosteriorFitter.run(problem)
        try fit.validate(requireComplete:true)
        let values=try fit.posterior.checkpoint!.particles.map { try prepared.physicalValues($0) }
        let posteriorMeans=(0..<3).map { col in values.reduce(0) { $0+$1[col]/Double(values.count) } }
        try check(abs(posteriorMeans[0]-100000)<25000 && abs(posteriorMeans[1]-0.2)<0.08,"end-to-end kinetic fit with inferred SD")
        try check(posteriorMeans[2]>0.002 && posteriorMeans[2]<0.03,"noise scale is inferred, not silently fixed")
        try rejects("mean-only information cannot qualify learned assay noise") { _ = try VivoTargetPosteriorSensitivity.evaluate(fit) }
        let predictive = try await VivoTargetPosteriorPredictor.predict(fit)
        try check(predictive.schemaVersion == 2 && predictive.cases.allSatisfy { $0.failures.isEmpty }, "assay posterior predictive path")
        try check(predictive.cases.flatMap(\.points).allSatisfy { $0.measurementSD == nil && $0.predictiveLower != nil }, "inferred noise propagates without inventing reported precision")
        try check(predictive.cases.allSatisfy { $0.jointLogObservationLikelihood?.isFinite == true }, "joint assay predictive likelihood")
        let correlatedAssays = problem.study.cases.map {
            VivoTargetAssayModel(caseIdentifier:$0.identifier,noise:.init(correlatedFraction:0.4,correlationTimeSeconds:3),
                evidence:.init(source:"synthetic",locator:"correlated full-path check"))
        }
        let correlatedProblem=VivoTargetPosteriorProblem(study:problem.study,bindings:problem.bindings,sampler:problem.sampler,
            numerics:problem.numerics,parallelEvaluations:problem.parallelEvaluations,assays:correlatedAssays)
        let correlatedFit=try await VivoTargetPosteriorFitter.run(correlatedProblem)
        try correlatedFit.validate(requireComplete:true)
        let correlatedReport=try await VivoTargetPosteriorPredictor.predict(correlatedFit)
        try check(correlatedReport.cases.allSatisfy{$0.failures.isEmpty},"correlated assay fit and prediction execute")
        let c0=correlatedReport.cases[0]
        let marginalSum=c0.points.compactMap(\.logPredictiveDensity).reduce(0,+)
        try check(abs(c0.jointLogObservationLikelihood!-marginalSum)>1e-4,"joint correlated score is not sum of marginal scores")
        let censored=VivoTargetAssayModel(caseIdentifier:problem.study.cases[0].identifier,
            support:[problem.study.cases[0].observations[0].identifier:.below(limit:0.15)],
            evidence:.init(source:"synthetic",locator:"censored full-path check"))
        let censoredProblem=VivoTargetPosteriorProblem(study:problem.study,bindings:problem.bindings,sampler:problem.sampler,
            numerics:problem.numerics,parallelEvaluations:problem.parallelEvaluations,assays:[censored])
        let censoredFit=try await VivoTargetPosteriorFitter.run(censoredProblem)
        try censoredFit.validate(requireComplete:true)
        let censoredReport=try await VivoTargetPosteriorPredictor.predict(censoredFit)
        let point=censoredReport.cases[0].points.first!
        try check(point.logObservationProbability != nil && point.logObservationProbability!<=0,"censored predictive event probability")
        try check(point.logPredictiveDensity==nil && point.predictiveIntervalContainsObservation==nil,"censoring never scores placeholder as exact observed value")
        let recordData=try VivoCanonicalJSON.encode(fit)
        if let outputDirectory = ProcessInfo.processInfo.environment["NUMIVIVO_CHECK_OUTPUT"] {
            try recordData.write(to: URL(fileURLWithPath: outputDirectory).appendingPathComponent("posterior-fit.json"))
        }
        struct Report: Encodable {
            let checksPassed: Int
            let checks: [String]
            let gaussianMeans: [[Double]]
            let priorPredictiveTrials: [SBCOutcome]
            let sbcInterpretation: String
            let fittedMeans: [Double]
            let applePackageBuilt: Bool
            let metalExecuted: Bool
        }
        let report=Report(checksPassed:names.count,checks:names,gaussianMeans:means,priorPredictiveTrials:sbc,
            sbcInterpretation:"24 fixed-seed prior-predictive trials against an analytic posterior. Ranks are descriptive: SMC particles are dependent. Not a general calibration certificate or biological validation.",fittedMeans:posteriorMeans,applePackageBuilt:false,metalExecuted:false)
        let dataOut=try VivoCanonicalJSON.encode(report)
        if CommandLine.arguments.count==2 { try dataOut.write(to: URL(fileURLWithPath:CommandLine.arguments[1])) }
        print("\(names.count) portable checks passed; prior-predictive 90% coverage \(sbc.filter(\.covered90).count)/\(sbc.count); kinetic means \(posteriorMeans)")
    }
}
