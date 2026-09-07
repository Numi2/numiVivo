import Foundation

@main @MainActor struct AdaptiveChemistryChecks {
    static func main() throws {
        let output = URL(fileURLWithPath: CommandLine.arguments.count > 1 ? CommandLine.arguments[1] : "/tmp/numivivo-adaptive-checks")
        try FileManager.default.createDirectory(at: output,withIntermediateDirectories: true)
        var checks: [String] = []
        func expect(_ condition: Bool,_ label: String) throws {
            guard condition else { throw VivoChemistryError.invalid("assertion failed: \(label)") }
            checks.append(label)
        }
        func expectThrows(_ label: String,_ body: () throws -> Void) throws {
            do { try body() } catch { checks.append(label); return }
            throw VivoChemistryError.invalid("expected rejection: \(label)")
        }
        func save<T: Encodable>(_ value: T,_ name: String) throws {
            let encoder = JSONEncoder(); encoder.outputFormatting = [.prettyPrinted,.sortedKeys]
            try encoder.encode(value).write(to: output.appendingPathComponent(name))
        }
        func mixture(seed: UInt64,count: Int) -> VivoMBARTargetRefinementRequest {
            let energies = [0.0,log(2.0),log(3.0),log(4.0)]
            var samples: [VivoMBARTargetSample] = [], random = VivoSplitMix64(state: seed)
            for origin in 0..<2 {
                let weights = energies.map { origin == 0 ? 1.0 : exp(-$0) }, total = weights.reduce(0,+)
                for draw in 0..<count {
                    var target = random.unit()*total, microstate = weights.count-1
                    for i in weights.indices { target -= weights[i]; if target <= 0 { microstate = i; break } }
                    let identifier = "sample-\(seed)-\(origin)-\(draw)"
                    samples.append(.init(identifier: identifier,originStateIndex: origin,independentBlockIdentifier: identifier,
                        coordinate: microstate < 2 ? 0 : 1,sampledReducedPotentials: [0,energies[microstate]],targetReducedPotentials: [0,energies[microstate]]))
                }
            }
            let source = VivoMBARConfiguration(minimumTargetEffectiveSamples: 10,bootstrapReplicates: 0)
            let cfg = VivoMBARTargetConfiguration(sourceSolver: source,bootstrapReplicates: 32,minimumBinEffectiveSamples: 10,
                maximumBinWeight: 0.05,maximumPrimitiveElements: 100_000_000)
            return .init(sampledStateIdentifiers: ["uniform","tilted"],targetIdentifiers: ["reference","refined"],
                matchingSampledStateIndices: [0,1],binEdges: [-0.5,0.5,1.5],coordinateUnit: "finite-state-group",samples: samples,configuration: cfg)
        }
        let request = mixture(seed: 19,count: 1024), result = try VivoMBARTargetRefinement.calculate(request)
        let exact = log(18.0/7.0), correction = result.corrections[1].deltaRelativeFreeEnergy!
        try expect(abs(correction-exact) < 0.06,"paired target bin correction agrees with enumerated partition functions")
        try expect(result.qualification == .sampledTargetsPassDiagnostics,"independently sampled targets pass conditional weight diagnostics")
        try expect(result.completedBootstrapReplicates == 32 && result.correctionCovariance != nil,"shared bootstrap retains full covariance")
        try expect(result.corrections[0].deltaRelativeFreeEnergy == 0,"reference-bin correction cancels exactly")
        try expect(abs(result.correctionCovariance![0][0]) < 1e-25,"shared bootstrap preserves exact covariance cancellation")
        try expect(result.chargedPrimitiveElements <= request.configuration.maximumPrimitiveElements,"target correction charges aggregate numerical work")
        try VivoMBARTargetRefinement.validate(result,request: request)
        try expect(true,"target profile validates by deterministic reconstruction")
        try save(result,"target-correction-native.json")
        var shifted = request
        for i in shifted.samples.indices {
            let offset = sin(Double(i))*17
            shifted.samples[i].sampledReducedPotentials = shifted.samples[i].sampledReducedPotentials.map { $0+offset }
            shifted.samples[i].targetReducedPotentials = shifted.samples[i].targetReducedPotentials.map { $0+offset }
        }
        shifted.configuration.bootstrapReplicates = 0
        let shiftResult = try VivoMBARTargetRefinement.calculate(shifted)
        try expect(abs(shiftResult.corrections[1].deltaRelativeFreeEnergy!-correction) < 2e-10,"per-configuration energy gauges cancel")
        var unsampled = request; unsampled.matchingSampledStateIndices[1] = nil
        let pilot = try VivoMBARTargetRefinement.calculate(unsampled)
        try expect(pilot.qualification == .exploratoryUnsampledTarget,"unsampled target is never accepted from ESS alone")
        var falselyMatched = request; falselyMatched.samples[0].targetReducedPotentials[1] += 0.01
        try expectThrows("forged sampled-target Hamiltonian column rejects") { _ = try VivoMBARTargetRefinement.calculate(falselyMatched) }
        var crossed = request
        crossed.samples[1024].independentBlockIdentifier = crossed.samples[0].independentBlockIdentifier
        try expectThrows("cross-origin blocks require joint exchange analysis") { _ = try VivoMBARTargetRefinement.calculate(crossed) }
        var tiny = request; tiny.configuration.maximumPrimitiveElements = 1
        try expectThrows("target correction budget exhaustion does not return acceptance") { _ = try VivoMBARTargetRefinement.calculate(tiny) }
        var one = request
        one.sampledStateIdentifiers = ["uniform"]; one.matchingSampledStateIndices = [0,nil]
        one.samples = request.samples.filter { $0.originStateIndex == 0 }.map {
            .init(identifier: $0.identifier,originStateIndex: 0,independentBlockIdentifier: $0.independentBlockIdentifier,
                  coordinate: $0.coordinate,sampledReducedPotentials: [0],targetReducedPotentials: $0.targetReducedPotentials)
        }
        try expect(try VivoMBARTargetRefinement.calculate(one).qualification == .exploratoryUnsampledTarget,"single-source MBAR limit supports exploratory target corrections")
        var empty = request
        for i in empty.samples.indices { empty.samples[i].coordinate = 0 }
        let emptyResult = try VivoMBARTargetRefinement.calculate(empty)
        try expect(emptyResult.qualification == .insufficientEvidence && emptyResult.corrections[1].deltaRelativeFreeEnergy == nil,
                   "unvisited bins remain missing rather than fabricated finite barriers")

        let q = try VivoQMMatrix(rows: 1,columns: 1,values: [-2.0])
        let sensitivity = try VivoLinearKineticSensitivity.calculate(generator: q,initialProbability: [1],
            observationTimesSeconds: [0,0.2,1],derivatives: [.init(identifier: "log-k",generatorDerivativePerSecond: q)])
        for point in sensitivity.observations {
            let t = point.timeSeconds, p = exp(-2*t), derivative = point.derivatives[0]
            try expect(abs(point.survivalProbability-p) < 2e-12,"one-state survival \(t)")
            try expect(abs(derivative.survivalProbability+2*t*p) < 2e-10,"analytic log-rate survival derivative \(t)")
            try expect(abs(derivative.hazardPerSecond!-2) < 2e-10,"analytic log-rate hazard derivative \(t)")
        }
        let generator = try VivoQMMatrix(rows: 2,columns: 2,values: [-0.8,0.5,0.2,-0.6])
        let b = try VivoQMMatrix(rows: 2,columns: 2,values: [-0.5,0.5,0,0])
        let times = [0.0,0.15,1.0,7.0]
        let exactSensitivity = try VivoLinearKineticSensitivity.calculate(generator: generator,initialProbability: [0.6,0.4],
            observationTimesSeconds: times,derivatives: [.init(identifier: "log-exchange-0-1",generatorDerivativePerSecond: b)])
        struct KineticExport: Encodable { let generator: VivoQMMatrix; let derivative: VivoQMMatrix; let initial: [Double]; let result: VivoLinearKineticSensitivityResult }
        try save(KineticExport(generator: generator,derivative: b,initial: [0.6,0.4],result: exactSensitivity),"kinetic-sensitivity-native.json")
        let epsilon = 1e-5
        let plus = try generator.adding(b.scaled(epsilon)), minus = try generator.adding(b.scaled(-epsilon))
        let fp = try VivoLinearKineticSensitivity.calculate(generator: plus,initialProbability: [0.6,0.4],observationTimesSeconds: times)
        let fm = try VivoLinearKineticSensitivity.calculate(generator: minus,initialProbability: [0.6,0.4],observationTimesSeconds: times)
        var derivativeError = 0.0
        for i in times.indices { for state in 0..<2 {
            derivativeError = max(derivativeError,abs((fp.observations[i].probabilityByState[state]-fm.observations[i].probabilityByState[state])/(2*epsilon)
                - exactSensitivity.observations[i].derivatives[0].probabilityByState[state]))
        } }
        try expect(derivativeError < 1e-8,"shared recurrence agrees with independent finite differences")
        let exchange = try VivoQMMatrix(rows: 2,columns: 2,values: [-3,3,2,-2])
        let conservative = try VivoLinearKineticSensitivity.calculate(generator: exchange,initialProbability: [1,0],observationTimesSeconds: [0,10],
            derivatives: [.init(identifier: "exchange",generatorDerivativePerSecond: exchange)])
        try expect(abs(conservative.observations[1].survivalProbability-1) < 1e-11,"closed generator conserves probability")
        try expect(abs(conservative.observations[1].derivatives[0].survivalProbability) < 1e-10,"closed-generator total sensitivity is zero")
        let zero = try VivoLinearKineticSensitivity.calculate(generator: VivoQMMatrix(1,1),initialProbability: [1],observationTimesSeconds: [2],
            derivatives: [.init(identifier: "zero-limit",parameterization: "dimensionless-parameter",generatorDerivativePerSecond: q)])
        try expect(zero.observations[0].derivatives[0].survivalProbability == -4,"zero-generator derivative has correct analytic limit")
        try expectThrows("negative transition generator rejects") {
            _ = try VivoLinearKineticSensitivity.calculate(generator: try .init(rows: 2,columns: 2,values: [-1,-1,0,-1]),initialProbability: [1,0],observationTimesSeconds: [1])
        }
        try expectThrows("joint kinetic work budget rejects") {
            _ = try VivoLinearKineticSensitivity.calculate(generator: generator,initialProbability: [1,0],observationTimesSeconds: [100],configuration: .init(maximumPrimitiveWork: 10))
        }
        let covariance = try VivoQMMatrix(rows: 2,columns: 2,values: [0.04,0.03,0.03,0.04])
        let projection = try VivoObservableCovariance.project(derivatives: [1,-1],covariance: covariance)
        try expect(abs(projection.variance-0.02) < 1e-12,"observable propagation retains off-diagonal covariance")
        try expectThrows("indefinite covariance rejects") {
            _ = try VivoObservableCovariance.project(derivatives: [1,1],covariance: try .init(rows: 2,columns: 2,values: [1,2,2,1]))
        }

        let context = String(repeating: "a",count: 64)
        let criterion = VivoRefinementCriterion(identifier: "rate",metricUnit: "log-rate-sd",observableFingerprint: context,maximumMetric: 0.1)
        let proposals: [VivoRefinementActionProposal] = [
            .init(identifier: "d",criterionIdentifier: "rate",role: .discovery,costClass: "same-native-backend",declaredWorkUnits: 10,initialEstimatedSeconds: 20,expectedMetricReduction: 0.3),
            .init(identifier: "c",criterionIdentifier: "rate",role: .confirmation,prerequisites: ["d"],costClass: "same-native-backend",declaredWorkUnits: 20,initialEstimatedSeconds: 40,expectedMetricReduction: 0)]
        func metric(_ value: Double,_ passed: Bool,_ id: String,_ sources: [String],observable: String? = nil) -> VivoRefinementMetricEvidence {
            .init(metricUnit: "log-rate-sd",observableFingerprint: observable ?? context,value: value,nativeChecksPassed: passed,
                  evidenceIdentifier: String(repeating: id,count: 64),executionSourceIdentifiers: sources,interpretation: "test fixture metric")
        }
        let first = try VivoAdaptiveRefinementPolicy.decide(criteria: [criterion],proposals: proposals,records: [],remainingDeclaredWorkUnits: 100)
        try expect(first.selectedActionIdentifier == "d" && !first.allCriteriaSatisfied,"discovery precedes confirmation")
        let observed = VivoRefinementActionRecord(actionIdentifier: "d",metric: metric(0.05,true,"1",["discovery-source"]),elapsedSeconds: 2,entirelyUncachedExecution: true,failure: nil)
        let second = try VivoAdaptiveRefinementPolicy.decide(criteria: [criterion],proposals: proposals,records: [observed],remainingDeclaredWorkUnits: 90)
        try expect(second.selectedActionIdentifier == "c" && second.rankedActions[0].estimatedSeconds == 4,"fresh measured cost calibrates subsequent priority")
        let cached = VivoRefinementActionRecord(actionIdentifier: "d",metric: observed.metric,elapsedSeconds: 0.01,entirelyUncachedExecution: false,failure: nil)
        let cachedDecision = try VivoAdaptiveRefinementPolicy.decide(criteria: [criterion],proposals: proposals,records: [cached],remainingDeclaredWorkUnits: 90)
        try expect(cachedDecision.rankedActions[0].estimatedSeconds == 40,"cache replay latency is not fresh-compute calibration")
        let confirmed = VivoRefinementActionRecord(actionIdentifier: "c",metric: metric(0.06,true,"2",["independent-confirmation"]),elapsedSeconds: 3,entirelyUncachedExecution: true,failure: nil)
        try expect(try VivoAdaptiveRefinementPolicy.decide(criteria: [criterion],proposals: proposals,records: [observed,confirmed],remainingDeclaredWorkUnits: 70).allCriteriaSatisfied,"native disjoint confirmation satisfies declared criterion")
        let transplanted = VivoRefinementActionRecord(actionIdentifier: "c",metric: metric(0.01,true,"3",["discovery-source"]),elapsedSeconds: 3,entirelyUncachedExecution: true,failure: nil)
        let rejected = try VivoAdaptiveRefinementPolicy.decide(criteria: [criterion],proposals: proposals,records: [observed,transplanted],remainingDeclaredWorkUnits: 70)
        try expect(rejected.criteria[0].confirmationFailed && rejected.rankedActions.isEmpty,"reused held-out source terminates without further tuning")
        let otherTarget = VivoRefinementActionRecord(actionIdentifier: "d",metric: metric(0.01,true,"4",["other"],observable: String(repeating: "b",count: 64)),elapsedSeconds: 1,entirelyUncachedExecution: true,failure: nil)
        try expectThrows("a metric from another observable cannot be transplanted") {
            _ = try VivoAdaptiveRefinementPolicy.decide(criteria: [criterion],proposals: proposals,records: [otherTarget],remainingDeclaredWorkUnits: 90)
        }
        try expectThrows("confirmation cannot precede its successful prerequisite") {
            _ = try VivoAdaptiveRefinementPolicy.decide(criteria: [criterion],proposals: proposals,
                records: [confirmed],remainingDeclaredWorkUnits: 100)
        }
        let tune = VivoRefinementActionProposal(identifier: "tune",criterionIdentifier: "rate",role: .discovery,
            prerequisites: ["c"],costClass: "same-native-backend",declaredWorkUnits: 10,initialEstimatedSeconds: 1,expectedMetricReduction: 1)
        let tuned = VivoRefinementActionRecord(actionIdentifier: "tune",metric: metric(0.02,true,"5",["tuned"]),elapsedSeconds: 1,entirelyUncachedExecution: true,failure: nil)
        try expectThrows("a confirmation cannot become discovery tuning evidence") {
            _ = try VivoAdaptiveRefinementPolicy.decide(criteria: [criterion],proposals: proposals+[tune],
                records: [observed,confirmed,tuned],remainingDeclaredWorkUnits: 100)
        }
        let tinyCovariance = try VivoQMMatrix(rows: 2,columns: 2,values: [1e-320,0,0,1e-320])
        let tiny = try VivoObservableCovariance.project(derivatives: [1,-1],covariance: tinyCovariance)
        try expect(tiny.standardDeviation.isFinite && tiny.standardDeviation > 0,"subnormal covariance normalization stays finite")
        try expect(tiny.variance == 2e-320,"subnormal covariance retains the represented variance")
        var loop = proposals; loop[0].prerequisites = ["c"]
        try expectThrows("cyclic adaptive actions reject") { try VivoAdaptiveRefinementPolicy.validate(criteria: [criterion],proposals: loop) }
        try expect(try VivoAdaptiveRefinementPolicy.decide(criteria: [criterion],proposals: proposals,records: [],remainingDeclaredWorkUnits: 1).selectedActionIdentifier == nil,
                   "admission work exhaustion never means criteria satisfied")
        struct Report: Encodable { let schema: String; let assertions: Int; let checks: [String]; let exactFiniteStateCorrection: Double; let nativeFiniteStateCorrection: Double; let kineticFiniteDifferenceMaximumError: Double }
        try save(Report(schema: "numivivo.org/adaptive-numerical-observations/v1",assertions: checks.count,checks: checks,
            exactFiniteStateCorrection: exact,nativeFiniteStateCorrection: correction,kineticFiniteDifferenceMaximumError: derivativeError),"results.json")
        print("PASS \(checks.count) adaptive chemistry assertions")
    }
}
