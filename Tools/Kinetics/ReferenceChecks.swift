import Foundation

@main struct ReferenceChecks {
    static func main() throws {
        var passed = 0
        func check(_ value: Bool, _ label: String) throws {
            guard value else { throw VivoKineticsError.numerical("regression: \(label)") }
            passed += 1
        }
        func near(_ a: Double, _ b: Double, _ label: String, tolerance: Double = 1e-9) throws {
            try check(abs(a - b) <= tolerance, "\(label): \(a) vs \(b)")
        }
        func rejects(_ label: String, _ body: () throws -> Void) throws {
            do { try body() } catch { passed += 1; return }
            throw VivoKineticsError.numerical("expected rejection: \(label)")
        }
        let evidence = VivoKineticEvidence(source: "synthetic mathematical fixture", locator: "ReferenceChecks.swift")
        let context = VivoKineticContext(compound: "synthetic-ligand", target: "synthetic-target", targetVariant: "reference",
            site: "one-site", chemicalState: "declared-synthetic-state", hostContext: "well-mixed-fixture",
            temperatureK: 300, pH: 7, ionicStrengthM: 0.1)
        func parameter(_ value: Double, _ unit: VivoKineticUnit) -> VivoKineticParameter {
            .init(value: value, unit: unit, origin: .assumed, evidence: evidence)
        }
        func model(a: Double = 1e5, b: Double = 0.2, k: Double = 0, d: Double = 0,
                   competitor: VivoKineticCompetitor? = nil) -> VivoCovalentKineticPack {
            .init(identifier: "synthetic", context: context, association: parameter(a, .perMolarSecond),
                  dissociation: parameter(b, .perSecond), inactivation: parameter(k, .perSecond),
                  targetTurnover: parameter(d, .perSecond), baselineTarget: parameter(1e-7, .molar),
                  competitor: competitor, maximumUnboundDrugM: 1e-3)
        }
        func experiment(_ m: VivoCovalentKineticPack, _ knots: [VivoExposureKnot],
                        _ times: [Double], initial: VivoTargetFractions = .init()) -> VivoTargetEngagementExperiment {
            .init(kinetics: m, exposure: .init(context: context, knots: knots, origin: .assumed, evidence: evidence),
                  initial: initial, sampleTimesSeconds: times)
        }
        func constant(_ m: VivoCovalentKineticPack, time: Double, concentration: Double = 1e-6,
                      initial: VivoTargetFractions = .init()) throws -> VivoTargetEngagementResult {
            try VivoTargetEngagementReference.run(experiment(m,
                [.init(timeSeconds: 0, unboundDrugM: concentration), .init(timeSeconds: time, unboundDrugM: concentration)],
                [0, time], initial: initial))
        }
        let reversible = try constant(model(), time: 10)
        try near(reversible.samples.last!.drugOccupancy, 0.1 / 0.3 * (1 - exp(-3)), "reversible analytic solution")
        try near(reversible.samples.first!.fractions.free, 1, "initial publication")
        try near(reversible.samples.last!.totalTargetM, 1e-7, "target conservation", tolerance: 1e-16)
        let covalent = try constant(model(a: 0, b: 0, k: 0.01), time: 100,
                                   initial: .init(free: 0, reversible: 1))
        try near(covalent.samples.last!.fractions.reversible, exp(-1), "covalent precursor decay")
        try near(covalent.samples.last!.fractions.covalent, 1 - exp(-1), "covalent product formation")
        let recovery = try constant(model(a: 0, b: 0, d: 0.01), time: 100, concentration: 0,
                                   initial: .init(free: 0, covalent: 1))
        try near(recovery.samples.last!.fractions.covalent, exp(-1), "target turnover removes covalent complex")
        try near(recovery.samples.last!.fractions.free, 1 - exp(-1), "turnover restores free target")
        let inert = try constant(model(a: 0, b: 0), time: 1e9)
        try near(inert.samples.last!.fractions.free, 1, "zero generator")
        let competitor = VivoKineticCompetitor(identifier: "synthetic-competitor", association: parameter(3e5, .perMolarSecond),
            dissociation: parameter(0.4, .perSecond), unboundConcentration: parameter(1e-6, .molar))
        let competition = try constant(model(competitor: competitor), time: 1000)
        let free: Double = 1.0 / (1.0 + 0.5 + 0.75)
        try near(competition.samples.last!.fractions.free, free, "competitive equilibrium free")
        try near(competition.samples.last!.fractions.reversible, free * 0.5, "competitive equilibrium drug")
        try near(competition.samples.last!.fractions.competitor, free * 0.75, "competitive equilibrium competitor")
        let pulse = experiment(model(), [.init(timeSeconds: 0, unboundDrugM: 1e-6),
            .init(timeSeconds: 10, unboundDrugM: 0), .init(timeSeconds: 20, unboundDrugM: 0)], [0, 10, 20])
        let washout = try VivoTargetEngagementReference.run(pulse)
        try near(washout.samples[1].drugOccupancy, reversible.samples.last!.drugOccupancy, "exposure discontinuity does not reset target")
        try near(washout.samples[2].drugOccupancy, washout.samples[1].drugOccupancy * exp(-2), "washout dissociation")
        try near(washout.samples[1].unboundDrugM, 0, "right-continuous exposure")
        let stiff = try constant(model(a: 1e11, b: 1e5), time: 1)
        try near(stiff.samples.last!.drugOccupancy, 0.5, "stiff relaxation", tolerance: 1e-8)
        let one = try constant(model(k: 0.02, d: 0.001), time: 100)
        let many = try VivoTargetEngagementReference.run(experiment(model(k: 0.02, d: 0.001),
            [.init(timeSeconds: 0, unboundDrugM: 1e-6), .init(timeSeconds: 100, unboundDrugM: 1e-6)],
            (0...100).map(Double.init)))
        try near(one.samples.last!.drugOccupancy, many.samples.last!.drugOccupancy, "observation partition invariance")
        try rejects("units") { try parameter(1, .molar).validate(unit: .perSecond, label: "rate") }
        try rejects("unbound domain") { _ = try constant(model(), time: 10, concentration: 0.1) }
        try rejects("unknown provenance") {
            try VivoKineticParameter(value: 1, unit: .perSecond, origin: .measured, evidence: evidence)
                .validate(unit: .perSecond, label: "measured")
        }
        try rejects("repeated exposure time") {
            try experiment(model(), [.init(timeSeconds: 0, unboundDrugM: 0), .init(timeSeconds: 0, unboundDrugM: 0)], [0]).validate()
        }
        try rejects("initial fractions") { try VivoTargetFractions(free: 1, covalent: 1).validate(hasCompetitor: false) }
        try rejects("missing competitor") { try VivoTargetFractions(free: 0, competitor: 1).validate(hasCompetitor: false) }
        try rejects("numerical budget") {
            _ = try VivoTargetEngagementReference.run(pulse, numerics: .init(maximumPropagations: 1))
        }
        let data = try JSONEncoder().encode(pulse)
        let decoded = try JSONDecoder().decode(VivoTargetEngagementExperiment.self, from: data)
        try check(decoded == pulse, "experiment round trip")
        try rejects("initial numerical tolerance also applies at time zero") {
            let e = experiment(model(), [.init(timeSeconds: 0, unboundDrugM: 0), .init(timeSeconds: 1, unboundDrugM: 0)],
                [0], initial: .init(free: 1 - 1e-11))
            _ = try VivoTargetEngagementReference.run(e, numerics: .init(fractionConservationTolerance: 1e-13))
        }
        let dummyFingerprint = String(repeating: "a", count: 64)
        let source = try VivoTargetEngagementProgramSource.make(pulse,
            kineticsFingerprint: dummyFingerprint, experimentFingerprint: dummyFingerprint)
        let json = try JSONSerialization.jsonObject(with: source) as! [String: Any]
        let spec = json["spec"] as! [String: Any]
        try check(spec["minimumFidelity"] as? String == "F1", "source fidelity")
        let species = spec["species"] as! [[String: Any]]
        try check(species.count == 4 && species.allSatisfy { $0["unit"] as? String == "fraction" }, "normalized target lowering")
        let params = spec["parameters"] as! [[String: Any]]
        try check(params.first { $0["id"] as? String == "kon" }?["unit"] as? String == "M^-1 s^-1", "association unit lowering")
        try rejects("unbound source identities") {
            _ = try VivoTargetEngagementProgramSource.make(pulse, kineticsFingerprint: "not-a-digest", experimentFingerprint: dummyFingerprint)
        }
        func rateRequest(value: Double = 80, unit: VivoMolarEnergyUnit = .kilojoulesPerMol,
                         kind: VivoBarrierQuantity = .activationGibbsFreeEnergy,
                         reference: VivoBarrierReferenceState = .preReactiveBoundComplex,
                         kappa: Double = 0.5, sd: Double? = 2) -> VivoTransitionStateRateRequest {
            .init(barrier: .init(context: context, quantity: kind, referenceState: reference, value: value,
                unit: unit, conditionalStandardDeviation: sd, method: "synthetic analytic fixture",
                samplingDescription: "no molecular sampling; synthetic test value", origin: .assumed, evidence: evidence),
                transmissionProbability: kappa, transmissionOrigin: .assumed, transmissionEvidence: evidence)
        }
        let rate = try VivoTransitionStateRateEstimator.estimate(rateRequest())
        let rt = 8.31446261815324 * 300
        try near(rate.ratePerSecond, 0.5 * (1.380649e-23 * 300 / 6.62607015e-34) * exp(-80000 / rt), "TST formula")
        try near(rate.conditionalLogRateStandardDeviation!, 2000 / rt, "conditional barrier uncertainty propagation")
        let joules = try VivoTransitionStateRateEstimator.estimate(rateRequest(value: 80000, unit: .joulesPerMol, sd: 2000))
        try near(joules.ratePerSecond, rate.ratePerSecond, "TST joule conversion")
        let kcal = try VivoTransitionStateRateEstimator.estimate(rateRequest(value: 80000 / 4184, unit: .kilocaloriesPerMol))
        try near(kcal.ratePerSecond, rate.ratePerSecond, "TST kcal conversion")
        let unknown = try VivoTransitionStateRateEstimator.estimate(rateRequest(sd: nil))
        try check(unknown.conditionalLogRateStandardDeviation == nil && unknown.containsAssumptions, "unknown uncertainty and assumptions retained")
        try rejects("electronic-only barrier") { _ = try VivoTransitionStateRateEstimator.estimate(rateRequest(kind: .electronicEnergyDifference)) }
        try rejects("wrong kinetic standard state") { _ = try VivoTransitionStateRateEstimator.estimate(rateRequest(reference: .separatedReactants)) }
        try rejects("invalid transmission") { _ = try VivoTransitionStateRateEstimator.estimate(rateRequest(kappa: 2)) }
        try rejects("rate underflow is not zero reactivity") { _ = try VivoTransitionStateRateEstimator.estimate(rateRequest(value: 1e8)) }
        let observed = VivoOccupancyObservation(identifier: "synthetic-point", timeSeconds: 10, observable: .drugOccupancy,
            value: 0.1 / 0.3 * (1 - exp(-3)), standardDeviation: 0.02, origin: .assumed, evidence: evidence)
        let study = VivoTargetEngagementStudy(identifier: "synthetic-study", cases: [
            .init(identifier: "synthetic-heldout", leakageGroup: "synthetic-condition-A", partition: .validation,
                  experiment: pulse, observations: [observed])])
        let evaluation = try VivoTargetEngagementStudyEvaluator.evaluate(study)
        try near(evaluation.partitions.first { $0.partition == .validation }!.rootMeanSquaredError!, 0, "study observation mapping")
        try check(evaluation.partitions.first { $0.partition == .test }!.rootMeanSquaredError == nil, "absent partition is not perfect prediction")
        try check(evaluation.partitions.first { $0.partition == .validation }!.measuredObservationCount == 0, "synthetic is not measured")
        try rejects("held-out group leakage") {
            try VivoTargetEngagementStudy(identifier: "leak", cases: study.cases + [
                .init(identifier: "leaked", leakageGroup: "synthetic-condition-A", partition: .calibration,
                      experiment: pulse, observations: [observed])]).validate()
        }
        let failedStudy = try VivoTargetEngagementStudyEvaluator.evaluate(study, numerics: .init(maximumPropagations: 1))
        try check(failedStudy.cases.count == 1 && failedStudy.cases[0].failure != nil, "failed study case retained")
        try check(failedStudy.partitions.first { $0.partition == .validation }!.rootMeanSquaredError == nil, "failed case suppresses aggregate score")
        print("\(passed) kinetic reference/source/chemistry/study checks passed; no Metal or biological validation claimed.")
        if CommandLine.arguments.count == 2 {
            let output = URL(fileURLWithPath: CommandLine.arguments[1], isDirectory: true)
            try FileManager.default.createDirectory(at: output, withIntermediateDirectories: true)
            let encoder = JSONEncoder(); encoder.outputFormatting = [.prettyPrinted, .sortedKeys, .withoutEscapingSlashes]
            let covalentPulse = experiment(model(k: 0.02, d: 0.001), [.init(timeSeconds: 0, unboundDrugM: 1e-6),
                .init(timeSeconds: 10, unboundDrugM: 0), .init(timeSeconds: 100, unboundDrugM: 0)], [0, 1, 2, 5, 10, 20, 40, 100])
            try encoder.encode(covalentPulse).write(to: output.appendingPathComponent("synthetic-pulse.json"))
            try encoder.encode(rateRequest()).write(to: output.appendingPathComponent("synthetic-rate-request.json"))
            try encoder.encode(study).write(to: output.appendingPathComponent("synthetic-study.json"))
        }
    }
}
