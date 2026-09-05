import Foundation
#if NUMIVIVO_PORTABLE_SEPARATE
@testable import NumiVivoPortable
#endif
struct DomainFailure: Error { let message: String }
@main struct DomainChecks {
    static func main() throws {
        var names: [String] = []
        func check(_ b: Bool, _ name: String) throws { guard b else { throw DomainFailure(message: name) }; names.append(name) }
        func reject(_ name: String, _ work: () throws -> Void) throws {
            do { try work() } catch { names.append(name); return }
            throw DomainFailure(message: "unexpected acceptance " + name)
        }
        let evidence = VivoKineticEvidence(source: "synthetic operator check", locator: "DomainChecks.swift")
        let context = VivoKineticContext(compound: "synthetic", target: "synthetic", targetVariant: "reference",
            site: "one", chemicalState: "specified", hostContext: "fixed-volume", temperatureK: 300,pH: 7,ionicStrengthM: 0.1)
        func p(_ x: Double, _ unit: VivoKineticUnit) -> VivoKineticParameter { .init(value: x,unit: unit,origin: .assumed,evidence:evidence) }
        func model(kon: Double = 1e6, koff: Double = 0.2, k: Double = 0.1, d: Double = 0.003,
                   clearance: Double = 0.01, metabolism: Double = 0.02, mclear: Double = 0.005) -> VivoFiniteDrugModel {
            .init(context: context, association:p(kon,.perMolarSecond),dissociation:p(koff,.perSecond),
                  inactivation:p(k,.perSecond),targetTurnover:p(d,.perSecond),baselineTarget:p(2e-6,.molar),
                  freeDrugClearance:p(clearance,.perSecond),metabolism:p(metabolism,.perSecond),
                  metaboliteClearance:p(mclear,.perSecond),maximumFreeDrugM:1e-3)
        }
        let initial=VivoFiniteDrugState(freeDrugM:1e-6,freeTargetM:2e-6)
        let output=try VivoFiniteDrugReactions.advance(initial,model:model(),durationSeconds:100,
            policy:.init(relativeTolerance:1e-8,absoluteToleranceM:1e-16))
        struct References:Decodable { let expectedStates:[String:[Double]] }
        let refs=try JSONDecoder().decode(References.self,from:Data(contentsOf:URL(fileURLWithPath:"Tools/Posterior/domain-reference.json")))
        let error=zip(output.candidate.values,refs.expectedStates["Radau"]!).map{abs($0-$1)}.max()!
        try check(error<2e-11,"finite pool agrees with independent Radau and DOP853")
        try check(output.maximumRelativeMaterialError<1e-10,"drug and target balance through clearance, synthesis and turnover")
        try check(output.candidate.drugBalanceM<=initial.drugBalanceM*(1+1e-10),"finite drug cannot be created by binding")
        try check(output.candidate.inactiveMetaboliteM>0 && output.candidate.eliminatedDrugEquivalentM>0,"metabolite and elimination are separate accounted pools")
        let bindOnly=model(koff:0,k:0,d:0,clearance:0,metabolism:0,mclear:0)
        let limited=try VivoFiniteDrugReactions.advance(initial,model:bindOnly,durationSeconds:100)
        try check(abs(limited.candidate.reversibleComplexM-1e-6)<1e-12,"stoichiometric ligand limitation")
        try check(limited.candidate.occupancy!<=0.5000000001,"occupancy cannot exceed finite ligand supply")
        let equal=try VivoFiniteDrugReactions.advance(.init(freeDrugM:1e-6,freeTargetM:1e-6),model:bindOnly,durationSeconds:10)
        try check(abs(equal.candidate.freeDrugM-1e-6/11)<1e-15,"equal reactant exact association limit")
        let synthesis=try VivoFiniteDrugReactions.advance(.init(freeDrugM:0,freeTargetM:0),model:model(kon:0,koff:0,k:0,clearance:0,metabolism:0,mclear:0),durationSeconds:100)
        try check(abs(synthesis.candidate.freeTargetM-2e-6*(1-exp(-0.3)))<1e-14,"target synthesis and removal exact affine solution")
        let unchanged=try VivoFiniteDrugReactions.advance(initial,model:model(),durationSeconds:0)
        try check(unchanged.initial==unchanged.candidate && unchanged.acceptedSubsteps==0,"zero duration leaves state unchanged")
        try reject("finite-pool budget fails explicitly") { _ = try VivoFiniteDrugReactions.advance(initial,model:model(),durationSeconds:100,policy:.init(maximumAttempts:1)) }
        try reject("invalid physical state rejected") { try VivoFiniteDrugState(freeDrugM:-1,freeTargetM:1).validate() }
        // Finite empirical mixture: identical predictions have zero information;
        // widely separated low-noise predictions reveal one of two components.
        var random=VivoSplitMix64(state:76431)
        let zero=try VivoTargetExperimentDesigner.information(means:[0,0],scales:[0.1,0.1],draws:1024,random:&random)
        try check(abs(zero.mean)<1e-12,"uninformative design has zero information")
        let separated=try VivoTargetExperimentDesigner.information(means:[0,1],scales:[0.01,0.01],draws:1024,random:&random)
        try check(abs(separated.mean-log(2))<1e-10,"separated design approaches log two information")
        let noisy=try VivoTargetExperimentDesigner.information(means:[0,1],scales:[10,10],draws:4096,random:&random)
        try check(noisy.mean<0.02,"large measurement noise reduces information")
        let finiteExample=try VivoCanonicalJSON.decode(VivoFiniteDrugExperiment.self,
            from:Data(contentsOf:URL(fileURLWithPath:"Examples/target-engagement/finite-drug.json")))
        let runRecord=try VivoFiniteDrugRunRecord.run(finiteExample)
        try check(runRecord.result.candidate==output.candidate,"finite-pool example and workflow agree with tested operator")
        let benchmark=try VivoCanonicalJSON.decode(VivoAggregateOccupancyBenchmark.self,
            from:Data(contentsOf:URL(fileURLWithPath:"Examples/btk-aggregate-benchmark/observations.json")))
        try benchmark.validate()
        try check(benchmark.observations.count==2 && benchmark.observations[0].patientCount==12,"published aggregate manifest shape")
        let syntheticEvidence=VivoKineticEvidence(source:"synthetic check only",locator:"all-one virtual cohorts, not a clinical prediction",
            sourceFingerprint:try VivoCanonicalJSON.fingerprint(Data("synthetic all-one population".utf8)).hex)
        let predictions=benchmark.observations.map { obs in
            VivoPredictedOccupancyCohort(observationIdentifier:obs.identifier,compound:obs.compound,target:obs.target,tissue:obs.tissue,
                visit:obs.visit,regimenLabel:obs.regimenLabel,cohortDraws:[[Double](repeating:1,count:obs.patientCount)],
                populationModelEvidence:syntheticEvidence,timingMappingEvidence:syntheticEvidence,predictionEvidence:syntheticEvidence)
        }
        let comparison=try VivoAggregateOccupancyEvaluator.compare(.init(benchmark:benchmark,predictions:predictions))
        try check(abs(comparison.comparisons[0].medianDifference-0.06)<1e-12,"aggregate median comparison does not invent individual observations")
        try check(abs(comparison.comparisons[1].thresholdFractionDifference-0.11)<1e-12,"reported rounded population percentage retained")
        try reject("missing benchmark cohort cannot be filtered away") { _ = try VivoAggregateOccupancyEvaluator.compare(.init(benchmark:benchmark,predictions:Array(predictions.prefix(1)))) }
        if let directory=ProcessInfo.processInfo.environment["NUMIVIVO_CHECK_OUTPUT"] {
            let posterior=try VivoCanonicalJSON.decode(VivoTargetPosteriorRecord.self,
                from:Data(contentsOf:URL(fileURLWithPath:directory).appendingPathComponent("posterior-fit.json")))
            let template=posterior.problem.study.cases[0]
            let candidates=[1.0,10.0,20.0].map { time in VivoTargetDesignCandidate(identifier:"sample-\(time)",
                templateCaseIdentifier:template.identifier,exposureKnots:template.experiment.exposure.knots,
                sampleTimeSeconds:time,observable:.drugOccupancy,reportedStandardDeviation:nil,evidence:evidence) }
            let loaded=try VivoCanonicalJSON.decode([VivoTargetDesignCandidate].self,
                from:Data(contentsOf:URL(fileURLWithPath:"Examples/target-engagement/design-candidates.json")))
            try check(loaded.count==4,"design example decodes")
            let first=try VivoTargetExperimentDesigner.rank(record:posterior,candidates:candidates)
            let reversed=try VivoTargetExperimentDesigner.rank(record:posterior,candidates:Array(candidates.reversed()))
            try check(first.ranking==reversed.ranking,"design ranking invariant to candidate ordering")
            try check(first.ranking.count==3 && first.ranking.allSatisfy{$0.conditionalMonteCarloStandardError.isFinite},"end-to-end kinetic experiment design")
            try VivoCanonicalJSON.encode(first).write(to:URL(fileURLWithPath:directory).appendingPathComponent("design-result.json"))
        }
        struct Report:Codable { let checksPassed:Int;let checks:[String];let finitePoolMaximumAbsoluteErrorM:Double;let finitePool:VivoFiniteDrugStep;let clinicalPredictionExecuted:Bool;let metalExecuted:Bool }
        let report=Report(checksPassed:names.count,checks:names,finitePoolMaximumAbsoluteErrorM:error,finitePool:output,clinicalPredictionExecuted:false,metalExecuted:false)
        if CommandLine.arguments.count==2 { try VivoCanonicalJSON.encode(report).write(to:URL(fileURLWithPath:CommandLine.arguments[1])) }
        print("\(names.count) domain checks passed; finite-pool max reference error \(error) M")
    }
}
