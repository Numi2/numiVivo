import Foundation
#if NUMIVIVO_PORTABLE_CHECKS
@testable import NumiVivoNumerics
#else
import Testing
@testable import NumiVivoKit
#endif

struct BarrierConvergenceRegressions {
    struct Report: Codable {
        let checks: [String]
        let observations: [String: Double]
        let scope: String
    }
    struct Check {
        var labels: [String] = []
        mutating func require(_ condition: Bool, _ label: String) throws {
            guard condition else { throw VivoChemistryError.invalid("barrier regression: \(label)") }
            labels.append(label)
        }
        mutating func rejects(_ label: String, _ f: () throws -> Void) throws {
            var rejected = false
            do { try f() } catch is VivoChemistryError { rejected = true } catch is VivoECCPathError { rejected = true }
            try require(rejected, label)
        }
    }
    static func modified<T: Codable>(_ value: T, _ edit: (inout [String: Any]) -> Void) throws -> T {
        var object = try JSONSerialization.jsonObject(with: JSONEncoder().encode(value)) as! [String: Any]
        edit(&object)
        return try JSONDecoder().decode(T.self, from: JSONSerialization.data(withJSONObject: object))
    }
    static func record(_ report: Report, _ name: String) throws {
        guard let folder = ProcessInfo.processInfo.environment["NUMIVIVO_TEST_ARTIFACTS"] else { return }
        let path = URL(fileURLWithPath: folder, isDirectory: true)
        try FileManager.default.createDirectory(at: path, withIntermediateDirectories: true)
        let encoder = JSONEncoder(); encoder.outputFormatting = [.prettyPrinted, .sortedKeys]
        try encoder.encode(report).write(to: path.appendingPathComponent("barrier-\(name).json"), options: .atomic)
    }
    static func physical() throws -> Report {
        var c = Check()
        let request = VivoBarrierBenchmarks.hydrogenExchange631G(), result = try VivoBarrierConvergence.run(request)
        try c.require(result.levels.count == 4 && result.referenceEnergiesHartree.count == 9, "all four levels and nine mapped geometries execute")
        try c.require(result.referenceDifferences.forwardHartree > 0 && result.referenceDifferences.reverseHartree > 0, "full-space electronic reference has a positive central profile difference")
        try c.require(result.levels.allSatisfy { $0.outcome.evaluation != nil }, "numerically converged levels are not confused with accurate barriers")
        try c.require(result.levels.dropLast().allSatisfy { $0.outcome.evaluation!.isGenuinelyReduced && !$0.outcome.evaluation!.meetsReferenceAccuracy }, "three genuinely reduced spaces fail the predeclared accuracy threshold")
        try c.require(!result.levels.last!.outcome.evaluation!.isGenuinelyReduced && result.levels.last!.outcome.evaluation!.meetsReferenceAccuracy, "full-space agreement is retained as a separate control")
        try c.require(result.assessment == .reducedAccuracyNotEstablished && result.acceptedReducedLevelIdentifier == nil, "full-space control does not qualify an inaccurate reduction")
        try c.require(result.referenceResiduals.allSatisfy { $0 <= 1.01*request.referenceResidualTolerance }, "full-space physical eigen residuals meet reference contract")
        try c.require(result.chargedPointEvaluations == 45, "reference and CAS point-solve accounting is exact")
        try c.require(result.transportMinimumSingularValues.allSatisfy { $0 >= request.minimumTransportSingularValue }, "actual cross-geometry subspaces pass without lowering overlap threshold")
        try c.require(result.adjacentReferenceOverlapsSquared.count == 8 && result.adjacentReferenceOverlapsSquared.allSatisfy { $0 >= request.minimumReferenceOverlapSquared }, "physical ground-state overlap is retained at every edge")
        try VivoBarrierConvergence.validate(result, request: request)
        try c.require(true, "negative scientific assessment passes full physical reconstruction")
        return .init(checks: c.labels, observations: ["referenceForwardHartree": result.referenceDifferences.forwardHartree,
            "fiveOrbitalBarrierErrorHartree": result.levels[2].outcome.evaluation!.maximumBarrierErrorHartree],
            scope: "H3/6-31G declared scan; unsuccessful truncated models, no saddle or rate certification")
    }
    static func tampering() throws -> Report {
        var c = Check()
        let request = VivoBarrierBenchmarks.hydrogenExchange631G(), result = try VivoBarrierConvergence.run(request)
        try c.rejects("forged acceptance flag rejected") {
            let x = try modified(result) { $0["assessment"] = "reducedAccuracyEstablished"; $0["acceptedReducedLevelIdentifier"] = "CAS-3e-5o" }
            try VivoBarrierConvergence.validate(x, request: request)
        }
        try c.rejects("forged full-space reference energy rejected") {
            let x = try modified(result) { $0["referenceEnergiesHartree"] = [Double](repeating: 99, count: 9) }
            try VivoBarrierConvergence.validate(x, request: request)
        }
        try c.rejects("missing unsuccessful approximation cannot be cherry-picked away") {
            let x = try modified(result) { var levels = $0["levels"] as! [Any]; levels.removeFirst(); $0["levels"] = levels }
            try VivoBarrierConvergence.validate(x, request: request)
        }
        try c.rejects("forged stage barrier rejected") {
            let x = try modified(result) { object in
                var levels = object["levels"] as! [[String: Any]]
                var outcome = levels[0]["outcome"] as! [String: Any]
                var evaluated = outcome["evaluated"] as! [String: Any]
                var value = evaluated["result"] as! [String: Any]
                value["maximumBarrierErrorHartree"] = 0; evaluated["result"] = value; outcome["evaluated"] = evaluated
                levels[0]["outcome"] = outcome; object["levels"] = levels
            }
            try VivoBarrierConvergence.validate(x, request: request)
        }
        try c.rejects("electronic report cannot acquire a saddle or rate label") {
            let x = try modified(result) { $0["meaning"] = "qualified kinetic rate" }
            try VivoBarrierConvergence.validate(x, request: request)
        }
        try c.rejects("incorrect work charge rejected") {
            let x = try modified(result) { $0["chargedPointEvaluations"] = 1 }
            try VivoBarrierConvergence.validate(x, request: request)
        }
        try c.rejects("forged physical reference overlap rejected") {
            let x = try modified(result) { $0["adjacentReferenceOverlapsSquared"] = [Double](repeating: 0, count: 8) }
            try VivoBarrierConvergence.validate(x, request: request)
        }
        return .init(checks: c.labels, observations: [:], scope: "stored report acceptance and scientific-meaning integrity")
    }
    static func contracts() throws -> Report {
        var c = Check(); let request = VivoBarrierBenchmarks.hydrogenExchange631G()
        try c.rejects("non-nested CAS levels rejected before execution") {
            let r = try modified(request) { object in
                var levels = object["levels"] as! [[String: Any]]; levels.swapAt(0,1); object["levels"] = levels
            }; try r.validate()
        }
        try c.rejects("unfunded campaign rejected before allocation") {
            let r = try modified(request) { $0["maximumPointEvaluations"] = 44 }; try r.validate()
        }
        try c.rejects("changed atom identity rejected") {
            let r = try modified(request) { object in
                var points = object["snapshots"] as! [[String: Any]], system = points[0]["system"] as! [String: Any]
                var nuclei = system["nuclei"] as! [[String: Any]]; nuclei[0]["structureAtomIndex"] = 123
                system["nuclei"] = nuclei; points[0]["system"] = system; object["snapshots"] = points
            }; try r.validate()
        }
        try c.rejects("barrier evaluation point must be internal") {
            let r = try modified(request) { $0["barrierPointIdentifier"] = "h3-0" }; try r.validate()
        }
        try c.rejects("reference tolerance cannot exceed relative-profile accuracy contract") {
            let r = try modified(request) { object in
                var acceptance = object["acceptance"] as! [String: Any]
                acceptance["maximumRelativeProfileErrorHartree"] = 1e-14; object["acceptance"] = acceptance
            }; try r.validate()
        }
        try c.rejects("empty transport group rejected") {
            let r = try modified(request) { $0["transportGroups"] = [[], [0,1,2,3,4,5]] }; try r.validate()
        }
        // Repeated identical geometries are only an algebraic zero-difference
        // fixture to exercise the positive decision path, NOT chemical evidence.
        let constant = try modified(request) { object in
            var points = object["snapshots"] as! [[String: Any]]
            let system = points[4]["system"]!
            for i in points.indices { points[i]["system"] = system }
            object["snapshots"] = points
        }
        let equal = try VivoBarrierConvergence.run(constant)
        try c.require(equal.assessment == .reducedAccuracyEstablished && equal.acceptedReducedLevelIdentifier == "CAS-3e-5o", "two stable reduced zero-difference levels can satisfy explicit numerical criteria")
        let oneReduced = try modified(constant) { object in
            let levels = object["levels"] as! [[String: Any]]; object["levels"] = Array(levels.suffix(2))
        }
        let onlyFull = try VivoBarrierConvergence.run(oneReduced)
        try c.require(onlyFull.assessment == .reducedAccuracyNotEstablished, "one reduced level plus full space is insufficient for the stability window")
        return .init(checks: c.labels, observations: [:], scope: "invalid requests and positive/negative decision controls; identical-geometry case is algebra only")
    }
    static func ensemble() throws -> Report {
        var c = Check()
        let plain = try VivoBarrierConvergence.run(VivoBarrierBenchmarks.hydrogenExchange631G())
        let request = VivoBarrierBenchmarks.hydrogenExchange631G(ensemble: true), result = try VivoBarrierConvergence.run(request)
        try c.require(result.ensembleOccupations?.count == 6 && result.ensembleRotation?.rows == 6, "common ensemble-density frame is recorded explicitly")
        try c.require(abs(result.ensembleOccupations!.reduce(0,+)-3) < 1e-8, "ensemble natural occupation trace preserves the three electrons")
        try c.require(zip(result.referenceEnergiesHartree,plain.referenceEnergiesHartree).allSatisfy { abs($0-$1) < 1e-9 }, "reference energies are invariant to the common ensemble rotation")
        let improved = result.levels[2].outcome.evaluation!, baseline = plain.levels[2].outcome.evaluation!
        try c.require(improved.maximumRelativeProfileErrorHartree < baseline.maximumRelativeProfileErrorHartree/100, "density-informed shared orbitals improve this declared five-orbital profile without fitting energies")
        try c.require(improved.meetsReferenceAccuracy && !result.levels[1].outcome.evaluation!.meetsReferenceAccuracy, "five orbitals meet the accuracy target but the adjacent four-orbital level does not")
        try c.require(result.assessment == .reducedAccuracyNotEstablished, "one accurate reduced space does not establish a stable reduced ladder")
        try VivoBarrierConvergence.validate(result,request:request)
        try c.require(true, "ensemble frame, occupations and accuracy claims reconstruct")
        try c.rejects("unnormalized ensemble density weights rejected") {
            let r = try modified(request) { object in
                var ensemble = object["ensembleOrbitals"] as! [String:Any]
                ensemble["pointWeights"] = [Double](repeating:1,count:9); object["ensembleOrbitals"] = ensemble
            }; try r.validate()
        }
        let nearDegenerate = try modified(request) { object in
            var ensemble = object["ensembleOrbitals"] as! [String:Any]
            ensemble["minimumOccupationBoundaryGap"] = 3.0; object["ensembleOrbitals"] = ensemble
        }
        let rejected = try VivoBarrierConvergence.run(nearDegenerate)
        try c.require(rejected.assessment == .incompleteLevelExecution && rejected.levels.last!.outcome.evaluation != nil,
            "occupation-gap failures remain visible even when full space succeeds")
        return .init(checks:c.labels,observations:["fiveOrbitalProfileErrorHartree":improved.maximumRelativeProfileErrorHartree,
            "fiveOrbitalBarrierErrorHartree":improved.maximumBarrierErrorHartree],
            scope:"reference-assisted H3/6-31G shared-density orbitals; accuracy improved, reduced-level stability not established")
    }
    static func ecc() throws -> Report {
        var c = Check(); let h3 = VivoBarrierBenchmarks.hydrogenExchange631G()
        let basis = VivoGaussianBasis(identifier: h3.basis.identifier,
            shells: h3.basis.shells.filter { $0.nucleusIndex < 2 }, source: h3.basis.source)
        let points = [1.4,1.5,1.6].enumerated().map { i,d in VivoMolecularPathSnapshot(identifier: "p\(i)", coordinate: d,
            system: .init(nuclei: [.init(atomicNumber: 1, positionBohr: .init(0,0,-d/2), structureAtomIndex: 0),
                                  .init(atomicNumber: 1, positionBohr: .init(0,0,d/2), structureAtomIndex: 1)], alphaElectrons: 1, betaElectrons: 1)) }
        let levels = (1...3).map { bath -> VivoBarrierLevel in
            let cfg = VivoECCDMETConfiguration(mode: .singleFragment,
                fragments: [.init(identifier: "F", orbitals: [0], maximumBathOrbitals: bath, clusterAlphaElectrons: 1, clusterBetaElectrons: 1)],
                bathSelection: .init(minimumBathOrbitals: bath))
            return .init(identifier: "bath-\(bath)", method: .eccPath(configuration: .init(embedding: cfg,
                transportGroups: [[0,1,2,3]], expectedBathOrbitals: [bath], pointWeights: [1.0/3,1.0/3,1.0/3], maximumPointEvaluations: 3)))
        }
        let request = VivoBarrierConvergenceRequest(identifier: "H2-6-31G-ECC-bath-closure-test", atomIdentifiers: ["left","right"],
            coordinateUnit: "Bohr", snapshots: points, basis: basis, anchorPointIdentifier: "p1", barrierPointIdentifier: "p1",
            transportGroups: [[0,1,2,3]], levels: levels, maximumPointEvaluations: 21)
        let result = try VivoBarrierConvergence.run(request)
        try c.require(result.levels.count == 3, "all ECC bath levels are retained")
        try c.require(result.levels.first!.outcome.evaluation == nil, "fractional inactive population rejects the reduced fixed-electron cluster")
        try c.require(result.levels.last!.outcome.evaluation?.meetsReferenceAccuracy == true, "full-bath ECC still executes after an earlier level fails")
        try c.require(result.assessment == .incompleteLevelExecution && result.acceptedReducedLevelIdentifier == nil, "later ECC success cannot conceal earlier numerical failure")
        try c.require(result.chargedPointEvaluations == 21, "ECC execution and reconstruction budgets are reserved even for failed levels")
        try VivoBarrierConvergence.validate(result, request: request)
        try c.require(true, "ECC failure and full-bath result reconstruct without a CAS substitute")
        return .init(checks: c.labels, observations: [:], scope: "three physical H2/6-31G ECC bath sizes with explicit particle-closure failures; not a barrier accuracy claim")
    }
}

#if !NUMIVIVO_PORTABLE_CHECKS
@Suite(.serialized) struct BarrierConvergenceTests {
    @Test func physicalBarrierAccuracy() throws { try BarrierConvergenceRegressions.record(BarrierConvergenceRegressions.physical(), "physical") }
    @Test func reportTamperRejection() throws { try BarrierConvergenceRegressions.record(BarrierConvergenceRegressions.tampering(), "tamper") }
    @Test func inputAndAssessmentContracts() throws { try BarrierConvergenceRegressions.record(BarrierConvergenceRegressions.contracts(), "contracts") }
    @Test func ensembleReferenceFrame() throws { try BarrierConvergenceRegressions.record(BarrierConvergenceRegressions.ensemble(), "ensemble") }
    @Test func failedECCLevelsAreRetained() throws { try BarrierConvergenceRegressions.record(BarrierConvergenceRegressions.ecc(), "ecc") }
}
#endif
