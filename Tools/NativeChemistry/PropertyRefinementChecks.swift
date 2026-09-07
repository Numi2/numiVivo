import Foundation
#if NUMIVIVO_SCOPED_NUMERICS
@testable import NumiVivoNumerics
#endif

/// Small deterministic execution checks, not a chemical-accuracy benchmark.
@main struct PropertyRefinementChecks {
    static func main() throws {
        var passed: [String] = []
        func require(_ condition: @autoclosure () -> Bool, _ name: String) throws {
            guard condition() else { throw VivoChemistryError.invalid("refinement regression: \(name)") }
            passed.append(name)
        }
        func rejects(_ name: String, _ body: () throws -> Void) throws {
            do { try body() } catch { passed.append(name); return }
            throw VivoChemistryError.invalid("refinement regression unexpectedly accepted: \(name)")
        }
        func changed<T: Codable>(_ value: T, _ change: (inout [String:Any]) throws -> Void) throws -> T {
            var json = try JSONSerialization.jsonObject(with: JSONEncoder().encode(value)) as! [String:Any]
            try change(&json)
            return try JSONDecoder().decode(T.self,from: JSONSerialization.data(withJSONObject: json,options: [.sortedKeys]))
        }
        let request = try VivoPropertyRefinementExamples.algebraicPath()
        let result = try VivoPropertyDirectedSpace.run(request)
        try require(result.sensitivityEstablishedWithinDeclaredPool,"direct CI establishes declared-pool sensitivity")
        try require(result.finalSpace.active == [0,2],"property-directed selection retains reaction-sensitive block")
        try require(result.rounds[0].selectedBlockIdentifier == "reaction-sensitive","first expansion is the measured reactive block")
        let common = result.rounds[0].trials.first { $0.blockIdentifier == "common-correction" }!
        try require(common.shift.maximumAbsoluteEnergyHartree > 0.01 && abs(common.shift.forwardHartree) < 1e-10,
                    "large absolute common-mode correction is not mistaken for barrier sensitivity")
        try require(result.confirmation?.chosen.assessedAdjacentEdges == 4,"confirmation checks every physical adjacent overlap")
        try require(result.confirmation?.candidateUnion.partition.active == [0,1,2],"joint candidate union is checked even when individual shifts are small")
        try require(result.pointEvaluations <= request.maximumPointEvaluations,"aggregate point budget respected")
        try VivoPropertyDirectedSpace.validate(result,request: request)
        try require(true,"bounded replay reconstructs result")
        let selectedRequest = try VivoPropertyRefinementExamples.algebraicPath(selectedCI: true)
        let selected = try VivoPropertyDirectedSpace.run(selectedRequest)
        try require(selected.sensitivityEstablishedWithinDeclaredPool && selected.finalSpace == result.finalSpace,"selected-CI and direct-CI refinement agree")
        for (a,b) in zip(result.confirmation!.chosen.points,selected.confirmation!.chosen.points) {
            try require(abs(a.energyHartree-b.energyHartree) < 1e-10,"selected/direct energy \(a.pointIdentifier)")
        }
        let capped = try changed(request) { $0["maximumPointEvaluations"] = 1 }
        let capResult = try VivoPropertyDirectedSpace.run(capped)
        try require(!capResult.sensitivityEstablishedWithinDeclaredPool && capResult.termination == .incompleteExecution,
                    "budget exhaustion does not establish low sensitivity")
        let unionCapped = try changed(request) { $0["maximumActiveOrbitals"] = 2 }
        let unionResult = try VivoPropertyDirectedSpace.run(unionCapped)
        try require(!unionResult.sensitivityEstablishedWithinDeclaredPool && unionResult.termination == .activeSpaceLimit,
                    "infeasible joint candidate check blocks acceptance")
        let holdout = try changed(request) { json in
            var points = json["points"] as! [[String:Any]]
            for i in [1,3] {
                var h = points[i]["hamiltonian"] as! [String:Any]
                var matrix = h["oneElectron"] as! [String:Any]
                var values = matrix["values"] as! [Double]
                values[1] = 0.8; values[3] = 0.8
                matrix["values"] = values; h["oneElectron"] = matrix; points[i]["hamiltonian"] = h
            }
            json["points"] = points; json["minimumStateOverlapSquared"] = 0.25
        }
        let heldoutResult = try VivoPropertyDirectedSpace.run(holdout)
        try require(heldoutResult.termination == .confirmationFailed,"held-out electronic sensitivity fails confirmation")
        try require(heldoutResult.finalSpace.active == [0,2],"failed holdout is not reused to tune the candidate space")
        let forged = try changed(result) { $0["finalSpace"] = ["active":[0],"doublyOccupiedCore":[],"frozenOrbitals":[]] }
        try rejects("forged active-space evidence rejects") { try VivoPropertyDirectedSpace.validate(forged,request: request) }
        let split = try changed(request) { $0["transportGroups"] = [[0],[1,2]] }
        try rejects("splitting a transported orbital group rejects") { try split.validate() }
        let overlap = try changed(request) { $0["confirmationPointIdentifiers"] = ["reactant"] }
        try rejects("overlapping discovery and holdout identities reject") { try overlap.validate() }
        let state = VivoCIState(orbitalCount: 3,alphaElectrons: 1,betaElectrons: 0,
                               determinants: [1,4],coefficients: [1/sqrt(2),1/sqrt(2)])
        let information = try VivoSelectiveOrbitalInformation.analyze(state,
            selection: .init(orbitals: [0,1],pairs: [.init(0,1)]))
        try require(abs(information.singles[0].entropyNats-log(2)) < 1e-12,"single-orbital entropy uses natural logarithms")
        try require(abs(information.pairs[0].mutualInformationNats-2*log(2)) < 1e-12,"mutual-information convention has no half factor")
        let boundary = try VivoSelectiveOrbitalInformation.boundary(information,inside: [0],outside: [1,2])
        try require(!boundary.complete && boundary.measuredPairs == 1 && boundary.possiblePairs == 2,"missing pair coverage is not reported as zero")
        let bad = try changed(information) { $0["completePairCoverage"] = true }
        try rejects("forged orbital-information coverage rejects") {
            try VivoSelectiveOrbitalInformation.validate(bad,state: state,selection: information.selection)
        }
        let multiRequest = try VivoPropertyRefinementExamples.multistateCancellation()
        let multi = try VivoPropertyDirectedSpace.run(multiRequest)
        try require(multi.sensitivityEstablishedWithinDeclaredPool && multi.finalSpace.active == [0,1,2,3],
                    "multistate refinement expands despite average-energy cancellation")
        let shifts = multi.rounds[0].trials[0].shift.rootShifts!
        try require(abs(shifts[0].forwardHartree + shifts[1].forwardHartree) < 1e-10 &&
                    abs(shifts[0].forwardHartree) > 0.01,"opposite state errors do not cancel the acceptance criterion")
        try require(multi.rounds[0].trials[0].shift.maximumStateGapShiftHartree! > 0.01,
                    "every pairwise state gap participates in qualification")
        var mixing = VivoQMMatrix(2,2)
        mixing[0,0] = 1/sqrt(2); mixing[0,1] = -1/sqrt(2)
        mixing[1,0] = 1/sqrt(2); mixing[1,1] = 1/sqrt(2)
        let grouped = VivoRefinementStatePolicy(labels: ["a","b"],groups: [[0,1]])
        let retained = try grouped.minimumRetainedOverlapSquared(mixing,threshold: 0.9)
        try require(retained > 0.999999,"degenerate state-subspace rotations preserve physical overlap")
        try rejects("separate state identities cannot silently mix") {
            _ = try VivoRefinementStatePolicy(labels: ["a","b"],groups: [[0],[1]])
                .minimumRetainedOverlapSquared(mixing,threshold: 0.9)
        }
        try rejects("near-degenerate roots across separate groups reject") {
            try VivoRefinementStatePolicy(labels: ["a","b"],groups: [[0],[1]]).validateSpectrum([0,1e-12])
        }
        let saRequest = try changed(multiRequest) { json in
            var cfg = VivoMultiStateCASSCFConfiguration(weights: [0.5,0.5],rootLabels: ["lower","upper"],
                optimization: .init(gradientTolerance: 1e-8,energyToleranceHartree: 1e-12),followRoots: false)
            cfg.davidson.residualTolerance = 1e-12
            let solver = VivoSpaceRefinementSolver.stateAveragedCASSCF(configuration: cfg,
                states: .init(labels: ["lower","upper"],groups: [[0],[1]]))
            json["solver"] = try JSONSerialization.jsonObject(with: JSONEncoder().encode(solver))
        }
        let sa = try VivoPropertyDirectedSpace.run(saRequest)
        try require(sa.sensitivityEstablishedWithinDeclaredPool,"state-averaged CASSCF uses shared refinement loop")
        try require(sa.confirmation!.chosen.points.allSatisfy { $0.optimizedOrbitalRotation != nil },
                    "optimized orbital frames survive the refinement result")
        try require(sa.reservedAuxiliaryWork > 0 && sa.hamiltonianOperatorApplications > 0,
                    "optimized refinement reports aggregate solver work")
        let molecularRequest = try VivoPropertyRefinementExamples.molecularHydrogenStretch()
        let preparation = try VivoMolecularSpacePreparation.prepare(molecularRequest)
        try require(preparation.refinementRequest.initialSpace.active == [0,1] &&
                    preparation.refinementRequest.candidateBlocks.isEmpty,"atomic projection produces a complete small molecular seed")
        let molecular = try VivoPropertyDirectedSpace.run(preparation.refinementRequest)
        try require(molecular.sensitivityEstablishedWithinDeclaredPool,"molecular preparation feeds the actual refinement solver")
        try require(preparation.refinementRequest.points[1].overlapWithPrevious != (try? VivoQMMatrix.identity(2)),
                    "cross-geometry molecular overlap is integrated rather than assigned identity")
        let atomH = VivoElement.from(symbol: "H")!
        let atoms = [VivoMolecularAtom(index: 0,name: "H-left",element: atomH),
                     VivoMolecularAtom(index: 1,name: "H-right",element: atomH)]
        let seeds = try VivoReactionAtomSeeds.derive(source: .mappedEndpoints(
            reactant: .init(identifier: "bonded",atoms: atoms,bonds: [.init(atomA: 0,atomB: 1)]),
            product: .init(identifier: "separated",atoms: atoms),additionalAtomIndices: []),
            atomIdentifiers: ["H-left","H-right"],system: molecularRequest.snapshots[0].system)
        try require(seeds.map(\.atomIndex) == [0,1],"mapped bond changes generate mandatory reactive atoms")
        let unknown = try changed(molecularRequest) { $0["targetShellIndices"] = [99] }
        try rejects("unknown atomic target shell rejects") { try unknown.validate() }
        if CommandLine.arguments.count == 2 {
            let directory = URL(fileURLWithPath: CommandLine.arguments[1])
            try FileManager.default.createDirectory(at: directory,withIntermediateDirectories: true)
            let encoder = JSONEncoder(); encoder.outputFormatting = [.sortedKeys,.prettyPrinted]
            try encoder.encode(request).write(to: directory.appendingPathComponent("algebraic.request.json"))
            try encoder.encode(result).write(to: directory.appendingPathComponent("algebraic.result.json"))
            try encoder.encode(selected).write(to: directory.appendingPathComponent("selected.result.json"))
            try encoder.encode(multi).write(to: directory.appendingPathComponent("multistate.result.json"))
            try encoder.encode(sa).write(to: directory.appendingPathComponent("optimized.result.json"))
            try encoder.encode(preparation).write(to: directory.appendingPathComponent("molecular.preparation.json"))
            try encoder.encode(molecular).write(to: directory.appendingPathComponent("molecular.result.json"))
            try encoder.encode(heldoutResult).write(to: directory.appendingPathComponent("holdout-rejection.result.json"))
            try JSONSerialization.data(withJSONObject: ["scope":"portable native refinement subset; algebraic and H2 preparation fixtures",
                "passed":passed,"assertions":passed.count],options: [.sortedKeys,.prettyPrinted])
                .write(to: directory.appendingPathComponent("results.json"))
        }
        print("\(passed.count) property-refinement assertions passed; algebraic and H2 preparation fixtures; not chemical qualification")
    }
}
