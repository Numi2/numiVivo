import Foundation
#if NUMIVIVO_SCOPED_NUMERICS
import NumiVivoNumerics
#else
import NumiVivoKit
#endif

@main struct ECCPathChecks {
    static func main() throws {
        guard CommandLine.arguments.count==2 else { throw VivoChemistryError.invalid("expected output directory") }
        let out=URL(fileURLWithPath:CommandLine.arguments[1],isDirectory:true)
        try FileManager.default.createDirectory(at:out,withIntermediateDirectories:true)
        var checks:[String]=[], observations:[String:Double]=[:]
        func require(_ value: Bool, _ name: String) throws {
            guard value else { throw VivoChemistryError.invalid("ECC path check: \(name)") }
            checks.append(name); print("PASS \(name)")
        }
        func rejects(_ name: String, _ operation: () throws -> Void) throws {
            var rejected=false
            do { try operation() } catch is VivoChemistryError { rejected=true } catch is VivoECCPathError { rejected=true }
            try require(rejected,name)
        }
        func write<T:Encodable>(_ value: T, _ name: String) throws {
            let encoder=JSONEncoder(); encoder.outputFormatting=[.sortedKeys,.prettyPrinted]
            try encoder.encode(value).write(to:out.appendingPathComponent(name),options:.atomic)
        }
        func mutate<T:Codable>(_ value: T, _ change: (inout [String:Any]) -> Void) throws -> T {
            var object=try JSONSerialization.jsonObject(with:JSONEncoder().encode(value)) as! [String:Any]
            change(&object)
            return try JSONDecoder().decode(T.self,from:JSONSerialization.data(withJSONObject:object))
        }
        let request=try VivoMolecularECCPath.hydrogenStretchTemplate()
        let ao=try request.snapshots.map { try VivoGaussianIntegralEngine.compute(system:$0.system,basis:request.basis,budget:request.budget) }
        let references=try ao.map { try VivoHartreeFock.solve(system:$0.sourceSystem,integrals:$0) }
        let points=try VivoMolecularECCPath.prepare(request,integrals:ao,references:references)
        let cfg=request.configuration
        let result=try VivoECCReactionPath.solve(points:points,configuration:cfg,initialSharedRotation:request.initialSharedRotation)
        try require(result.converged,"five physical H2 geometries reach shared QIO stationarity")
        try require(result.iterations.count>2 && result.iterations.last!.objective<result.iterations.first!.objective-1e-3,
                    "shared rotation performs nontrivial objective-reducing updates")
        try require(result.pointResults.allSatisfy { $0.configuration.localityGroups.isEmpty && $0.orbitalGradient.isEmpty },
                    "independent pointwise orbital optimization is disabled")
        try require(result.continuity.count==4,"every adjacent geometry has a continuity record")
        try require(result.continuity.allSatisfy { $0.referenceStateOverlapSquared>=cfg.minimumReferenceOverlapSquared },
                    "physical determinant overlaps preserve reference continuity")
        try require(result.pointResults.allSatisfy { abs($0.frame.electronCountDefect)<=cfg.embedding.electronTolerance },
                    "electron closure holds at every geometry")
        try VivoECCReactionPath.validate(result,points:points,configuration:cfg)
        try require(true,"final profile reconstructs and passes independent shared-gradient validation")
        let exact=try points.map { try VivoDirectCI.solve($0.hamiltonian).roots[0].energyHartree }
        let error=zip(exact,result.pointResults).map { abs($0-$1.energyHartree) }.max()!
        try require(error<1e-8,"all full-bath ECC energies match separate native FCI calculations")
        observations["maximumNativeFCIErrorHartree"]=error
        observations["objectiveBefore"]=result.iterations.first!.objective
        observations["objectiveAfter"]=result.objective
        observations["sharedGradientNorm"]=result.sharedGradient.reduce(0.0) { hypot($0,$1) }
        observations["pointEvaluations"]=Double(result.pointEvaluations)
        // Alternate signs in raw SCF columns must not affect the physical result.
        let gauges=try points.indices.map { i in try VivoQMMatrix(rows:2,columns:2,values:[1,0,0,i%2==0 ? 1:-1]) }
        let rephased=try points.indices.map { i in
            VivoECCPathPoint(identifier:points[i].identifier,hamiltonian:try points[i].hamiltonian.rotated(by:gauges[i]),
                overlapWithPrevious:i==0 ? nil : try gauges[i-1].transposed.multiplied(by:points[i].overlapWithPrevious!).multiplied(by:gauges[i]))
        }
        let gaugeResult=try VivoECCReactionPath.solve(points:rephased,configuration:cfg,initialSharedRotation:request.initialSharedRotation)
        try require(gaugeResult.converged && abs(gaugeResult.objective-result.objective)<1e-8,"SCF phase changes preserve the shared objective")
        try require(zip(gaugeResult.relativeElectronicEnergiesHartree,result.relativeElectronicEnergiesHartree).allSatisfy { abs($0-$1)<1e-8 },
                    "SCF phase changes preserve the entire relative-energy profile")
        try rejects("non-normalized geometry weights are rejected") {
            var c=cfg; c.pointWeights[0]+=0.1; _=try VivoECCReactionPath.solve(points:points,configuration:c)
        }
        try rejects("missing physical cross overlap is rejected") {
            var p=points; p[1] = .init(identifier:p[1].identifier,hamiltonian:p[1].hamiltonian)
            _=try VivoECCReactionPath.solve(points:p,configuration:cfg)
        }
        try rejects("an overlap larger than an orthonormal contraction is rejected") {
            var p=points; p[1] = .init(identifier:p[1].identifier,hamiltonian:p[1].hamiltonian,
                overlapWithPrevious:try .init(rows:2,columns:2,values:[2,0,0,2]))
            _=try VivoECCReactionPath.solve(points:p,configuration:cfg)
        }
        try rejects("lost cross-geometry orbital subspace is rejected") {
            var p=points; p[1] = .init(identifier:p[1].identifier,hamiltonian:p[1].hamiltonian,
                overlapWithPrevious:try .init(rows:2,columns:2,values:[0.01,0,0,0.01]))
            _=try VivoECCReactionPath.solve(points:p,configuration:cfg)
        }
        try rejects("a different constructed bath size is not silently accepted") {
            var c=cfg; c.expectedBathOrbitals=[0]; _=try VivoECCReactionPath.solve(points:points,configuration:c)
        }
        try rejects("same dimension does not bypass physical subspace continuity") {
            var c=cfg; c.minimumPathSubspaceSingularValue=1
            _=try VivoECCReactionPath.solve(points:points,configuration:c)
        }
        var limited=cfg; limited.maximumPointEvaluations=points.count
        let unfinished=try VivoECCReactionPath.solve(points:points,configuration:limited,initialSharedRotation:request.initialSharedRotation)
        try require(!unfinished.converged && unfinished.termination=="point-evaluation-limit","aggregate point-evaluation exhaustion is not convergence")
        try rejects("unfinished profile cannot pass output validation") {
            try VivoECCReactionPath.validate(unfinished,points:points,configuration:limited)
        }
        try rejects("forged relative energy is rejected") {
            let tampered=try mutate(result) { $0["relativeElectronicEnergiesHartree"]=[Double](repeating:99,count:points.count) }
            try VivoECCReactionPath.validate(tampered,points:points,configuration:cfg)
        }
        try rejects("forged objective is rejected") {
            let tampered=try mutate(result) { $0["objective"]=99 }
            try VivoECCReactionPath.validate(tampered,points:points,configuration:cfg)
        }
        try rejects("pointwise convergence does not authorize a different shared frame") {
            let altered=try JSONSerialization.jsonObject(with:JSONEncoder().encode(request.initialSharedRotation!))
            let tampered=try mutate(result) { $0["sharedRotation"]=altered; $0["sharedGradient"]=[0.0] }
            try VivoECCReactionPath.validate(tampered,points:points,configuration:cfg)
        }
        try rejects("nuclear mapping changes are rejected before preparation") {
            let changed=try mutate(request) { object in
                var snapshots=object["snapshots"] as! [[String:Any]]
                var system=snapshots[1]["system"] as! [String:Any]
                var nuclei=system["nuclei"] as! [[String:Any]]
                nuclei[0]["structureAtomIndex"]=23; system["nuclei"]=nuclei; snapshots[1]["system"]=system; object["snapshots"]=snapshots
            }; try changed.validate()
        }
        try rejects("unordered reaction coordinates are rejected") {
            let changed=try mutate(request) { object in
                var snapshots=object["snapshots"] as! [[String:Any]]; snapshots[1]["coordinate"]=snapshots[0]["coordinate"]; object["snapshots"]=snapshots
            }; try changed.validate()
        }
        try rejects("changing electron population is rejected") {
            let changed=try mutate(request) { object in
                var snapshots=object["snapshots"] as! [[String:Any]]
                var system=snapshots[1]["system"] as! [String:Any]; system["alphaElectrons"]=0; snapshots[1]["system"]=system; object["snapshots"]=snapshots
            }; try changed.validate()
        }
        try rejects("changing frozen embedding charges is rejected") {
            let changed=try mutate(request) { object in
                var snapshots=object["snapshots"] as! [[String:Any]]
                var system=snapshots[1]["system"] as! [String:Any]
                system["pointCharges"]=[["chargeE":0.1,"positionBohr":[5.0,0.0,0.0]]]; snapshots[1]["system"]=system; object["snapshots"]=snapshots
            }; try changed.validate()
        }
        // Same basis and same bath dimension, but a changed ground state.
        func twoLevel(_ sign: Double) throws -> VivoEmbeddedHamiltonian {
            .init(orbitalIdentifiers:["a","b"],alphaElectrons:1,betaElectrons:1,
                oneElectron:try .init(rows:2,columns:2,values:[-sign,0,0,sign]),twoElectron:[Double](repeating:0,count:16),
                constantEnergyHartree:0,energyReference:"synthetic level crossing")
        }
        var crossing=cfg; crossing.pointWeights=[0.5,0.5]; crossing.embedding.localityGroups=[]
        try rejects("reference-state discontinuity is rejected despite identical bath dimensions") {
            _=try VivoECCReactionPath.solve(points:[.init(identifier:"left",hamiltonian:twoLevel(1)),
                .init(identifier:"right",hamiltonian:twoLevel(-1),overlapWithPrevious:.identity(2))],configuration:crossing)
        }
        try write(request,"request.json"); try write(points,"points.json"); try write(result,"result.json")
        try write(gaugeResult,"rephased-result.json")
        struct Report: Codable { let checks:[String]; let observations:[String:Double]; let passed:Bool; let scope:String }
        try write(Report(checks:checks,observations:observations,passed:true,
            scope:"five-geometry H2/STO-3G full-bath shared ECC limit and explicit rejection tests; no activation barrier or paper reaction"),"checks.json")
        print("PASS \(checks.count) ECC path checks")
    }
}
