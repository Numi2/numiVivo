import Foundation
import CryptoKit
import Testing
@testable import NumiVivoKit

@Suite(.serialized) struct ScientificClosureTests {
    private func record<T:Encodable>(_ value:T,_ name:String) throws {
        guard let folder=ProcessInfo.processInfo.environment["NUMIVIVO_TEST_ARTIFACTS"] else {return}
        let root=URL(fileURLWithPath:folder,isDirectory:true)
        try FileManager.default.createDirectory(at:root,withIntermediateDirectories:true)
        try VivoCanonicalJSON.encode(value).write(to:root.appendingPathComponent(name+".json"),options:.atomic)
    }
    private func rejects(_ f:() throws -> Void) {
        do {try f();Issue.record("expected explicit rejection") } catch {}
    }
    private func mutate<T:Codable>(_ value:T,_ edit:(inout [String:Any])->Void) throws -> T {
        var o=try JSONSerialization.jsonObject(with:JSONEncoder().encode(value)) as! [String:Any]
        edit(&o);return try JSONDecoder().decode(T.self,from:JSONSerialization.data(withJSONObject:o))
    }
    @Test func sharedResidualHierarchyIsActuallyReduced() throws {
        let request=VivoResidualBarrierCampaign.template(),result=try VivoResidualBarrierCampaign.run(request)
        #expect(result.baseline.assessment == .reducedAccuracyNotEstablished)
        #expect(result.reducedAccuracyEstablished && result.acceptedRound==6)
        #expect(result.levels.map(\.variationalDimension)==[9,18,27,36,45,54,63])
        #expect(result.levels.allSatisfy { $0.fullSectorDimension==90 && $0.stateSpaceIsReduced })
        #expect(result.levels.suffix(2).allSatisfy { $0.meetsReferenceAccuracy })
        let successiveTolerance=request.baseline.acceptance.maximumSuccessiveChangeHartree
        #expect(result.levels.suffix(2).map(\.variationalDimension)==[54,63])
        #expect(result.levels[2].meetsReferenceAccuracy)
        #expect(result.levels[3].meetsReferenceAccuracy)
        #expect((result.levels[2].maximumSuccessiveChangeHartree ?? .infinity)>successiveTolerance)
        #expect((result.levels[3].maximumSuccessiveChangeHartree ?? .infinity)<=successiveTolerance)
        #expect(result.levels.last!.maximumProfileErrorHartree<1e-7)
        #expect((result.levels.last!.maximumSuccessiveChangeHartree ?? .infinity)<successiveTolerance)
        #expect(result.levels.allSatisfy { $0.points.allSatisfy { $0.state.orbitalCount==6 } })
        #expect(result.hamiltonianOperatorApplications<=request.baseline.budget.maximumOperatorApplications)
        for (a,b) in zip(result.levels,result.levels.dropFirst()) {
            #expect(b.variationalDimension>a.variationalDimension)
            #expect(zip(a.points,b.points).allSatisfy { $1.energyHartree <= $0.energyHartree+1e-10 })
        }
        try VivoResidualBarrierCampaign.validate(result,request:request)
        try record(result,"shared-residual-barrier")
    }
    @Test func residualFailureAndForgeryCannotBecomeSuccess() throws {
        let original=VivoResidualBarrierCampaign.template()
        let short=VivoResidualBarrierRequest(baseline:original.baseline,seed:original.seed,maximumRefinementRounds:1,space:original.space)
        let result=try VivoResidualBarrierCampaign.run(short)
        #expect(!result.reducedAccuracyEstablished && result.acceptedRound==nil)
        let bad=try mutate(result) { $0["reducedAccuracyEstablished"]=true;$0["acceptedRound"]=1 }
        rejects { try VivoResidualBarrierCampaign.validate(bad,request:short) }
        let tooSmall=VivoResidualBarrierRequest(baseline:original.baseline,seed:original.seed,maximumRefinementRounds:2,
            space:.init(maximumDimension:10))
        rejects { _=try VivoResidualBarrierCampaign.run(tooSmall) }
        let tinyBudget=try mutate(original) { o in
            var b=o["baseline"] as! [String:Any],budget=b["budget"] as! [String:Any]
            budget["maximumOperatorApplications"]=1;b["budget"]=budget;o["baseline"]=b
        }
        rejects { _=try VivoResidualBarrierCampaign.run(tinyBudget) }
    }
    @Test func overlappingProjectorsHaveOnePhysicalDensity() throws {
        let request=VivoVariationalEmbedding.template(),result=try VivoVariationalEmbedding.run(request)
        #expect(result.variationalDimension==14 && result.fullSectorDimension==16)
        #expect(result.redundantSeedColumns==4)
        #expect(abs(result.occupations.reduce(0,+)-2)<1e-10)
        #expect(result.occupations.allSatisfy { $0 >= -1e-10 && $0 <= 2+1e-10 })
        #expect(abs(result.state.coefficients.reduce(0) { $0+$1*$1 }-1)<1e-10)
        #expect(result.selfConsistentProjectedResidualHartree<request.stationarityToleranceHartree)
        // This symmetric H2 projector can contain the exact ground state even
        // while omitting two determinants. Roundoff residual components need
        // not have a prescribed ordering; use a genuinely truncated polar case.
        #expect(result.externalResidualHartree.isFinite && result.externalResidualHartree>=0)
        #expect(result.externalResidualHartree<1e-8)
        let root=URL(fileURLWithPath:#filePath).deletingLastPathComponent().deletingLastPathComponent().deletingLastPathComponent()
        let input=try VivoCanonicalJSON.decode(VivoReactionCalculationRequest.self,
            from:Data(contentsOf:root.appendingPathComponent("Examples/reaction/global-lih-overlap.json")))
        guard case .globalEmbedding(let polarRequest)=input.calculation else {throw VivoChemistryError.invalid("polar test fixture method")}
        let polar=try VivoVariationalEmbedding.run(polarRequest)
        #expect(polar.variationalDimension==14 && polar.fullSectorDimension==225)
        #expect(polar.externalResidualHartree>1e-4)
        #expect(polar.selfConsistentProjectedResidualHartree<=polarRequest.stationarityToleranceHartree)
        try record(polar,"global-polar-overlap-solvent")
        var fragments=request.fragments
        fragments.append(.init(identifier:"duplicate-projector",partition:fragments[0].partition))
        let duplicate=try VivoVariationalEmbedding.run(.init(molecule:request.molecule,fragments:fragments,space:request.space))
        #expect(duplicate.variationalDimension==result.variationalDimension)
        #expect(abs(duplicate.energyHartree-result.energyHartree)<1e-11)
        #expect(try duplicate.densityAO.adding(result.densityAO,scale:-1).frobeniusNorm<1e-10)
        let ao=try VivoGaussianIntegralEngine.compute(system:request.molecule.system,basis:request.molecule.basis)
        let h=try VivoEmbeddedHamiltonian.fromAO(ao,coefficients:result.coefficients,alphaElectrons:1,betaElectrons:1,
            orbitalIdentifiers:(0..<ao.count).map { "check-\($0)" },energyReference:"physical gas energy")
        let rdms=try VivoCIDensityMatrices.compute(result.state)
        #expect(abs(try rdms.energy(of:h)-result.gasEnergyHartree)<1e-9)
        #expect(abs(result.gasEnergyHartree+result.equilibriumField.polarizationEnergyHartree-result.energyHartree)<1e-12)
        try VivoVariationalEmbedding.validate(result,request:request)
        let forged=try mutate(result) { o in
            var matrix=o["densityAO"] as! [String:Any]
            matrix["values"]=(matrix["values"] as! [Double]).map {2*$0};o["densityAO"]=matrix
        }
        rejects { try VivoVariationalEmbedding.validate(forged,request:request) }
        try record(result,"global-overlap-solvent")
    }
    @Test func rankLossAndOrbitalFramesAreExplicit() throws {
        let request=VivoVariationalEmbedding.template()
        let full=VivoVariationalEmbeddingRequest(molecule:request.molecule,
            fragments:[.init(identifier:"full",partition:.init(active:[0,1,2,3]))],space:request.space)
        let result=try VivoVariationalEmbedding.run(full),fci=try VivoCorrelatedSolvation.solve(request.molecule)
        #expect(result.variationalDimension==16)
        #expect(abs(result.energyHartree-fci.energyHartree)<1e-9)
        #expect(try result.densityAO.adding(fci.totalDensityAO,scale:-1).frobeniusNorm<1e-8)
        let invalid=VivoVariationalEmbeddingRequest(molecule:request.molecule,
            fragments:[.init(identifier:"bad",partition:.init(active:[0,1]),orbitalRotation:try .identity(3))])
        rejects { _=try VivoVariationalEmbedding.run(invalid) }
        rejects { _=try VivoVariationalEmbedding.run(.init(molecule:request.molecule,fragments:request.fragments,space:.init(maximumDimension:1))) }
        let wrongCore=VivoVariationalEmbeddingRequest(molecule:request.molecule,
            fragments:[.init(identifier:"bad-core",partition:.init(doublyOccupiedCore:[0],active:[0,1]))])
        rejects { _=try VivoVariationalEmbedding.run(wrongCore) }
        let partial=try VivoVariationalEmbedding.run(request)
        var generator=VivoQMMatrix(4,4)
        generator[0,3]=0.27;generator[3,0] = -0.27;generator[1,2]=0.19;generator[2,1] = -0.19
        let rotation=try VivoQMDenseAlgebra.orbitalRotation(generator:generator)
        let rotatedMolecule=VivoCorrelatedSolventRequest(system:request.molecule.system,basis:request.molecule.basis,
            solvent:request.molecule.solvent,configuration:request.molecule.configuration,
            coefficients:try partial.coefficients.multiplied(by:rotation),budget:request.molecule.budget)
        let transformed=request.fragments.map {VivoFockFragment(identifier:$0.identifier,partition:$0.partition,orbitalRotation:rotation.transposed)}
        let covariant=try VivoVariationalEmbedding.run(.init(molecule:rotatedMolecule,fragments:transformed,space:request.space))
        #expect(covariant.variationalDimension==partial.variationalDimension)
        #expect(abs(covariant.energyHartree-partial.energyHartree)<1e-9)
        #expect(try covariant.densityAO.adding(partial.densityAO,scale:-1).frobeniusNorm<1e-8)
        try record(result,"global-full-space-control")
    }
    @Test func validatedECCFramesCanEnterGlobalFunctional() throws {
        let system=VivoElectronicSystem(nuclei:[.init(atomicNumber:1,positionBohr:.init(0,0,-0.7)),
            .init(atomicNumber:1,positionBohr:.init(0,0,0.7))],alphaElectrons:1,betaElectrons:1)
        let basis=VivoGaussianBasis.hydrogenSTO3G(nucleusIndices:[0,1]),ao=try VivoGaussianIntegralEngine.compute(system:system,basis:basis)
        let eig=try VivoQMDenseAlgebra.symmetricEigen(ao.overlap)
        var c=eig.vectors;for i in 0..<2 {for j in 0..<2 {c[i,j]/=sqrt(eig.values[j])}}
        let h=try VivoEmbeddedHamiltonian.fromAO(ao,coefficients:c,alphaElectrons:1,betaElectrons:1,
            orbitalIdentifiers:["a","b"],energyReference:"physical electronic energy; AO scalar once; equilibrium solvent polarization accounted separately")
        let ecc=try VivoECCDMET.solve(h,configuration:.init(mode:.singleFragment,fragments:[.init(identifier:"F",orbitals:[0],maximumBathOrbitals:1,clusterAlphaElectrons:1,clusterBetaElectrons:1)],bathSelection:.init(minimumBathOrbitals:1)))
        let molecular=VivoCorrelatedSolventRequest(system:system,basis:basis,solvent:.init(dielectricConstant:4,angularPoints:50),coefficients:c)
        let request=try VivoVariationalEmbedding.fromECC(molecule:molecular,sourceHamiltonian:h,result:ecc,inactiveDoublyOccupiedColumns:[[]])
        let result=try VivoVariationalEmbedding.run(request),reference=try VivoCorrelatedSolvation.solve(molecular)
        #expect(abs(result.energyHartree-reference.energyHartree)<1e-9)
        rejects { _=try VivoVariationalEmbedding.fromECC(molecule:molecular,sourceHamiltonian:h,result:ecc,inactiveDoublyOccupiedColumns:[[1]]) }
        #expect(result.method != ecc.frame.energyConvention)
        try record(result,"ecc-global-import")
    }
    @Test func paperInputAbsenceAndIntegrityAreNotInvented() throws {
        let missing=VivoReproductionPackage(target:.acrylamideMethanethiolate,resultIdentifier:"table-to-reproduce")
        let report=try VivoReproductionPreflight.inspect(missing)
        #expect(report.missingRoles.count==5 && !report.suppliedInputPackageConsistent && !report.reproductionExecuted)
        let btk=try VivoReproductionPreflight.inspect(.init(target:.btkSnapshot,resultIdentifier:"snapshot-to-reproduce"))
        #expect(btk.missingRoles.count==8 && !btk.sourceAuthenticityVerified)
        let data=Data("{}".utf8)
        let malformed=VivoReproductionAsset(role:.precomplex,origin:.authorSupplied,sourceIdentifier:"unverified-source",
            sha256:String(repeating:"0",count:64),payload:data)
        rejects { _=try VivoReproductionPreflight.inspect(.init(target:.acrylamideMethanethiolate,resultIdentifier:"table",assets:[malformed])) }
        let req=VivoVariationalEmbedding.template().molecule
        let geometry=VivoReproductionGeometry(atomIdentifiers:["H0","H1"],system:req.system)
        let bytes=try JSONEncoder().encode(geometry),hash=SHA256.hash(data:bytes).map {String(format:"%02x",$0)}.joined()
        let independent=VivoReproductionAsset(role:.precomplex,origin:.independentReconstruction,sourceIdentifier:"synthetic-conformance-not-paper",sha256:hash,payload:bytes)
        let rejectedClaim=try VivoReproductionPreflight.inspect(.init(target:.acrylamideMethanethiolate,resultIdentifier:"table",assets:[independent]))
        #expect(rejectedClaim.blockers.count==1 && !rejectedClaim.suppliedInputPackageConsistent)
        #expect(rejectedClaim.verifiedContentDigests["precomplex"]==hash)
        try record(report,"missing-michael-inputs");try record(btk,"missing-btk-inputs")
    }
}
