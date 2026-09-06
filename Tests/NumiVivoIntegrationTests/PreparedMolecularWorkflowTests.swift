import Foundation
import Testing
@testable import NumiVivoKit

/// Synthetic contract fixtures, not physical force fields or molecular benchmarks.
@Suite(.serialized) struct PreparedMolecularWorkflowTests {
    private func rejects(_ work:() throws -> Void) {
        do { try work();Issue.record("Expected explicit molecular preparation/sampling rejection") }
        catch { }
    }
    private func fingerprint(_ text:String) throws -> VivoFingerprint {
        try VivoCanonicalJSON.fingerprint(Data(text.utf8))
    }
    private func water() -> VivoMolecularStructure {
        let atoms:[VivoMolecularAtom]=[
            .init(index:0,name:"O",element:.init(atomicNumber:8,symbol:"O"),residueIndex:0,alternateLocation:"A"),
            .init(index:1,name:"H1",element:.init(atomicNumber:1,symbol:"H"),residueIndex:0,alternateLocation:"A"),
            .init(index:2,name:"H2",element:.init(atomicNumber:1,symbol:"H"),residueIndex:0,alternateLocation:"A")]
        let p:[VivoVector3D]=[.zero,.init(0.096,0,0),.init(-0.025,0.093,0)]
        return .init(identifier:"synthetic-water-preparation",atoms:atoms,
            bonds:[.init(atomA:0,atomB:1),.init(atomA:0,atomB:2)],
            residues:[.init(index:0,name:"HOH",atomIndices:[0,1,2])],conformers:[
                .init(identifier:"first",positionsNM:p,potentialEnergyKJPerMol:42),
                .init(identifier:"translated",positionsNM:p.map{$0 + .init(1,2,3)},potentialEnergyKJPerMol:43)])
    }
    @Test func protonationPreservesMapsAndEditsEveryConformer() throws {
        let source=water()
        let request=VivoMolecularPreparationRequest(structure:source,microstateIdentifier:"H3O+",
            protonationSourceIdentifier:"explicit synthetic test",pH:1,expectedFormalCharge:1,
            addHydrogens:[.init(identifier:"added-proton",name:"H3",parent:0,angleReference:1,planeReference:2,
                lengthNM:0.1,angleRadians:1.9,dihedralRadians:1.2)],chargeEdits:[.init(sourceAtom:0,formalCharge:1)],
            residueEdits:[.init(sourceResidue:0,name:"H3O")])
        let result=try VivoMolecularPreparation.prepare(request)
        #expect(result.sourceToPrepared == [0,1,2])
        #expect(result.preparedToSource == [0,1,2,nil])
        #expect(result.hydrogenIndices["added-proton"] == 3)
        #expect(result.structure.atoms[0].formalCharge == 1)
        #expect(result.structure.conformers.allSatisfy{$0.positionsNM.count==4 && $0.potentialEnergyKJPerMol==nil})
        #expect(result.requiresCoordinateRelaxation)
        #expect(result.structure.atoms[0].alternateLocation == "A")
        let d=result.structure.conformers[1].positionsNM[3]-result.structure.conformers[0].positionsNM[3]
        #expect((d - .init(1,2,3)).norm < 1e-12)
        #expect(source.conformers[0].potentialEnergyKJPerMol == 42)
        rejects { var invalid=request;invalid.expectedFormalCharge=0;_ = try VivoMolecularPreparation.prepare(invalid) }
        rejects { var invalid=request;invalid.removeHydrogens=[2];_ = try VivoMolecularPreparation.prepare(invalid) }
    }
    @Test func invertedTetrahedralConformerDoesNotPassPreparation() throws {
        let vectors:[VivoVector3D]=[.init(0.1,0.1,0.1),.init(-0.1,-0.1,0.1),.init(-0.1,0.1,-0.1),.init(0.1,-0.1,-0.1)]
        let positive=(vectors[0]-vectors[3]).dot((vectors[1]-vectors[3]).cross(vectors[2]-vectors[3]))>0
        let atoms: [VivoMolecularAtom] = (0..<5).map { i in
            let element = VivoElement(atomicNumber: i == 0 ? UInt16(6) : UInt16(1), symbol: i == 0 ? "C" : "H")
            return VivoMolecularAtom(index: UInt32(i), name: "atom-\(i)", element: element)
        }
        let structure=VivoMolecularStructure(identifier:"ordered-neighbor-fixture",atoms:atoms,
            bonds:(1..<5).map{.init(atomA:0,atomB:UInt32($0))},conformers:[.init(positionsNM:[.zero]+vectors)])
        var request=VivoMolecularPreparationRequest(structure:structure,microstateIdentifier:"explicit",
            protonationSourceIdentifier:"synthetic tetrahedron",pH:7,expectedFormalCharge:0,
            stereochemistry:[.tetrahedral(center:.source(0),orderedNeighbors:(1..<5).map{.source(UInt32($0))},positive:positive)])
        _ = try VivoMolecularPreparation.prepare(request)
        let mirror = ([VivoVector3D.zero] + vectors).map { VivoVector3D(-$0.x, $0.y, $0.z) }
        request.structure.conformers.append(VivoMolecularConformer(identifier: "mirror", positionsNM: mirror))
        rejects { _ = try VivoMolecularPreparation.prepare(request) }
    }
    @Test func triclinicMinimumImageAgreesWithExplicitLatticeEnumeration() throws {
        let cell=VivoPeriodicCell(a:.init(2,0,0),b:.init(1.8,0.5,0),c:.init(0,0,2))
        let d=VivoVector3D(2.28,0.3,0.1),actual=try cell.minimumImage(d)
        var reference=Double.infinity
        for i in -4...4 { for j in -4...4 { for k in -4...4 {
            reference=min(reference,(d-cell.a*Double(i)-cell.b*Double(j)-cell.c*Double(k)).squaredNorm)
        } } }
        #expect(abs(actual.squaredNorm-reference)<1e-13)
        let rebuilt=try VivoPeriodicCell.crystallographic(lengthsNM:.init(2,2,2),anglesDegrees:.init(90,90,60))
        #expect(rebuilt.isValid)
        try rebuilt.validateSingleImageRadius(0.2)
        rejects { try rebuilt.validateSingleImageRadius(1.1) }
        rejects { _ = try VivoPeriodicCell.crystallographic(lengthsNM:.init(2,2,2),anglesDegrees:.init(10,10,170)) }
    }
    private func improperFixture() -> (VivoMolecularStructure,VivoForceFieldLibrary) {
        let names=["C","H1","H2","H3"]
        let atoms=names.enumerated().map{VivoMolecularAtom(index:UInt32($0.offset),name:$0.element,
            element:.init(atomicNumber:$0.offset==0 ? 6:1,symbol:$0.offset==0 ? "C":"H"),residueIndex:0)}
        let structure=VivoMolecularStructure(identifier:"synthetic-improper",atoms:atoms,
            bonds:(1..<4).map{.init(atomA:0,atomB:UInt32($0))},residues:[.init(index:0,name:"X",atomIndices:[0,1,2,3])],
            conformers:[.init(positionsNM:[.zero,.init(0.1,0,0),.init(0,0.1,0),.init(-0.1,-0.1,0.01)])])
        let template=VivoForceFieldResidueTemplate(residueNames:["X"],
            atoms:names.map{.init(name:$0,typeIdentifier:$0=="C" ? "C":"H",chargeE:0)},
            bonds:names.dropFirst().map{.init(atomA:"C",atomB:$0)},
            impropers:[.init(orderedAtoms:["H1","C","H2","H3"],center:"C")])
        let library=VivoForceFieldLibrary(identifier:"synthetic-unit-test-only",version:"1",
            atomTypes:[.init(identifier:"C",elementAtomicNumber:6,massDa:12,sigmaNM:0.3,epsilonKJPerMol:0.1),
                       .init(identifier:"H",elementAtomicNumber:1,massDa:1,sigmaNM:0.1,epsilonKJPerMol:0.01)],
            residueTemplates:[template],bondParameters:[.init(typeA:"C",typeB:"H",lengthNM:0.1,forceConstant:100)],
            angleParameters:[.init(typeA:"H",typeB:"C",typeC:"H",angleRadians:2.0,forceConstant:20)],
            torsionParameters:[.init(typeA:"H",typeB:"C",typeC:"H",typeD:"H",periodicity:2,
                                    phaseRadians:Double.pi,barrierKJPerMol:0.5,improper:true)],
            provenance:["purpose":"synthetic contract fixture; not a physical force field"])
        return (structure,library)
    }
    @Test func orderedNativeImpropersReachTheExistingClassicalSystem() throws {
        let (structure,library)=improperFixture()
        let assignment=try VivoResidueTemplateAssigner.assignPrepared(structure:structure,library:library)
        let result=try VivoForceFieldCompiler.compile(structure:structure,library:library,assignment:assignment)
        #expect(result.system.torsions.count==1)
        let term=try #require(result.system.torsions.first)
        #expect(term.improper)
        #expect([term.a,term.b,term.c,term.d] == [1,0,2,3])
        rejects { var legacy=library;legacy.residueTemplates[0].impropers=nil
            _ = try VivoForceFieldCompiler.compile(structure:structure,library:legacy) }
        rejects { var wrong=structure;wrong.bonds[0].order = .double
            _ = try VivoResidueTemplateAssigner.assignPrepared(structure:wrong,library:library) }
    }
    @Test func torsionMatcherRejectsCompetingEqualSpecificityFamilies() throws {
        var (_,library)=improperFixture()
        library.torsionParameters=[
            .init(typeA:"*",typeB:"C",typeC:"H",typeD:"H",periodicity:2,phaseRadians:0,barrierKJPerMol:1,improper:true),
            .init(typeA:"H",typeB:"C",typeC:"H",typeD:"*",periodicity:3,phaseRadians:0,barrierKJPerMol:1,improper:true)]
        rejects { _ = try VivoForceFieldCompiler.selectedTorsions(["H","C","H","H"],improper:true,library:library) }
    }
    @Test func constantReplicaTracesCannotCertifyMolecularSampling() throws {
        let fp=try fingerprint("synthetic sampling structure")
        let observable=VivoMolecularObservable(identifier:"distance",kind:.distance(atomA:0,atomB:1),maximumMeanStandardError:0.01)
        var replicas:[VivoMolecularReplicaSeries]=[]
        for i in 0..<4 {
            var md=VivoMDConfiguration();md.randomSeed=UInt64(i+1)
            replicas.append(.init(identifier:"replica-\(i)",sourceFingerprint:try fingerprint("synthetic replica \(i)"),configuration:md,
                steps:(1...256).map{UInt64($0)},timesPS:(1...256).map{Double($0)*0.1},valuesByObservable:[Array(repeating:1,count:256)]))
        }
        let request=VivoMolecularSamplingRequest(structureFingerprint:fp,systemFingerprint:fp,contextIdentifier:"synthetic",
            observables:[observable],replicas:replicas)
        let result=try VivoMolecularSampling.analyze(request)
        #expect(!result.converged)
        #expect(!result.issues.isEmpty)
        #expect(result.diagnostics[0].bulkEffectiveSamples==nil)
        _ = try VivoCanonicalJSON.encode(result)
        try VivoMolecularSampling.validate(result)
        rejects { var duplicate=request;duplicate.replicas[1].configuration.randomSeed=1
            _ = try VivoMolecularSampling.analyze(duplicate) }
    }
    @Test func rectangularRIFactorsMatchTheFullTransform() throws {
        let values=try VivoQMMatrix(rows:6,columns:2,values:[1,0.2,0.4,0.1,0.8,0.3,0.5,-0.1,0.2,0.6,1.2,0.7])
        let factors=VivoCoulombFactors(orbitalCount:3,values:values,method:"synthetic positive-semidefinite factors")
        let left=try VivoQMMatrix(rows:3,columns:1,values:[1,0,0])
        let right=try VivoQMMatrix(rows:3,columns:2,values:[0,0,1,0,0,1])
        let panel=try factors.orbitalPairPanel(left:left,right:right)
        for a in 0..<2 { for k in 0..<2 {
            #expect(abs(panel[a,k]-values[VivoCoulombFactors.pair(0,a+1),k])<1e-12)
        } }
    }
}
