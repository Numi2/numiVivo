import Foundation
import Testing
@testable import NumiVivoKit

@Suite(.serialized) struct MDBondedForceConformanceTests {
    /// OpenMM 8.6 Reference, independently checked against central energy
    /// differences at 1e-5 nm. Generator: Tools/Benchmarks/prepare_bonded_references.py.
    private let forces:[[VivoVector3D]] = [
        [.init(0,-11.295975472359858,11.295975472359858),.init(-7.790327911972316,20.25485257112802,-4.67419674718339),
         .init(17.138721406339094,-15.191139428346016,-19.086303384332172),.init(-9.348393494366778,6.2322623295778525,12.464524659155705)],
        [.init(0,-14.713812418395007,14.713812418395007),.init(-10.147456840272419,26.383387784708287,-6.088474104163451),
         .init(22.32440504859932,-19.787540838531214,-24.861269258667424),.init(-12.176948208326902,8.117965472217934,16.235930944435868)],
        [.init(0,11.295975472359862,-11.295975472359862),.init(7.7903279119723186,-20.254852571128026,4.674196747183389),
         .init(-17.1387214063391,15.191139428346021,19.086303384332183),.init(9.348393494366782,-6.232262329577855,-12.46452465915571)]
    ]
    private func system(_ count:Int) throws -> VivoClassicalSystem {
        .init(identifier:"bonded-reference",structureFingerprint:try VivoCanonicalJSON.fingerprint(Data("openmm-8.6-bonded-reference".utf8)),
            particles:(0..<count).map { .init(index:UInt32($0),atomIndex:UInt32($0),typeIdentifier:"C",massDa:12,chargeE:0,sigmaNM:0,epsilonKJPerMol:0) })
    }
    private func configuration(iterations:UInt32=128) -> VivoMDConfiguration {
        .init(timeStepPS:0.001,cutoffNM:10,electrostatics:.cutoff,ensemble:.nve,thermostat:.none,
            targetTemperatureK:nil,frictionPerPS:nil,maximumConstraintIterations:iterations,neighborListEnabled:false)
    }
    @Test(arguments:[0,1,2]) func torsionForcesAgreeWithIndependentReference(_ index:Int) async throws {
        var system=try system(4)
        system.torsions=[.init(a:0,b:1,c:2,d:3,periodicity:3,phaseRadians:[0,0.37,Double.pi][index],barrierKJPerMol:2.2)]
        let geometry=try VivoMDCandidateGeometry(particlePositionsNM:[.init(0,0.25,0.125),.init(0.25,0.125,0),.init(0.5,0.25,0.125),.init(0.75,0.125,0.375)],periodicCell:nil)
        let ref=VivoMDBenchmarkReference(identifier:"torsion",geometry:geometry,
            energyKJPerMol:[3.7738660711048126,3.111487633408774,0.626133928895188][index],forcesKJPerMolNM:forces[index])
        let request=VivoMDBenchmarkRequest(identifier:"torsion",system:system,configuration:configuration(),references:[ref],referenceProvenance:["engine":"OpenMM 8.6 Reference"])
        let report=try await VivoMDBenchmark.run(request)
        #expect(report.outcome == .passed)
        #expect(report.comparisons[0].forceNormalizedRMS<1e-5)
    }
    @Test func explicitProjectionPreservesClockAndFailedProjectionPreservesState() async throws {
        var model=try system(3)
        model.constraints=[.init(a:0,b:1,distanceNM:0.125),.init(a:1,b:2,distanceNM:0.125),.init(a:0,b:2,distanceNM:0.1875)]
        let positions:[VivoVector3D]=[.zero,.init(0.15625,0,0),.init(0.15625,0.15625,0)]
        let initial=VivoClassicalInitialState(systemFingerprint:try model.fingerprint(),positionsNM:positions)
        let good=try await VivoMDMetalRuntime.make(system:model,initialState:initial,configuration:configuration())
        let before=try await good.checkpoint(),after=try await good.projectConstraints()
        #expect(after.acceptedStep==before.acceptedStep && after.timePS==before.timePS)
        #expect(after.positionsNM != before.positionsNM)
        for constraint in model.constraints {
            let distance=(after.positionsNM[Int(constraint.a)]-after.positionsNM[Int(constraint.b)]).norm
            #expect(abs(distance-constraint.distanceNM)<2e-7)
        }
        let bad=try await VivoMDMetalRuntime.make(system:model,initialState:initial,configuration:configuration(iterations:1))
        let checkpoint=try await bad.checkpoint()
        await #expect(throws:Error.self) { try await bad.projectConstraints() }
        let retained=try await bad.checkpoint()
        #expect(retained==checkpoint)
        var old=after;old.numericalContract="numivivo.org/md-metal-numerics/v3"
        #expect(throws:Error.self) { try old.validate(particleCount:3) }
    }
}
