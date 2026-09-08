import Foundation
import Testing
@testable import NumiVivoKit

struct MDBenchmarkTests {
    private func fixture() throws -> (VivoMDBenchmarkReference,VivoMDHamiltonianEvaluation) {
        let id=try VivoCanonicalJSON.fingerprint(Data("benchmark-test".utf8))
        let geometry=try VivoMDCandidateGeometry(particlePositionsNM:[.zero,.init(1,0,0)],periodicCell:nil)
        let forces:[VivoVector3D]=[.init(2,0,0),.init(-2,0,0)]
        return (.init(identifier:"pair",geometry:geometry,energyKJPerMol:1,forcesKJPerMolNM:forces),
                .init(systemFingerprint:id,configurationFingerprint:id,evaluatedGeometry:geometry,
                    energyKJPerMol:1,physicalParticleForcesKJPerMolNM:forces,normalizedForceResidual:0))
    }
    @Test func exactReferencePassesAndIndependentErrorsFail() throws {
        let (reference,value)=try fixture()
        #expect(try VivoMDBenchmarkComparison.compare(value,reference:reference,limits:.init()).outcome == .passed)
        var changed=reference; changed.energyKJPerMol += 0.1
        #expect(try VivoMDBenchmarkComparison.compare(value,reference:changed,limits:.init()).outcome == .failed)
        changed=reference;changed.forcesKJPerMolNM[0].x += 1
        #expect(try VivoMDBenchmarkComparison.compare(value,reference:changed,limits:.init()).outcome == .failed)
        changed=reference;changed.geometry=try .init(particlePositionsNM:[.zero,.init(2,0,0)],periodicCell:nil)
        #expect(try VivoMDBenchmarkComparison.compare(value,reference:changed,limits:.init()).outcome == .failed)
    }
    @Test func NonfiniteAndMismatchedReferenceCannotPass() throws {
        let (reference,value)=try fixture()
        var changed=reference;changed.energyKJPerMol = .nan
        #expect(throws:Error.self) { try VivoMDBenchmarkComparison.compare(value,reference:changed,limits:.init()) }
        changed=reference;changed.forcesKJPerMolNM=[]
        #expect(throws:Error.self) { try VivoMDBenchmarkComparison.compare(value,reference:changed,limits:.init()) }
        var limits=VivoMDBenchmarkLimits();limits.forceNormalizedRMS = .infinity
        #expect(throws:Error.self) { try VivoMDBenchmarkComparison.compare(value,reference:reference,limits:limits) }
    }
    @Test func unsupportedCellPublishesNoPassingObservation() async throws {
        let (ref,_)=try fixture(),id=try VivoCanonicalJSON.fingerprint(Data("cell".utf8))
        let system=VivoClassicalSystem(identifier:"skew",structureFingerprint:id,particles:(0..<2).map {
            .init(index:UInt32($0),atomIndex:UInt32($0),typeIdentifier:"A",massDa:12,chargeE:0,sigmaNM:0,epsilonKJPerMol:0)
        })
        var ref=ref
        ref.geometry=try .init(particlePositionsNM:[.zero,.init(1,0,0)],
            periodicCell:.init(a:.init(4,0,0),b:.init(1,4,0),c:.init(0,0,4)))
        let request=VivoMDBenchmarkRequest(identifier:"skew",system:system,configuration:.init(),references:[ref],referenceProvenance:["fixture":"synthetic"])
        let report=try await VivoMDBenchmark.run(request)
        #expect(report.outcome == .unsupported)
        #expect(report.deviceName == nil && report.comparisons.isEmpty && report.dynamics == nil)
    }
}
