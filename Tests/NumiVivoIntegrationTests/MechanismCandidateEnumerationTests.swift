import Foundation
import Testing
@testable import NumiVivoKit

@Suite(.serialized) struct MechanismCandidateEnumerationTests {
    private func structure() -> VivoMolecularStructure {
        let carbon = VivoElement(atomicNumber: 6, symbol: "C")
        let atoms = (0..<4).map { VivoMolecularAtom(index: UInt32($0), name: "C\($0)", element: carbon) }
        return .init(identifier: "mechanism-fixture", atoms: atoms,
            bonds: [.init(atomA:0,atomB:1,order:.single),.init(atomA:2,atomB:3,order:.single)],
            conformers:[.init(positionsNM:[.zero,.init(0.15,0,0),.init(0.30,0,0),.init(0.45,0,0)])])
    }

    @Test func enumeratesMappedFormBreakHypothesisAndCoordinateSeed() throws {
        let request = VivoMechanismCandidateRequest(identifier:"mapped",
            structure:structure(), editableBonds:[
                .init(atomA:1,atomB:2,allowedTargets:[.single]),
                .init(atomA:2,atomB:3,allowedTargets:[.absent])
            ], maximumEditsPerCandidate:2)
        let result = try VivoMechanismCandidateEnumeration.enumerate(request)
        #expect(result.candidates.count == 3)
        let coupled = try #require(result.candidates.first { $0.changes.count == 2 })
        #expect(coupled.changes.contains { $0.atomA == 1 && $0.atomB == 2 && $0.from == .absent && $0.to == .single })
        #expect(coupled.changes.contains { $0.atomA == 2 && $0.atomB == 3 && $0.from == .single && $0.to == .absent })
        let seed = try #require(coupled.suggestedCoordinate)
        #expect(seed.kind == .distanceDifference)
        #expect(seed.atomIndices == [1,2,2,3])
        #expect(coupled.interpretation.contains("graph hypothesis only"))
    }

    @Test func explicitValenceLimitRejectsOverbondedCandidateButKeepsCoupledExchange() throws {
        let request = VivoMechanismCandidateRequest(identifier:"valence",
            structure:structure(), editableBonds:[
                .init(atomA:1,atomB:2,allowedTargets:[.single]),
                .init(atomA:2,atomB:3,allowedTargets:[.absent])
            ], valenceLimits:[.init(atomIndex:2,maximumBondOrderSum:1)], maximumEditsPerCandidate:2)
        let result = try VivoMechanismCandidateEnumeration.enumerate(request)
        #expect(result.rejectedByValence == 1)
        #expect(result.candidates.count == 2)
        #expect(result.candidates.contains { $0.changes.count == 2 })
    }

    @Test func unknownEditableBondOrderAndNoOpTargetsFailClosed() throws {
        var unknown = structure(); unknown.bonds[0].order = .unknown
        let unknownRequest = VivoMechanismCandidateRequest(identifier:"unknown",structure:unknown,
            editableBonds:[.init(atomA:0,atomB:1,allowedTargets:[.absent])])
        #expect(throws:(any Error).self) { _ = try VivoMechanismCandidateEnumeration.enumerate(unknownRequest) }
        let noOp = VivoMechanismCandidateRequest(identifier:"noop",structure:structure(),
            editableBonds:[.init(atomA:0,atomB:1,allowedTargets:[.single])])
        #expect(throws:(any Error).self) { _ = try VivoMechanismCandidateEnumeration.enumerate(noOp) }
    }
}
