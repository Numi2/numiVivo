import Foundation
import Testing
@testable import NumiVivoKit

@Suite struct SingleCellProgramTests {
    func definition(_ members: [VivoSingleCellProgramMember],coverage: Double = 1) -> VivoSingleCellProgramDefinition {
        .init(id: "signed",organism: "human",featureNamespace: "fixture-id",sourceURI: "urn:numivivo:fixture",sourceVersion: "1",
              sourceDescription: "Numerical fixture only",members: members,minimumWeightCoverage: coverage)
    }
    func accumulator(_ options: VivoSingleCellProgramOptions,organism: String = "human") throws -> VivoSingleCellProgramAccumulator {
        let features=["A","B","C"].map { VivoOmicsFeature(id: $0,name: $0) }
        let cells=(0..<3).map { VivoOmicsCell(barcode: String($0),sampleID: "s") }
        var quality: [VivoCellQuality] = []
        for row in 0..<3 {
            let total: UInt64 = row == 2 ? 0 : 10
            let count: Int = row == 2 ? 0 : 1
            quality.append(VivoCellQuality(sampleID: "s",barcode: String(row),totalCounts: total,
                detectedFeatures: count,mitochondrialCounts: 0,mitochondrialFeatureCount: 0,mitochondrialFraction: nil))
        }
        return try .init(features: features,cells: cells,samples: [.init(id: "s",biologicalReplicateID: "d",condition: "c",batchID: "b",organism: organism)],
                         quality: quality,normalizationTarget: 10000,options: options)
    }
    @Test func signedScoresSeparateAbsentExpressionFromEmptyLibraries() throws {
        let d=definition([.init(featureID: "A",weight: 2),.init(featureID: "B",weight: -1)])
        let a=try accumulator(.init(definitions: [d]))
        try a.add(row: 0,feature: 0,logValue: 3);try a.add(row: 0,feature: 1,logValue: 1.5)
        try a.add(row: 1,feature: 2,logValue: 8)
        let r=try a.finish()
        #expect(r.scores==[[1.5],[0],[nil]])
        #expect(r.detectedMembers==[[2],[0],[0]] && r.updates==2)
        #expect(r.programs[0].featureIndices==[0,1] && r.programs[0].weightCoverage==1)
        #expect(throws: (any Error).self) { try a.add(row: 2,feature: 0,logValue: 1) }
    }
    @Test func missingFeaturesRequireExplicitCoverageAndChangeTheFingerprint() throws {
        #expect(throws: (any Error).self) { try accumulator(.init(definitions: [definition([
            .init(featureID: "A"),.init(featureID: "missing",weight: 1e-100)])])) }
        let members: [VivoSingleCellProgramMember]=[.init(featureID: "A"),.init(featureID: "missing")]
        #expect(throws: (any Error).self) { try accumulator(.init(definitions: [definition(members)])) }
        let a=try accumulator(.init(definitions: [definition(members,coverage: 0.5)]))
        try a.add(row: 0,feature: 0,logValue: 2)
        let r=try a.finish();#expect(r.scores[0][0]==2)
        #expect(r.programs[0].missingFeatureIDs==["missing"] && r.programs[0].weightCoverage==0.5)
        let b=try accumulator(.init(definitions: [definition([.init(featureID: "A")])]))
        #expect(r.programs[0].definitionFingerprint != b.programs[0].definitionFingerprint)
        #expect(throws: (any Error).self) { try accumulator(.init(definitions: [definition([.init(featureID: "missing")],coverage: 0.01)])) }
    }
    @Test func sourceOrganismBudgetsAndInvalidDefinitionsFailClosed() throws {
        let d=definition([.init(featureID: "A")]);var options=VivoSingleCellProgramOptions(definitions: [d])
        #expect(throws: (any Error).self) { try accumulator(options,organism: "mouse") }
        options.maximumScoreValues=2
        #expect(throws: (any Error).self) { try accumulator(options) }
        options.maximumScoreValues=3;options.maximumUpdates=1
        let a=try accumulator(options);try a.add(row: 0,feature: 0,logValue: 1)
        #expect(throws: (any Error).self) { try a.add(row: 1,feature: 0,logValue: 1) }
        let invalidMembers: [[VivoSingleCellProgramMember]] = [[.init(featureID: "A",weight: 0)],[.init(featureID: "A",weight: .nan)],
            [.init(featureID: "A"),.init(featureID: "A")]]
        for members in invalidMembers {
            #expect(throws: (any Error).self) { try definition(members).validate() }
        }
        #expect(throws: (any Error).self) { try JSONDecoder().decode(VivoSingleCellProgramMember.self,from: Data("{\"featureID\":\"A\",\"weight\":1,\"alias\":\"B\"}".utf8)) }
    }
    @Test func finiteWeightScalingAndMemberOrderDoNotChangeScores() throws {
        let a=try accumulator(.init(definitions: [definition([.init(featureID: "A",weight: Double.greatestFiniteMagnitude),.init(featureID: "B",weight: -Double.greatestFiniteMagnitude)])]))
        let b=try accumulator(.init(definitions: [definition([.init(featureID: "B",weight: -1),.init(featureID: "A",weight: 1)])]))
        for x in [a,b] { try x.add(row: 0,feature: 0,logValue: 4);try x.add(row: 0,feature: 1,logValue: 2) }
        #expect(try a.finish().scores==b.finish().scores)
    }
}
