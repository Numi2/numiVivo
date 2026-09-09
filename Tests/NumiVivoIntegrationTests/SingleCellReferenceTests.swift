import Foundation
import Testing
@testable import NumiVivoKit

@Suite struct SingleCellReferenceTests {
    @Test func exactReferenceNeighborsPreserveDistanceAndIndexTies() throws {
        let result=try VivoSingleCellReference.nearest([0,0],reference: [[1,0],[-1,0],[0,2],[0,0],[0,-1]],neighbors: 3)
        #expect(result.map(\.index)==[3,0,1])
        #expect(result.map(\.squaredDistance)==[0,1,1])
        #expect(throws: (any Error).self) { try VivoSingleCellReference.nearest([Double.greatestFiniteMagnitude],reference: [[-Double.greatestFiniteMagnitude]],neighbors: 1) }
    }
    @Test func projectionCentersAreOptInAndReproduceSparseTrainingScores() throws {
        let dataset=try VivoSingleCellExamples.pairedCounts()
        let processed=try VivoSingleCellProcessing.run(dataset,policy: .init(),normalizationTarget: 10_000)
        var options=VivoSingleCellReductionOptions();options.components=2;options.highlyVariableFeatures=4;options.maximumBasis=32
        let original=try VivoSingleCellReduction.run(processed,options: options)
        #expect(original.projectionCenters==nil)
        options.retainProjectionCenters=true
        let frozen=try VivoSingleCellReduction.run(processed,options: options)
        #expect(frozen.scores==original.scores)
        let centers=try #require(frozen.projectionCenters)
        #expect(centers==frozen.selectedFeatureIndices.map { processed.features[$0].meanLogNormalized! })
        let selected=Dictionary(uniqueKeysWithValues: frozen.selectedFeatureIndices.enumerated().map { ($0.element,$0.offset) })
        let normal=processed.normalized
        for row in frozen.cells.indices {
            var score=[Double](repeating: 0,count: 2)
            for j in centers.indices { for c in 0..<2 { score[c]-=centers[j]*frozen.loadings[j][c] } }
            for k in normal.rowOffsets[row]..<normal.rowOffsets[row+1] {
                if let j=selected[normal.featureIndices[k]] { for c in 0..<2 { score[c]+=normal.values[k]*frozen.loadings[j][c] } }
            }
            for c in 0..<2 { #expect(abs(score[c]-frozen.scores[row][c])<1e-10) }
        }
    }
}
