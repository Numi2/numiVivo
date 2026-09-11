import Foundation
import Testing
@testable import NumiVivoKit

@Suite struct SingleCellCountObservationTests {
    @Test func zeroCountsRetainProperPriorUncertainty() throws {
        let p=try VivoCountObservation.posterior(counts: [0,0,0],libraryCounts: [1000,2000,4000],cellDispersion: 0,gammaPriorShape: 0.5,gammaPriorRatePerCPM: 0.001)
        // Exact Gamma-Poisson conjugacy, independently known analytic moments.
        #expect(abs(p.meanCPM-62.5)<1e-7)
        #expect(abs(p.varianceCPM-7812.5)<1e-5)
        #expect(p.lowerLog1pCPM>0 && p.upperLog1pCPM>p.lowerLog1pCPM)
        #expect(p.maximumScaledMomentTailBound<1e-12)
    }
    @Test func poissonCountsAggregateButNBCellsCannotBeCollapsed() throws {
        func run(_ counts: [UInt64],_ depths: [UInt64],_ phi: Double) throws -> VivoCountObservationPosterior {
            try VivoCountObservation.posterior(counts: counts,libraryCounts: depths,cellDispersion: phi,gammaPriorShape: 1,gammaPriorRatePerCPM: 0.001)
        }
        let p=try run([5,10],[1000,2000],0),aggregate=try run([15],[3000],0)
        #expect(abs(p.meanCPM-4000)<1e-6 && abs(p.varianceCPM-1_000_000)<1e-3)
        #expect(abs(p.meanLog1pCPM-aggregate.meanLog1pCPM)<1e-10)
        let nb=try run([5,10],[1000,2000],1),wrong=try run([15],[3000],1)
        #expect(abs(nb.varianceLog1pCPM-wrong.varianceLog1pCPM)>0.01)
    }
    @Test func plannedCellsAffectOverdispersionAtFixedTotalDepth() throws {
        let p=try VivoCountObservation.posterior(counts: [0,3,1],libraryCounts: [1000,2000,3000],cellDispersion: 2,gammaPriorShape: 0.5,gammaPriorRatePerCPM: 0.001)
        let a=try VivoCountObservation.predictiveMoments(p,plannedLibraryCounts: [2000]),b=try VivoCountObservation.predictiveMoments(p,plannedLibraryCounts: [1000,1000])
        #expect(a.meanGeneCounts==b.meanGeneCounts && a.latentRateVariance==b.latentRateVariance)
        #expect(abs(a.conditionalCellOverdispersionVariance-2*b.conditionalCellOverdispersionVariance)<1e-10)
        #expect(a.totalGeneCountVariance>b.totalGeneCountVariance)
        #expect(try VivoCountObservation.posterior(counts: [0,3,1],libraryCounts: [1000,2000,3000],cellDispersion: 2,gammaPriorShape: 0.5,gammaPriorRatePerCPM: 0.001)==p)
    }
    @Test func concentratedPosteriorDoesNotLoseItsVariance() throws {
        let p=try VivoCountObservation.posterior(counts: [100_000_000],libraryCounts: [1_000_000_000],cellDispersion: 0,gammaPriorShape: 1,gammaPriorRatePerCPM: 0.001)
        let shape=100_000_001.0,rate=1000.001
        #expect(abs(p.meanCPM-shape/rate)/(shape/rate)<1e-9)
        #expect(abs(p.varianceCPM-shape/(rate*rate))/(shape/(rate*rate))<1e-8)
    }
    @Test func invalidCountsPriorsAndPlansFail() throws {
        let invalid: [([UInt64],[UInt64],Double,Double,Double)] = [([],[],0.0,1.0,1.0),([2],[1],0,1,1),([0],[0],0,1,1),([0],[1],-1,1,1),([0],[1],0,0,1),([0],[1],0,1,0)]
        for (counts,depths,phi,a,b) in invalid {
            #expect(throws: (any Error).self) { try VivoCountObservation.posterior(counts: counts,libraryCounts: depths,cellDispersion: phi,gammaPriorShape: a,gammaPriorRatePerCPM: b) }
        }
    }
}
