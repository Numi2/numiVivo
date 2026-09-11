import Foundation
import Testing
@testable import NumiVivoKit

@Suite struct SingleCellJointCountTests {
    let source=try! VivoFingerprint(bytes: Array(repeating: 11,count: 32))
    func stratum(_ y: UInt64,_ n: Int=1) -> VivoCountDepthStratum {
        .init(libraryCounts: [1000],cellsPerLibrary: [n],geneCountsPerLibrary: [y*UInt64(n)])
    }
    func fit(_ pairs: [VivoJointCountPair],phi: Double=0,iterations: Int=20_000) throws -> VivoJointCountModel {
        try VivoJointCountResponse.fit(pairs: pairs,featureID: "g",controlConditionID: "C",treatedConditionID: "T",
            controlCellDispersion: phi,treatedCellDispersion: phi,trainingSource: source,
            plan: .init(gridPointsPerAxis: 17,maximumIterations: iterations))
    }
    @Test func binnedCountLikelihoodMatchesIndividualCellRateTerms() throws {
        let s=VivoCountDepthStratum(libraryCounts: [100,200],cellsPerLibrary: [1,2],geneCountsPerLibrary: [0,5])
        for phi in [0.0,0.7,100.0] {
            let l=try VivoCountRateLikelihood(s,cellDispersion: phi)
            let r=1234.0,m=l.maximumLikelihoodRateCPM
            var expected=0.0
            for (y,e) in [(0.0,0.0001),(2.0,0.0002),(3.0,0.0002)] {
                if phi==0 { expected+=y*log(r/m)-e*(r-m) }
                else { expected+=y*log(r/m)-(y+1/phi)*(log1p(phi*e*r)-log1p(phi*e*m)) }
            }
            #expect(abs(try l.logRelativeLikelihood(rateCPM: r)-expected)<1e-10)
            #expect(abs(try l.logRelativeLikelihood(rateCPM: m))<1e-12)
            #expect(try l.logRelativeLikelihood(rateCPM: 0)==(-Double.infinity))
        }
    }
    @Test func exactZeroRateAndInvalidDepthBinsAreExplicit() throws {
        let l=try VivoCountRateLikelihood(stratum(0,4),cellDispersion: 1)
        #expect(l.maximumLikelihoodRateCPM==0)
        #expect(try l.logRelativeLikelihood(rateCPM: 0)==0)
        #expect(abs(try l.logRelativeLikelihood(rateCPM: 1000)+4*log(2))<1e-12)
        #expect(throws: (any Error).self) { try VivoCountRateLikelihood(.init(libraryCounts: [100,100],cellsPerLibrary: [1,1],geneCountsPerLibrary: [0,0]),cellDispersion: 1) }
        #expect(throws: (any Error).self) { try VivoCountRateLikelihood(stratum(1001),cellDispersion: 1) }
        #expect(throws: (any Error).self) { try VivoCountObservation.predictiveMoments(meanCPM: 0,varianceCPM: 1,cellDispersion: 0,plannedLibraryCounts: [100]) }
    }
    @Test func negativeJointAssociationPredictsFromControlOnly() throws {
        let pairs=zip([UInt64(100),300,900],[UInt64(900),300,100]).enumerated().map {
            VivoJointCountPair(donorID: "d\($0.offset)",control: stratum($0.element.0,100),treated: stratum($0.element.1,100))
        }
        let m=try fit(pairs)
        #expect(m.status=="convergedFiniteGridLikelihood")
        #expect(m.meanLogLikelihoodGap<=1e-8)
        #expect(m.moments.covarianceCPM2<0)
        #expect(abs(m.relativeLogLikelihood+3*log(3))<1e-10)
        #expect(try fit(pairs.reversed())==m)
        let p=try VivoJointCountResponse.predict(control: stratum(100,100),featureID: "g",queryDonorID: "new",controlConditionID: "C",querySource: source,model: m,plannedTreatedLibraryCounts: Array(repeating: 1000,count: 20))
        #expect(p.moments.treatedMeanCPM>800000)
        #expect(p.moments.controlMeanCPM<200000)
        #expect(p.plannedTreatedCountMoments.meanGeneCounts>16000)
        #expect(!p.degenerateTreatedRateDistribution)
        #expect(p.underflowedPosteriorComponents>0)
        #expect(throws: (any Error).self) { try VivoJointCountResponse.predict(control: stratum(1),featureID: "g",queryDonorID: "d0",controlConditionID: "C",querySource: source,model: m,plannedTreatedLibraryCounts: [1000]) }
    }
    @Test func iterationExhaustionIsNotAnAdmittedPredictionModel() throws {
        let pairs=zip([UInt64(0),1,2,4,9],[UInt64(6),3,1,2,0]).enumerated().map {
            VivoJointCountPair(donorID: "d\($0.offset)",control: stratum($0.element.0),treated: stratum($0.element.1))
        }
        let m=try fit(pairs,phi: 0.2,iterations: 1)
        #expect(m.status != "convergedFiniteGridLikelihood")
        #expect(throws: (any Error).self) { try VivoJointCountResponse.predict(control: stratum(1),featureID: "g",queryDonorID: "new",controlConditionID: "C",querySource: source,model: m,plannedTreatedLibraryCounts: [1000]) }
    }
    @Test func fittedPointMassIsNotHiddenAsBiologicalCertainty() throws {
        let pairs=(0..<3).map { VivoJointCountPair(donorID: "d\($0)",control: stratum(0),treated: stratum(0)) }
        let m=try fit(pairs)
        let p=try VivoJointCountResponse.predict(control: stratum(0),featureID: "g",queryDonorID: "new",controlConditionID: "C",querySource: source,model: m,plannedTreatedLibraryCounts: [1000])
        #expect(p.degenerateTreatedRateDistribution)
        #expect(p.plannedTreatedCountMoments.totalGeneCountVariance==0)
        #expect(throws: (any Error).self) { try VivoJointCountResponse.predict(control: stratum(1),featureID: "g",queryDonorID: "new",controlConditionID: "C",querySource: source,model: m,plannedTreatedLibraryCounts: [1000]) }
    }
    @Test func constantTreatedCoordinateRemainsExactlyConstantAfterConditioning() throws {
        let pairs=[UInt64(100),300,900].enumerated().map { VivoJointCountPair(donorID: "d\($0.offset)",control: stratum($0.element,100),treated: stratum(3,100)) }
        let m=try fit(pairs)
        for y in stride(from: 50,through: 350,by: 17) {
            let p=try VivoJointCountResponse.predict(control: stratum(UInt64(y)),featureID: "g",queryDonorID: "new",controlConditionID: "C",querySource: source,model: m,plannedTreatedLibraryCounts: [1000])
            #expect(p.degenerateTreatedRateDistribution)
            #expect(!p.numericallyDegenerateTreatedRateMoments)
            #expect(p.moments.treatedVarianceCPM2==0)
            #expect(p.moments.treatedMeanCPM==3000)
        }
    }
}
