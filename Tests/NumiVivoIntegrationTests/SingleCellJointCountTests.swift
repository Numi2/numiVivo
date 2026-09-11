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
    @Test func concaveCoordinateTangentsBoundCountLikelihood() throws {
        for phi in [0.0,0.2,100.0] { for y: UInt64 in [0,7] {
            let l=try VivoCountRateLikelihood(.init(libraryCounts: [100,2000],cellsPerLibrary: [2,3],geneCountsPerLibrary: [y,y*2]),cellDispersion: phi)
            let scale=l.concaveCoordinateScaleCPM
            for x in [0.1,1.0,5.0,10.0] {
                let r=scale*expm1(x),ell=try l.logRelativeLikelihood(rateCPM: r),slope=try l.coordinateSlope(rateCPM: r,scaleCPM: scale)
                for xx in [0.0,0.05,0.2,2.0,7.0,11.0] {
                    #expect(try l.logRelativeLikelihood(rateCPM: scale*expm1(xx))<=ell+slope*(xx-x)+1e-8)
                }
                let h=1e-5,d=(try l.logRelativeLikelihood(rateCPM: scale*expm1(x+h))-l.logRelativeLikelihood(rateCPM: scale*expm1(x-h)))/(2*h)
                #expect(abs(d-slope)<1e-6*max(1,abs(slope)))
            }
        } }
    }
    @Test func adaptiveIdenticalDonorsHaveAnExplicitContinuousBound() throws {
        let pairs=(0..<3).map { VivoJointCountPair(donorID: "d\($0)",control: stratum(0),treated: stratum(3)) }
        let m=try VivoAdaptiveJointCountResponse.fit(pairs: pairs,featureID: "g",controlConditionID: "C",treatedConditionID: "T",controlCellDispersion: 0.2,treatedCellDispersion: 0.2,trainingSource: source)
        #expect(m.status=="boundedContinuousLikelihood")
        #expect(m.certificate!.boxes.count==1)
        #expect(m.certificate!.maximumMeanDirectionalUpperBound<=1+1e-6)
        let p=try VivoAdaptiveJointCountResponse.predict(control: stratum(0),featureID: "g",queryDonorID: "new",controlConditionID: "C",querySource: source,model: m,plannedTreatedLibraryCounts: [1000])
        #expect(p.degenerateTreatedRateDistribution)
        #expect(abs(p.moments.treatedMeanCPM-3000)<1e-8)
    }
    @Test func adaptiveBudgetExhaustionDoesNotAdmitPrediction() throws {
        let pairs=zip([UInt64(0),1,2,4,9],[UInt64(6),3,1,2,0]).enumerated().map {
            VivoJointCountPair(donorID: "d\($0.offset)",control: stratum($0.element.0),treated: stratum($0.element.1))
        }
        let m=try VivoAdaptiveJointCountResponse.fit(pairs: pairs,featureID: "g",controlConditionID: "C",treatedConditionID: "T",controlCellDispersion: 0.2,treatedCellDispersion: 0.2,trainingSource: source,plan: .init(initialGridPointsPerAxis: 3,maximumSupportAdditions: 0,maximumOracleLeaves: 1))
        #expect(m.status != "boundedContinuousLikelihood")
        #expect(m.certificate != nil)
        #expect(throws: (any Error).self) { try VivoAdaptiveJointCountResponse.predict(control: stratum(1),featureID: "g",queryDonorID: "new",controlConditionID: "C",querySource: source,model: m,plannedTreatedLibraryCounts: [1000]) }
    }

    @Test func numericalAllowanceDoesNotTrapTheRealDNAJC8SupportSearch() throws {
        struct Group: Decodable { let donorID: String;let conditionID: String;let libraryCounts: [UInt64];let cellsPerLibrary: [Int] }
        struct Header: Decodable { let groups: [Group] }
        struct Counts: Decodable { let bins: [Int];let counts: [UInt64] }
        struct Gene: Decodable { let featureID: String;let controlCellDispersion: Double;let treatedCellDispersion: Double;let groups: [Counts] }
        struct Fixture: Decodable { let header: Header;let gene: Gene }
        let root=URL(fileURLWithPath: #filePath).deletingLastPathComponent().deletingLastPathComponent().deletingLastPathComponent()
        let file=root.appendingPathComponent("Tools/Omics/CountObservation/Joint/Adaptive/Full/DNAJC8-regression.json")
        let fixture=try JSONDecoder().decode(Fixture.self,from: Data(contentsOf: file))
        let groups=fixture.header.groups,gene=fixture.gene
        let strata=groups.indices.map { i -> VivoCountDepthStratum in
            var y=Array(repeating: UInt64(0),count: groups[i].libraryCounts.count)
            for j in gene.groups[i].bins.indices { y[gene.groups[i].bins[j]]=gene.groups[i].counts[j] }
            return .init(libraryCounts: groups[i].libraryCounts,cellsPerLibrary: groups[i].cellsPerLibrary,geneCountsPerLibrary: y)
        }
        let pairs=Set(groups.map(\.donorID)).sorted().map { donor in
            let c=groups.firstIndex { $0.donorID==donor && $0.conditionID=="control" }!
            let t=groups.firstIndex { $0.donorID==donor && $0.conditionID=="IFNB" }!
            return VivoJointCountPair(donorID: donor,control: strata[c],treated: strata[t])
        }
        let model=try VivoAdaptiveJointCountResponse.fit(pairs: pairs,featureID: gene.featureID,controlConditionID: "control",treatedConditionID: "IFNB",controlCellDispersion: gene.controlCellDispersion,treatedCellDispersion: gene.treatedCellDispersion,trainingSource: source)
        #expect(model.status=="boundedContinuousLikelihood")
        #expect(model.certificate!.maximumMeanDirectionalUpperBound<=1+model.plan.meanLogLikelihoodGapTolerance)
    }

}
