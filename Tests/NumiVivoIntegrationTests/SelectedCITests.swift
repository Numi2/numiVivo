import Foundation
import Testing
@testable import NumiVivoKit

@Suite(.serialized) struct SelectedCITests {
    private func h2Hamiltonian() throws -> VivoEmbeddedHamiltonian {
        let source=VivoBarrierBenchmarks.hydrogenExchange631G()
        let basis=VivoGaussianBasis(identifier:"H2-6-31G-selected-ci",
            shells:source.basis.shells.filter{$0.nucleusIndex<2},source:source.basis.source)
        let system=VivoElectronicSystem(nuclei:[.init(atomicNumber:1,positionBohr:.init(0,0,-0.7)),
            .init(atomicNumber:1,positionBohr:.init(0,0,0.7))],alphaElectrons:1,betaElectrons:1)
        let ao=try VivoGaussianIntegralEngine.compute(system:system,basis:basis)
        let hf=try VivoHartreeFock.solve(system:system,integrals:ao,configuration:.init(reference:.restricted,energyToleranceHartree:1e-12,densityTolerance:1e-10))
        return try VivoEmbeddedHamiltonian.fromAO(ao,coefficients:hf.alphaCoefficients,
            alphaElectrons:1,betaElectrons:1,orbitalIdentifiers:(0..<ao.count).map{"h2-mo-\($0)"},
            energyReference:"H2/6-31G RHF orbital frame for selected-CI regression")
    }
    private func algebraic(_ diagonal: [Double], _ couplings: [(Int,Int,Double)] = [],
                           alpha: Int = 1, beta: Int = 0) -> VivoEmbeddedHamiltonian {
        let n=diagonal.count;var one=VivoQMMatrix(n,n)
        for i in 0..<n {one[i,i]=diagonal[i]}
        for (i,j,v) in couplings {one[i,j]=v;one[j,i]=v}
        return .init(orbitalIdentifiers:(0..<n).map{"algebraic-\($0)"},alphaElectrons:alpha,betaElectrons:beta,
            oneElectron:one,twoElectron:[Double](repeating:0,count:n*n*n*n),constantEnergyHartree:0,
            energyReference:"explicit algebraic regression; not molecular accuracy evidence")
    }

    @Test func adaptiveSelectionApproachesFullCIWithoutClaimingPT2ErrorBound() throws {
        let h=try h2Hamiltonian(),fci=try VivoDirectCI.solve(h,configuration:.init(residualTolerance:1e-12)).roots[0]
        let selected=try VivoSelectedCI.solve(h,configuration:.init(maximumIterations:16,maximumDeterminants:16,
            selectionBatchSize:4,minimumSelectionContributionHartree:0,pt2ToleranceHartree:1e-10,
            eigenResidualTolerance:1e-11,minimumDenominatorHartree:1e-6,maximumDavidsonSubspace:16))
        #expect(selected.converged)
        #expect(selected.fullSectorDimension==16)
        #expect(selected.selectedDeterminantCount<=16)
        #expect(selected.iterations.count>=1)
        #expect(selected.iterations.map(\.determinantCount)==selected.iterations.map(\.determinantCount).sorted())
        #expect(abs(selected.variationalEnergyHartree-fci.energyHartree)<1e-8)
        let correction=try #require(selected.pt2CorrectionHartree)
        #expect(abs(correction)<=1e-10)
        #expect(selected.eigenResidual<=1e-11)
        #expect(selected.method.contains("PT2 is a diagnostic remainder"))
        #expect(selected.diagnostics.fullResidualNorm<=selected.configuration.effectiveFullResidualTolerance)
        try selected.state.validate()
        try VivoSelectedCI.validate(selected,hamiltonian:h)
    }

    @Test func determinantCapProducesExplicitUnconvergedRemainder() throws {
        let h=try h2Hamiltonian(),fci=try VivoDirectCI.solve(h).roots[0]
        let limited=try VivoSelectedCI.solve(h,configuration:.init(maximumIterations:4,maximumDeterminants:1,
            selectionBatchSize:1,minimumSelectionContributionHartree:0,pt2ToleranceHartree:1e-12,
            eigenResidualTolerance:1e-10,minimumDenominatorHartree:1e-6,maximumDavidsonSubspace:8))
        #expect(!limited.converged)
        #expect(limited.termination == .determinantLimit)
        #expect(limited.selectedDeterminantCount==1)
        #expect(limited.variationalEnergyHartree>=fci.energyHartree-1e-10)
        let correction=try #require(limited.pt2CorrectionHartree)
        #expect(abs(correction)>1e-12)
        #expect(limited.iterations.last?.externalCandidateCount ?? 0 > 0)
        try VivoSelectedCI.validate(limited,hamiltonian:h)
    }

    @Test func oppositeDenominatorsCannotCancelIntoFalseConvergence() throws {
        let h=algebraic([0,-1,1],[(0,1,0.1),(0,2,0.1)])
        let result=try VivoSelectedCI.solve(h,configuration:.init(maximumDeterminants:1,selectionBatchSize:1))
        #expect(!result.converged)
        #expect(result.diagnostics.intruderCount==1)
        #expect(result.pt2CorrectionHartree==nil)
        #expect(result.pt2CorrectedEnergyHartree==nil)
        #expect(abs(result.diagnostics.externalResidualNorm-sqrt(0.02))<1e-12)
        try VivoSelectedCI.validate(result,hamiltonian:h)
    }

    @Test func nearZeroDenominatorIsNotSilentlyRegularized() throws {
        let h=algebraic([0,0],[(0,1,0.01)])
        let result=try VivoSelectedCI.solve(h,configuration:.init(maximumDeterminants:1,selectionBatchSize:1))
        #expect(result.diagnostics.intruderCount==1)
        #expect(result.pt2CorrectionHartree==nil)
        #expect(!result.converged)
        let expanded=try VivoSelectedCI.solve(h,configuration:.init(maximumDeterminants:2,selectionBatchSize:1))
        #expect(expanded.converged)
        #expect(abs(expanded.variationalEnergyHartree+0.01)<1e-12)
    }

    @Test func completeResidualDoesNotDropTinyCoefficients() throws {
        let h=algebraic([0,1,2],[(1,2,1)])
        let state=VivoCIState(orbitalCount:3,alphaElectrons:1,betaElectrons:0,
            determinants:[1,4],coefficients:[1,1e-16])
        let checked=try VivoSelectedCI.diagnose(h,state:state,energyHartree:0)
        #expect(checked.externalResidualNorm>0)
        #expect(abs(checked.externalResidualNorm-1e-16)<1e-30)
    }

    @Test func connectedResidualMatchesEnumeratedFullSector() throws {
        let h=algebraic([0,1,2],[(0,1,0.1),(0,2,0.2),(1,2,0.3)])
        let state=VivoCIState(orbitalCount:3,alphaElectrons:1,betaElectrons:0,
            determinants:[4,1],coefficients:[0.6,0.8])
        let sparse=try VivoSelectedCI.diagnose(h,state:state,energyHartree:0.7)
        let dense=try VivoDirectCI.residualNorm(hamiltonian:h,state:state,energyHartree:0.7)
        #expect(abs(sparse.fullResidualNorm-dense)<1e-12)
    }

    @Test func externalCapacityThrowsInsteadOfPublishingTruncatedResidual() throws {
        let h=algebraic([0,1,2],[(0,1,0.1),(0,2,0.2)])
        #expect(throws:VivoChemistryError.self) {
            try VivoSelectedCI.solve(h,configuration:.init(maximumDeterminants:1,selectionBatchSize:1,
                maximumExternalDeterminants:1))
        }
    }

    @Test func selectedSpaceDoesNotEnumerateLargeFullSector() throws {
        let h=algebraic((1...10).map(Double.init),alpha:5,beta:5)
        #expect(throws:VivoChemistryError.self) {try VivoDirectCI.solve(h)}
        let result=try VivoSelectedCI.solve(h,configuration:.init(maximumDeterminants:1,selectionBatchSize:1))
        #expect(result.converged)
        #expect(result.fullSectorDimension==63_504)
        #expect(result.selectedDeterminantCount==1)
        #expect(abs(result.variationalEnergyHartree-30)<1e-12)
    }

    @Test func duplicateSeedsAndForgedConvergenceReject() throws {
        let h=algebraic([0,1],[(0,1,0.1)])
        #expect(throws:VivoChemistryError.self) {try VivoSelectedCI.solve(h,initialDeterminants:[1,1])}
        let result=try VivoSelectedCI.solve(h,configuration:.init(maximumDeterminants:1,selectionBatchSize:1))
        var object=try #require(JSONSerialization.jsonObject(with:JSONEncoder().encode(result)) as? [String:Any])
        object["converged"]=true;object["termination"]="residualConverged"
        let forged=try JSONDecoder().decode(VivoSelectedCIResult.self,from:JSONSerialization.data(withJSONObject:object))
        #expect(throws:VivoChemistryError.self) {try VivoSelectedCI.validate(forged,hamiltonian:h)}
    }
}
