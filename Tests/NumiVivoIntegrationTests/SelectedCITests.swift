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
        #expect(abs(selected.pt2CorrectionHartree)<=1e-10)
        #expect(selected.eigenResidual<=1e-11)
        #expect(selected.method.contains("PT2 is a diagnostic remainder"))
        try selected.state.validate()
    }

    @Test func determinantCapProducesExplicitUnconvergedRemainder() throws {
        let h=try h2Hamiltonian(),fci=try VivoDirectCI.solve(h).roots[0]
        let limited=try VivoSelectedCI.solve(h,configuration:.init(maximumIterations:4,maximumDeterminants:1,
            selectionBatchSize:1,minimumSelectionContributionHartree:0,pt2ToleranceHartree:1e-12,
            eigenResidualTolerance:1e-10,minimumDenominatorHartree:1e-6,maximumDavidsonSubspace:8))
        #expect(!limited.converged)
        #expect(limited.selectedDeterminantCount==1)
        #expect(limited.variationalEnergyHartree>=fci.energyHartree-1e-10)
        #expect(abs(limited.pt2CorrectionHartree)>1e-12)
        #expect(limited.iterations.last?.externalCandidateCount ?? 0 > 0)
    }
}
