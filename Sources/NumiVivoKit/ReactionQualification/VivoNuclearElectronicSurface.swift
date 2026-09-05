import Foundation

public enum VivoNuclearSolver: String, Codable, Sendable { case hartreeFock, fullCI, anchoredECC }
public struct VivoAnchoredECCFrame: Codable, Sendable, Equatable {
    public var anchorSystem: VivoElectronicSystem
    /// For a point from a shared path: transportRotations[i] * sharedRotation.
    public var anchorRotation: VivoQMMatrix
    public var transportGroups: [[Int]]
    public var expectedBathOrbitals: [Int]
    public var minimumTransportSingularValue: Double
    public var embedding: VivoECCDMETConfiguration
    public init(anchorSystem: VivoElectronicSystem, anchorRotation: VivoQMMatrix, transportGroups: [[Int]],
                expectedBathOrbitals: [Int], embedding: VivoECCDMETConfiguration, minimumTransportSingularValue: Double = 0.8) {
        self.anchorSystem=anchorSystem;self.anchorRotation=anchorRotation;self.transportGroups=transportGroups
        self.expectedBathOrbitals=expectedBathOrbitals;self.embedding=embedding;self.minimumTransportSingularValue=minimumTransportSingularValue
    }
}
public struct VivoNuclearElectronicModel: Codable, Sendable, Equatable {
    public var system: VivoElectronicSystem
    public var basis: VivoGaussianBasis
    public var solver: VivoNuclearSolver
    public var scf: VivoSCFConfiguration
    /// HF-equilibrated smooth polarization. For correlated solvers the field is
    /// frozen at the reference density at each geometry, not at correlated D.
    public var solvent: VivoSmoothCPCMConfiguration?
    public var eccFrame: VivoAnchoredECCFrame?
    public var budget: VivoChemistryBudget
    public init(system: VivoElectronicSystem, basis: VivoGaussianBasis, solver: VivoNuclearSolver,
                scf: VivoSCFConfiguration = .init(), solvent: VivoSmoothCPCMConfiguration? = nil,
                eccFrame: VivoAnchoredECCFrame? = nil, budget: VivoChemistryBudget = .init()) {
        self.system=system;self.basis=basis;self.solver=solver;self.scf=scf;self.solvent=solvent;self.eccFrame=eccFrame;self.budget=budget
    }
    public func validate() throws {
        try system.validate();try basis.validate(nucleusCount:system.nuclei.count);try budget.validate();try scf.validate()
        guard (solver == .anchoredECC)==(eccFrame != nil) else {throw VivoChemistryError.invalid("nuclear ECC frame/method mismatch")}
        if let frame=eccFrame {
            try frame.anchorSystem.validate()
            guard frame.embedding.localityGroups.isEmpty,scf.reference == .restricted,
                  frame.anchorSystem.nuclei.map(\.atomicNumber)==system.nuclei.map(\.atomicNumber),
                  frame.anchorSystem.nuclei.map(\.structureAtomIndex)==system.nuclei.map(\.structureAtomIndex),
                  frame.anchorSystem.alphaElectrons==system.alphaElectrons,frame.anchorSystem.betaElectrons==system.betaElectrons,
                  frame.anchorSystem.pointCharges==system.pointCharges,
                  frame.minimumTransportSingularValue.isFinite,frame.minimumTransportSingularValue>0,frame.minimumTransportSingularValue<=1 else {
                throw VivoChemistryError.invalid("anchored nuclear ECC requires fixed orbital frame and unchanged system identity")
            }
        }
    }
}
public struct VivoNuclearDifferenceConfiguration: Codable, Sendable, Equatable {
    public var gradientStepBohr: Double
    public var hessianStepBohr: Double
    public var maximumGradientDifference: Double
    public var maximumHessianDifference: Double
    public var maximumEnergyEvaluations: Int
    public init(gradientStepBohr: Double = 2e-4, hessianStepBohr: Double = 1e-3,
                maximumGradientDifference: Double = 2e-6, maximumHessianDifference: Double = 1e-4,
                maximumEnergyEvaluations: Int = 100000) {
        self.gradientStepBohr=gradientStepBohr;self.hessianStepBohr=hessianStepBohr
        self.maximumGradientDifference=maximumGradientDifference;self.maximumHessianDifference=maximumHessianDifference
        self.maximumEnergyEvaluations=maximumEnergyEvaluations
    }
    public func validate() throws {
        guard gradientStepBohr.isFinite,(1e-6...1e-2).contains(gradientStepBohr),hessianStepBohr.isFinite,
              (1e-5...1e-1).contains(hessianStepBohr),maximumGradientDifference.isFinite,maximumGradientDifference>0,
              maximumHessianDifference.isFinite,maximumHessianDifference>0,(1...10000000).contains(maximumEnergyEvaluations) else {
            throw VivoChemistryError.invalid("nuclear finite-difference steps, agreement bounds or aggregate work")
        }
    }
}
public struct VivoReferencePolarizedHamiltonian: Sendable {
    public let hamiltonian: VivoEmbeddedHamiltonian
    public let coefficients: VivoQMMatrix
    public let reference: VivoHartreeFockResult
    public let solvent: VivoSmoothCPCMResult?
    public let frozenFieldConstantHartree: Double
}
public enum VivoReferencePolarization {
    /// At fixed reference surface charge q, physical energy is
    /// <H_gas + V(q)> + [Gpol(Dref) - Tr(Dref V(q))]. The second term retains
    /// nuclear-surface + surface-self energy; simply adding Gpol double counts V.
    public static func prepare(system: VivoElectronicSystem, ao: VivoAOIntegrals,
                               scf: VivoSCFConfiguration, solvent: VivoSmoothCPCMConfiguration?,
                               coefficients selected: VivoQMMatrix? = nil,
                               budget: VivoChemistryBudget = .init()) throws -> VivoReferencePolarizedHamiltonian {
        let reference:VivoHartreeFockResult,pcm:VivoSmoothCPCMResult?
        if let solvent {
            let result=try VivoSmoothSolventSCF.solve(system:system,integrals:ao,solvent:solvent,scf:scf,budget:budget)
            reference=result.scf;pcm=result.solvent
        } else {reference=try VivoHartreeFock.solve(system:system,integrals:ao,configuration:scf,budget:budget);pcm=nil}
        let c=selected ?? reference.alphaCoefficients
        let gas=try VivoEmbeddedHamiltonian.fromAO(ao,coefficients:c,alphaElectrons:system.alphaElectrons,
            betaElectrons:system.betaElectrons,orbitalIdentifiers:(0..<ao.count).map{"nuclear-orbital-\($0)"},
            energyReference:solvent == nil ? "gas electronic energy; AO scalar once" : "HF-reference-polarized smooth C-PCM; frozen-field nuclear and self term included once",budget:budget)
        var h=gas.oneElectron,constant=0.0
        if let pcm {
            h=try h.adding(pcm.reactionPotentialMatrix.congruence(c))
            let d=try reference.alphaDensity.adding(reference.betaDensity)
            constant=pcm.polarizationEnergyHartree-zip(d.values,pcm.reactionPotentialMatrix.values).reduce(0.0){$0+$1.0*$1.1}
        }
        let result=VivoEmbeddedHamiltonian(orbitalIdentifiers:gas.orbitalIdentifiers,alphaElectrons:gas.alphaElectrons,betaElectrons:gas.betaElectrons,
            oneElectron:h,twoElectron:gas.twoElectron,constantEnergyHartree:gas.constantEnergyHartree+constant,
            energyReference:gas.energyReference,provenance:["polarization":pcm == nil ? "none":"HF-reference frozen surface charge; not correlated self-consistent PCM"])
        try result.validate(budget:budget)
        return .init(hamiltonian:result,coefficients:c,reference:reference,solvent:pcm,frozenFieldConstantHartree:constant)
    }
}
/// One serial calculation-local evaluator; mutable work counters are not shared
/// across tasks. The energy cache stores only bounded scalar results.
public final class VivoNuclearElectronicSurface {
    public let model: VivoNuclearElectronicModel
    public let differences: VivoNuclearDifferenceConfiguration
    public private(set) var energyEvaluations=0
    private var cache:[[UInt64]:Double]=[:]
    private var anchorCoefficients:VivoQMMatrix?
    public init(model: VivoNuclearElectronicModel, differences: VivoNuclearDifferenceConfiguration = .init()) throws {
        try model.validate();try differences.validate();self.model=model;self.differences=differences
        if let frame=model.eccFrame {
            let ao=try VivoGaussianIntegralEngine.compute(system:frame.anchorSystem,basis:model.basis,budget:model.budget)
            let base=try VivoReferencePolarization.prepare(system:frame.anchorSystem,ao:ao,scf:model.scf,solvent:model.solvent,budget:model.budget)
            energyEvaluations=1
            anchorCoefficients=base.coefficients
            let rotatedAnchor=try base.coefficients.multiplied(by:frame.anchorRotation)
            let metric=try ao.overlap.congruence(rotatedAnchor)
            guard try metric.adding(.identity(ao.count),scale:-1).frobeniusNorm<1e-8 else {throw VivoChemistryError.invalid("anchored ECC rotation is not metric-orthonormal")}
        }
    }
    public func energy(_ positions: [SIMD3<Double>]) throws -> Double {
        guard positions.count==model.system.nuclei.count,positions.allSatisfy(vivoQMFinite) else {throw VivoChemistryError.invalid("nuclear energy coordinate shape")}
        let key=positions.flatMap{[$0.x.bitPattern,$0.y.bitPattern,$0.z.bitPattern]}
        if let value=cache[key] {return value}
        guard energyEvaluations<differences.maximumEnergyEvaluations else {throw VivoChemistryError.resourceLimit("aggregate nuclear electronic solve budget")};energyEvaluations+=1
        var system=model.system
        for i in positions.indices {system.nuclei[i].positionBohr=positions[i]}
        let ao=try VivoGaussianIntegralEngine.compute(system:system,basis:model.basis,budget:model.budget)
        let value:Double
        if model.solver == .hartreeFock {
            if let solvent=model.solvent {value=try VivoSmoothSolventSCF.solve(system:system,integrals:ao,solvent:solvent,scf:model.scf,budget:model.budget).scf.energyHartree}
            else {value=try VivoHartreeFock.solve(system:system,integrals:ao,configuration:model.scf,budget:model.budget).energyHartree}
        } else if model.solver == .fullCI && model.solvent == nil {
            let eig=try VivoQMDenseAlgebra.symmetricEigen(ao.overlap)
            guard eig.values[0]>=model.scf.minimumOverlapEigenvalue else {throw VivoChemistryError.invalid("nuclear AO rank loss")}
            var c=eig.vectors
            for i in 0..<ao.count {for j in 0..<ao.count {c[i,j]/=sqrt(eig.values[j])}}
            let h=try VivoEmbeddedHamiltonian.fromAO(ao,coefficients:c,alphaElectrons:system.alphaElectrons,betaElectrons:system.betaElectrons,
                orbitalIdentifiers:(0..<ao.count).map{"ao-orthogonal-\($0)"},energyReference:"gas FCI; AO scalar once",budget:model.budget)
            value=try VivoDirectCI.solve(h,configuration:.init(residualTolerance:1e-12),budget:model.budget).roots[0].energyHartree
        } else {
            let prepared=try VivoReferencePolarization.prepare(system:system,ao:ao,scf:model.scf,solvent:model.solvent,budget:model.budget)
            if let frame=model.eccFrame,let anchor=anchorCoefficients {
                let cross=try VivoGaussianIntegralEngine.crossOverlap(leftSystem:frame.anchorSystem,leftBasis:model.basis,rightSystem:system,rightBasis:model.basis,budget:model.budget)
                let aligned=try VivoOrbitalPathTransport.align(referenceCoefficients:anchor,currentCoefficients:prepared.coefficients,
                    crossAOOverlap:cross,groups:frame.transportGroups,minimumSingularValue:frame.minimumTransportSingularValue)
                let anchored=try aligned.rotation.multiplied(by:frame.anchorRotation)
                let h=try prepared.hamiltonian.rotated(by:anchored,budget:model.budget)
                let result=try VivoECCDMET.solve(h,configuration:frame.embedding,budget:model.budget)
                guard result.converged,result.frame.clusters.map({$0.coefficients.columns-$0.fragment.orbitals.count})==frame.expectedBathOrbitals else {
                    throw VivoChemistryError.convergence("nuclear ECC inner convergence or active-dimension change")
                }
                value=result.energyHartree
            } else {
                value=try VivoDirectCI.solve(prepared.hamiltonian,configuration:.init(residualTolerance:1e-12),budget:model.budget).roots[0].energyHartree
            }
        }
        guard value.isFinite else {throw VivoChemistryError.convergence("nonfinite nuclear energy")}
        if cache.count>=8192 {cache.removeAll(keepingCapacity:true)}
        cache[key]=value;return value
    }
    public func gradient(_ positions: [SIMD3<Double>], step: Double? = nil) throws -> VivoGeometryEvaluation {
        let delta=step ?? differences.gradientStepBohr,central=try energy(positions)
        guard delta.isFinite,delta>0 else {throw VivoChemistryError.invalid("nuclear derivative step")}
        var g=[SIMD3<Double>](repeating:.zero,count:positions.count)
        for i in positions.indices {for a in 0..<3 {
            var plus=positions,minus=positions;plus[i][a]+=delta;minus[i][a]-=delta
            g[i][a]=(try energy(plus)-energy(minus))/(2*delta)
        }}
        return .init(energyHartree:central,gradientHartreePerBohr:g)
    }
    public func hessian(_ positions: [SIMD3<Double>], step: Double? = nil) throws -> VivoQMMatrix {
        let delta=step ?? differences.hessianStepBohr,n=positions.count*3,e=try energy(positions)
        guard delta.isFinite,delta>0 else {throw VivoChemistryError.invalid("nuclear Hessian displacement")}
        _=try model.budget.elements([n,n],simultaneousArrays:4)
        var h=VivoQMMatrix(n,n)
        for i in 0..<n {
            var plus=positions,minus=positions;plus[i/3][i%3]+=delta;minus[i/3][i%3]-=delta
            h[i,i]=(try energy(plus)-2*e+energy(minus))/(delta*delta)
            for j in 0..<i {
                var pp=plus,pm=plus,mp=minus,mm=minus
                pp[j/3][j%3]+=delta;pm[j/3][j%3]-=delta;mp[j/3][j%3]+=delta;mm[j/3][j%3]-=delta
                let v=(try energy(pp)-energy(pm)-energy(mp)+energy(mm))/(4*delta*delta)
                h[i,j]=v;h[j,i]=v
            }
        }
        return h
    }
}
