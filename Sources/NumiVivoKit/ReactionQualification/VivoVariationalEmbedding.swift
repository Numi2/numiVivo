import Foundation

public struct VivoVariationalEmbeddingRequest: Codable, Sendable, Equatable {
    public static let schema="numivivo.org/global-variational-embedding/v1"
    public let schema:String
    public let molecule:VivoCorrelatedSolventRequest
    public let fragments:[VivoFockFragment]
    public let space:VivoVariationalSpaceConfiguration
    /// Fixed gas-phase residual enrichment, before solvent iterations. The
    /// variational projector does not move with the density or applied field.
    public let gasResidualRounds:Int
    public let stationarityToleranceHartree:Double
    public init(molecule:VivoCorrelatedSolventRequest,fragments:[VivoFockFragment],
                space:VivoVariationalSpaceConfiguration = .init(),gasResidualRounds:Int=0,
                stationarityToleranceHartree:Double=1e-8) {
        schema=Self.schema;self.molecule=molecule;self.fragments=fragments;self.space=space
        self.gasResidualRounds=gasResidualRounds;self.stationarityToleranceHartree=stationarityToleranceHartree
    }
    public func validate() throws {
        try molecule.validate();try space.validate()
        guard schema==Self.schema,molecule.configuration.partition==nil,(1...32).contains(fragments.count),
              Set(fragments.map(\.identifier)).count==fragments.count,(0...16).contains(gasResidualRounds),
              stationarityToleranceHartree.isFinite,stationarityToleranceHartree>0,stationarityToleranceHartree<=1e-5 else {
            throw VivoChemistryError.invalid("global variational embedding schema, fragment identities or stationarity contract")
        }
    }
}
public struct VivoVariationalEmbeddingResult: Codable, Sendable, Equatable {
    public static let schema="numivivo.org/global-variational-embedding-result/v1"
    public let schema:String
    public let request:VivoVariationalEmbeddingRequest
    public let coefficients:VivoQMMatrix
    public let globalSubspace:VivoQMMatrix
    public let redundantSeedColumns:Int
    public let state:VivoCIState
    public let densityAO:VivoQMMatrix
    public let inputDensityAO:VivoQMMatrix
    public let occupations:[Double]
    public let variationalDimension:Int
    public let fullSectorDimension:Int
    public let gasEnergyHartree:Double
    public let energyHartree:Double
    public let equilibriumField:VivoSmoothCPCMResult
    public let selfConsistentProjectedResidualHartree:Double
    public let externalResidualHartree:Double
    public let history:[VivoCorrelatedSolventIteration]
    public let hamiltonianOperatorApplications:Int
    public let method:String
}

/// A well-defined global variational closure of overlapping orbital projectors.
/// Solvent responds to one normalized global CI state, with all interference
/// terms. This is NOT the nonvariational democratic ECC-DMET energy functional,
/// nor an assertion that a finite-subspace approximation is chemically accurate.
public enum VivoVariationalEmbedding {
    public static let method="coherent union of global CAS/Fock subspaces + variational CI + equilibrium smooth C-PCM; globally N-representable CI density; not democratic ECC-DMET or an accuracy/rate certificate; full determinant-sector vectors retained"
    private struct Workspace {
        let c:VivoQMMatrix
        let gas:VivoEmbeddedHamiltonian
        let pcm:VivoSmoothCPCMOperator
        init(_ r:VivoCorrelatedSolventRequest) throws {
            try r.validate()
            let ao=try VivoGaussianIntegralEngine.compute(system:r.system,basis:r.basis,budget:r.budget)
            if let c=r.coefficients { self.c=c }
            else {
                let eig=try VivoQMDenseAlgebra.symmetricEigen(ao.overlap)
                guard eig.values.first!>=1e-8 else { throw VivoChemistryError.invalid("global embedding AO rank loss") }
                var c=eig.vectors
                for i in 0..<ao.count { for j in 0..<ao.count { c[i,j]/=sqrt(eig.values[j]) } }
                self.c=c
            }
            guard c.rows==ao.count,c.columns==ao.count,c.values.allSatisfy(\.isFinite),
                  try ao.overlap.congruence(c).adding(.identity(ao.count),scale:-1).frobeniusNorm<1e-8 else {
                throw VivoChemistryError.invalid("global embedding molecular coefficient frame")
            }
            gas=try VivoEmbeddedHamiltonian.fromAO(ao,coefficients:c,alphaElectrons:r.system.alphaElectrons,
                betaElectrons:r.system.betaElectrons,orbitalIdentifiers:(0..<ao.count).map { "equilibrium-orbital-\($0)" },
                energyReference:"physical electronic energy; AO scalar once; equilibrium solvent polarization accounted separately",budget:r.budget)
            pcm=try .init(system:r.system,basis:r.basis,configuration:r.solvent,budget:r.budget)
        }
    }
    private static func effective(_ gas:VivoEmbeddedHamiltonian,_ potential:VivoQMMatrix) throws -> VivoEmbeddedHamiltonian {
        .init(orbitalIdentifiers:gas.orbitalIdentifiers,alphaElectrons:gas.alphaElectrons,betaElectrons:gas.betaElectrons,
            oneElectron:try gas.oneElectron.adding(potential),twoElectron:gas.twoElectron,
            constantEnergyHartree:gas.constantEnergyHartree,energyReference:gas.energyReference)
    }
    public static func run(_ request:VivoVariationalEmbeddingRequest) throws -> VivoVariationalEmbeddingResult {
        try request.validate()
        let r=request.molecule,w=try Workspace(r),cfg=r.configuration
        var space=try VivoFockSpace(hamiltonian:w.gas,configuration:request.space,budget:r.budget),work=0
        for fragment in request.fragments { try space.add(fragment,hamiltonian:w.gas) }
        let redundant=space.dependentColumns
        var initial=try space.solve(w.gas,work:&work)
        if request.gasResidualRounds>0 { for _ in 0..<request.gasResidualRounds {
            let added=try space.enrich(hamiltonians:[w.gas],states:[initial],work:&work)
            if added==0 { break }
            initial=try space.solve(w.gas,work:&work)
        } }
        var p=try VivoCIOneParticleDensity.spatial(initial.state,budget:r.budget)
        var previous:Double?,history:[VivoCorrelatedSolventIteration]=[]
        func ao(_ p:VivoQMMatrix) throws -> VivoQMMatrix { try w.c.multiplied(by:p).multiplied(by:w.c.transposed) }
        for iteration in 1...cfg.maximumIterations {
            let input=try ao(p),oldField=try w.pcm.evaluate(totalDensity:input)
            let v=try oldField.reactionPotentialMatrix.congruence(w.c)
            let solved=try space.solve(effective(w.gas,v),work:&work)
            let next=try VivoCIOneParticleDensity.spatial(solved.state,budget:r.budget),density=try ao(next)
            let field=try w.pcm.evaluate(totalDensity:density)
            let dError=try next.adding(p,scale:-1).frobeniusNorm
            let potential=try field.reactionPotentialMatrix.congruence(w.c)
            let vError=try potential.adding(v,scale:-1).frobeniusNorm
            let gas=solved.energyHartree-zip(next.values,v.values).reduce(0.0) { $0+$1.0*$1.1 }
            let energy=gas+field.polarizationEnergyHartree,delta=previous.map { abs($0-energy) }
            guard energy.isFinite else { throw VivoChemistryError.convergence("global embedding nonfinite energy") }
            history.append(.init(iteration:iteration,energyHartree:energy,densityResidual:dError,
                potentialResidualHartree:vError,energyChangeHartree:delta))
            // Inspect stationarity in the actual returned correlated field,
            // not merely the preceding field used by the eigensolver.
            let physical=try effective(w.gas,potential)
            let action=try VivoDirectHamiltonian(physical,determinants:space.determinants,budget:r.budget)
            let image=try action.apply(solved.state.coefficients,work:&work)
            let expectation=VivoFockSpace.dot(image,solved.state.coefficients)
            let residual=zip(image,solved.state.coefficients).map { $0-expectation*$1 }
            let inside=VivoFockSpace.norm(space.columns.map { VivoFockSpace.dot($0,residual) })
            let outside=VivoFockSpace.norm(space.removeProjection(residual))
            if let delta,delta<=cfg.energyToleranceHartree,dError<=cfg.densityTolerance,
               vError<=cfg.potentialToleranceHartree,inside<=request.stationarityToleranceHartree {
                let closure=try space.solve(physical,work:&work)
                let closureD=try VivoCIOneParticleDensity.spatial(closure.state,budget:r.budget)
                guard try closureD.adding(next,scale:-1).frobeniusNorm<=2*cfg.densityTolerance else {
                    throw VivoChemistryError.convergence("global fragment/solvent equilibrium closure")
                }
                let occupations=try VivoQMDenseAlgebra.symmetricEigen(next,tolerance:1e-13).values
                return .init(schema:VivoVariationalEmbeddingResult.schema,request:request,coefficients:w.c,
                    globalSubspace:space.matrix,redundantSeedColumns:redundant,state:solved.state,densityAO:density,inputDensityAO:input,
                    occupations:occupations,variationalDimension:space.dimension,fullSectorDimension:space.determinants.count,
                    gasEnergyHartree:gas,energyHartree:energy,equilibriumField:field,
                    selfConsistentProjectedResidualHartree:inside,externalResidualHartree:outside,
                    history:history,hamiltonianOperatorApplications:work,method:method)
            }
            previous=energy;p=try next.scaled(1-cfg.damping).adding(p,scale:cfg.damping)
        }
        throw VivoChemistryError.convergence("global variational embedding and solvent failed joint stationarity")
    }
    public static func validate(_ result:VivoVariationalEmbeddingResult,request:VivoVariationalEmbeddingRequest) throws {
        guard result.schema==VivoVariationalEmbeddingResult.schema,result.request==request,result.method==method else {
            throw VivoChemistryError.invalid("global embedding request or functional binding")
        }
        let rebuilt=try run(request)
        guard rebuilt==result else { throw VivoChemistryError.invalid("global embedding wavefunction, density, field or energy reconstruction") }
    }
    /// Import validated ECC fragment+bath orbital projectors. Occupied inactive
    /// columns are explicit physical assumptions; fractional environment RDMs
    /// are NOT rounded to integer occupations. All returned states live in the
    /// full global electron sector. The newly solved variational energy does
    /// not reuse or relabel the input ECC democratic energy.
    public static func fromECC(molecule:VivoCorrelatedSolventRequest,sourceHamiltonian:VivoEmbeddedHamiltonian,
                               result:VivoECCDMETResult,inactiveDoublyOccupiedColumns:[[Int]],
                               space:VivoVariationalSpaceConfiguration = .init()) throws -> VivoVariationalEmbeddingRequest {
        let w=try Workspace(molecule)
        guard inactiveDoublyOccupiedColumns.count==result.frame.clusters.count,
              sourceHamiltonian.orbitalCount==w.gas.orbitalCount,
              sourceHamiltonian.alphaElectrons==w.gas.alphaElectrons,sourceHamiltonian.betaElectrons==w.gas.betaElectrons,
              abs(sourceHamiltonian.constantEnergyHartree-w.gas.constantEnergyHartree)<1e-10,
              try sourceHamiltonian.oneElectron.adding(w.gas.oneElectron,scale:-1).frobeniusNorm<1e-10,
              sourceHamiltonian.twoElectron.count==w.gas.twoElectron.count,
              zip(sourceHamiltonian.twoElectron,w.gas.twoElectron).allSatisfy({abs($0-$1)<1e-10}) else {
            throw VivoChemistryError.invalid("ECC projector import does not match the molecular orbital Hamiltonian")
        }
        try VivoECCDMET.validate(result,hamiltonian:sourceHamiltonian,configuration:result.configuration,budget:molecule.budget)
        var fragments:[VivoFockFragment]=[]
        for (i,cluster) in result.frame.clusters.enumerated() {
            let size=cluster.coefficients.columns,core=inactiveDoublyOccupiedColumns[i]
            guard core.allSatisfy({$0>=size && $0<w.gas.orbitalCount}),Set(core).count==core.count,
                  w.gas.alphaElectrons-core.count==cluster.fragment.clusterAlphaElectrons,
                  w.gas.betaElectrons-core.count==cluster.fragment.clusterBetaElectrons else {
                throw VivoChemistryError.invalid("ECC inactive occupation assumptions are absent, fractional or inconsistent with the declared cluster population")
            }
            let rotation=try result.orbitalRotation.multiplied(by:cluster.completeFrame)
            fragments.append(.init(identifier:cluster.fragment.identifier,partition:.init(doublyOccupiedCore:core,active:Array(0..<size)),orbitalRotation:rotation))
        }
        let request=VivoVariationalEmbeddingRequest(molecule:molecule,fragments:fragments,space:space)
        try request.validate();return request
    }
    /// Two genuinely overlapping three-orbital projectors in H2/6-31G (16
    /// full-sector determinants). Neither projector alone is the full space.
    public static func template() -> VivoVariationalEmbeddingRequest {
        let whole=VivoBarrierBenchmarks.hydrogenExchange631G(),basis=VivoGaussianBasis(identifier:whole.basis.identifier,
            shells:whole.basis.shells.filter { $0.nucleusIndex<2 },source:whole.basis.source)
        let molecule=VivoCorrelatedSolventRequest(system:.init(nuclei:[.init(atomicNumber:1,positionBohr:.init(0,0,-0.7)),
            .init(atomicNumber:1,positionBohr:.init(0,0,0.7))],alphaElectrons:1,betaElectrons:1),basis:basis,
            solvent:.init(dielectricConstant:4,angularPoints:50))
        return .init(molecule:molecule,fragments:[.init(identifier:"left-subspace",partition:.init(active:[0,1,2])),
            .init(identifier:"right-subspace",partition:.init(active:[1,2,3]))],space:.init(maximumDimension:16))
    }
}
