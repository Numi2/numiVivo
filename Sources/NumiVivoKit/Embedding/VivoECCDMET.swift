import Foundation

public enum VivoECCDMETMode: String, Codable, Sendable {
    case singleFragment
    case selfConsistentPartition
}
public struct VivoECCDMETConfiguration: Codable, Sendable, Equatable {
    public var mode: VivoECCDMETMode
    public var fragments: [VivoECCFragment]
    public var correlationOperators: [VivoECCCorrelationOperator]
    public var matching: VivoECCSelfConsistencyConfiguration
    public var bathSelection: VivoECCBathSelection
    public var qioWeights: VivoQIOWeights
    /// Allowed orbital rotations are block diagonal in these explicitly declared
    /// locality classes. Empty means every orbital is fixed (not unrestricted).
    public var localityGroups: [[Int]]
    public var maximumOrbitalIterations: Int
    public var maximumFrameEvaluations: Int
    public var orbitalGradientTolerance: Double
    public var orbitalDifferenceStep: Double
    public var maximumOrbitalStep: Double
    public var energyToleranceHartree: Double
    public var electronTolerance: Double
    public init(mode:VivoECCDMETMode,fragments:[VivoECCFragment],correlationOperators:[VivoECCCorrelationOperator]=[],
                matching:VivoECCSelfConsistencyConfiguration = .init(referenceMethod:.fci),
                bathSelection:VivoECCBathSelection = .init(),qioWeights:VivoQIOWeights = .init(),localityGroups:[[Int]]=[],
                maximumOrbitalIterations:Int=40,maximumFrameEvaluations:Int=2000,orbitalGradientTolerance:Double=1e-5,
                orbitalDifferenceStep:Double=1e-4,maximumOrbitalStep:Double=0.15,energyToleranceHartree:Double=1e-8,
                electronTolerance:Double=1e-7) {
        self.mode=mode;self.fragments=fragments;self.correlationOperators=correlationOperators;self.matching=matching
        self.bathSelection=bathSelection;self.qioWeights=qioWeights;self.localityGroups=localityGroups
        self.maximumOrbitalIterations=maximumOrbitalIterations;self.maximumFrameEvaluations=maximumFrameEvaluations
        self.orbitalGradientTolerance=orbitalGradientTolerance;self.orbitalDifferenceStep=orbitalDifferenceStep
        self.maximumOrbitalStep=maximumOrbitalStep;self.energyToleranceHartree=energyToleranceHartree;self.electronTolerance=electronTolerance
    }
    public func validate(hamiltonian h:VivoEmbeddedHamiltonian,budget:VivoChemistryBudget = .init()) throws {
        try h.validate(budget:budget);try bathSelection.validate();try qioWeights.validate();try matching.numberMatching.validate()
        let n=h.orbitalCount,ids=fragments.map(\.identifier),indices=fragments.flatMap(\.orbitals),groups=localityGroups.flatMap{$0}
        guard (1...16).contains(fragments.count),Set(ids).count==ids.count,ids.allSatisfy({!$0.isEmpty}),
              Set(indices).count==indices.count,indices.allSatisfy({$0>=0 && $0<n}),
              fragments.allSatisfy({!$0.orbitals.isEmpty && $0.maximumBathOrbitals>=0 && $0.maximumBathOrbitals<=n-$0.orbitals.count &&
                $0.clusterAlphaElectrons>=0 && $0.clusterAlphaElectrons==$0.clusterBetaElectrons}),
              h.alphaElectrons==h.betaElectrons,(1...500).contains(maximumOrbitalIterations),(1...100000).contains(maximumFrameEvaluations),
              [orbitalGradientTolerance,orbitalDifferenceStep,maximumOrbitalStep,energyToleranceHartree,electronTolerance].allSatisfy({$0.isFinite && $0>0}),
              maximumOrbitalStep<=0.5,orbitalDifferenceStep<=0.01,
              (localityGroups.isEmpty || (groups.count==n && Set(groups)==Set(0..<n) && localityGroups.allSatisfy({!$0.isEmpty}))),
              matching.bathDiscardedWeight.isFinite,matching.bathDiscardedWeight>=0 else {
            throw VivoChemistryError.invalid("ECC-DMET frame, fragment, locality or convergence settings")
        }
        switch mode {
        case .singleFragment:
            guard fragments.count==1,correlationOperators.isEmpty else {throw VivoChemistryError.invalid("single-fragment ECC-DMET has no partition property-matching potential")}
        case .selfConsistentPartition:
            guard indices.count==n,!correlationOperators.isEmpty else {throw VivoChemistryError.invalid("self-consistent ECC-DMET requires a complete disjoint partition and explicit correlation operators")}
        }
    }
}
public struct VivoECCDMETIteration: Codable, Sendable, Equatable {
    public let iteration: Int
    public let leakageObjective: Double
    public let orbitalGradientNorm: Double
    public let momentResidualNorm: Double
    public let electronDefect: Double
    public let energyHartree: Double
    public let acceptedStep: Double
}
public struct VivoECCDMETResult: Codable, Sendable, Equatable {
    public static let schema="numivivo.org/ecc-dmet-result/v1"
    public let schema: String
    public let configuration: VivoECCDMETConfiguration
    public let orbitalRotation: VivoQMMatrix
    public let orbitalGradient: [Double]
    public let correlationPotentialHartree: [Double]
    public let matching: VivoECCSelfConsistencyResult?
    public let frame: VivoECCFrameResult
    public let iterations: [VivoECCDMETIteration]
    public let frameEvaluations: Int
    public let converged: Bool
    public let termination: String
    public var energyHartree: Double { frame.energyHartree }
    public let scope: String
}

/// Integrated correlated reference -> ranked/SVD bath -> locality-constrained QIO
/// -> impurity -> physical energy. Partition mode nests one-/two-body moment
/// fitting and global particle feedback inside every trial orbital frame.
/// This is a finite-basis algorithm; success is not chemical validation.
public enum VivoECCDMET {
    private struct Evaluation {
        let frame:VivoECCFrameResult
        let matching:VivoECCSelfConsistencyResult?
        var objective:Double {frame.clusters.reduce(0.0){$0+$1.qio.objective}/Double(frame.clusters.count)}
        var potential:[Double] {matching?.correlationPotentialHartree ?? []}
        var momentNorm:Double {matching?.momentResiduals.reduce(0.0,{hypot($0,$1)}) ?? 0}
    }
    private static func evaluate(_ h:VivoEmbeddedHamiltonian,configuration cfg:VivoECCDMETConfiguration,
                                 rotation:VivoQMMatrix,potential:[Double]?,budget:VivoChemistryBudget) throws -> Evaluation {
        let physical=try h.rotated(by:rotation,budget:budget)
        switch cfg.mode {
        case .singleFragment:
            let ref:VivoCIResult
            if cfg.matching.referenceMethod == .fci {ref=try VivoDirectCI.solve(physical,budget:budget).roots[0]}
            else {ref=try VivoConfigurationInteraction.solve(physical,method:.cisd,budget:budget)}
            let frame=try VivoECCClusterConstruction.evaluate(physical,reference:ref,fragments:cfg.fragments,
                partitionMatching:nil,selection:cfg.bathSelection,discardedWeight:cfg.matching.bathDiscardedWeight,
                qioWeights:cfg.qioWeights,budget:budget)
            return .init(frame:frame,matching:nil)
        case .selfConsistentPartition:
            let result=try VivoECCSelfConsistency.solve(physical,fragments:cfg.fragments,operators:cfg.correlationOperators,
                initialPotentialHartree:potential,configuration:cfg.matching,bathSelection:cfg.bathSelection,qioWeights:cfg.qioWeights,budget:budget)
            guard result.converged,let frame=result.frame else {throw VivoChemistryError.convergence("ECC-DMET correlation feedback: \(result.termination)")}
            return .init(frame:frame,matching:result)
        }
    }
    public static func solve(_ h:VivoEmbeddedHamiltonian,configuration cfg:VivoECCDMETConfiguration,
                             initialRotation:VivoQMMatrix?=nil,initialPotentialHartree:[Double]?=nil,
                             budget:VivoChemistryBudget = .init()) throws -> VivoECCDMETResult {
        try cfg.validate(hamiltonian:h,budget:budget)
        let n=h.orbitalCount,groups=cfg.localityGroups.isEmpty ? (0..<n).map{[$0]}:cfg.localityGroups
        let pairs=groups.flatMap{g in g.indices.flatMap{i in ((i+1)..<g.count).map{(g[i],g[$0])}}}
        var rotation=try initialRotation ?? VivoQMMatrix.identity(n)
        try validateRotation(rotation,groups:groups,n:n)
        var calls=0
        func trial(_ u:VivoQMMatrix,_ potential:[Double]?) throws -> Evaluation {
            guard calls<cfg.maximumFrameEvaluations else {throw VivoChemistryError.resourceLimit("ECC-DMET aggregate orbital-frame evaluation budget")}
            calls+=1;return try evaluate(h,configuration:cfg,rotation:u,potential:potential,budget:budget)
        }
        func delta(_ vector:[Double]) throws -> VivoQMMatrix {
            var generator=VivoQMMatrix(n,n)
            for (i,pair) in pairs.enumerated(){generator[pair.0,pair.1]=vector[i];generator[pair.1,pair.0] = -vector[i]}
            return try VivoQMDenseAlgebra.orbitalRotation(generator:generator)
        }
        var current=try trial(rotation,initialPotentialHartree),gradient:[Double]=[],history:[VivoECCDMETIteration]=[]
        var previousEnergy:Double?,acceptedStep=0.0,termination="orbital-iteration-limit",converged=false
        let budgetPerGradient=2*pairs.count
        for iteration in 0...cfg.maximumOrbitalIterations {
            if calls+budgetPerGradient>cfg.maximumFrameEvaluations {termination="frame-evaluation-limit";break}
            gradient=[Double](repeating:0,count:pairs.count)
            let dimensions=current.frame.clusters.map{$0.coefficients.columns}
            for i in pairs.indices {
                var step=[Double](repeating:0,count:pairs.count);step[i]=cfg.orbitalDifferenceStep
                let plus=try trial(rotation.multiplied(by:delta(step)),current.potential.isEmpty ? nil:current.potential)
                step[i] = -cfg.orbitalDifferenceStep
                let minus=try trial(rotation.multiplied(by:delta(step)),current.potential.isEmpty ? nil:current.potential)
                guard plus.frame.clusters.map({$0.coefficients.columns})==dimensions,minus.frame.clusters.map({$0.coefficients.columns})==dimensions else {
                    throw VivoChemistryError.convergence("ECC-DMET bath rank changes across an orbital derivative; enlarge the shared bath before optimization")
                }
                gradient[i]=(plus.objective-minus.objective)/(2*cfg.orbitalDifferenceStep)
            }
            let norm=gradient.reduce(0.0){hypot($0,$1)}
            guard norm.isFinite else {throw VivoChemistryError.convergence("ECC-DMET nonfinite orbital derivative")}
            history.append(.init(iteration:iteration,leakageObjective:current.objective,orbitalGradientNorm:norm,
                momentResidualNorm:current.momentNorm,electronDefect:current.frame.electronCountDefect,
                energyHartree:current.frame.energyHartree,acceptedStep:acceptedStep))
            if norm<=cfg.orbitalGradientTolerance {
                if abs(current.frame.electronCountDefect)>cfg.electronTolerance {termination="inactive-particle-closure-failed";break}
                if previousEnergy == nil || abs(current.frame.energyHartree-previousEnergy!)<=cfg.energyToleranceHartree {
                    converged=true;termination=pairs.isEmpty ? "converged-fixed-orbital-contract":"converged-orbital-moment-and-number-contract";break
                }
                previousEnergy=current.frame.energyHartree;continue
            }
            if iteration==cfg.maximumOrbitalIterations{break}
            let bound=min(1,cfg.maximumOrbitalStep/max(norm,1e-300)),direction=gradient.map{-$0*bound}
            let slope = -bound*norm*norm
            var scale=1.0,accepted=false
            for _ in 0..<20 {
                if calls>=cfg.maximumFrameEvaluations{termination="frame-evaluation-limit";break}
                let u=try rotation.multiplied(by:delta(direction.map{scale*$0}))
                let value=try trial(u,current.potential.isEmpty ? nil:current.potential)
                if value.frame.clusters.map({$0.coefficients.columns})==dimensions,value.objective<=current.objective+1e-4*scale*slope {
                    previousEnergy=current.frame.energyHartree;current=value;rotation=u;acceptedStep=scale*bound*norm;accepted=true;break
                }
                scale*=0.5
            }
            if !accepted {if termination != "frame-evaluation-limit"{termination="orbital-line-search-stalled"};break}
        }
        return .init(schema:VivoECCDMETResult.schema,configuration:cfg,orbitalRotation:rotation,orbitalGradient:gradient,
            correlationPotentialHartree:current.potential,matching:current.matching,frame:current.frame,iterations:history,
            frameEvaluations:calls,converged:converged,termination:termination,
            scope:"rank-aware 1-RDM bath with explicit information supplements; correlated CI reference and FCI impurities; locality-constrained QIO; partition one/two-moment and particle feedback; reference-restored physical electronic energy, not Gibbs energy")
    }
    private static func validateRotation(_ u:VivoQMMatrix,groups:[[Int]],n:Int) throws {
        guard u.rows==n,u.columns==n,u.values.allSatisfy(\.isFinite),
              try u.transposed.multiplied(by:u).adding(.identity(n),scale:-1).frobeniusNorm<1e-8 else {throw VivoChemistryError.invalid("ECC orbital rotation")}
        var labels=[Int](repeating:0,count:n)
        for (i,g) in groups.enumerated(){for p in g{labels[p]=i}}
        for p in 0..<n{for q in 0..<n where labels[p] != labels[q] && abs(u[p,q])>1e-10{throw VivoChemistryError.invalid("ECC rotation violates locality classes")}}
    }
    /// Rebuild the final correlated reference, baths, moment/particle feedback,
    /// impurity states and physical energies from the task inputs. Claimed output
    /// flags or residuals are not sufficient evidence for a cache hit.
    public static func validate(_ result:VivoECCDMETResult,hamiltonian h:VivoEmbeddedHamiltonian,
                                configuration cfg:VivoECCDMETConfiguration,budget:VivoChemistryBudget = .init()) throws {
        try cfg.validate(hamiltonian:h,budget:budget)
        let n=h.orbitalCount,groups=cfg.localityGroups.isEmpty ? (0..<n).map{[$0]}:cfg.localityGroups
        let count=groups.reduce(0){$0+$1.count*($1.count-1)/2}
        try validateRotation(result.orbitalRotation,groups:groups,n:n)
        guard result.schema==VivoECCDMETResult.schema,result.configuration==cfg,result.converged,
              result.orbitalGradient.count==count,result.orbitalGradient.allSatisfy(\.isFinite),
              result.orbitalGradient.reduce(0.0,{hypot($0,$1)})<=cfg.orbitalGradientTolerance,
              result.frame.energyHartree.isFinite,abs(result.frame.electronCountDefect)<=cfg.electronTolerance,
              result.correlationPotentialHartree.count==(cfg.mode == .singleFragment ? 0:cfg.correlationOperators.count),
              result.correlationPotentialHartree.allSatisfy({$0.isFinite && abs($0)<=cfg.matching.maximumPotentialHartree}) else {
            throw VivoChemistryError.invalid("ECC-DMET output convergence, particle or method binding")
        }
        let rebuilt=try evaluate(h,configuration:cfg,rotation:result.orbitalRotation,
            potential:result.correlationPotentialHartree.isEmpty ? nil:result.correlationPotentialHartree,budget:budget)
        guard abs(rebuilt.frame.energyHartree-result.energyHartree)<max(1e-8,10*cfg.energyToleranceHartree),
              abs(rebuilt.frame.electronCountDefect)<=cfg.electronTolerance,rebuilt.momentNorm<=cfg.matching.momentTolerance,
              rebuilt.frame.clusters.count==result.frame.clusters.count else {throw VivoChemistryError.invalid("ECC rebuilt energy, moments or particle count differ")}
        guard result.frame.states.count==rebuilt.frame.states.count,
              result.frame.energyContributions.count==rebuilt.frame.energyContributions.count,
              result.frame.energyConvention==rebuilt.frame.energyConvention,
              abs(result.frame.referencePhysicalEnergyHartree-rebuilt.frame.referencePhysicalEnergyHartree)<1e-8,
              abs(result.frame.reconstructedElectronCount-rebuilt.frame.reconstructedElectronCount)<cfg.electronTolerance,
              (cfg.mode == .singleFragment ? result.matching==nil && result.frame.numberMatching==nil : result.matching != nil && result.frame.numberMatching != nil),
              !result.iterations.isEmpty,result.iterations.count<=cfg.maximumOrbitalIterations+1,
              result.frameEvaluations>0,result.frameEvaluations<=cfg.maximumFrameEvaluations else {
            throw VivoChemistryError.invalid("ECC output energy decomposition, mode or iteration trace")
        }
        if let matching=rebuilt.matching {
            // A valid stored potential already satisfies the moments. Do not
            // silently repair a bad cached potential by iterating it to a root.
            guard matching.iterations.count==1,
                  zip(matching.correlationPotentialHartree,result.correlationPotentialHartree).allSatisfy({abs($0-$1)<1e-10}) else {
                throw VivoChemistryError.invalid("ECC stored correlation potential does not itself satisfy the moment equations")
            }
        }
        func sameState(_ a:VivoCIState,_ b:VivoCIState) throws {
            try a.validate(budget:budget);try b.validate(budget:budget)
            guard a.orbitalCount==b.orbitalCount,a.alphaElectrons==b.alphaElectrons,a.betaElectrons==b.betaElectrons else {
                throw VivoChemistryError.invalid("ECC output CI sector")
            }
            let index=Dictionary(uniqueKeysWithValues:a.determinants.enumerated().map{($0.element,$0.offset)})
            let overlap=b.determinants.enumerated().reduce(0.0){sum,item in
                sum+(index[item.element].map{a.coefficients[$0]*b.coefficients[item.offset]} ?? 0)
            }
            guard abs(abs(overlap)-1)<1e-7 else {throw VivoChemistryError.invalid("ECC output state differs from recomputed root")}
        }
        func near(_ a:[Double],_ b:[Double],tolerance:Double=1e-7) -> Bool {
            a.count==b.count && b.allSatisfy(\.isFinite) && zip(a,b).allSatisfy{abs($0-$1)<=tolerance}
        }
        func nearMatrix(_ a:VivoQMMatrix,_ b:VivoQMMatrix) -> Bool {
            a.rows==b.rows && a.columns==b.columns && near(a.values,b.values)
        }
        guard result.frame.reference.method==rebuilt.frame.reference.method,
              abs(result.frame.reference.energyHartree-rebuilt.frame.reference.energyHartree)<1e-8 else {
            throw VivoChemistryError.invalid("ECC auxiliary reference energy/method")
        }
        if let original=rebuilt.matching,let stored=result.matching {
            guard stored.converged,stored.numberMatching.converged,
                  near(original.momentResiduals,stored.momentResiduals),near(original.referenceMoments,stored.referenceMoments),
                  near(original.fragmentMoments,stored.fragmentMoments),near(original.correlationPotentialHartree,stored.correlationPotentialHartree),
                  abs(original.numberMatching.chemicalPotentialHartree-stored.numberMatching.chemicalPotentialHartree)<1e-7 else {
                throw VivoChemistryError.invalid("ECC cached moment/chemical-potential diagnostics differ")
            }
        }
        try sameState(rebuilt.frame.reference.state,result.frame.reference.state)
        for i in rebuilt.frame.clusters.indices {
            let a=rebuilt.frame.clusters[i],b=result.frame.clusters[i]
            guard a.fragment==b.fragment,a.coefficients.rows==b.coefficients.rows,a.coefficients.columns==b.coefficients.columns,
                  try a.coefficients.adding(b.coefficients,scale:-1).frobeniusNorm<1e-7,
                  a.spectralBathRank==b.spectralBathRank,a.supplementalBathRank==b.supplementalBathRank,
                  b.qio.objective.isFinite,abs(a.qio.objective-b.qio.objective)<1e-7,
                  abs(a.omittedCandidateCouplingSquared-b.omittedCandidateCouplingSquared)<1e-9,
                  abs(a.discardedCouplingSquared-b.discardedCouplingSquared)<1e-9 else {
                throw VivoChemistryError.invalid("ECC output bath frame, ranks or leakage differs from the recomputed correlated bath")
            }
            try b.physicalHamiltonian.validate(budget:budget)
            guard nearMatrix(a.completeFrame,b.completeFrame),nearMatrix(a.bareOneElectron,b.bareOneElectron),
                  nearMatrix(a.environmentPotential,b.environmentPotential),nearMatrix(a.physicalHamiltonian.oneElectron,b.physicalHamiltonian.oneElectron),
                  near(a.physicalHamiltonian.twoElectron,b.physicalHamiltonian.twoElectron),b.physicalHamiltonian.constantEnergyHartree==0,
                  a.physicalHamiltonian.alphaElectrons==b.physicalHamiltonian.alphaElectrons,a.physicalHamiltonian.betaElectrons==b.physicalHamiltonian.betaElectrons,
                  nearMatrix(a.referenceDensityInFrame.one,b.referenceDensityInFrame.one),near(a.referenceDensityInFrame.two,b.referenceDensityInFrame.two),
                  near(a.qio.entropyNats,b.qio.entropyNats),nearMatrix(a.qio.mutualInformation,b.qio.mutualInformation),
                  abs(rebuilt.frame.states[i].fragmentPopulation-result.frame.states[i].fragmentPopulation)<1e-7 else {
                throw VivoChemistryError.invalid("ECC cached projected Hamiltonian, densities or information values differ")
            }
            try sameState(a.referenceStateInFrame,b.referenceStateInFrame)
            try sameState(rebuilt.frame.states[i].biasedCI.state,result.frame.states[i].biasedCI.state)
            let original=rebuilt.frame.energyContributions[i],stored=result.frame.energyContributions[i]
            guard original.fragmentIdentifier==stored.fragmentIdentifier,
                  abs(original.referenceContributionHartree-stored.referenceContributionHartree)<1e-8,
                  abs(original.impurityContributionHartree-stored.impurityContributionHartree)<1e-8,
                  abs(rebuilt.frame.states[i].physicalClusterEnergyHartree-result.frame.states[i].physicalClusterEnergyHartree)<1e-8 else {
                throw VivoChemistryError.invalid("ECC output physical energy contributions")
            }
        }
        // Independently recompute the orbital derivative. A claimed small
        // gradient does not authenticate a stationary QIO orbital frame.
        let pairs=groups.flatMap{g in g.indices.flatMap{i in ((i+1)..<g.count).map{(g[i],g[$0])}}}
        var trueGradient:[Double]=[]
        for pair in pairs {
            var generator=VivoQMMatrix(n,n);generator[pair.0,pair.1]=cfg.orbitalDifferenceStep;generator[pair.1,pair.0] = -cfg.orbitalDifferenceStep
            let plusRotation=try result.orbitalRotation.multiplied(by:VivoQMDenseAlgebra.orbitalRotation(generator:generator))
            let minusRotation=try result.orbitalRotation.multiplied(by:VivoQMDenseAlgebra.orbitalRotation(generator:generator.scaled(-1)))
            let potential=result.correlationPotentialHartree.isEmpty ? nil:result.correlationPotentialHartree
            let plus=try evaluate(h,configuration:cfg,rotation:plusRotation,potential:potential,budget:budget)
            let minus=try evaluate(h,configuration:cfg,rotation:minusRotation,potential:potential,budget:budget)
            guard plus.frame.clusters.map({$0.coefficients.columns})==rebuilt.frame.clusters.map({$0.coefficients.columns}),
                  minus.frame.clusters.map({$0.coefficients.columns})==rebuilt.frame.clusters.map({$0.coefficients.columns}) else {
                throw VivoChemistryError.invalid("ECC cached orbital derivative crosses a bath-rank boundary")
            }
            trueGradient.append((plus.objective-minus.objective)/(2*cfg.orbitalDifferenceStep))
        }
        guard trueGradient.allSatisfy(\.isFinite),trueGradient.reduce(0.0,{hypot($0,$1)})<=1.01*cfg.orbitalGradientTolerance,
              zip(trueGradient,result.orbitalGradient).reduce(0.0,{hypot($0,$1.0-$1.1)})<=max(1e-8,0.1*cfg.orbitalGradientTolerance) else {
            throw VivoChemistryError.invalid("ECC cached QIO derivative is not stationary")
        }
    }
}
