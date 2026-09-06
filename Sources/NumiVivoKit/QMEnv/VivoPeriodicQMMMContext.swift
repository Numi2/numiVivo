import Foundation

public struct VivoPeriodicQMMMConfiguration: Codable, Sendable, Equatable {
    public var ewald: VivoPeriodicElectrostaticConfiguration
    public var exactNearRadiusBohr: Double
    public var exactNearSwitchOffBohr: Double
    /// A model-admission center-separation bound, NOT a proof of negligible exchange.
    public var minimumQMImageSeparationBohr: Double
    public var polarization: VivoInducedDipoleConfiguration?
    public var reciprocalMesh: VivoMultipoleMeshConfiguration?
    public var scf: VivoSCFConfiguration
    public init(ewald: VivoPeriodicElectrostaticConfiguration = .init(), exactNearRadiusBohr: Double = 8,
                exactNearSwitchOffBohr: Double = 12,minimumQMImageSeparationBohr: Double = 8,
                scf: VivoSCFConfiguration = .init(), reciprocalMesh: VivoMultipoleMeshConfiguration? = nil,polarization: VivoInducedDipoleConfiguration? = nil) {
        self.ewald = ewald; self.exactNearRadiusBohr = exactNearRadiusBohr; self.exactNearSwitchOffBohr = exactNearSwitchOffBohr
        self.minimumQMImageSeparationBohr = minimumQMImageSeparationBohr; self.scf = scf; self.reciprocalMesh = reciprocalMesh; self.polarization = polarization
    }
    public func validate() throws {
        try ewald.validate(); try scf.validate(); try reciprocalMesh?.validate();try polarization?.validate()
        if reciprocalMesh != nil,scf.commutatorTolerance < 1e-6 {
            throw VivoChemistryError.unsupported("FP32 multipolar PME requires an explicit SCF commutator tolerance >=1e-6; use FP64 Ewald for tighter qualification")
        }
        guard exactNearRadiusBohr.isFinite,exactNearRadiusBohr > 0,exactNearSwitchOffBohr.isFinite,
              exactNearSwitchOffBohr > exactNearRadiusBohr,minimumQMImageSeparationBohr.isFinite,minimumQMImageSeparationBohr > 0,
              ewald.chargeConvention == .requireNeutral else {
            throw VivoChemistryError.invalid("periodic QM/MM near-field, image-separation or neutral-cell convention")
        }
    }
}
struct VivoPeriodicEmbeddingEvaluation {
    let energy: Double
    let fock: VivoQMMatrix
    let nuclearGradients: [VivoVector3D]
    let chargeGradients: [VivoVector3D]
    let affineStrain: VivoQMMatrix
    /// Derivative with respect to an induced dipole added at each MM center.
    let inducedDipoleDerivatives: [VivoVector3D]
    var polarization: VivoInducedDipoleResult? = nil
}

/// Exact Gaussian near-field / distributed quadrupolar far-field Hamiltonian.
/// The isolated primary QM electrostatics and the embedding-MM self channel are
/// subtracted from the periodic multipolar functional before near corrections.
/// Thus retained physical MM PME belongs to MD, not to this electronic term.
struct VivoPeriodicQMMMContext {
    struct NearPair {
        let nucleus: Int
        let charge: Int
        let imagePosition: SIMD3<Double>
        let delta: VivoVector3D
        let switching: Double
        let switchDerivative: Double
        let exactNuclearCoefficients: [Double]
        let exactAO: [VivoQMMatrix]
    }
    let integrals: VivoAOIntegrals
    let source: VivoElectronicSystem
    let cell: VivoPeriodicCell
    let configuration: VivoPeriodicQMMMConfiguration
    let momentOperators: VivoQMMultipoleOperators
    let mmSources: [VivoCartesianMultipole]
    let mmEwald: VivoPeriodicElectrostaticResult
    let nearPairs: [NearPair]
    let reciprocalOperator: VivoReciprocalElectrostaticOperator?
    let budget: VivoChemistryBudget

    init(source: VivoElectronicSystem,integrals: VivoAOIntegrals,cell: VivoPeriodicCell,
         configuration cfg: VivoPeriodicQMMMConfiguration,budget: VivoChemistryBudget,
         reciprocalOperator: VivoReciprocalElectrostaticOperator? = nil) throws {
        try cfg.validate();try source.validate();try budget.validate()
        guard cfg.reciprocalMesh == reciprocalOperator?.configuration else {
            throw VivoChemistryError.invalid("periodic reciprocal operator must match the fingerprinted mesh profile")
        }
        guard cell.isValid,integrals.sourceSystem.nuclei == source.nuclei,integrals.sourceSystem.pointCharges.isEmpty,
              integrals.sourceSystem.alphaElectrons == source.alphaElectrons,integrals.sourceSystem.betaElectrons == source.betaElectrons else {
            throw VivoChemistryError.invalid("periodic QM/MM source, cell or isolated integral identity")
        }
        let lattice = [cell.a,cell.b,cell.c].map { $0/VivoAtomicUnits.bohrInNM }
        let det = lattice[0].dot(lattice[1].cross(lattice[2]))
        let reciprocal = [lattice[1].cross(lattice[2])/det,lattice[2].cross(lattice[0])/det,lattice[0].cross(lattice[1])/det]
        let center = source.nuclei.map(\.positionBohr).reduce(SIMD3<Double>.zero,+)/Double(source.nuclei.count)
        let radius = source.nuclei.map { vivoQMNorm($0.positionBohr-center) }.max() ?? 0
        let minimumHeight = reciprocal.map { 1/$0.norm }.min() ?? 0
        guard minimumHeight-2*radius >= cfg.minimumQMImageSeparationBohr else {
            throw VivoChemistryError.unsupported("localized QM cluster violates its declared image-separation model bound")
        }
        let physicalCharge = source.nuclei.reduce(0.0) { $0+Double($1.atomicNumber) }-Double(source.alphaElectrons+source.betaElectrons)
            + source.pointCharges.reduce(0) { $0+$1.chargeE }
        guard abs(physicalCharge) <= cfg.ewald.neutralToleranceE else { throw VivoChemistryError.invalid("periodic QM/MM complete cell is not neutral") }
        let operators = try VivoQMMultipoleOperators(integrals: integrals,budget: budget)
        let mm = source.pointCharges.map { VivoCartesianMultipole(positionBohr: .init($0.positionBohr.x,$0.positionBohr.y,$0.positionBohr.z),chargeE: $0.chargeE) }
        var mmConfiguration = cfg.ewald;mmConfiguration.chargeConvention = .uniformNeutralizingBackground
        let mmValue = try VivoPeriodicElectrostatics.evaluate(sources: mm,cell: cell,configuration: mmConfiguration,reciprocalOperator: reciprocalOperator,budget: budget)
        let n = integrals.count
        var near: [NearPair] = [],work = 0
        let componentCount = cfg.polarization == nil ? 1 : 4
        for nucleus in source.nuclei.indices { for charge in source.pointCharges.indices where source.pointCharges[charge].chargeE != 0 || componentCount > 1 {
            let rn = source.nuclei[nucleus].positionBohr,rq = source.pointCharges[charge].positionBohr
            let raw = VivoVector3D(rn.x-rq.x,rn.y-rq.y,rn.z-rq.z)
            var reduced = raw
            let fractions = reciprocal.map { $0.dot(raw).rounded() }
            for axis in 0..<3 { reduced = reduced-lattice[axis]*fractions[axis] }
            var ranges: [ClosedRange<Int>] = [],candidateCount = 1
            for axis in 0..<3 {
                let c = reciprocal[axis].dot(reduced),bound = reciprocal[axis].norm*cfg.exactNearSwitchOffBohr
                let low = ceil(c-bound-1e-12),high = floor(c+bound+1e-12)
                guard low.isFinite,high.isFinite,low >= -4096,high <= 4096 else { throw VivoChemistryError.resourceLimit("QM/MM near-image range") }
                if low > high { candidateCount = 0;break }
                let product = candidateCount.multipliedReportingOverflow(by: Int(high-low)+1)
                guard !product.overflow,product.partialValue <= cfg.ewald.maximumImagePairs-work else { throw VivoChemistryError.resourceLimit("QM/MM exact near-image capacity") }
                candidateCount = product.partialValue;ranges.append(Int(low)...Int(high))
            }
            if candidateCount == 0 { continue };work += candidateCount
            for x in ranges[0] { for y in ranges[1] { for z in ranges[2] {
                let delta = reduced-lattice[0]*Double(x)-lattice[1]*Double(y)-lattice[2]*Double(z),r = delta.norm
                if r >= cfg.exactNearSwitchOffBohr { continue }
                guard r > 1e-8 else { throw VivoChemistryError.invalid("embedding charge overlaps QM center") }
                let image = rn-SIMD3<Double>(delta.x,delta.y,delta.z)
                let blend = Self.switching(r,on: cfg.exactNearRadiusBohr,off: cfg.exactNearSwitchOffBohr)
                var operators: [VivoQMMatrix] = []
                for component in 0..<componentCount {
                    var moment = [Double](repeating:0,count:4);moment[component]=1
                    var ao = VivoQMMatrix(n,n)
                    for p in 0..<n { for q in 0...p {
                        let a = integrals.orbitals[p],b = integrals.orbitals[q],ownership = VivoQMMultipoleOperators.ownership(nucleus,a,b)
                        if ownership == 0 { continue }
                        let ra = source.nuclei[a.nucleusIndex].positionBohr,rb = source.nuclei[b.nucleusIndex].positionBohr
                        let value = -ownership*VivoGaussianDerivative.contracted(a,b,ra: ra,rb: rb) { ea,la,aa,eb,lb,bb in
                            VivoGaussianDerivative.multipolarPotential(ea,la,aa,eb,lb,bb,center:image,moments:moment)
                        }
                        ao[p,q]=value;ao[q,p]=value
                    } }
                    operators.append(ao)
                }
                let bytes = (near.count+1).multipliedReportingOverflow(by: componentCount*n*n)
                guard !bytes.overflow else { throw VivoChemistryError.resourceLimit("near-field operator count") }
                _ = try budget.elements([bytes.partialValue],simultaneousArrays: 2)
                var nucleusMoment = [Double](repeating:0,count:10);nucleusMoment[0]=Double(source.nuclei[nucleus].atomicNumber)
                let nuclear = VivoPeriodicElectrostatics.pair(nucleusMoment,Array(repeating:0,count:10),tensor:try .screened(delta,alpha:0))
                near.append(.init(nucleus:nucleus,charge:charge,imagePosition:image,delta:delta,switching:blend.value,
                    switchDerivative:blend.derivative,exactNuclearCoefficients:Array(nuclear.right.prefix(componentCount)),exactAO:operators))
            } } }
        } }
        self.integrals = integrals;self.source = source;self.cell = cell;configuration = cfg;momentOperators = operators
        mmSources = mm;mmEwald = mmValue;nearPairs = near;self.budget = budget; self.reciprocalOperator = reciprocalOperator
    }
    static func switching(_ r: Double,on: Double,off: Double) -> (value: Double,derivative: Double) {
        if r <= on { return (1,0) };if r >= off { return (0,0) }
        let x = (r-on)/(off-on),x2 = x*x,x3 = x2*x
        return (1-10*x3+15*x3*x-6*x3*x2,(-30*x2+60*x3-30*x3*x)/(off-on))
    }
    func evaluate(density: VivoQMMatrix,derivatives: Bool,inducedDipoles: [VivoVector3D]? = nil) throws -> VivoPeriodicEmbeddingEvaluation {
        guard inducedDipoles == nil || (configuration.polarization != nil && inducedDipoles!.count == mmSources.count
            && inducedDipoles!.allSatisfy(\.isFinite)) else { throw VivoChemistryError.invalid("induced dipole source shape or missing near operators") }
        var activeMM = mmSources
        if let inducedDipoles { for i in activeMM.indices { activeMM[i].dipoleEBohr = inducedDipoles[i] } }

        let qm = try momentOperators.sources(density: density),nq = qm.count
        let periodic = try VivoPeriodicElectrostatics.evaluate(sources: qm+activeMM,cell: cell,configuration: configuration.ewald,reciprocalOperator: reciprocalOperator,budget: budget)
        var energy = periodic.energyHartree-mmEwald.energyHartree
        var lambda = Array(periodic.momentDerivatives.prefix(nq)),gn = periodic.forcesHartreePerBohr.prefix(nq).map { $0 * -1 }
        var gm = (0..<mmSources.count).map { periodic.forcesHartreePerBohr[nq+$0] * -1+mmEwald.forcesHartreePerBohr[$0] }
        var dipoleDerivatives = Array(periodic.momentDerivatives.dropFirst(nq)).map { VivoVector3D($0[1],$0[2],$0[3]) }
        var exactNearAO = VivoQMMatrix(integrals.count,integrals.count)
        var strain = try periodic.affineStrainDerivativeHartree.adding(mmEwald.affineStrainDerivativeHartree,scale: -1)
        func addStrain(_ gradient: VivoVector3D,_ displacement: VivoVector3D,scale: Double) {
            let g = [gradient.x,gradient.y,gradient.z],r = [displacement.x,displacement.y,displacement.z]
            for a in 0..<3 { for b in 0..<3 { strain[a,b] += scale*g[a]*r[b] } }
        }
        // Isolated primary QM Coulomb interactions already belong to the Gaussian
        // electronic Hamiltonian. Subtract their multipole approximation exactly once.
        for i in 0..<nq { for j in 0..<i {
            let d = qm[i].positionBohr-qm[j].positionBohr
            let pair = VivoPeriodicElectrostatics.pair(qm[i].moments,qm[j].moments,tensor: try .screened(d,alpha: 0))
            energy -= pair.energy
            for c in 0..<10 { lambda[i][c] -= pair.left[c];lambda[j][c] -= pair.right[c] }
            gn[i] = gn[i]-pair.gradient;gn[j] = gn[j]+pair.gradient
            addStrain(pair.gradient,d,scale: -1)
        } }
        for item in nearPairs {
            let i = item.nucleus,j = item.charge
            let pair = VivoPeriodicElectrostatics.pair(qm[i].moments,activeMM[j].moments,tensor: try .screened(item.delta,alpha: 0))
            let exactCoefficients = item.exactAO.indices.map { c in
                item.exactNuclearCoefficients[c]+zip(density.values,item.exactAO[c].values).reduce(0.0) { $0+$1.0*$1.1 }
            }
            let exact = zip(exactCoefficients,activeMM[j].moments).reduce(0) { $0+$1.0*$1.1 }
            for c in item.exactAO.indices { for slot in exactNearAO.values.indices {
                exactNearAO.values[slot] += item.switching*activeMM[j].moments[c]*item.exactAO[c].values[slot]
            } }
            if exactCoefficients.count == 4 {
                dipoleDerivatives[j] = dipoleDerivatives[j]+VivoVector3D(exactCoefficients[1]-pair.right[1],exactCoefficients[2]-pair.right[2],exactCoefficients[3]-pair.right[3])*item.switching
            }
            energy += item.switching*(exact-pair.energy)
            for c in 0..<10 { lambda[i][c] -= item.switching*pair.left[c] }
            if derivatives {
                let r = item.delta.norm
                let switchGradient = item.delta*((exact-pair.energy)*item.switchDerivative/r)
                let multipoleGradient = pair.gradient * -item.switching+switchGradient
                gn[i] = gn[i]+multipoleGradient;gm[j] = gm[j]-multipoleGradient
                addStrain(multipoleGradient,item.delta,scale: 1)
                var nuclearMoment = [Double](repeating:0,count:10);nuclearMoment[0]=Double(source.nuclei[i].atomicNumber)
                let nuclear = try VivoPeriodicElectrostatics.pair(nuclearMoment,activeMM[j].moments,tensor: .screened(item.delta,alpha:0)).gradient*item.switching
                gn[i] = gn[i]+nuclear;gm[j] = gm[j]-nuclear
                addStrain(nuclear,item.delta,scale: 1)
                for p in 0..<integrals.count { for q in 0..<integrals.count {
                    let a = integrals.orbitals[p],b = integrals.orbitals[q],ia = a.nucleusIndex,ib = b.nucleusIndex
                    let weight = -item.switching*density[p,q]*VivoQMMultipoleOperators.ownership(i,a,b)
                    if weight == 0 { continue }
                    let ra = source.nuclei[ia].positionBohr,rb = source.nuclei[ib].positionBohr
                    let integral: VivoGaussianDerivative.PrimitivePair = { ea,la,aa,eb,lb,bb in
                        VivoGaussianDerivative.multipolarPotential(ea,la,aa,eb,lb,bb,center:item.imagePosition,moments:Array(activeMM[j].moments.prefix(4)))
                    }
                    for axis in 0..<3 {
                        let ga = weight*VivoGaussianDerivative.pair(a,b,ra: ra,rb: rb,axis: axis,onFirst: true,operator: integral)
                        let gb = weight*VivoGaussianDerivative.pair(a,b,ra: ra,rb: rb,axis: axis,onFirst: false,operator: integral)
                        VivoGaussianDerivative.add(&gn,center: ia,axis: axis,value: ga)
                        VivoGaussianDerivative.add(&gn,center: ib,axis: axis,value: gb)
                        VivoGaussianDerivative.add(&gm,center: j,axis: axis,value: -ga-gb)
                        // Use the actual image coordinate, not the primary-cell
                        // charge coordinate, for the affine lattice derivative.
                        for axisB in 0..<3 { strain[axis,axisB] += ga*(ra[axisB]-item.imagePosition[axisB])+gb*(rb[axisB]-item.imagePosition[axisB]) }
                    }
                } }
            }
        }
        let potential = try momentOperators.fock(momentDerivatives: lambda).adding(exactNearAO)
        if derivatives {
            let response = try momentOperators.geometryGradient(density: density,momentDerivatives: lambda,budget: budget)
            for i in 0..<nq { gn[i] = gn[i]+response[i];addStrain(response[i],qm[i].positionBohr,scale: 1) }
        }
        guard energy.isFinite,potential.values.allSatisfy(\.isFinite),gn.allSatisfy(\.isFinite),gm.allSatisfy(\.isFinite),strain.values.allSatisfy(\.isFinite) else {
            throw VivoChemistryError.convergence("nonfinite periodic electronic embedding")
        }
        return .init(energy: energy,fock: potential,nuclearGradients: gn,chargeGradients: gm,affineStrain: strain,inducedDipoleDerivatives: dipoleDerivatives)
    }
}
