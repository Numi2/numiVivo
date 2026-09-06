import Foundation

struct VivoInducedDipoleContext {
    struct Frame {
        let center: Int
        let zReference: Int?
        let xReference: Int?
        let axes: [VivoVector3D]
        let inversePrincipal: [Double]
        let zDisplacement: VivoVector3D
        let xDisplacement: VivoVector3D
        let projectedXLength: Double
        let inverseTensor: VivoQMMatrix
    }
    struct Pair {
        let first: Int
        let second: Int
        let displacement: VivoVector3D
        let field: Double
        let induction: Double
        let screeningPerBohr3: Double?
    }
    struct Correction {
        var energy=0.0
        var selfEnergy=0.0
        var dampingEnergy=0.0
        var channelEnergy=0.0
        var gradients: [VivoVector3D]
        var dipoleDerivatives: [VivoVector3D]
        var strain=VivoQMMatrix(3,3)
    }
    let configuration: VivoInducedDipoleConfiguration
    let chargeIndices: [Int]
    let frames: [Frame]
    let pairs: [Pair]
    let embedding: [VivoCartesianMultipole]
    let permanent: [VivoCartesianMultipole]
    let embeddingReference: VivoPeriodicElectrostaticResult
    let permanentReference: VivoPeriodicElectrostaticResult
    let differentChargeChannels: Bool
    let factor: VivoQMPositiveDefiniteFactor
    let cell: VivoPeriodicCell
    let ewald: VivoPeriodicElectrostaticConfiguration
    let reciprocalOperator: VivoReciprocalElectrostaticOperator?
    let budget: VivoChemistryBudget

    init(pointCharges: [VivoQMPointCharge],cell: VivoPeriodicCell,configuration cfg: VivoInducedDipoleConfiguration,
         ewald: VivoPeriodicElectrostaticConfiguration,reciprocalOperator: VivoReciprocalElectrostaticOperator?,
         physicalPermanentChargesE: [Double]? = nil,embeddingReference supplied: VivoPeriodicElectrostaticResult? = nil,
         budget: VivoChemistryBudget) throws {
        try cfg.validate();try budget.validate();try ewald.validate()
        if reciprocalOperator != nil,cfg.residualToleranceHartreePerEBohr<5e-6 {
            throw VivoChemistryError.unsupported("FP32 reciprocal polarization requires declared response tolerance >=5e-6; use direct FP64 for tighter response")
        }
        guard cell.isValid,pointCharges.allSatisfy({ $0.classicalParticleIndex != nil && vivoQMFinite($0.positionBohr) && $0.chargeE.isFinite }),
              Set(pointCharges.compactMap(\.classicalParticleIndex)).count==pointCharges.count,
              physicalPermanentChargesE == nil || (physicalPermanentChargesE!.count==pointCharges.count && physicalPermanentChargesE!.allSatisfy(\.isFinite)) else {
            throw VivoChemistryError.invalid("polarization requires unique mapped MM centers and a finite physical charge channel")
        }
        let index=Dictionary(uniqueKeysWithValues:pointCharges.enumerated().map { ($0.element.classicalParticleIndex!,$0.offset) })
        func resolve(_ id: UInt32) throws -> Int {
            guard let i=index[id] else { throw VivoChemistryError.invalid("polarizable site/frame/pair references absent MM particle \(id)") };return i
        }
        let mm=pointCharges.map { VivoCartesianMultipole(positionBohr:.init($0.positionBohr.x,$0.positionBohr.y,$0.positionBohr.z),chargeE:$0.chargeE) }
        var phys=mm
        if let physicalPermanentChargesE { for i in phys.indices { phys[i].chargeE=physicalPermanentChargesE[i] } }
        guard abs(phys.reduce(0) { $0+$1.chargeE }-mm.reduce(0) { $0+$1.chargeE })<1e-8 else {
            throw VivoChemistryError.invalid("physical and embedding charge channels must conserve total MM charge")
        }
        var evalCfg=ewald;evalCfg.chargeConvention = .uniformNeutralizingBackground
        let mmReference=try supplied ?? VivoPeriodicElectrostatics.evaluate(sources:mm,cell:cell,configuration:evalCfg,reciprocalOperator:reciprocalOperator,budget:budget)
        let different=phys.map(\.chargeE) != mm.map(\.chargeE)
        let physReference=try different ? VivoPeriodicElectrostatics.evaluate(sources:phys,cell:cell,configuration:evalCfg,reciprocalOperator:reciprocalOperator,budget:budget) : mmReference
        func displacement(_ a: Int,_ b: Int) throws -> VivoVector3D {
            try VivoMDPreparationGeometry.minimumImage((mm[a].positionBohr-mm[b].positionBohr)*VivoAtomicUnits.bohrInNM,cell:cell)/VivoAtomicUnits.bohrInNM
        }
        var frames: [Frame]=[],indices: [Int]=[]
        let volumeUnit=pow(VivoAtomicUnits.bohrInNM,3)
        for site in cfg.sites {
            let center=try resolve(site.particleIndex),alpha=site.principalPolarizabilitiesNM3/volumeUnit
            let inverse=[1/alpha.x,1/alpha.y,1/alpha.z]
            var axes=[VivoVector3D(1,0,0),.init(0,1,0),.init(0,0,1)]
            var zr: Int?,xr: Int?,dz=VivoVector3D.zero,dx=VivoVector3D.zero,projectionLength=1.0
            if let z=site.zParticle,let x=site.xParticle {
                zr=try resolve(z);xr=try resolve(x);dz=try displacement(zr!,center);dx=try displacement(xr!,center)
                guard dz.norm>1e-8 else { throw VivoChemistryError.invalid("degenerate polarizability z frame") }
                let uz=dz/dz.norm,projected=dx-uz*uz.dot(dx)
                projectionLength=projected.norm
                guard projectionLength>1e-8 else { throw VivoChemistryError.invalid("collinear polarizability frame") }
                let ux=projected/projectionLength;axes=[ux,uz.cross(ux),uz]
            }
            var tensor=VivoQMMatrix(3,3)
            for axis in 0..<3 {
                let e=[axes[axis].x,axes[axis].y,axes[axis].z]
                for a in 0..<3 { for b in 0..<3 { tensor[a,b]+=inverse[axis]*e[a]*e[b] } }
            }
            indices.append(center)
            frames.append(.init(center:center,zReference:zr,xReference:xr,axes:axes,inversePrincipal:inverse,
                zDisplacement:dz,xDisplacement:dx,projectedXLength:projectionLength,inverseTensor:tensor))
        }
        let pairs=try cfg.pairs.map { pair -> Pair in
            let a=try resolve(pair.firstParticle),b=try resolve(pair.secondParticle),d=try displacement(a,b)
            guard d.norm>1e-8 else { throw VivoChemistryError.invalid("overlapping polarization pair centers") }
            return .init(first:a,second:b,displacement:d,field:pair.permanentFieldScale,induction:pair.mutualInductionScale,
                         screeningPerBohr3:pair.screeningPerNM3.map { $0*volumeUnit })
        }
        let dimensions=3*frames.count
        _ = try budget.elements([dimensions,dimensions],simultaneousArrays:4)
        var matrix=VivoQMMatrix(dimensions,dimensions)
        let pure=indices.map { VivoCartesianMultipole(positionBohr:mm[$0].positionBohr,chargeE:0) }
        let modeCount=ewald.reciprocalHalfWidths.reduce(1) { $0*(2*$1+1) }
        let cost=dimensions.multipliedReportingOverflow(by:pure.count*modeCount)
        guard !cost.overflow,cost.partialValue<=budget.maximumOperatorApplications else {
            throw VivoChemistryError.resourceLimit("bounded polarization response-matrix work; no silent site truncation")
        }
        for column in 0..<dimensions {
            var vector=pure
            let axis=column%3,site=column/3
            vector[site].dipoleEBohr = axis==0 ? .init(1,0,0) : (axis==1 ? .init(0,1,0):.init(0,0,1))
            let result=try VivoPeriodicElectrostatics.evaluate(sources:vector,cell:cell,configuration:evalCfg,reciprocalOperator:reciprocalOperator,budget:budget)
            for row in 0..<dimensions { matrix[row,column]=result.momentDerivatives[row/3][1+row%3] }
        }
        for site in frames.indices { for a in 0..<3 { for b in 0..<3 {
            matrix[3*site+a,3*site+b]+=frames[site].inverseTensor[a,b]
        } } }
        let polarIndex=Dictionary(uniqueKeysWithValues:indices.enumerated().map { ($0.element,$0.offset) })
        for pair in pairs {
            guard let i=polarIndex[pair.first],let j=polarIndex[pair.second] else { continue }
            let c=Self.coefficients(pair),r=[pair.displacement.x,pair.displacement.y,pair.displacement.z]
            for a in 0..<3 { for b in 0..<3 {
                let value=(a==b ? c.b:0)+c.c*r[a]*r[b]
                matrix[3*i+a,3*j+b]+=value;matrix[3*j+b,3*i+a]+=value
            } }
        }
        configuration=cfg;chargeIndices=indices;self.frames=frames;self.pairs=pairs;embedding=mm;permanent=phys
        embeddingReference=mmReference;permanentReference=physReference;differentChargeChannels=different
        factor=try VivoQMPositiveDefiniteFactor(matrix,pivotFloor:cfg.minimumResponsePivot,symmetryTolerance:reciprocalOperator == nil ? 1e-9:5e-6)
        self.cell=cell;self.ewald=evalCfg;self.reciprocalOperator=reciprocalOperator;self.budget=budget
    }
    private static func coefficients(_ pair: Pair) -> (a:Double,b:Double,c:Double,da:Double,db:Double,dc:Double) {
        let r=pair.displacement.norm,r2=r*r,r3=r2*r,r5=r3*r2
        var f3=1.0,f5=1.0,d3=0.0,d5=0.0
        if let screening=pair.screeningPerBohr3 {
            let v=screening*r3,expv=exp(-v),dv=3*screening*r2
            if v < 50 { f3 = -expm1(-v);f5=f3-v*expv;d3=expv*dv;d5=v*expv*dv }
        }
        let a=(pair.field*f3-1)/r3,b=(pair.induction*f3-1)/r3,c = -3*(pair.induction*f5-1)/r5
        return (a,b,c,pair.field*d3/r3-3*a/r,pair.induction*d3/r3-3*b/r,-3*pair.induction*d5/r5-5*c/r)
    }
    func response(linearDipoleDerivatives: [VivoVector3D]) throws -> [VivoVector3D] {
        guard linearDipoleDerivatives.count==embedding.count,linearDipoleDerivatives.allSatisfy(\.isFinite) else {
            throw VivoChemistryError.invalid("polarization driving field shape")
        }
        var rhs: [Double]=[]
        for i in chargeIndices { let g=linearDipoleDerivatives[i];rhs += [-g.x,-g.y,-g.z] }
        let solution=try factor.solve(rhs)
        var dipoles=[VivoVector3D](repeating:.zero,count:embedding.count)
        for (slot,index) in chargeIndices.enumerated() { dipoles[index] = .init(solution[3*slot],solution[3*slot+1],solution[3*slot+2]) }
        return dipoles
    }
    func correction(dipoles: [VivoVector3D],derivatives: Bool) throws -> Correction {
        guard dipoles.count==embedding.count,dipoles.allSatisfy(\.isFinite) else { throw VivoChemistryError.invalid("polarization moment dimensions") }
        var out=Correction(gradients:Array(repeating:.zero,count:embedding.count),dipoleDerivatives:Array(repeating:.zero,count:embedding.count))
        func strain(_ gradient: VivoVector3D,_ d: VivoVector3D) {
            let g=[gradient.x,gradient.y,gradient.z],r=[d.x,d.y,d.z]
            for a in 0..<3 { for b in 0..<3 { out.strain[a,b]+=g[a]*r[b] } }
        }
        if differentChargeChannels {
            var phys=permanent,embed=embedding
            for i in dipoles.indices { phys[i].dipoleEBohr=dipoles[i];embed[i].dipoleEBohr=dipoles[i] }
            let p=try VivoPeriodicElectrostatics.evaluate(sources:phys,cell:cell,configuration:ewald,reciprocalOperator:reciprocalOperator,budget:budget)
            let e=try VivoPeriodicElectrostatics.evaluate(sources:embed,cell:cell,configuration:ewald,reciprocalOperator:reciprocalOperator,budget:budget)
            out.channelEnergy=p.energyHartree-e.energyHartree-permanentReference.energyHartree+embeddingReference.energyHartree
            out.strain=try p.affineStrainDerivativeHartree.adding(e.affineStrainDerivativeHartree,scale:-1)
                .adding(permanentReference.affineStrainDerivativeHartree,scale:-1).adding(embeddingReference.affineStrainDerivativeHartree)
            for i in dipoles.indices {
                out.gradients[i]=p.forcesHartreePerBohr[i] * -1+e.forcesHartreePerBohr[i]+permanentReference.forcesHartreePerBohr[i]-embeddingReference.forcesHartreePerBohr[i]
                let v=zip(p.momentDerivatives[i],e.momentDerivatives[i]).map(-)
                out.dipoleDerivatives[i] = .init(v[1],v[2],v[3])
            }
        }
        for frame in frames {
            let mu=dipoles[frame.center],components=frame.axes.map { $0.dot(mu) }
            let local=zip(components,frame.inversePrincipal).map(*)
            out.selfEnergy+=0.5*zip(components,local).reduce(0) { $0+$1.0*$1.1 }
            out.dipoleDerivatives[frame.center]=out.dipoleDerivatives[frame.center]+zip(frame.axes,local).reduce(VivoVector3D.zero) { $0+$1.0*$1.1 }
            if derivatives,let zi=frame.zReference,let xi=frame.xReference {
                let x=frame.axes[0],z=frame.axes[2],gy=mu*local[1]
                let gx=mu*local[0]+gy.cross(z)
                var gz=mu*local[2]+x.cross(gy)
                let gt=(gx-x*x.dot(gx))/frame.projectedXLength
                let gdx=gt-z*z.dot(gt)
                gz=gz-gt*z.dot(frame.xDisplacement)-frame.xDisplacement*z.dot(gt)
                let gdz=(gz-z*z.dot(gz))/frame.zDisplacement.norm
                out.gradients[xi]=out.gradients[xi]+gdx;out.gradients[zi]=out.gradients[zi]+gdz
                out.gradients[frame.center]=out.gradients[frame.center]-gdx-gdz
                strain(gdx,frame.xDisplacement);strain(gdz,frame.zDisplacement)
            }
        }
        for pair in pairs {
            let i=pair.first,j=pair.second,d=pair.displacement,mi=dipoles[i],mj=dipoles[j],c=Self.coefficients(pair)
            let vector=mj*permanent[i].chargeE-mi*permanent[j].chargeE
            let cross=vector.dot(d),dot=mi.dot(mj),ri=mi.dot(d),rj=mj.dot(d),r=d.norm
            out.dampingEnergy += cross*c.a+dot*c.b+ri*rj*c.c
            out.dipoleDerivatives[i]=out.dipoleDerivatives[i]-d*(permanent[j].chargeE*c.a)+mj*c.b+d*(rj*c.c)
            out.dipoleDerivatives[j]=out.dipoleDerivatives[j]+d*(permanent[i].chargeE*c.a)+mi*c.b+d*(ri*c.c)
            if derivatives {
                let gradient=vector*c.a+d*((cross*c.da+dot*c.db+ri*rj*c.dc)/r)+(mi*rj+mj*ri)*c.c
                out.gradients[i]=out.gradients[i]+gradient;out.gradients[j]=out.gradients[j]-gradient;strain(gradient,d)
            }
        }
        out.energy=out.selfEnergy+out.dampingEnergy+out.channelEnergy
        guard out.energy.isFinite,out.gradients.allSatisfy(\.isFinite),out.dipoleDerivatives.allSatisfy(\.isFinite),out.strain.values.allSatisfy(\.isFinite) else {
            throw VivoChemistryError.convergence("nonfinite polarization functional")
        }
        return out
    }
    func summary(dipoles: [VivoVector3D],derivatives: [VivoVector3D],correction: Correction) throws -> VivoInducedDipoleResult {
        let residual=chargeIndices.map { derivatives[$0].norm }.max() ?? 0
        guard residual.isFinite,residual<=configuration.residualToleranceHartreePerEBohr else {
            throw VivoChemistryError.convergence("mutual polarization did not converge to its declared response residual: \(residual)")
        }
        return .init(model:configuration.model,moments:zip(configuration.sites,chargeIndices).map { .init(particleIndex:$0.0.particleIndex,dipoleEBohr:dipoles[$0.1]) },
            maximumResponseResidual:residual,minimumResponsePivot:factor.minimumPivot,selfEnergyHartree:correction.selfEnergy,
            dampingCorrectionEnergyHartree:correction.dampingEnergy,physicalChargeChannelCorrectionHartree:correction.channelEnergy,
            physicalPermanentChargesE:permanent.map(\.chargeE))
    }
    /// Eliminate the stationary dipoles for each AO density. This is mutually
    /// self-consistent mean-field induction, NOT a direct-reaction-field model.
    func evaluate(context: VivoPeriodicQMMMContext,density: VivoQMMatrix,derivatives: Bool) throws -> VivoPeriodicEmbeddingEvaluation {
        let zero=[VivoVector3D](repeating:.zero,count:embedding.count)
        let initial=try context.evaluate(density:density,derivatives:false,inducedDipoles:zero)
        let initialCorrection=try correction(dipoles:zero,derivatives:false)
        let dipoles=try response(linearDipoleDerivatives:zip(initial.inducedDipoleDerivatives,initialCorrection.dipoleDerivatives).map(+))
        let evaluated=try context.evaluate(density:density,derivatives:derivatives,inducedDipoles:dipoles)
        let correction=try correction(dipoles:dipoles,derivatives:derivatives)
        let gradient=zip(evaluated.inducedDipoleDerivatives,correction.dipoleDerivatives).map(+)
        let state=try summary(dipoles:dipoles,derivatives:gradient,correction:correction)
        return .init(energy:evaluated.energy+correction.energy,fock:evaluated.fock,nuclearGradients:evaluated.nuclearGradients,
            chargeGradients:zip(evaluated.chargeGradients,correction.gradients).map(+),affineStrain:try evaluated.affineStrain.adding(correction.strain),
            inducedDipoleDerivatives:gradient,polarization:state)
    }
}

public struct VivoPolarizationEnergyForces: Codable, Sendable, Equatable {
    public let energyHartree: Double
    public let particleIndices: [UInt32]
    public let centerForcesHartreePerBohr: [VivoVector3D]
    public let affineStrainDerivativeHartree: VivoQMMatrix
    public let response: VivoInducedDipoleResult
}
public enum VivoInducedDipoles {
    /// Classical mutual polarization correction to a permanent-charge Ewald
    /// Hamiltonian. Permanent self electrostatics remain owned by the MD engine.
    public static func evaluate(pointCharges: [VivoQMPointCharge],cell: VivoPeriodicCell,
                                configuration: VivoInducedDipoleConfiguration,ewald: VivoPeriodicElectrostaticConfiguration,
                                reciprocalOperator: VivoReciprocalElectrostaticOperator? = nil,
                                budget: VivoChemistryBudget = .init()) throws -> VivoPolarizationEnergyForces {
        let context=try VivoInducedDipoleContext(pointCharges:pointCharges,cell:cell,configuration:configuration,ewald:ewald,
            reciprocalOperator:reciprocalOperator,budget:budget)
        let zero=[VivoVector3D](repeating:.zero,count:pointCharges.count)
        let correction0=try context.correction(dipoles:zero,derivatives:false)
        let linear=context.permanentReference.momentDerivatives.map { VivoVector3D($0[1],$0[2],$0[3]) }
        let dipoles=try context.response(linearDipoleDerivatives:zip(linear,correction0.dipoleDerivatives).map(+))
        var sources=context.permanent
        for i in sources.indices { sources[i].dipoleEBohr=dipoles[i] }
        let periodic=try VivoPeriodicElectrostatics.evaluate(sources:sources,cell:cell,configuration:context.ewald,reciprocalOperator:reciprocalOperator,budget:budget)
        let correction=try context.correction(dipoles:dipoles,derivatives:true)
        let derivatives=zip(periodic.momentDerivatives,correction.dipoleDerivatives).map { VivoVector3D($0.0[1],$0.0[2],$0.0[3])+$0.1 }
        let response=try context.summary(dipoles:dipoles,derivatives:derivatives,correction:correction)
        let forces=sources.indices.map { periodic.forcesHartreePerBohr[$0]-context.permanentReference.forcesHartreePerBohr[$0]-correction.gradients[$0] }
        return .init(energyHartree:periodic.energyHartree-context.permanentReference.energyHartree+correction.energy,
            particleIndices:pointCharges.map { $0.classicalParticleIndex! },centerForcesHartreePerBohr:forces,
            affineStrainDerivativeHartree:try periodic.affineStrainDerivativeHartree.adding(context.permanentReference.affineStrainDerivativeHartree,scale:-1).adding(correction.strain),response:response)
    }
}
