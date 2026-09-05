import Foundation

public struct VivoTensorCCSDConfiguration:Codable,Sendable,Equatable {
    public var maximumIterations:Int
    public var residualTolerance:Double
    public var energyTolerance:Double
    public var diisHistory:Int
    public init(maximumIterations:Int=100,residualTolerance:Double=1e-8,energyTolerance:Double=1e-10,diisHistory:Int=7) {
        self.maximumIterations=maximumIterations;self.residualTolerance=residualTolerance;self.energyTolerance=energyTolerance;self.diisHistory=diisHistory
    }
}
public struct VivoTensorCCSDResult:Codable,Sendable,Equatable {
    public let converged:Bool
    public let energyHartree:Double
    public let referenceEnergyHartree:Double
    public let singles:VivoQMMatrix
    public let doubles:[Double]
    public let occupiedSpinModes:[Int]
    public let virtualSpinModes:[Int]
    public let residual:Double
    public let iterations:Int
    public let contractionOperations:Int
}
private struct CCTensor {
    let dims:[Int];var x:[Double]
    init(_ dims:[Int],_ x:[Double]?=nil) { self.dims=dims;self.x=x ?? .init(repeating:0,count:dims.reduce(1,*)) }
    func perm(_ axes:[Int])->Self {
        var out=Self(axes.map { dims[$0] })
        for i in x.indices {
            var rest=i,index=[Int](repeating:0,count:dims.count)
            for d in dims.indices.reversed() { index[d]=rest%dims[d];rest/=dims[d] }
            var k=0;for d in axes { k=k*dims[d]+index[d] };out.x[k]=x[i]
        };return out
    }
    func add(_ b:Self,_ scale:Double=1)->Self { precondition(dims==b.dims);return Self(dims,zip(x,b.x).map { $0+scale*$1 }) }
    func scale(_ s:Double)->Self { Self(dims,x.map { $0*s }) }
}
/// Conventional spin-orbital CCSD tensor equations, not determinant-space
/// exponentiation. O(N^6) contractions and O(N^4) tensors, including vvvv.
/// FP64 Einstein contractions are packed into the common Accelerate GEMM path.
/// This is not a local-correlation approximation or a low-memory DF-CCSD method.
public enum VivoTensorCCSD {
    public static func solve(_ h:VivoEmbeddedHamiltonian,configuration cfg:VivoTensorCCSDConfiguration = .init(),
                             budget:VivoChemistryBudget = .init()) throws -> VivoTensorCCSDResult {
        try h.validate(budget:budget)
        guard (1...10000).contains(cfg.maximumIterations),(2...20).contains(cfg.diisHistory),cfg.residualTolerance.isFinite,cfg.residualTolerance>0,
              cfg.energyTolerance.isFinite,cfg.energyTolerance>0 else { throw VivoChemistryError.invalid("tensor CCSD configuration") }
        let n=h.orbitalCount,occ=(0..<h.alphaElectrons).map { 2*$0 }+(0..<h.betaElectrons).map { 2*$0+1 }
        let vir=(h.alphaElectrons..<n).map { 2*$0 }+(h.betaElectrons..<n).map { 2*$0+1 },o=occ.count,v=vir.count
        _ = try budget.elements([max(1,o,v),max(1,o,v),max(1,o,v),max(1,o,v)],simultaneousArrays:32+2*cfg.diisHistory)
        func g(_ p:Int,_ q:Int,_ r:Int,_ s:Int)->Double {
            var value=0.0
            if p%2==r%2 && q%2==s%2 { value+=h.eri(p/2,r/2,q/2,s/2) }
            if p%2==s%2 && q%2==r%2 { value-=h.eri(p/2,s/2,q/2,r/2) };return value
        }
        let modes=occ+vir
        var f=VivoQMMatrix(o+v,o+v),reference=h.constantEnergyHartree
        for (i,p) in modes.enumerated() { for (j,q) in modes.enumerated() {
            f[i,j]=(p%2==q%2 ? h.oneElectron[p/2,q/2]:0)+occ.reduce(0.0) { $0+g(p,$1,q,$1) }
        } }
        for p in occ { reference+=h.oneElectron[p/2,p/2];for q in occ { reference+=0.5*g(p,q,p,q) } }
        if o==0 || v==0 { return .init(converged:true,energyHartree:reference,referenceEnergyHartree:reference,
            singles:try .init(rows:o,columns:v,values:[]),doubles:[],occupiedSpinModes:occ,virtualSpinModes:vir,residual:0,iterations:0,contractionOperations:0) }
        func block(_ labels:String)->CCTensor {
            let groups=Array(labels).map { $0=="o" ? occ:vir },dims=groups.map(\.count);var b=CCTensor(dims)
            for i in groups[0].indices { for j in groups[1].indices { for k in groups[2].indices { for l in groups[3].indices {
                b.x[((i*dims[1]+j)*dims[2]+k)*dims[3]+l]=g(groups[0][i],groups[1][j],groups[2][k],groups[3][l])
            } } } };return b
        }
        func fblock(_ rows:Range<Int>,_ columns:Range<Int>)->CCTensor {
            CCTensor([rows.count,columns.count],rows.flatMap { i in columns.map { f[i,$0] } })
        }
        let oooo=block("oooo"),ooov=block("ooov"),oovv=block("oovv"),ovov=block("ovov"),ovvv=block("ovvv"),vvvv=block("vvvv")
        let foo=fblock(0..<o,0..<o),fov=fblock(0..<o,o..<(o+v)),fvv=fblock(o..<(o+v),o..<(o+v))
        var operations=0
        func e(_ a:CCTensor,_ an:String,_ b:CCTensor,_ bn:String,_ on:String) throws -> CCTensor {
            let al=Array(an),bl=Array(bn),ol=Array(on),shared=al.filter { bl.contains($0) },left=al.filter { !shared.contains($0) },right=bl.filter { !shared.contains($0) }
            guard al.count==a.dims.count,bl.count==b.dims.count,Set(al).count==al.count,Set(bl).count==bl.count,
                  Set(ol)==Set(left+right),shared.allSatisfy({!ol.contains($0)}) else { throw VivoChemistryError.invalid("CC tensor contraction labels") }
            for label in shared { guard a.dims[al.firstIndex(of:label)!]==b.dims[bl.firstIndex(of:label)!] else { throw VivoChemistryError.invalid("CC contraction shape") } }
            let ap=a.perm((left+shared).map { al.firstIndex(of:$0)! }),bp=b.perm((shared+right).map { bl.firstIndex(of:$0)! })
            let m=left.reduce(1) { $0*a.dims[al.firstIndex(of:$1)!] },k=shared.reduce(1) { $0*a.dims[al.firstIndex(of:$1)!] },nn=right.reduce(1) { $0*b.dims[bl.firstIndex(of:$1)!] }
            let (count,overflow)=(m*k).multipliedReportingOverflow(by:nn)
            guard !overflow,count<=budget.maximumOperatorApplications-operations else { throw VivoChemistryError.resourceLimit("tensor CCSD aggregate work") };operations+=count
            let product=try VivoQMMatrix(rows:m,columns:k,values:ap.x).multiplied(by:VivoQMMatrix(rows:k,columns:nn,values:bp.x))
            let labels=left+right,dims=left.map { a.dims[al.firstIndex(of:$0)!] }+right.map { b.dims[bl.firstIndex(of:$0)!] }
            return CCTensor(dims,product.values).perm(ol.map { labels.firstIndex(of:$0)! })
        }
        var gap1=[Double](repeating:0,count:o*v),gap2=[Double](repeating:0,count:o*o*v*v)
        for i in 0..<o { for a in 0..<v {
            gap1[i*v+a]=f[i,i]-f[o+a,o+a]
            guard gap1[i*v+a] < -1e-8 else { throw VivoChemistryError.convergence("tensor CCSD inverted/near-zero gap") }
            for j in 0..<o { for b in 0..<v { gap2[((i*o+j)*v+a)*v+b]=f[i,i]+f[j,j]-f[o+a,o+a]-f[o+b,o+b] } }
        } }
        var t1=CCTensor([o,v]),t2=CCTensor([o,o,v,v],zip(oovv.x,gap2).map { $0/$1 })
        var history:[([Double],[Double])]=[],previous:Double?,residual=Double.infinity,energy=reference
        for iteration in 1...cfg.maximumIterations {
            let product=try e(t1,"ia",t1,"jb","ijab"),antisym=product.add(product.perm([0,1,3,2]),-1)
            let tau=t2.add(antisym),halfTau=t2.add(antisym,0.5)
            energy=reference+zip(t1.x,fov.x).reduce(0.0) { $0+$1.0*$1.1 }
                + 0.25*zip(t2.add(product,2).x,oovv.x).reduce(0.0) { $0+$1.0*$1.1 }
            var Fvv=try fvv.add(e(fov,"me",t1,"ma","ae"),-0.5).add(e(t1,"mf",ovvv.perm([1,0,3,2]),"amef","ae"))
                .add(e(halfTau,"mnaf",oovv,"mnef","ae"),-0.5)
            var Foo=try foo.add(e(fov,"me",t1,"ie","mi"),0.5).add(e(t1,"ne",ooov,"mnie","mi"))
                .add(e(halfTau,"inef",oovv,"mnef","mi"),0.5)
            let Fov=try fov.add(e(t1,"nf",oovv,"mnef","me"))
            let tmpO=try e(t1,"je",ooov,"mnie","mnij")
            let Woooo=try oooo.add(tmpO).add(tmpO.perm([0,1,3,2]),-1).add(e(tau,"ijef",oovv,"mnef","mnij"),0.25)
            let tmpV=try e(t1,"ma",ovvv,"mbef","abef")
            let Wvvvv=try vvvv.add(tmpV,-1).add(tmpV.perm([1,0,2,3])).add(e(tau,"mnab",oovv,"mnef","abef"),0.25)
            let t1pair=try e(t1,"jf",t1,"nb","jnfb")
            let Wovvo=try ovov.perm([0,1,3,2]).scale(-1).add(e(t1,"jf",ovvv,"mbef","mbej"))
                .add(e(t1,"nb",ooov.perm([0,1,3,2]),"mnej","mbej"))
                .add(e(t2,"jnfb",oovv,"mnef","mbej"),-0.5).add(e(t1pair,"jnfb",oovv,"mnef","mbej"),-1)
            for a in 0..<v { Fvv.x[a*v+a]-=f[o+a,o+a] };for i in 0..<o { Foo.x[i*o+i]-=f[i,i] }
            let s1=try fov.add(e(t1,"ie",Fvv,"ae","ia")).add(e(t1,"ma",Foo,"mi","ia"),-1)
                .add(e(t2,"imae",Fov,"me","ia")).add(e(t1,"nf",ovov,"naif","ia"),-1)
                .add(e(t2,"imef",ovvv,"maef","ia"),-0.5).add(e(t2,"mnae",ooov,"mnie","ia"),-0.5)
            let tmp1=try e(t2,"ijae",Fvv.add(e(t1,"mb",Fov,"me","be"),-0.5),"be","ijab")
            let tmp2=try e(t2,"imab",Foo.add(e(t1,"je",Fov,"me","mj"),0.5),"mj","ijab")
            var s2=try oovv.add(tmp1).add(tmp1.perm([0,1,3,2]),-1).add(tmp2,-1).add(tmp2.perm([1,0,2,3]))
                .add(e(tau,"mnab",Woooo,"mnij","ijab"),0.5).add(e(tau,"ijef",Wvvvv,"abef","ijab"),0.5)
            let mixed1=try e(t1,"ie",ovov,"mbje","mbji")
            var mixed=try e(t2,"imae",Wovvo,"mbej","ijab").add(e(t1,"ma",mixed1,"mbji","ijab"))
            mixed=mixed.add(mixed.perm([1,0,2,3]),-1);mixed=mixed.add(mixed.perm([0,1,3,2]),-1);s2=s2.add(mixed)
            let tmp3=try e(t1,"ie",ovvv,"jeba","ijab"),tmp4=try e(t1,"ma",ooov,"ijmb","ijab")
            s2=s2.add(tmp3).add(tmp3.perm([1,0,2,3]),-1).add(tmp4,-1).add(tmp4.perm([0,1,3,2]))
            let old=t1.x+t2.x,gaps=gap1+gap2,numerator=s1.x+s2.x
            let errors=zip(numerator,zip(gaps,old)).map { $0-$1.0*$1.1 }
            residual=errors.map(abs).max() ?? 0
            guard residual.isFinite,energy.isFinite else { throw VivoChemistryError.convergence("tensor CCSD nonfinite amplitudes/energy") }
            if residual<=cfg.residualTolerance,let previous,abs(energy-previous)<=cfg.energyTolerance {
                return .init(converged:true,energyHartree:energy,referenceEnergyHartree:reference,
                    singles:try .init(rows:o,columns:v,values:t1.x),doubles:t2.x,occupiedSpinModes:occ,virtualSpinModes:vir,
                    residual:residual,iterations:iteration,contractionOperations:operations)
            }
            if iteration==cfg.maximumIterations { break }
            var next=zip(numerator,gaps).map { $0/$1 }
            history.append((old,errors));if history.count>cfg.diisHistory { history.removeFirst() }
            if history.count>=2 {
                let m=history.count;var matrix=VivoQMMatrix(m+1,m+1),rhs=[Double](repeating:0,count:m+1);rhs[m] = -1
                let scale=max(1e-30,history.map { $0.1.reduce(0) { $0+$1*$1 } }.max() ?? 1)
                for i in 0..<m { for j in 0..<m { matrix[i,j]=zip(history[i].1,history[j].1).reduce(0) { $0+$1.0*$1.1 }/scale };matrix[i,m] = -1;matrix[m,i] = -1 }
                if let weights=try? VivoQMDenseAlgebra.solve(matrix,rhs:rhs),weights.allSatisfy({abs($0)<1e5}) {
                    next=[Double](repeating:0,count:old.count)
                    for i in 0..<m { for k in next.indices { next[k]+=weights[i]*(history[i].0[k]+history[i].1[k]/gaps[k]) } }
                }
            }
            guard next.allSatisfy(\.isFinite) else { throw VivoChemistryError.convergence("tensor CCSD update overflow") }
            previous=energy;t1=CCTensor([o,v],Array(next.prefix(o*v)));t2=CCTensor([o,o,v,v],Array(next.dropFirst(o*v)))
        }
        return .init(converged:false,energyHartree:energy,referenceEnergyHartree:reference,singles:try .init(rows:o,columns:v,values:t1.x),doubles:t2.x,
            occupiedSpinModes:occ,virtualSpinModes:vir,residual:residual,iterations:cfg.maximumIterations,contractionOperations:operations)
    }
}
