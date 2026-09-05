import Foundation

public struct VivoDavidsonConfiguration:Codable,Sendable,Equatable {
    public var roots:Int
    public var maximumIterations:Int
    public var maximumSubspace:Int
    public var residualTolerance:Double
    public init(roots:Int=1,maximumIterations:Int=150,maximumSubspace:Int=40,residualTolerance:Double=1e-9) {
        self.roots=roots;self.maximumIterations=maximumIterations;self.maximumSubspace=maximumSubspace;self.residualTolerance=residualTolerance
    }
}
public struct VivoCISpectrum:Codable,Sendable,Equatable {
    public let roots:[VivoCIResult]
    public let iterations:Int
    public let matrixVectorProducts:Int
    public let operatorApplications:Int
}
/// Connected Slater-Condon actions. No determinant-square Hamiltonian is stored.
/// Full CI still has combinatorial state-space growth; matrix-free is not DMRG.
public enum VivoDirectCI {
    static func determinants(n:Int,na:Int,nb:Int,budget:VivoChemistryBudget) throws -> [UInt64] {
        guard (1...31).contains(n),(0...n).contains(na),(0...n).contains(nb) else { throw VivoChemistryError.invalid("CI electron sector") }
        func count(_ k:Int) throws -> Int {
            let k=min(k,n-k);var c=1
            if k>0 { for i in 1...k {
                let next=c.multipliedReportingOverflow(by:n-k+i)
                guard !next.overflow else { throw VivoChemistryError.resourceLimit("CI count overflow") }
                c=next.partialValue/i;guard c<=budget.maximumDeterminants else { throw VivoChemistryError.resourceLimit("CI determinant capacity") }
            } };return c
        }
        let ca=try count(na),cb=try count(nb),size=ca.multipliedReportingOverflow(by:cb)
        guard !size.overflow,size.partialValue<=budget.maximumDeterminants else { throw VivoChemistryError.resourceLimit("CI determinant capacity") }
        func combinations(_ population:Int,_ spin:Int)->[UInt64] {
            var values:[UInt64]=[]
            func visit(_ first:Int,_ remaining:Int,_ bits:UInt64) {
                if remaining==0 { values.append(bits);return };if first>n-remaining { return }
                for p in first...(n-remaining) { visit(p+1,remaining-1,bits | (UInt64(1)<<(2*p+spin))) }
            };visit(0,population,0);return values
        }
        let a=combinations(na,0),b=combinations(nb,1);return a.flatMap { x in b.map { x | $0 } }.sorted()
    }
    public static func solve(_ h:VivoEmbeddedHamiltonian,configuration cfg:VivoDavidsonConfiguration = .init(),
                             budget:VivoChemistryBudget = .init()) throws -> VivoCISpectrum {
        try h.validate(budget:budget);try budget.validate()
        guard (1...32).contains(cfg.roots),cfg.maximumIterations>0,cfg.maximumIterations<=10000,
              cfg.maximumSubspace>=2*cfg.roots+2,cfg.maximumSubspace<=1024,
              cfg.residualTolerance.isFinite,cfg.residualTolerance>0 else { throw VivoChemistryError.invalid("Davidson settings") }
        let n=h.orbitalCount,dets=try determinants(n:n,na:h.alphaElectrons,nb:h.betaElectrons,budget:budget)
        let d=dets.count,roots=cfg.roots,capacity=min(d,cfg.maximumSubspace)
        guard roots<=d else { throw VivoChemistryError.invalid("more roots than determinants") }
        _ = try budget.elements([d,capacity],simultaneousArrays:8)
        let action=try VivoDirectHamiltonian(h,determinants:dets,budget:budget)
        let diagonal=action.diagonal
        var work=0,products=0
        func apply(_ vector:[Double]) throws -> [Double] {
            products+=1
            return try action.apply(vector,work:&work)
        }
        func dot(_ a:[Double],_ b:[Double])->Double { zip(a,b).reduce(0.0) { $0+$1.0*$1.1 } }
        func orthogonal(_ a:[Double],_ space:[[Double]])->[Double]? {
            var q=a
            for _ in 0..<2 { for b in space { let s=dot(q,b);for k in q.indices { q[k]-=s*b[k] } } }
            let norm=q.reduce(0.0) { hypot($0,$1) };return norm>1e-12 && norm.isFinite ? q.map { $0/norm }:nil
        }
        var vectors:[[Double]]=[],images:[[Double]]=[]
        let order=diagonal.indices.sorted { diagonal[$0]==diagonal[$1] ? $0<$1:diagonal[$0]<diagonal[$1] }
        for i in order.prefix(min(d,roots+2)) {
            var v=[Double](repeating:0,count:d);v[i]=1;vectors.append(v);images.append(try apply(v))
        }
        if vectors.count<capacity {
            let seed=(0..<d).map { sin(Double($0+1)*1.6180339887498948)+cos(Double($0+1)*0.731) }
            if let v=orthogonal(seed,vectors) { vectors.append(v);images.append(try apply(v)) }
        }
        for iteration in 1...cfg.maximumIterations {
            let m=vectors.count;var projected=VivoQMMatrix(m,m)
            for i in 0..<m { for j in 0...i {
                let a=dot(vectors[i],images[j]),b=dot(vectors[j],images[i])
                guard abs(a-b)<1e-8*max(1,abs(a),abs(b)) else { throw VivoChemistryError.invalid("direct CI Hermiticity") }
                projected[i,j]=0.5*(a+b);projected[j,i]=projected[i,j]
            } }
            let eig=try VivoQMDenseAlgebra.symmetricEigen(projected,tolerance:1e-14,maximumSweeps:150)
            var ritz:[[Double]]=[],aritz:[[Double]]=[],residuals:[[Double]]=[],norms:[Double]=[]
            for root in 0..<roots {
                var v=[Double](repeating:0,count:d),av=v
                for j in 0..<m { for k in 0..<d { v[k]+=eig.vectors[j,root]*vectors[j][k];av[k]+=eig.vectors[j,root]*images[j][k] } }
                let residual=(0..<d).map { av[$0]-eig.values[root]*v[$0] }
                ritz.append(v);aritz.append(av);residuals.append(residual);norms.append(residual.reduce(0.0) { hypot($0,$1) })
            }
            if norms.allSatisfy({$0<=cfg.residualTolerance}) {
                var results:[VivoCIResult]=[]
                for root in 0..<roots {
                    var v=ritz[root]
                    if let k=v.indices.max(by:{abs(v[$0])<abs(v[$1])}),v[k]<0 { v=v.map { -$0 } }
                    let exact=try apply(v),error=(0..<d).map { exact[$0]-eig.values[root]*v[$0] }.reduce(0.0) { hypot($0,$1) }
                    guard error<=1.01*cfg.residualTolerance else { throw VivoChemistryError.convergence("Davidson true residual") }
                    let state=VivoCIState(orbitalCount:n,alphaElectrons:h.alphaElectrons,betaElectrons:h.betaElectrons,determinants:dets,coefficients:v)
                    try state.validate(budget:budget)
                    results.append(.init(method:.fci,energyHartree:eig.values[root]+h.constantEnergyHartree,state:state,eigenResidualNorm:error,
                        nextStateGapHartree:root+1<roots ? eig.values[root+1]-eig.values[root]:nil))
                };return .init(roots:results,iterations:iteration,matrixVectorProducts:products,operatorApplications:work)
            }
            var corrections:[[Double]]=[]
            for root in 0..<roots where norms[root]>cfg.residualTolerance {
                let candidate=(0..<d).map { k -> Double in
                    let gap=eig.values[root]-diagonal[k];return residuals[root][k]/((gap<0 ? -1.0:1.0)*max(abs(gap),1e-4))
                }
                if let q=orthogonal(candidate,vectors+corrections) { corrections.append(q) }
            }
            guard !corrections.isEmpty else { throw VivoChemistryError.convergence("Davidson stagnation") }
            if vectors.count+corrections.count>capacity { vectors=ritz;images=aritz }
            for candidate in corrections {
                if vectors.count>=capacity { break }
                if let q=orthogonal(candidate,vectors) { vectors.append(q);images.append(try apply(q)) }
            }
        };throw VivoChemistryError.convergence("Davidson iteration limit")
    }
}
