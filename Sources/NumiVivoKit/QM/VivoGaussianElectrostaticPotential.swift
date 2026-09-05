import Foundation

private struct VivoPotentialHermiteKey: Hashable { let i:Int,j:Int,t:Int }
private struct VivoPotentialCoulombKey: Hashable { let t:Int,u:Int,v:Int,n:Int }

public enum VivoGaussianElectrostaticPotential {
    private static func boys(_ n:Int,_ t:Double)->Double {
        if t < 60 {
            var term=1/Double(2*n+1),sum=term
            for k in 1...512 { term *= 2*t/Double(2*n+2*k+1);sum += term;if abs(term)<abs(sum)*2e-16 { break } }
            return exp(-t)*sum
        }
        let exponential=exp(-t);var f=0.5*sqrt(Double.pi/t)*erf(sqrt(t))
        if n>0 { for m in 0..<n { f=(Double(2*m+1)*f-exponential)/(2*t) } };return f
    }
    private static func hermite(_ i:Int,_ j:Int,_ t:Int,_ q:Double,_ a:Double,_ b:Double,
                                _ cache:inout [VivoPotentialHermiteKey:Double])->Double {
        if i<0||j<0||t<0||t>i+j{return 0};let key=VivoPotentialHermiteKey(i:i,j:j,t:t);if let v=cache[key]{return v};let p=a+b;let value:Double
        if i==0&&j==0 { value=exp(-a*b/p*q*q) }
        else if j==0 { value=hermite(i-1,j,t-1,q,a,b,&cache)/(2*p)-b*q/p*hermite(i-1,j,t,q,a,b,&cache)+Double(t+1)*hermite(i-1,j,t+1,q,a,b,&cache) }
        else { value=hermite(i,j-1,t-1,q,a,b,&cache)/(2*p)+a*q/p*hermite(i,j-1,t,q,a,b,&cache)+Double(t+1)*hermite(i,j-1,t+1,q,a,b,&cache) }
        cache[key]=value;return value
    }
    private static func coefficients(_ i:Int,_ j:Int,_ q:Double,_ a:Double,_ b:Double)->[Double] {
        var cache:[VivoPotentialHermiteKey:Double]=[:];return (0...(i+j)).map{hermite(i,j,$0,q,a,b,&cache)}
    }
    private static func coulomb(_ t:Int,_ u:Int,_ v:Int,_ n:Int,_ p:Double,_ pc:SIMD3<Double>,
                                _ cache:inout [VivoPotentialCoulombKey:Double])->Double {
        if t<0||u<0||v<0{return 0};let key=VivoPotentialCoulombKey(t:t,u:u,v:v,n:n);if let x=cache[key]{return x};let value:Double
        if t==0&&u==0&&v==0 { value=pow(-2*p,Double(n))*boys(n,p*vivoQMDot(pc,pc)) }
        else if t>0 { value=(t>1 ? Double(t-1)*coulomb(t-2,u,v,n+1,p,pc,&cache):0)+pc.x*coulomb(t-1,u,v,n+1,p,pc,&cache) }
        else if u>0 { value=(u>1 ? Double(u-1)*coulomb(t,u-2,v,n+1,p,pc,&cache):0)+pc.y*coulomb(t,u-1,v,n+1,p,pc,&cache) }
        else { value=(v>1 ? Double(v-1)*coulomb(t,u,v-2,n+1,p,pc,&cache):0)+pc.z*coulomb(t,u,v-1,n+1,p,pc,&cache) }
        cache[key]=value;return value
    }
    private static func primitive(_ a:Double,_ la:[Int],_ ra:SIMD3<Double>,_ b:Double,_ lb:[Int],_ rb:SIMD3<Double>,
                                  _ point:SIMD3<Double>)->Double {
        let p=a+b,pc=(a*ra+b*rb)/p-point,e=(0..<3).map{coefficients(la[$0],lb[$0],ra[$0]-rb[$0],a,b)}
        var cache:[VivoPotentialCoulombKey:Double]=[:],value=0.0
        for t in e[0].indices { for u in e[1].indices { for v in e[2].indices {
            value += e[0][t]*e[1][u]*e[2][v]*coulomb(t,u,v,0,p,pc,&cache)
        } } }
        return 2*Double.pi/p*value
    }
    public static func matrix(integrals ao:VivoAOIntegrals,at point:SIMD3<Double>,
                              budget:VivoChemistryBudget = .init())throws->VivoQMMatrix {
        try ao.validate(budget:budget);guard vivoQMFinite(point) else{throw VivoChemistryError.invalid("electrostatic-potential point")}
        let n=ao.count;_ = try budget.elements([n,n],simultaneousArrays:4);var result=VivoQMMatrix(n,n)
        for mu in 0..<n { for nu in 0...mu {
            let a=ao.orbitals[mu],b=ao.orbitals[nu],ra=ao.sourceSystem.nuclei[a.nucleusIndex].positionBohr,rb=ao.sourceSystem.nuclei[b.nucleusIndex].positionBohr;var value=0.0
            for i in a.weights.indices { for j in b.weights.indices {
                value += a.weights[i]*b.weights[j]*primitive(a.primitiveExponents[i],a.angular,ra,b.primitiveExponents[j],b.angular,rb,point)
            } }
            result[mu,nu]=value;result[nu,mu]=value
        } };return result
    }
    public static func solutePotential(system:VivoElectronicSystem,integrals ao:VivoAOIntegrals,totalDensity:VivoQMMatrix,
                                       at point:SIMD3<Double>,budget:VivoChemistryBudget = .init())throws->Double {
        guard totalDensity.rows==ao.count,totalDensity.columns==ao.count,totalDensity.values.allSatisfy(\.isFinite),system==ao.sourceSystem else {
            throw VivoChemistryError.invalid("electrostatic potential density/source")
        }
        let matrix=try self.matrix(integrals:ao,at:point,budget:budget);var value=0.0
        for nucleus in system.nuclei { let r=vivoQMNorm(point-nucleus.positionBohr);guard r>1e-10 else{throw VivoChemistryError.invalid("potential point on nucleus")};value += Double(nucleus.atomicNumber)/r }
        for charge in system.pointCharges { let r=vivoQMNorm(point-charge.positionBohr);guard r>1e-10 else{throw VivoChemistryError.invalid("potential point on external charge")};value += charge.chargeE/r }
        value -= zip(totalDensity.values,matrix.values).reduce(0){$0+$1.0*$1.1}
        guard value.isFinite else{throw VivoChemistryError.convergence("nonfinite solute potential")};return value
    }
}
