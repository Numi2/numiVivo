import Foundation

/// Adaptive sum over integer counts using cubic Hermite blocks. The error
/// bound is a derivative remainder for a discrete sum, not a Gamma surrogate.
enum VivoOmicsNBAdaptiveMoments {
    struct Node {
        let count: UInt64, probability: Double, deviance: Double
        let logSlope: Double, devianceSlope: Double
        var values: [Double] { [probability,probability*deviance,probability*deviance*deviance] }
        var slopes: [Double] { [probability*logSlope,
            probability*(logSlope*deviance+devianceSlope),
            probability*(logSlope*deviance*deviance+2*deviance*devianceSlope)] }
    }
    private struct Sum {
        var value = 0.0, correction = 0.0
        mutating func add(_ x: Double) {
            let y = x-correction, next = value+y
            correction = (next-value)-y; value = next
        }
    }
    // Stirling correction for log Gamma(z), with the next omitted term below
    // ordinary FP64 rounding for z >= 16. Small z uses the library log Gamma.
    static func stirlingCorrection(_ z: Double) -> Double {
        if z < 16 { return lgamma(z)-(z-0.5)*log(z)+z-0.5*log(2*Double.pi) }
        let t = 1/(z*z)
        return (1.0/12+t*(-1.0/360+t*(1.0/1260+t*(-1.0/1680+t*(1.0/1188+t*(-691.0/360360+t/156))))))/z
    }
    static func logProbability(count x: Double, mean mu: Double, size r: Double, deviance: Double) -> Double {
        if x == 0 { return -r*log1p(mu/r) }
        return (-deviance/2+stirlingCorrection(x+r)-stirlingCorrection(x)-stirlingCorrection(r) - 0.5*(log(2*Double.pi)+log(x)+log1p(x/r)))
    }
    private static func digamma(_ value: Double) -> Double {
        var z = value, result = 0.0
        while z < 16 { result -= 1/z; z += 1 }
        let t = 1/(z*z)
        return result+log(z)-0.5/z+t*(-1.0/12+t*(1.0/120+t*(-1.0/252+t*(1.0/240+t*(-1.0/132+t*691/32760)))))
    }
    static func logProbabilitySlope(count x: Double, mean mu: Double, size r: Double) -> Double {
        let z = x+1, delta = r-1
        if min(z,x+r) < 16 { return digamma(x+r)-digamma(z)-log1p(r/mu) }
        let logRatio = log1p(delta/z)
        func difference(_ n: Int) -> Double { pow(z,-Double(n))*expm1(-Double(n)*logRatio) }
        return (logRatio-0.5*difference(1)-difference(2)/12+difference(4)/120-difference(6)/252 + difference(8)/240-difference(10)/132+691*difference(12)/32760-log1p(r/mu))
    }
    static func polynomialSum(left a: Double, right b: Double, leftSlope da: Double, rightSlope db: Double, span n: Double) -> Double {
        (n+1)*(a+b)/2+(n*n-1)*(da-db)/12
    }
    static func fourthDerivativeBounds(left a: Node, right b: Node, size r: Double, mean mu: Double, modeProbability: Double) -> [Double] {
        let x = Double(a.count), z = x+min(1,r)
        let l1 = max(abs(a.logSlope),abs(b.logSlope))+1e-14
        let l2 = abs(r-1)*(1/(z*z)+2/(z*z*z))
        let l3 = abs(r-1)*(2/pow(z,3)+6/pow(z,4))
        let l4 = abs(r-1)*(6/pow(z,4)+24/pow(z,5))
        let p: Double
        if a.logSlope < -1e-14 && b.logSlope < -1e-14 { p = a.probability*(1+1e-13) }
        else if a.logSlope > 1e-14 && b.logSlope > 1e-14 { p = b.probability*(1+1e-13) }
        else if r <= 1 { p = max(a.probability,b.probability)*(1+1e-13) }
        else {
            // A continuous maximum is within 1/2 of an integer; the discrete
            // modal mass bounds that integer. Bound the intervening log slope.
            let near = x+0.5
            let slope = log1p(r/mu)+abs(r-1)*(1/near+1/(near*near))
            p = modeProbability*exp(0.5*slope)*(1+1e-13)
        }
        let pd = [p,p*l1,p*(l1*l1+l2),p*(l1*l1*l1+3*l1*l2+l3),
            p*(pow(l1,4)+6*l1*l1*l2+3*l2*l2+4*l1*l3+l4)]
        let d = max(a.deviance,b.deviance), d1 = max(abs(a.devianceSlope),abs(b.devianceSlope))
        let d2 = 2/(x*(1+x/r))
        let d3 = -2*expm1(-2*log1p(r/x))/(x*x)
        let d4 = -4*expm1(-3*log1p(r/x))/(x*x*x)
        let gd = [d,d1,d2,d3,d4]
        let gd2 = [d*d,2*d*d1,2*d1*d1+2*d*d2,6*d1*d2+2*d*d3,6*d2*d2+8*d1*d3+2*d*d4]
        let coefficients = [1.0,4,6,4,1]
        let first = (0...4).reduce(0.0) { $0+coefficients[$1]*pd[$1]*gd[4-$1] }
        let second = (0...4).reduce(0.0) { $0+coefficients[$1]*pd[$1]*gd2[4-$1] }
        return [pd[4],first,second]
    }
    static func evaluate(mean mu: Double, dispersion: Double, relativeTolerance: Double, maximumTerms: Int) throws -> VivoOmicsNBDevianceMoments {
        let r = 1/dispersion, q = mu/(mu+r), budget = relativeTolerance/128
        var evaluations = 0, blocks = 0
        func node(_ count: UInt64) throws -> Node {
            guard evaluations < maximumTerms else { throw VivoOmicsStatisticsError.invalid("adaptive deviance moments exhausted \(maximumTerms) endpoint evaluations") }
            evaluations += 1
            let x = Double(count), d = try VivoOmicsNegativeBinomial.unitDeviance(count: count,mean: mu,dispersion: dispersion)
            let p = exp(logProbability(count: x,mean: mu,size: r,deviance: d))
            let slope = x == 0 ? 0 : 2*log1p(((x-mu)/mu)/(1+x/r))
            guard p.isFinite, p >= 0, slope.isFinite else { throw VivoOmicsStatisticsError.invalid("adaptive moment node exceeds finite arithmetic") }
            return .init(count: count,probability: p,deviance: d,
                logSlope: logProbabilitySlope(count: x,mean: mu,size: r),devianceSlope: slope)
        }
        let mode = UInt64(floor(max(0,mu*(1-dispersion)))), modeNode = try node(mode)
        let zeroD = try VivoOmicsNegativeBinomial.unitDeviance(count: 0,mean: mu,dispersion: dispersion)
        func tail(_ v: Node, left: Bool) -> [Double] {
            let x = Double(v.count)
            if left {
                if v.count == 0 { return [0,0,0] }
                let rho = x/(q*(x+r-1))
                if rho >= 1 { return [.infinity,.infinity,.infinity] }
                let mass = v.probability*rho/(1-rho)
                return [mass,mass*zeroD,mass*zeroD*zeroD]
            }
            let rho = max(q,q*(x+r)/(x+1))
            if rho >= 1 { return [.infinity,.infinity,.infinity] }
            let t = 1/(1-rho), s0 = rho*t, s1 = s0*t, s2 = s1*(1+rho)*t
            let slope = 2*log1p(r/mu), b = slope*x+zeroD
            return [v.probability*s0,v.probability*(b*s0+slope*s1),
                v.probability*(b*b*s0+2*b*slope*s1+slope*slope*s2)]
        }
        func small(_ v: [Double]) -> Bool { v.allSatisfy { $0.isFinite && $0 <= budget } }
        var low = modeNode, left = tail(modeNode,left: true)
        while !small(left) { low = try node(low.count/2); left = tail(low,left: true) }
        var high = try node(max(mode+1,UInt64(ceil(mu+4*sqrt(mu+dispersion*mu*mu)))))
        var right = tail(high,left: false)
        while !small(right) {
            guard high.count < 4_503_599_627_370_000 else { throw VivoOmicsStatisticsError.invalid("adaptive moment tail exceeds exact count domain") }
            high = try node(high.count*2+1); right = tail(high,left: false)
        }
        var sums = [Sum(),Sum(),Sum()], errors = [Sum(),Sum(),Sum()]
        let firstCount = low.count, lastCount = high.count, covered = lastCount-firstCount+1
        func accept(_ values: [Double], _ bounds: [Double]) {
            for i in 0...2 { sums[i].add(values[i]); errors[i].add(bounds[i]) }
            blocks += 1
        }
        if low.count < 32 {
            for k in low.count...min(31,high.count) { accept(try node(k).values,[0,0,0]) }
            if high.count >= 32 { low = try node(32) }
        }
        var stack: [(Node,Node)] = high.count >= 32 ? [(low,high)] : []
        while let (a,b) = stack.popLast() {
            if a.count == b.count { accept(a.values,[0,0,0]); continue }
            if b.count-a.count == 1 { accept(zip(a.values,b.values).map(+),[0,0,0]); continue }
            let n = Double(b.count-a.count), va = a.values, vb = b.values, da = a.slopes, db = b.slopes
            let values = (0...2).map { polynomialSum(left: va[$0],right: vb[$0],leftSlope: da[$0],rightSlope: db[$0],span: n) }
            let bound = fourthDerivativeBounds(left: a,right: b,size: r,mean: mu,modeProbability: modeNode.probability)
                .map { $0*(pow(n,5)-n)/720 }
            let allowance = budget*(n+1)/Double(covered)
            if zip(values,bound).allSatisfy({ $0.isFinite && $0 >= 0 && $1.isFinite && $1 <= allowance }) { accept(values,bound) }
            else {
                let middle = a.count+(b.count-a.count)/2
                let m = try node(middle), next = try node(middle+1)
                stack.append((next,b)); stack.append((a,m))
            }
        }
        let mass = sums[0].value, m = sums[1].value/mass, rawSecond = sums[2].value/mass
        let variance = rawSecond-m*m, lowerMass = mass-errors[0].value
        let total = (0...2).map { errors[$0].value+left[$0]+right[$0] }
        let meanSumError = (errors[1].value+abs(m)*errors[0].value)/lowerMass
        let secondSumError = (errors[2].value+abs(rawSecond)*errors[0].value)/lowerMass
        let varianceSumError = secondSumError+(2*abs(m)+meanSumError)*meanSumError
        let meanError = (total[1]+abs(m)*total[0])/lowerMass
        let secondError = (total[2]+abs(rawSecond)*total[0])/lowerMass
        let varianceError = secondError+(2*abs(m)+meanError)*meanError
        guard mass.isFinite, lowerMass > 0, m.isFinite, m > 0, variance.isFinite, variance > 0,
              total[0]/lowerMass <= relativeTolerance, meanError <= relativeTolerance*m,
              varianceError <= relativeTolerance*variance else {
            throw VivoOmicsStatisticsError.invalid("adaptive deviance moments could not certify combined error")
        }
        let scale = 2*m/variance, df = scale*m
        guard scale.isFinite, df.isFinite else { throw VivoOmicsStatisticsError.invalid("adaptive moment scale exceeds finite arithmetic") }
        var result = VivoOmicsNBDevianceMoments(mean: m,variance: variance,devianceScale: scale,degreesOfFreedom: df,
            firstCount: firstCount,lastCount: lastCount,evaluatedCounts: evaluations,
            omittedProbabilityBound: (left[0]+right[0])/lowerMass,
            meanTruncationBound: max(0,meanError-meanSumError),varianceTruncationBound: max(0,varianceError-varianceSumError))
        result.summation = .init(method: .adaptive,coveredCounts: covered,acceptedBlocks: blocks,
            probabilityErrorBound: errors[0].value,firstRawMomentErrorBound: errors[1].value,
            secondRawMomentErrorBound: errors[2].value,meanErrorBound: meanSumError,varianceErrorBound: varianceSumError)
        return result
    }
}
