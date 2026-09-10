import Foundation

public struct VivoOmicsFLogTails: Codable, Sendable, Equatable {
    public let lower: Double
    public let upper: Double
}

/// Positive-shape special functions used by unequal-DF QL moderation.
/// Bounds and convergence failures are explicit; log tails retain tiny probabilities.
public enum VivoOmicsQLSpecialFunctions {
    static func validateShape(_ x: Double) throws {
        guard x.isFinite, (1e-8...1e6).contains(x) else {
            throw VivoOmicsStatisticsError.invalid("QL special-function shape outside [1e-8,1e6]")
        }
    }
    public static func logMinusDigamma(_ x: Double) throws -> Double {
        try validateShape(x)
        var z = x, sum = 0.0
        while z < 16 { sum += 1/z-log1p(1/z); z += 1 }
        let inverse = 1/z, squared = inverse*inverse
        let coefficients: [Double] = [1.0/12,-1.0/120,1.0/252,-1.0/240,1.0/132,-691.0/32760,1.0/12,-3617.0/8160]
        var power = squared
        sum += inverse/2
        for coefficient in coefficients { sum += coefficient*power; power *= squared }
        return sum
    }
    public static func trigamma(_ x: Double) throws -> Double {
        try validateShape(x)
        var z = x, sum = 0.0
        while z < 16 { sum += 1/(z*z); z += 1 }
        let inverse = 1/z, squared = inverse*inverse
        let coefficients: [Double] = [1.0/6,-1.0/30,1.0/42,-1.0/30,5.0/66,-691.0/2730,7.0/6,-3617.0/510]
        sum += inverse+squared/2
        var power = inverse*squared
        for coefficient in coefficients { sum += coefficient*power; power *= squared }
        return sum
    }
    static func stirlingCorrection(_ x: Double) -> Double {
        let inverse = 1/x, squared = inverse*inverse
        let coefficients: [Double] = [1.0/12,-1.0/360,1.0/1260,-1.0/1680,1.0/1188,-691.0/360360,1.0/156,-3617.0/122400]
        var power = inverse, sum = 0.0
        for coefficient in coefficients { sum += coefficient*power; power *= squared }
        return sum
    }
    public static func logGammaIncrement(base: Double, increment: Double) throws -> Double {
        try validateShape(base); try validateShape(increment)
        var x = base, shift = 0.0
        while x < 16 { shift -= log1p(increment/x); x += 1 }
        let main = shift+(x+increment-0.5)*log1p(increment/x)+increment*(log(x)-1)
        return main+stirlingCorrection(x+increment)-stirlingCorrection(x)
    }
    static func logBeta(_ a: Double, _ b: Double) throws -> Double {
        let small = min(a,b), large = max(a,b)
        return lgamma(small)-(try logGammaIncrement(base: large,increment: small))
    }
    static func betaFraction(_ a: Double, _ b: Double, _ x: Double, maximumIterations: Int) throws -> Double {
        func nonzero(_ value: Double) -> Double {
            abs(value) < 1e-300 ? (value < 0 ? -1e-300 : 1e-300) : value
        }
        var c = 1.0, d = 1/nonzero(1-(a+b)*x/(a+1)), h = d
        for k in 1...maximumIterations {
            let m = Double(k), twice = 2*m
            let first = m*(b-m)*x/((a-1+twice)*(a+twice))
            d = 1/nonzero(1+first*d); c = nonzero(1+first/c); h *= d*c
            let second = -(a+m)*(a+b+m)*x/((a+twice)*(a+1+twice))
            d = 1/nonzero(1+second*d); c = nonzero(1+second/c)
            let change = d*c; h *= change
            guard h.isFinite, h > 0 else { throw VivoOmicsStatisticsError.invalid("QL beta fraction lost positive finite arithmetic") }
            if abs(change-1) <= 2e-14 { return h }
        }
        throw VivoOmicsStatisticsError.invalid("QL beta fraction exhausted \(maximumIterations) iterations")
    }
    static func softplus(_ x: Double) -> Double { max(0,x)+log1p(exp(-abs(x))) }
    static func logComplement(_ logProbability: Double) -> Double {
        logProbability > -log(2) ? log(-expm1(logProbability)) : log1p(-exp(logProbability))
    }
    public static func logFTails(logStatistic: Double, numeratorDF: Double, denominatorDF: Double,
        maximumIterations: Int = 2048) throws -> VivoOmicsFLogTails {
        let a = numeratorDF/2, b = denominatorDF/2
        try validateShape(a); try validateShape(b)
        guard logStatistic.isFinite, (-1000...1000).contains(logStatistic), (1...8192).contains(maximumIterations) else {
            throw VivoOmicsStatisticsError.invalid("QL F-tail statistic or iteration bound")
        }
        let odds = log(numeratorDF)+logStatistic-log(denominatorDF)
        let logX = -softplus(-odds), logY = -softplus(odds)
        let lowerDirect = odds < log(a+1)-log(b+1)
        let first = lowerDirect ? a : b, second = lowerDirect ? b : a
        let coordinate = lowerDirect ? logX : logY
        let fraction = try betaFraction(first,second,exp(coordinate),maximumIterations: maximumIterations)
        let direct = a*logX+b*logY-(try logBeta(a,b))+log(fraction)-log(first)
        guard direct.isFinite, direct <= 1e-12 else {
            throw VivoOmicsStatisticsError.invalid("QL log beta probability outside its finite domain")
        }
        let bounded = min(0,direct), complement = logComplement(bounded)
        guard complement.isFinite else { throw VivoOmicsStatisticsError.invalid("QL complementary log tail lost finite precision") }
        return lowerDirect ? .init(lower: bounded,upper: complement) : .init(lower: complement,upper: bounded)
    }
}
