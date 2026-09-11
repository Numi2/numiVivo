import Foundation

/// A conditional, single-feature NB2 observation posterior. Rates are CPM and
/// exposures are full-RNA UMIs / 1e6. Cells, not donor aggregate rows, are the
/// independent likelihood observations. Dispersion is fixed, not inferred here.
public struct VivoCountObservationPosterior: Codable, Sendable, Equatable {
    public let method: String
    public let observedCells: Int
    public let geneCounts: UInt64
    public let libraryCounts: UInt64
    public let cellDispersion: Double
    public let gammaPriorShape: Double
    public let gammaPriorRatePerCPM: Double
    public let nominalCoverage: Double
    public let meanCPM: Double
    public let varianceCPM: Double
    public let meanLog1pCPM: Double
    public let varianceLog1pCPM: Double
    public let lowerLog1pCPM: Double
    public let upperLog1pCPM: Double
    public let quadraturePanels: Int
    public let maximumRefinementDifference: Double
    /// Conservative analytic tail bounds, normalized by retained mass, for
    /// mass and first/second moments of rate divided by the log-density mode.
    public let maximumScaledMomentTailBound: Double
}

public struct VivoCountObservationPredictiveMoments: Codable, Sendable, Equatable {
    public let plannedCells: Int
    public let plannedLibraryCounts: UInt64
    public let meanGeneCounts: Double
    public let conditionalPoissonVariance: Double
    public let conditionalCellOverdispersionVariance: Double
    public let latentRateVariance: Double
    public let totalGeneCountVariance: Double
}

public enum VivoCountObservation {
    /// Independent NB2 cells: E[Y_i|r] = L_i*r/1e6,
    /// Var[Y_i|r] = E[Y_i|r] + dispersion*E[Y_i|r]^2.
    /// r has an explicitly supplied proper Gamma(shape, rate) prior in CPM.
    /// The result is conditional on this prior, dispersion and independence;
    /// it is not calibrated biological uncertainty or a perturbation model.
    public static func posterior(counts: [UInt64], libraryCounts: [UInt64],
        cellDispersion: Double, gammaPriorShape: Double, gammaPriorRatePerCPM: Double,
        coverage: Double = 0.95) throws -> VivoCountObservationPosterior {
        guard !counts.isEmpty, counts.count == libraryCounts.count, counts.count <= 1_000_000,
              cellDispersion.isFinite, cellDispersion == 0 || (1e-8...100).contains(cellDispersion),
              gammaPriorShape.isFinite, (0.1...1e6).contains(gammaPriorShape),
              gammaPriorRatePerCPM.isFinite, (1e-12...1e12).contains(gammaPriorRatePerCPM),
              coverage.isFinite, (0.5...0.99).contains(coverage) else {
            throw VivoOmicsError.invalid("count observation dimensions, dispersion, proper Gamma prior or coverage")
        }
        var total: UInt64 = 0, libraries: UInt64 = 0
        for i in counts.indices {
            guard libraryCounts[i] > 0, libraryCounts[i] <= 1_000_000_000,
                  counts[i] <= libraryCounts[i] else {
                throw VivoOmicsError.invalid("count observation requires exact per-cell gene counts and full positive RNA libraries")
            }
            total = try vivoOmicsSum(total,counts[i]); libraries = try vivoOmicsSum(libraries,libraryCounts[i])
        }
        guard libraries <= 1_000_000_000_000 else { throw VivoOmicsError.limit("count observation total library budget") }
        let a = gammaPriorShape, b = gammaPriorRatePerCPM, phi = cellDispersion
        let exposures = libraryCounts.map { Double($0)/1e6 }
        let sumY = Double(total), shape = sumY+a
        let exposure = Double(libraries)/1e6
        let logOffsets = phi == 0 ? [] : exposures.map { log(phi*$0) }
        let weights = phi == 0 ? [] : counts.map { Double($0)+1/phi }
        func score(_ x: Double) -> Double {
            if phi == 0 { return shape-(b+exposure)*exp(x) }
            var value = a-b*exp(x)
            // y - (y+1/phi)*q, written to retain accuracy near q == 1.
            for i in counts.indices {
                let q = sigmoid(logOffsets[i]+x)
                value += Double(counts[i])*(1-q)-q/phi
            }
            return value
        }
        var lo = log(shape/(b+exposure)), hi = lo
        for _ in 0..<128 {
            if score(lo) >= 0 { break }; lo -= 2
        }
        for _ in 0..<128 {
            if score(hi) <= 0 { break }; hi += 2
        }
        guard score(lo) >= 0, score(hi) <= 0 else { throw VivoOmicsError.invalid("count posterior mode bracket") }
        for _ in 0..<100 {
            let mid = (lo+hi)/2
            if score(mid) > 0 { lo=mid } else { hi=mid }
        }
        let mode = (lo+hi)/2, modeRate = exp(mode), modeScore = score(mode)
        let q = logOffsets.map { sigmoid($0+mode) }
        let priorExposure = (phi == 0 ? b+exposure : b)*modeRate
        var curvature = priorExposure
        if phi > 0 { for i in counts.indices { curvature += weights[i]*q[i]*(1-q[i]) } }
        guard curvature.isFinite, curvature > 0 else { throw VivoOmicsError.invalid("count posterior curvature") }
        func centered(_ delta: Double) -> Double {
            var value = modeScore*delta-priorExposure*expm1Remainder(delta)
            if phi > 0 {
                for i in counts.indices {
                    value -= weights[i]*softplusRemainder(logOffsets[i]+mode,q[i],delta)
                }
            }
            return value
        }
        var left = -1/sqrt(curvature), right = -left
        func tail(_ d: Double, rightSide: Bool) -> Double {
            let slope = score(mode+d)
            if rightSide {
                guard slope < -2 else { return .infinity }
                return exp(centered(d)+2*max(0,d))/(-slope-2)
            }
            guard slope > 0 else { return .infinity }
            return exp(centered(d))/slope
        }
        let scale = 1/sqrt(curvature)
        for _ in 0..<256 {
            if tail(left,rightSide: false) < 1e-15*scale { break }; left *= 1.4
        }
        for _ in 0..<256 {
            if tail(right,rightSide: true) < 1e-15*scale { break }; right *= 1.4
        }
        guard left.isFinite,right.isFinite,mode+left > -700,mode+right < 700,
              tail(left,rightSide: false) < 1e-15*scale,
              tail(right,rightSide: true) < 1e-15*scale else {
            throw VivoOmicsError.invalid("count posterior tail bounds")
        }
        // Integrate centered moments, avoiding cancellation for narrow posteriors.
        let logCenter = softplus(mode), target = (1-coverage)/2
        func values(_ d: Double) -> [Double] {
            let density = exp(centered(d)), rateDelta = expm1(d)
            let logDelta = softplusDifference(mode,d)
            return [density,density*rateDelta,density*rateDelta*rateDelta,
                    density*logDelta,density*logDelta*logDelta]
        }
        func integrate(_ low: Double,_ high: Double) -> [Double] {
            let mid=(low+high)/2, half=(high-low)/2
            var sums=[Double](repeating: 0,count: 5)
            for k in nodes.indices {
                for sign in [-1.0,1.0] {
                    let v=values(mid+sign*half*nodes[k])
                    for j in sums.indices { sums[j] += half*weights16[k]*v[j] }
                }
            }
            return sums
        }
        var previous: [Double]?, panels=8, refinement=Double.infinity
        while panels <= 2048 {
            try Task.checkCancellation()
            let h=(right-left)/Double(panels)
            var sums=[Double](repeating: 0,count: 5), masses: [Double]=[]
            for i in 0..<panels {
                let v=integrate(left+Double(i)*h,left+Double(i+1)*h)
                masses.append(v[0]);for j in sums.indices { sums[j] += v[j] }
            }
            let z=sums[0]
            guard z.isFinite,z>0,sums.allSatisfy(\.isFinite) else { throw VivoOmicsError.invalid("count posterior quadrature arithmetic") }
            let rateDelta=sums[1]/z, logDelta=sums[3]/z
            let varianceRate=sums[2]/z-rateDelta*rateDelta
            let varianceLog=sums[4]/z-logDelta*logDelta
            func quantile(_ probability: Double) -> Double {
                let wanted=z*probability
                var cumulative=0.0
                for i in 0..<panels {
                    if cumulative+masses[i] >= wanted {
                        let start=left+Double(i)*h
                        var low=start, high=start+h
                        for _ in 0..<45 {
                            let mid=(low+high)/2
                            if cumulative+integrate(start,mid)[0] < wanted { low=mid } else { high=mid }
                        }
                        return softplus(mode+(low+high)/2)
                    }
                    cumulative += masses[i]
                }
                return softplus(mode+right)
            }
            let lower=quantile(target),upper=quantile(1-target)
            // Compare actual means, not centered corrections that can be
            // analytically zero (e.g. Gamma-Poisson rate mean at this mode).
            let current=[z,1+rateDelta,varianceRate,logCenter+logDelta,varianceLog,lower,upper]
            if let previous {
                refinement=zip(current,previous).map { abs($0-$1)/max(1e-10,abs($0),abs($1)) }.max()!
                if refinement < 1e-8 {
                    guard varianceRate>0,varianceLog>0,lower<upper else { throw VivoOmicsError.invalid("count posterior nonpositive uncertainty") }
                    let tailBound=(tail(left,rightSide: false)+tail(right,rightSide: true))/z
                    return .init(method: "independent-cell-nb2-gamma-rate-posterior-quadrature-v1",
                        observedCells: counts.count,geneCounts: total,libraryCounts: libraries,cellDispersion: phi,
                        gammaPriorShape: a,gammaPriorRatePerCPM: b,nominalCoverage: coverage,
                        meanCPM: modeRate*(1+rateDelta),varianceCPM: modeRate*modeRate*varianceRate,
                        meanLog1pCPM: logCenter+logDelta,varianceLog1pCPM: varianceLog,
                        lowerLog1pCPM: lower,upperLog1pCPM: upper,quadraturePanels: panels,
                        maximumRefinementDifference: refinement,maximumScaledMomentTailBound: tailBound)
                }
            }
            previous=current; panels *= 2
        }
        throw VivoOmicsError.invalid("count posterior quadrature did not converge")
    }

    /// Predict a new independent sample of the SAME latent condition. Planned
    /// per-cell RNA depths are required; no treated outcomes or inferred N.
    /// These are count moments, not an interval, perturbation response or
    /// distribution of log1p-normalized counts.
    public static func predictiveMoments(_ posterior: VivoCountObservationPosterior,
        plannedLibraryCounts: [UInt64]) throws -> VivoCountObservationPredictiveMoments {
        guard !plannedLibraryCounts.isEmpty,plannedLibraryCounts.count<=1_000_000,
              plannedLibraryCounts.allSatisfy({ $0>0 && $0<=1_000_000_000 }),
              posterior.meanCPM.isFinite,posterior.meanCPM>0,
              posterior.varianceCPM.isFinite,posterior.varianceCPM>0,
              posterior.cellDispersion.isFinite,posterior.cellDispersion==0 || (1e-8...100).contains(posterior.cellDispersion) else {
            throw VivoOmicsError.invalid("count predictive moments or planned per-cell libraries")
        }
        var total: UInt64=0,squared=0.0
        for library in plannedLibraryCounts { total=try vivoOmicsSum(total,library);let e=Double(library)/1e6;squared += e*e }
        guard total<=1_000_000_000_000 else { throw VivoOmicsError.limit("planned count observation total library budget") }
        let e=Double(total)/1e6,mean=e*posterior.meanCPM
        let over=posterior.cellDispersion*squared*(posterior.varianceCPM+posterior.meanCPM*posterior.meanCPM)
        let latent=e*e*posterior.varianceCPM,totalVariance=mean+over+latent
        guard totalVariance.isFinite else { throw VivoOmicsError.invalid("count predictive variance overflow") }
        return .init(plannedCells: plannedLibraryCounts.count,plannedLibraryCounts: total,meanGeneCounts: mean,
            conditionalPoissonVariance: mean,conditionalCellOverdispersionVariance: over,
            latentRateVariance: latent,totalGeneCountVariance: totalVariance)
    }
    private static func sigmoid(_ x: Double) -> Double { x>=0 ? 1/(1+exp(-x)) : exp(x)/(1+exp(x)) }
    private static func softplus(_ x: Double) -> Double { max(0,x)+log1p(exp(-abs(x))) }
    private static func softplusDifference(_ x: Double,_ d: Double) -> Double {
        abs(d)<0.5 ? log1p(sigmoid(x)*expm1(d)) : softplus(x+d)-softplus(x)
    }
    private static func expm1Remainder(_ d: Double) -> Double {
        if abs(d)>=0.01 { return expm1(d)-d }
        return d*d*(0.5+d*(1/6.0+d*(1/24.0+d*(1/120.0+d*(1/720.0+d/5040)))))
    }
    private static func softplusRemainder(_ x: Double,_ q: Double,_ d: Double) -> Double {
        if abs(d)>=0.001 { return softplusDifference(x,d)-q*d }
        let q2=q*q,q3=q2*q,q4=q3*q
        return q*(1-q)*d*d*(0.5+d*((1-2*q)/6+d*((1-6*q+6*q2)/24+d*((1-14*q+36*q2-24*q3)/120+d*(1-30*q+150*q2-240*q3+120*q4)/720))))
    }
    private static let nodes=[0.09501250983763744,0.2816035507792589,0.4580167776572274,0.6178762444026438,0.755404408355003,0.8656312023878318,0.9445750230732326,0.9894009349916499]
    private static let weights16=[0.1894506104550685,0.1826034150449236,0.16915651939500254,0.14959598881657673,0.12462897125553387,0.09515851168249278,0.062253523938647894,0.027152459411754096]
}
