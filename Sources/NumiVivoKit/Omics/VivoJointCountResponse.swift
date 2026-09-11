import Foundation

public struct VivoJointCountPair: Codable, Sendable, Equatable {
    public let donorID: String
    public let control: VivoCountDepthStratum
    public let treated: VivoCountDepthStratum
    public init(donorID: String,control: VivoCountDepthStratum,treated: VivoCountDepthStratum) {
        self.donorID=donorID;self.control=control;self.treated=treated
    }
}
public struct VivoJointCountPlan: Codable, Sendable, Equatable {
    public let gridPointsPerAxis: Int
    public let maximumIterations: Int
    public let meanLogLikelihoodGapTolerance: Double
    public init(gridPointsPerAxis: Int=33,maximumIterations: Int=20_000,meanLogLikelihoodGapTolerance: Double=1e-8) {
        self.gridPointsPerAxis=gridPointsPerAxis;self.maximumIterations=maximumIterations;self.meanLogLikelihoodGapTolerance=meanLogLikelihoodGapTolerance
    }
}
public struct VivoJointRateMoments: Codable, Sendable, Equatable {
    public let controlMeanCPM: Double
    public let treatedMeanCPM: Double
    public let controlVarianceCPM2: Double
    public let treatedVarianceCPM2: Double
    public let covarianceCPM2: Double
    public let responseMeanCPM: Double
    public let responseVarianceCPM2: Double
    public let controlMeanLog1pCPM: Double
    public let treatedMeanLog1pCPM: Double
    public let log1pResponseMean: Double
    public let log1pResponseVariance: Double
}
public struct VivoJointCountModel: Codable, Sendable, Equatable {
    public let method: String
    public let featureID: String
    public let trainingSource: VivoFingerprint
    public let trainingDonorIDs: [String]
    public let controlConditionID: String
    public let treatedConditionID: String
    public let controlCellDispersion: Double
    public let treatedCellDispersion: Double
    public let controlCells: [Int]
    public let treatedCells: [Int]
    public let controlDonorMLERatesCPM: [Double]
    public let treatedDonorMLERatesCPM: [Double]
    public let controlRatesCPM: [Double]
    public let treatedRatesCPM: [Double]
    /// Row-major control by treated support; includes explicit zero weights.
    public let probabilities: [Double]
    public let fittedDonorRelativeLikelihoods: [Double]
    public let relativeLogLikelihood: Double
    public let meanLogLikelihoodGap: Double
    public let iterations: Int
    public let status: String
    public let plan: VivoJointCountPlan
    public let moments: VivoJointRateMoments
    public let qualification: String
}
public struct VivoJointCountPrediction: Codable, Sendable, Equatable {
    public let querySource: VivoFingerprint
    public let queryDonorID: String
    public let featureID: String
    public let treatedConditionID: String
    public let probabilities: [Double]
    public let moments: VivoJointRateMoments
    public let plannedTreatedCountMoments: VivoCountObservationPredictiveMoments
    public let degenerateTreatedRateDistribution: Bool
    public let underflowedPosteriorComponents: Int
    public let numericallyDegenerateTreatedRateMoments: Bool
    public let qualification: String
}

public enum VivoJointCountResponse {
    /// Nonparametric ML on declared finite support, conditional on fixed cell
    /// dispersion. This estimates a joint rate distribution, not a covariance
    /// obtained by subtracting independent marginal noise estimates.
    public static func fit(pairs: [VivoJointCountPair],featureID: String,controlConditionID: String,
        treatedConditionID: String,controlCellDispersion: Double,treatedCellDispersion: Double,
        trainingSource: VivoFingerprint,plan: VivoJointCountPlan) throws -> VivoJointCountModel {
        guard (3...64).contains(pairs.count),Set(pairs.map(\.donorID)).count==pairs.count,
              pairs.allSatisfy({ vivoOmicsID($0.donorID) }),vivoOmicsID(featureID),vivoOmicsID(controlConditionID),
              vivoOmicsID(treatedConditionID),controlConditionID != treatedConditionID,
              (3...129).contains(plan.gridPointsPerAxis),(1...200_000).contains(plan.maximumIterations),
              (1e-12...1e-4).contains(plan.meanLogLikelihoodGapTolerance) else { throw VivoOmicsError.invalid("joint count identities or numerical plan") }
        let ordered=pairs.sorted { $0.donorID<$1.donorID }
        let control=try ordered.map { try VivoCountRateLikelihood($0.control,cellDispersion: controlCellDispersion) }
        let treated=try ordered.map { try VivoCountRateLikelihood($0.treated,cellDispersion: treatedCellDispersion) }
        let cm=control.map(\.maximumLikelihoodRateCPM),tm=treated.map(\.maximumLikelihoodRateCPM)
        func grid(_ rates: [Double]) -> [Double] {
            let lo=rates.min()!,hi=rates.max()!
            if lo==hi { return [lo] }
            let a=log1p(lo),span=log1p(hi)-a
            var result=Set(rates)
            for i in 1..<plan.gridPointsPerAxis-1 { result.insert(expm1(a+span*Double(i)/Double(plan.gridPointsPerAxis-1))) }
            return result.sorted()
        }
        return try fitRateGrid(pairs: ordered,featureID: featureID,controlConditionID: controlConditionID,
            treatedConditionID: treatedConditionID,controlCellDispersion: controlCellDispersion,
            treatedCellDispersion: treatedCellDispersion,trainingSource: trainingSource,plan: plan,
            control: control,treated: treated,controlRates: grid(cm),treatedRates: grid(tm),initial: nil)
    }

    static func fitRateGrid(pairs ordered: [VivoJointCountPair],featureID: String,controlConditionID: String,
        treatedConditionID: String,controlCellDispersion: Double,treatedCellDispersion: Double,
        trainingSource: VivoFingerprint,plan: VivoJointCountPlan,control: [VivoCountRateLikelihood],
        treated: [VivoCountRateLikelihood],controlRates cr: [Double],treatedRates tr: [Double],
        initial: VivoJointCountModel?,preferLargestExchangeGain: Bool=false) throws -> VivoJointCountModel {
        let cm=control.map(\.maximumLikelihoodRateCPM),tm=treated.map(\.maximumLikelihoodRateCPM)
        let k=cr.count*tr.count,n=ordered.count
        let cl=try control.map { l in try cr.map { try l.logRelativeLikelihood(rateCPM: $0) } }
        let tl=try treated.map { l in try tr.map { try l.logRelativeLikelihood(rateCPM: $0) } }
        var a=Array(repeating: Array(repeating: 0.0,count: k),count: n),w=Array(repeating: 0.0,count: k)
        for d in 0..<n {
            var best=0
            for i in cr.indices { for j in tr.indices {
                let index=i*tr.count+j;a[d][index]=exp(cl[d][i]+tl[d][j])
                if a[d][index]>a[d][best] { best=index }
            } }
            guard a[d][best]>0 else { throw VivoOmicsError.invalid("joint count support misses donor likelihood") }
            w[best]+=1/Double(n)
        }
        if let initial {
            w=Array(repeating: 0,count: k)
            for i in initial.controlRatesCPM.indices { for j in initial.treatedRatesCPM.indices {
                guard let x=cr.firstIndex(of: initial.controlRatesCPM[i]),let y=tr.firstIndex(of: initial.treatedRatesCPM[j]) else {
                    throw VivoOmicsError.invalid("adaptive count support lost prior coordinates")
                }
                w[x*tr.count+y]=initial.probabilities[i*initial.treatedRatesCPM.count+j]
            } }
        }
        var s=Array(repeating: 0.0,count: n),gradient=Array(repeating: 0.0,count: k),iteration=0,gap=Double.infinity,status="iterationLimit"
        while true {
            try Task.checkCancellation()
            let mass=w.reduce(0,+);guard mass.isFinite,mass>0 else { throw VivoOmicsError.invalid("joint count probability mass") }
            for j in w.indices { w[j]/=mass }
            s=Array(repeating: 0,count: n)
            for j in w.indices where w[j]>0 { for d in 0..<n { s[d]+=w[j]*a[d][j] } }
            guard s.allSatisfy({ $0.isFinite && $0>0 }) else { throw VivoOmicsError.invalid("joint count zero donor likelihood") }
            gradient=Array(repeating: 0,count: k)
            for d in 0..<n { for j in w.indices { gradient[j]+=a[d][j]/s[d] } }
            let best=gradient.indices.max { gradient[$0]<gradient[$1] }!
            gap=max(0,gradient[best]/Double(n)-1)
            if gap<=plan.meanLogLikelihoodGapTolerance { status="convergedFiniteGridLikelihood";break }
            if iteration>=plan.maximumIterations { break }
            func exchange(_ worst: Int) -> (step: Double,gain: Double) {
                let delta=(0..<n).map { a[$0][best]-a[$0][worst] },maximum=w[worst]
                func derivative(_ step: Double) -> Double {
                    var value=0.0
                    for d in 0..<n {
                        let den=s[d]+step*delta[d]
                        if den<=0 { return -.infinity }
                        value+=delta[d]/den
                    }
                    return value
                }
                var step=maximum
                if derivative(maximum)<0 {
                    var lo=0.0,hi=maximum
                    for _ in 0..<70 { let mid=(lo+hi)/2;if derivative(mid)>0 { lo=mid } else { hi=mid } }
                    step=(lo+hi)/2
                }
                var gain=0.0
                for d in 0..<n { gain+=log1p(step*delta[d]/s[d]) }
                return (step,gain)
            }
            var worst=w.indices.filter { w[$0]>0 }.min { gradient[$0]<gradient[$1] }!
            var move=exchange(worst)
            if preferLargestExchangeGain {
                for candidate in w.indices where w[candidate]>0 && candidate != best {
                    let option=exchange(candidate)
                    if option.gain>move.gain { worst=candidate;move=option }
                }
            }
            let step=move.step
            guard best != worst,step>0,w[best]+step>w[best] else { status="stalledBeforeCertificate";break }
            w[worst]-=step;w[best]+=step;iteration+=1
        }
        return .init(method: "paired-cell-nb2-finite-support-npmle-v1",featureID: featureID,trainingSource: trainingSource,
            trainingDonorIDs: ordered.map(\.donorID),controlConditionID: controlConditionID,treatedConditionID: treatedConditionID,
            controlCellDispersion: controlCellDispersion,treatedCellDispersion: treatedCellDispersion,
            controlCells: control.map(\.cells),treatedCells: treated.map(\.cells),controlDonorMLERatesCPM: cm,treatedDonorMLERatesCPM: tm,
            controlRatesCPM: cr,treatedRatesCPM: tr,probabilities: w,fittedDonorRelativeLikelihoods: s,
            relativeLogLikelihood: s.reduce(0) { $0+log($1) },meanLogLikelihoodGap: gap,iterations: iteration,status: status,plan: plan,
            moments: moments(control: cr,treated: tr,weights: w),
            qualification: "Conditional finite-support empirical likelihood estimate from paired donors and independent NB2 cells. Cell dispersion is fixed; mixing-distribution and dispersion estimation uncertainty are not integrated. The convex certificate covers only the declared finite grid, not continuous-support optimality, unique mixing weights, calibrated biological uncertainty or full RNA composition. Zero atoms and degenerate fitted distributions are explicit. Rate support is selected from training donor likelihood modes; new-donor extrapolation beyond that support is not learned.")
    }

    public static func predict(control: VivoCountDepthStratum,featureID: String,queryDonorID: String,
        controlConditionID: String,querySource: VivoFingerprint,model: VivoJointCountModel,
        plannedTreatedLibraryCounts: [UInt64]) throws -> VivoJointCountPrediction {
        guard model.method=="paired-cell-nb2-finite-support-npmle-v1",model.status=="convergedFiniteGridLikelihood",
              featureID==model.featureID,vivoOmicsID(queryDonorID),!model.trainingDonorIDs.contains(queryDonorID),
              controlConditionID==model.controlConditionID,!model.controlRatesCPM.isEmpty,!model.treatedRatesCPM.isEmpty,
              model.probabilities.count==model.controlRatesCPM.count*model.treatedRatesCPM.count,
              model.probabilities.allSatisfy({ $0.isFinite && $0>=0 }),abs(model.probabilities.reduce(0,+)-1)<1e-10 else {
            throw VivoOmicsError.invalid("joint count query identity, probability model or unconverged fit")
        }
        let l=try VivoCountRateLikelihood(control,cellDispersion: model.controlCellDispersion)
        let ell=try model.controlRatesCPM.map { try l.logRelativeLikelihood(rateCPM: $0) }
        let logWeights=model.probabilities.indices.map { i in model.probabilities[i]>0 ? log(model.probabilities[i])+ell[i/model.treatedRatesCPM.count] : -Double.infinity }
        let maximum=logWeights.max()!
        guard maximum.isFinite else { throw VivoOmicsError.invalid("query has zero probability under fitted control rate support") }
        var weights=logWeights.map { exp($0-maximum) };let total=weights.reduce(0,+)
        for i in weights.indices { weights[i]/=total }
        // Mathematical support is determined before exponentiation. A positive
        // but negligible probability can underflow without being a point mass.
        let supported=logWeights.indices.filter { logWeights[$0].isFinite }
        let firstT=model.treatedRatesCPM[supported[0]%model.treatedRatesCPM.count]
        let degenerate=supported.allSatisfy { model.treatedRatesCPM[$0%model.treatedRatesCPM.count]==firstT }
        let underflowed=supported.filter { weights[$0]==0 }.count
        let m=moments(control: model.controlRatesCPM,treated: model.treatedRatesCPM,weights: weights)
        let sampling=try VivoCountObservation.predictiveMoments(meanCPM: m.treatedMeanCPM,varianceCPM: m.treatedVarianceCPM2,
            cellDispersion: model.treatedCellDispersion,plannedLibraryCounts: plannedTreatedLibraryCounts)
        return .init(querySource: querySource,queryDonorID: queryDonorID,featureID: featureID,treatedConditionID: model.treatedConditionID,
            probabilities: weights,moments: m,plannedTreatedCountMoments: sampling,degenerateTreatedRateDistribution: degenerate,
            underflowedPosteriorComponents: underflowed,numericallyDegenerateTreatedRateMoments: !degenerate && m.treatedVarianceCPM2==0,
            qualification: "Control-only update of a frozen empirical joint rate distribution. No query-treated outcomes enter this prediction. Latent rate moments and planned treated-count moments are different endpoints. Fitted point-mass degeneracy is retained; it is not proof of zero biological or parameter uncertainty. This is not a calibrated future-count interval, independent biological validation or a repaired prior pseudobulk interval.")
    }

    private static func moments(control c: [Double],treated t: [Double],weights w: [Double]) -> VivoJointRateMoments {
        var mc=0.0,mt=0.0,lc=0.0,lt=0.0
        for i in w.indices where w[i]>0 { let x=c[i/t.count],y=t[i%t.count];mc+=w[i]*x;mt+=w[i]*y;lc+=w[i]*log1p(x);lt+=w[i]*log1p(y) }
        let active=w.indices.filter { w[$0]>0 },firstC=c[active[0]/t.count],firstT=t[active[0]%t.count]
        if active.allSatisfy({ c[$0/t.count]==firstC }) { mc=firstC;lc=log1p(firstC) }
        if active.allSatisfy({ t[$0%t.count]==firstT }) { mt=firstT;lt=log1p(firstT) }
        var vc=0.0,vt=0.0,cov=0.0,vr=0.0,vl=0.0
        for i in w.indices where w[i]>0 {
            let x=c[i/t.count],y=t[i%t.count],dx=x-mc,dy=y-mt
            vc+=w[i]*dx*dx;vt+=w[i]*dy*dy;cov+=w[i]*dx*dy
            vr+=w[i]*pow((y-x)-(mt-mc),2);vl+=w[i]*pow((log1p(y)-log1p(x))-(lt-lc),2)
        }
        return .init(controlMeanCPM: mc,treatedMeanCPM: mt,controlVarianceCPM2: vc,treatedVarianceCPM2: vt,covarianceCPM2: cov,
            responseMeanCPM: mt-mc,responseVarianceCPM2: vr,controlMeanLog1pCPM: lc,treatedMeanLog1pCPM: lt,log1pResponseMean: lt-lc,log1pResponseVariance: vl)
    }
}
