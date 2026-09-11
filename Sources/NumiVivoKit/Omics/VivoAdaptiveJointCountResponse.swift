import Foundation

public struct VivoAdaptiveJointCountPlan: Codable, Sendable, Equatable {
    public let initialGridPointsPerAxis: Int
    public let maximumSupportAdditions: Int
    public let maximumOracleLeaves: Int
    public let maximumWeightIterations: Int
    public let meanLogLikelihoodGapTolerance: Double
    public init(initialGridPointsPerAxis: Int=9,maximumSupportAdditions: Int=64,
        maximumOracleLeaves: Int=4096,maximumWeightIterations: Int=200_000,meanLogLikelihoodGapTolerance: Double=1e-6) {
        self.initialGridPointsPerAxis=initialGridPointsPerAxis;self.maximumSupportAdditions=maximumSupportAdditions
        self.maximumOracleLeaves=maximumOracleLeaves;self.maximumWeightIterations=maximumWeightIterations;self.meanLogLikelihoodGapTolerance=meanLogLikelihoodGapTolerance
    }
}
public struct VivoJointSupportBox: Codable, Sendable, Equatable {
    /// Binary subdivision path from the full training-MLE rate rectangle.
    public let path: String
    public let controlLow: Double
    public let controlHigh: Double
    public let treatedLow: Double
    public let treatedHigh: Double
    public let meanDirectionalUpperBound: Double
}
public struct VivoJointSupportCertificate: Codable, Sendable, Equatable {
    public let status: String
    public let controlScaleCPM: Double
    public let treatedScaleCPM: Double
    public let maximumMeanDirectionalLowerBound: Double
    public let maximumMeanDirectionalUpperBound: Double
    public let bestControlRateCPM: Double
    public let bestTreatedRateCPM: Double
    public let boxes: [VivoJointSupportBox]
    public let likelihoodEvaluations: Int
    public let qualification: String
}
public struct VivoAdaptiveJointCountStep: Codable, Sendable, Equatable {
    public let supportAdditions: Int
    public let controlSupportCount: Int
    public let treatedSupportCount: Int
    public let relativeLogLikelihood: Double
    public let finiteGridGap: Double
    public let oracleStatus: String
    public let directionalLowerBound: Double
    public let directionalUpperBound: Double
}
public struct VivoAdaptiveJointCountModel: Codable, Sendable, Equatable {
    public let status: String
    public let plan: VivoAdaptiveJointCountPlan
    public let model: VivoJointCountModel
    public let certificate: VivoJointSupportCertificate?
    public let steps: [VivoAdaptiveJointCountStep]
    public let qualification: String
}

public enum VivoAdaptiveJointCountResponse {
    public static func fit(pairs: [VivoJointCountPair],featureID: String,controlConditionID: String,
        treatedConditionID: String,controlCellDispersion: Double,treatedCellDispersion: Double,
        trainingSource: VivoFingerprint,plan: VivoAdaptiveJointCountPlan = .init()) throws -> VivoAdaptiveJointCountModel {
        guard (3...65).contains(plan.initialGridPointsPerAxis),(0...128).contains(plan.maximumSupportAdditions),
              (1...65_536).contains(plan.maximumOracleLeaves),(1...200_000).contains(plan.maximumWeightIterations),(1e-7...1e-3).contains(plan.meanLogLikelihoodGapTolerance) else {
            throw VivoOmicsError.invalid("adaptive joint count numerical plan")
        }
        let ordered=pairs.sorted { $0.donorID<$1.donorID }
        let finite=VivoJointCountPlan(gridPointsPerAxis: plan.initialGridPointsPerAxis,maximumIterations: plan.maximumWeightIterations)
        var m=try VivoJointCountResponse.fit(pairs: ordered,featureID: featureID,controlConditionID: controlConditionID,
            treatedConditionID: treatedConditionID,controlCellDispersion: controlCellDispersion,
            treatedCellDispersion: treatedCellDispersion,trainingSource: trainingSource,plan: finite)
        let c=try ordered.map { try VivoCountRateLikelihood($0.control,cellDispersion: controlCellDispersion) }
        let t=try ordered.map { try VivoCountRateLikelihood($0.treated,cellDispersion: treatedCellDispersion) }
        var steps=[VivoAdaptiveJointCountStep](),certificate: VivoJointSupportCertificate?,status="finiteSupportFitNotConverged"
        for addition in 0...plan.maximumSupportAdditions {
            try Task.checkCancellation()
            guard m.status=="convergedFiniteGridLikelihood" else { certificate=nil;status="finiteSupportFitNotConverged";break }
            let cert=try search(control: c,treated: t,model: m,plan: plan);certificate=cert
            steps.append(.init(supportAdditions: addition,controlSupportCount: m.controlRatesCPM.count,
                treatedSupportCount: m.treatedRatesCPM.count,relativeLogLikelihood: m.relativeLogLikelihood,
                finiteGridGap: m.meanLogLikelihoodGap,oracleStatus: cert.status,
                directionalLowerBound: cert.maximumMeanDirectionalLowerBound,directionalUpperBound: cert.maximumMeanDirectionalUpperBound))
            if cert.status=="boundedContinuousLikelihood" { status=cert.status;break }
            if cert.status != "improvingSupportPoint" { status=cert.status;break }
            if addition==plan.maximumSupportAdditions { status="supportAdditionLimit";break }
            let cr=Array(Set(m.controlRatesCPM+[cert.bestControlRateCPM])).sorted()
            let tr=Array(Set(m.treatedRatesCPM+[cert.bestTreatedRateCPM])).sorted()
            guard cr.count>m.controlRatesCPM.count || tr.count>m.treatedRatesCPM.count else { status="supportPointAlreadyPresent";break }
            m=try VivoJointCountResponse.fitRateGrid(pairs: ordered,featureID: featureID,controlConditionID: controlConditionID,
                treatedConditionID: treatedConditionID,controlCellDispersion: controlCellDispersion,
                treatedCellDispersion: treatedCellDispersion,trainingSource: trainingSource,plan: finite,
                control: c,treated: t,controlRates: cr,treatedRates: tr,initial: m,preferLargestExchangeGain: true)
        }
        return .init(status: status,plan: plan,model: m,certificate: certificate,steps: steps,
            qualification: "Adaptive nonnegative paired-rate likelihood, conditional on fixed cell dispersions. The analytic compact-domain upper bound is evaluated in FP64 with an explicit numerical allowance, not directed-rounding interval arithmetic. A small likelihood gap does not establish uniquely identified mixing weights, stable posterior moments, parameter uncertainty or biological prediction calibration. Original finite-grid refinement failures remain separate evidence.")
    }

    public static func predict(control: VivoCountDepthStratum,featureID: String,queryDonorID: String,
        controlConditionID: String,querySource: VivoFingerprint,model: VivoAdaptiveJointCountModel,
        plannedTreatedLibraryCounts: [UInt64]) throws -> VivoJointCountPrediction {
        guard model.status=="boundedContinuousLikelihood",let cert=model.certificate,
              cert.maximumMeanDirectionalUpperBound<=1+model.plan.meanLogLikelihoodGapTolerance else {
            throw VivoOmicsError.invalid("adaptive joint support not qualified for prediction")
        }
        return try VivoJointCountResponse.predict(control: control,featureID: featureID,queryDonorID: queryDonorID,
            controlConditionID: controlConditionID,querySource: querySource,model: model.model,
            plannedTreatedLibraryCounts: plannedTreatedLibraryCounts)
    }

    private struct AxisValue { let rates: [Double];let ell: [Double];let slope: [Double] }
    private static func search(control c: [VivoCountRateLikelihood],treated t: [VivoCountRateLikelihood],
        model m: VivoJointCountModel,plan: VivoAdaptiveJointCountPlan) throws -> VivoJointSupportCertificate {
        let cs=c.map(\.concaveCoordinateScaleCPM).min()!,ts=t.map(\.concaveCoordinateScaleCPM).min()!
        let cm=c.map(\.maximumLikelihoodRateCPM),tm=t.map(\.maximumLikelihoodRateCPM)
        let cl=log1p(cm.min()!/cs),ch=log1p(cm.max()!/cs),tl=log1p(tm.min()!/ts),th=log1p(tm.max()!/ts)
        var cc=[Double:AxisValue](),tc=[Double:AxisValue](),evaluations=0
        func axis(_ x: Double,_ control: Bool) throws -> AxisValue {
            if let v=control ? cc[x] : tc[x] { return v }
            let scale=control ? cs : ts,likes=control ? c : t,r=scale*expm1(x)
            let v=AxisValue(rates: Array(repeating: r,count: likes.count),ell: try likes.map { try $0.logRelativeLikelihood(rateCPM: r) },
                slope: try likes.map { try $0.coordinateSlope(rateCPM: r,scaleCPM: scale) })
            if control { cc[x]=v } else { tc[x]=v };evaluations+=likes.count;return v
        }
        let inv=m.fittedDonorRelativeLikelihoods.map { 1/($0*Double(c.count)) }
        var lower=0.0,bestC=cm[0],bestT=tm[0]
        func point(_ x: Double,_ y: Double) throws {
            let a=try axis(x,true),b=try axis(y,false)
            let score=c.indices.reduce(0.0) { $0+exp(a.ell[$1]+b.ell[$1])*inv[$1] }
            if score>lower { lower=score;bestC=a.rates[0];bestT=b.rates[0] }
        }
        func box(_ path: String,_ x0: Double,_ x1: Double,_ y0: Double,_ y1: Double) throws -> VivoJointSupportBox {
            let x=(x0+x1)/2,y=(y0+y1)/2
            let a0=try axis(x0,true),a1=try axis(x1,true),b0=try axis(y0,false),b1=try axis(y1,false)
            let am=try axis(x,true),bm=try axis(y,false)
            var individual=0.0
            for d in c.indices {
                let a=cm[d]>=a0.rates[d] && cm[d]<=a1.rates[d] ? 0 : max(a0.ell[d],a1.ell[d])
                let b=tm[d]>=b0.rates[d] && tm[d]<=b1.rates[d] ? 0 : max(b0.ell[d],b1.ell[d])
                individual+=exp(a+b)*inv[d]
            }
            var tangent=0.0
            for xx in [x0,x1] { for yy in [y0,y1] {
                var bound=0.0
                for d in c.indices {
                    // A zero-width axis has no tangent displacement, including
                    // positive-count likelihoods impossible at an all-zero rate.
                    let dx=xx==x ? 0 : am.slope[d]*(xx-x),dy=yy==y ? 0 : bm.slope[d]*(yy-y)
                    bound+=exp(am.ell[d]+bm.ell[d]+dx+dy)*inv[d]
                }
                if bound.isNaN { tangent = .infinity } else { tangent=max(tangent,bound) }
                try point(xx,yy)
            } }
            try point(x,y)
            let raw=min(individual,tangent)
            // This allowance is recorded, independently checked and is not a
            // substitute for directed rounding or uncertainty in fitted data.
            let upper=raw+1e-9*(1+raw)
            guard upper.isFinite,upper>=0 else { throw VivoOmicsError.invalid("joint support oracle upper bound") }
            return .init(path: path,controlLow: x0,controlHigh: x1,treatedLow: y0,treatedHigh: y1,meanDirectionalUpperBound: upper)
        }
        var leaves=[try box("",cl,ch,tl,th)],status="oracleLeafLimit"
        while true {
            try Task.checkCancellation()
            let j=leaves.indices.max { leaves[$0].meanDirectionalUpperBound<leaves[$1].meanDirectionalUpperBound }!
            let upper=leaves[j].meanDirectionalUpperBound
            if upper<=1+plan.meanLogLikelihoodGapTolerance { status="boundedContinuousLikelihood";break }
            // Reserve headroom for the recorded floating-point upper-bound
            // allowance. Waiting until the witness itself exceeds the final
            // target can subdivide forever when lower < target < lower+allowance.
            if lower>1+0.9*plan.meanLogLikelihoodGapTolerance { status="improvingSupportPoint";break }
            if leaves.count>=plan.maximumOracleLeaves { break }
            let b=leaves.remove(at: j)
            let splitC=(ch>cl) && ((th==tl) || (b.controlHigh-b.controlLow)/(ch-cl)>=(b.treatedHigh-b.treatedLow)/(th-tl))
            if splitC {
                let mid=(b.controlLow+b.controlHigh)/2
                guard mid>b.controlLow,mid<b.controlHigh else { leaves.append(b);status="oracleResolutionLimit";break }
                leaves.append(try box(b.path+"0",b.controlLow,mid,b.treatedLow,b.treatedHigh))
                leaves.append(try box(b.path+"1",mid,b.controlHigh,b.treatedLow,b.treatedHigh))
            } else {
                let mid=(b.treatedLow+b.treatedHigh)/2
                guard mid>b.treatedLow,mid<b.treatedHigh else { leaves.append(b);status="oracleResolutionLimit";break }
                leaves.append(try box(b.path+"0",b.controlLow,b.controlHigh,b.treatedLow,mid))
                leaves.append(try box(b.path+"1",b.controlLow,b.controlHigh,mid,b.treatedHigh))
            }
        }
        return .init(status: status,controlScaleCPM: cs,treatedScaleCPM: ts,maximumMeanDirectionalLowerBound: lower,
            maximumMeanDirectionalUpperBound: leaves.map(\.meanDirectionalUpperBound).max()!,bestControlRateCPM: bestC,
            bestTreatedRateCPM: bestT,boxes: leaves.sorted { $0.path<$1.path },likelihoodEvaluations: evaluations,
            qualification: "The binary leaves cover the compact product of training donor rate-MLE ranges, in x=log1p(r/scale) coordinates. Outside this rectangle moving inward improves every donor likelihood. Each upper bound is the minimum of the sum of separate donor maxima and the corner maximum of exponential log-likelihood tangent planes, plus 1e-9*(1+bound). FP64 evaluation is not a formal interval-arithmetic certificate. Lower bounds are evaluated support witnesses; a leaf or support budget is not convergence.")
    }
}
