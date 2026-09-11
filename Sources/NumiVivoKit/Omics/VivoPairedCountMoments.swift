import Foundation

/// A donor-condition stratum of complete RNA cell moments. The condition
/// samples must contain distinct, independent cells for the noise correction.
public struct VivoPairedCountGroup: Codable, Sendable, Equatable {
    public let donorID: String
    public let conditionID: String
    public let moments: VivoCellCountMoments
    public init(donorID: String, conditionID: String, moments: VivoCellCountMoments) {
        self.donorID=donorID;self.conditionID=conditionID;self.moments=moments
    }
}

public struct VivoPairedCountFeature: Codable, Sendable, Equatable {
    public let featureID: String
    public let controlCalibrationStatus: String
    public let treatedCalibrationStatus: String
    public let controlMeanCPM: Double
    public let treatedMeanCPM: Double
    public let observedControlVarianceCPM2: Double
    public let observedTreatedVarianceCPM2: Double
    public let observedCovarianceCPM2: Double
    public let observedResponseVarianceCPM2: Double
    public let controlMeasurementVarianceCPM2: Double?
    public let treatedMeasurementVarianceCPM2: Double?
    public let rawLatentControlVarianceCPM2: Double?
    public let rawLatentTreatedVarianceCPM2: Double?
    public let rawLatentResponseVarianceCPM2: Double?
    public let rawMinimumEigenvalueCPM2: Double?
    public let marginalClippedMinimumEigenvalueCPM2: Double?
    public let rawCovarianceStatus: String
    public let marginalClippedCovarianceStatus: String
    public let meanCellRateResponseCPM: Double
    public let meanPseudobulkResponseCPM: Double
    public let meanLog1pCellRateResponse: Double
    public let meanLog1pPseudobulkResponse: Double
    public let maximumAbsoluteEndpointDifferenceCPM: Double
    public let maximumAbsoluteLog1pEndpointDifference: Double
}

public struct VivoPairedCountReport: Codable, Sendable, Equatable {
    public let method: String
    public let trainingSource: VivoFingerprint
    public let controlConditionID: String
    public let treatedConditionID: String
    public let donorIDs: [String]
    public let controlCells: [Int]
    public let treatedCells: [Int]
    public let controlLibraryCounts: [UInt64]
    public let treatedLibraryCounts: [UInt64]
    public let covariancePSDRelativeTolerance: Double
    public let features: [VivoPairedCountFeature]
    public let qualification: String
}

public enum VivoPairedCountMoments {
    /// Diagnose whether separate marginal noise corrections describe a valid
    /// joint donor-rate covariance. An indefinite result is retained, never
    /// silently projected or used as predictive uncertainty. All donors have
    /// equal weight. Pairing is by explicit identity, never array position.
    public static func analyze(groups: [VivoPairedCountGroup],featureIDs: [String],
        controlConditionID: String,treatedConditionID: String,trainingSource: VivoFingerprint) throws -> VivoPairedCountReport {
        guard vivoOmicsID(controlConditionID),vivoOmicsID(treatedConditionID),controlConditionID != treatedConditionID,
              groups.count>=6,groups.count<=128,
              groups.allSatisfy({ vivoOmicsID($0.donorID) && ($0.conditionID==controlConditionID || $0.conditionID==treatedConditionID) }) else {
            throw VivoOmicsError.invalid("paired count conditions or groups")
        }
        let controls=groups.filter { $0.conditionID==controlConditionID }
        let treated=groups.filter { $0.conditionID==treatedConditionID }
        let ids=controls.map(\.donorID).sorted()
        guard Set(ids).count==ids.count,treated.count==ids.count,Set(treated.map(\.donorID))==Set(ids) else {
            throw VivoOmicsError.invalid("paired count duplicate or missing donor-condition stratum")
        }
        let c=ids.map { id in controls.first { $0.donorID==id }!.moments }
        let t=ids.map { id in treated.first { $0.donorID==id }!.moments }
        let cm=try VivoCountObservationCalibration.fit(c,featureIDs: featureIDs,donorIDs: ids,conditionID: controlConditionID,trainingSource: trainingSource)
        let tm=try VivoCountObservationCalibration.fit(t,featureIDs: featureIDs,donorIDs: ids,conditionID: treatedConditionID,trainingSource: trainingSource)
        let n=Double(ids.count),tolerance=1e-10
        // Eigenvalues are evaluated on a scaled symmetric 2x2 matrix. This
        // avoids overflow and does not turn a numerical tolerance into a prior.
        func eigen(_ a: Double,_ b: Double,_ off: Double) -> (Double,String) {
            let scale=max(abs(a),abs(b),abs(off))
            if scale==0 { return (0,"positiveSemidefinite") }
            let x=a/scale,y=b/scale,z=off/scale
            let minimum=0.5*(x+y-hypot(x-y,2*z))
            return (minimum*scale,minimum < -tolerance ? "indefinite" : "positiveSemidefinite")
        }
        var features: [VivoPairedCountFeature]=[]
        for j in featureIDs.indices {
            if j%128==0 { try Task.checkCancellation() }
            let cf=cm.features[j],tf=tm.features[j],meanC=cf.meanDonorRateCPM,meanT=tf.meanDonorRateCPM
            var covariance=0.0,responseVariance=0.0,pseudobulkResponse=0.0,cellLogResponse=0.0,bulkLogResponse=0.0,maxDifference=0.0,maxLogDifference=0.0
            for i in ids.indices {
                let x=c[i].meanCPM[j],y=t[i].meanCPM[j]
                let bx=Double(c[i].counts[j])/Double(c[i].libraryCounts)*1e6
                let by=Double(t[i].counts[j])/Double(t[i].libraryCounts)*1e6
                covariance+=(x-meanC)*(y-meanT)/(n-1)
                responseVariance+=pow((y-x)-(meanT-meanC),2)/(n-1)
                pseudobulkResponse+=(by-bx)/n
                cellLogResponse+=(log1p(y)-log1p(x))/n
                bulkLogResponse+=(log1p(by)-log1p(bx))/n
                maxDifference=max(maxDifference,abs(x-bx),abs(y-by))
                maxLogDifference=max(maxLogDifference,abs(log1p(x)-log1p(bx)),abs(log1p(y)-log1p(by)))
            }
            var rawC: Double?,rawT: Double?,rawResponse: Double?,rawEigen: Double?,clippedEigen: Double?
            var rawStatus="unavailableCellDispersion",clippedStatus=rawStatus
            // Marginal measurement moments remain usable when only the Gamma
            // prior is out of range; this analysis does not require that prior.
            if let vc=cf.meanDonorMeasurementVarianceCPM,let vt=tf.meanDonorMeasurementVarianceCPM {
                rawC=cf.observedDonorRateVarianceCPM-vc;rawT=tf.observedDonorRateVarianceCPM-vt
                rawResponse=responseVariance-vc-vt
                (rawEigen,rawStatus)=eigen(rawC!,rawT!,covariance)
                (clippedEigen,clippedStatus)=eigen(max(0,rawC!),max(0,rawT!),covariance)
            }
            features.append(.init(featureID: featureIDs[j],controlCalibrationStatus: cf.status,treatedCalibrationStatus: tf.status,
                controlMeanCPM: meanC,treatedMeanCPM: meanT,observedControlVarianceCPM2: cf.observedDonorRateVarianceCPM,
                observedTreatedVarianceCPM2: tf.observedDonorRateVarianceCPM,observedCovarianceCPM2: covariance,
                observedResponseVarianceCPM2: responseVariance,controlMeasurementVarianceCPM2: cf.meanDonorMeasurementVarianceCPM,
                treatedMeasurementVarianceCPM2: tf.meanDonorMeasurementVarianceCPM,rawLatentControlVarianceCPM2: rawC,
                rawLatentTreatedVarianceCPM2: rawT,rawLatentResponseVarianceCPM2: rawResponse,
                rawMinimumEigenvalueCPM2: rawEigen,marginalClippedMinimumEigenvalueCPM2: clippedEigen,
                rawCovarianceStatus: rawStatus,marginalClippedCovarianceStatus: clippedStatus,
                meanCellRateResponseCPM: meanT-meanC,meanPseudobulkResponseCPM: pseudobulkResponse,
                meanLog1pCellRateResponse: cellLogResponse,meanLog1pPseudobulkResponse: bulkLogResponse,
                maximumAbsoluteEndpointDifferenceCPM: maxDifference,maximumAbsoluteLog1pEndpointDifference: maxLogDifference))
        }
        return .init(method: "paired-donor-count-moment-compatibility-v1",trainingSource: trainingSource,
            controlConditionID: controlConditionID,treatedConditionID: treatedConditionID,donorIDs: ids,
            controlCells: c.map(\.cells),treatedCells: t.map(\.cells),controlLibraryCounts: c.map(\.libraryCounts),treatedLibraryCounts: t.map(\.libraryCounts),
            covariancePSDRelativeTolerance: tolerance,features: features,
            qualification: "Paired training-donor moment diagnosis, not a fitted joint likelihood or treatment prediction. Independent control/treated cell measurement errors contribute only diagonal noise. Cell dispersion is a conditional estimate; parameter uncertainty is not integrated. Raw noise-corrected and separately clipped marginal covariance failures are retained. Mean per-cell normalized rates and RNA-weighted pseudobulk endpoints are distinct. No PSD projection, biological calibration, causal effect, joint RNA composition or future-count interval is established.")
    }
}
