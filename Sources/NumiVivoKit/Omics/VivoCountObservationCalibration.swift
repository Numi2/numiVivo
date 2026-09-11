import Foundation

/// Sufficient moments of normalized per-cell rates; missing sparse entries are
/// measured zeros. Full-library normalization precedes any feature projection.
public struct VivoCellCountMoments: Codable, Sendable, Equatable {
    public let cells: Int
    public let libraryCounts: UInt64
    public let counts: [UInt64]
    public let positiveCells: [Int]
    public let meanCPM: [Double]
    public let sampleVarianceCPM: [Double]
    /// Mean(Y_i / e_i^2), estimating the average conditional Poisson variance.
    public let meanPoissonVarianceCPM: [Double]
    /// Mean(z_i*z_j) over distinct cells, an unbiased estimator of r^2.
    public let distinctCellRateProductCPM2: [Double]
}

/// One sparse cell at a time, with no dense cells-by-genes allocation.
/// A cell must include ALL RNA features, in strictly increasing index order.
public final class VivoCellCountMomentAccumulator {
    private let features: Int
    private var cells=0, libraries: UInt64=0
    private var counts: [UInt64], positives: [Int]
    private var means: [Double], m2: [Double], sums: [Double], products: [Double], shot: [Double]
    public init(featureCount: Int) throws {
        guard (1...100_000).contains(featureCount) else { throw VivoOmicsError.limit("cell moment feature count") }
        features=featureCount; counts=Array(repeating: 0,count: features);positives=Array(repeating: 0,count: features)
        means=Array(repeating: 0,count: features);m2=means;sums=means;products=means;shot=means
    }
    public func appendCell(featureIndices: [Int], counts values: [UInt64]) throws {
        guard cells<10_000_000,featureIndices.count==values.count,!values.isEmpty else {
            throw VivoOmicsError.invalid("cell moment dimensions, empty cell or cell limit")
        }
        var library: UInt64=0,previous = -1
        for k in values.indices {
            guard featureIndices[k]>previous,featureIndices[k]<features,values[k]>0 else {
                throw VivoOmicsError.invalid("cell moments require canonical positive sparse counts")
            }
            library=try vivoOmicsSum(library,values[k]);previous=featureIndices[k]
        }
        guard library<=1_000_000_000 else { throw VivoOmicsError.limit("cell moment RNA library budget") }
        // Validate every checked sum before mutation: a rejected cell is absent.
        let newLibraries=try vivoOmicsSum(libraries,library)
        for k in values.indices { _=try vivoOmicsSum(counts[featureIndices[k]],values[k]) }
        let exposure=Double(library)/1e6
        for k in values.indices {
            let j=featureIndices[k],z=Double(values[k])/exposure
            counts[j]+=values[k];positives[j]+=1
            let delta=z-means[j];means[j]+=delta/Double(positives[j]);m2[j]+=delta*(z-means[j])
            products[j]+=sums[j]*z;sums[j]+=z;shot[j]+=z/exposure
        }
        cells+=1;libraries=newLibraries
    }
    public func finish() throws -> VivoCellCountMoments {
        guard cells>=2 else { throw VivoOmicsError.invalid("cell moments require at least two observed cells") }
        let n=Double(cells)
        var mean=means,variance=means,poisson=means,pairs=means
        for j in 0..<features {
            mean[j]=sums[j]/n
            // Merge the positive-observation Welford state with implicit zeros.
            let zeros=n-Double(positives[j])
            variance[j]=(m2[j]+means[j]*means[j]*Double(positives[j])*zeros/n)/(n-1)
            poisson[j]=shot[j]/n;pairs[j]=2*products[j]/(n*(n-1))
        }
        return .init(cells: cells,libraryCounts: libraries,counts: counts,positiveCells: positives,
            meanCPM: mean,sampleVarianceCPM: variance,meanPoissonVarianceCPM: poisson,distinctCellRateProductCPM2: pairs)
    }
}

public struct VivoCountObservationCalibrationFeature: Codable, Sendable, Equatable {
    public let featureID: String
    public let status: String
    public let trainingCounts: UInt64
    public let positiveDonors: Int
    public let cellDispersionNumerator: Double
    public let cellDispersionDenominator: Double
    public let rawCellDispersion: Double?
    public let cellDispersion: Double?
    public let poissonBoundary: Bool
    public let meanDonorRateCPM: Double
    public let observedDonorRateVarianceCPM: Double
    public let meanDonorMeasurementVarianceCPM: Double?
    public let latentDonorRateVarianceCPM: Double?
    public let latentVarianceBoundary: Bool
    public let newDonorRateVarianceCPM: Double?
    public let gammaPriorShape: Double?
    public let gammaPriorRatePerCPM: Double?
}
public struct VivoCountObservationCalibrationModel: Codable, Sendable, Equatable {
    public let method: String
    public let trainingSource: VivoFingerprint
    public let conditionID: String
    public let trainingDonorIDs: [String]
    public let trainingCellCounts: [Int]
    public let features: [VivoCountObservationCalibrationFeature]
    public let qualification: String
}

public enum VivoCountObservationCalibration {
    /// Method-of-moments calibration, preserving donor intercepts. This is not
    /// the existing donor-level NB GLM's Cox-Reid profile likelihood estimator.
    public static func fit(_ groups: [VivoCellCountMoments],featureIDs: [String],donorIDs: [String],
        conditionID: String,trainingSource: VivoFingerprint) throws -> VivoCountObservationCalibrationModel {
        guard (3...64).contains(groups.count),groups.count==donorIDs.count,
              Set(donorIDs).count==donorIDs.count,donorIDs.allSatisfy(vivoOmicsID),vivoOmicsID(conditionID),
              !featureIDs.isEmpty,featureIDs.count<=100_000,Set(featureIDs).count==featureIDs.count,
              featureIDs.allSatisfy(vivoOmicsID) else { throw VivoOmicsError.invalid("cell calibration identities or independent donor count") }
        let m=featureIDs.count,n=Double(groups.count)
        for g in groups {
            guard g.cells>=2,g.libraryCounts>0,g.counts.count==m,g.positiveCells.count==m,
                  g.positiveCells.allSatisfy({ (0...g.cells).contains($0) }),
                  [g.meanCPM,g.sampleVarianceCPM,g.meanPoissonVarianceCPM,g.distinctCellRateProductCPM2].allSatisfy({ $0.count==m && $0.allSatisfy { $0.isFinite && $0>=0 } }) else {
                throw VivoOmicsError.invalid("cell calibration moments")
            }
        }
        var results: [VivoCountObservationCalibrationFeature]=[]
        for j in 0..<m {
            if j%128==0 { try Task.checkCancellation() }
            var numerator=0.0,denominator=0.0,mean=0.0,counts: UInt64=0,positive=0
            for g in groups {
                let df=Double(g.cells-1)
                numerator+=df*(g.sampleVarianceCPM[j]-g.meanPoissonVarianceCPM[j])
                denominator+=df*g.distinctCellRateProductCPM2[j]
                mean+=g.meanCPM[j]/n;counts=try vivoOmicsSum(counts,g.counts[j]);if g.counts[j]>0 { positive+=1 }
            }
            let observed=groups.reduce(0.0) { $0+pow($1.meanCPM[j]-mean,2) }/(n-1)
            var raw: Double?,phi: Double?,measurement: Double?,latent: Double?,future: Double?,shape: Double?,rate: Double?
            var status="insufficientWithinDonorCountPairs",poisson=false,boundary=false
            if denominator>0 {
                raw=numerator/denominator
                let value=max(0,raw!)
                if value.isFinite && (value==0 || (1e-8...100).contains(value)) {
                    phi=value;poisson=raw!<=0
                    measurement=groups.reduce(0.0) { $0+($1.meanPoissonVarianceCPM[j]+value*$1.distinctCellRateProductCPM2[j])/Double($1.cells) }/n
                    let rawLatent=observed-measurement!
                    latent=max(0,rawLatent);boundary=rawLatent<=0
                    // New-donor variance includes estimated uncertainty in the
                    // equally weighted training mean, even at zero latent tau.
                    future=latent!*(1+1/n)+measurement!/n
                    if mean>0,future!>0 {
                        let a=mean*mean/future!,b=mean/future!
                        if a.isFinite,b.isFinite,(0.1...1e6).contains(a),(1e-12...1e12).contains(b) {
                            shape=a;rate=b;status="availableConditionalMomentCalibration"
                        } else { status="gammaPriorOutsideObservationKernelDomain" }
                    } else { status="unidentifiedPositiveRatePrior" }
                } else { status="cellDispersionOutsideObservationKernelDomain" }
            }
            results.append(.init(featureID: featureIDs[j],status: status,trainingCounts: counts,positiveDonors: positive,
                cellDispersionNumerator: numerator,cellDispersionDenominator: denominator,rawCellDispersion: raw,
                cellDispersion: phi,poissonBoundary: poisson,meanDonorRateCPM: mean,observedDonorRateVarianceCPM: observed,
                meanDonorMeasurementVarianceCPM: measurement,latentDonorRateVarianceCPM: latent,
                latentVarianceBoundary: boundary,newDonorRateVarianceCPM: future,gammaPriorShape: shape,gammaPriorRatePerCPM: rate))
        }
        return .init(method: "donor-stratified-cell-nb2-moments-gamma-new-donor-rate-v1",trainingSource: trainingSource,
            conditionID: conditionID,trainingDonorIDs: donorIDs,trainingCellCounts: groups.map(\.cells),features: results,
            qualification: "Conditional empirical moment calibration from independent cells within donors and exchangeable donors. Full RNA depths are fixed offsets; cell rate is assumed independent of depth. Gene-specific dispersion is shared across donors within a condition. Hyperparameter estimation uncertainty is not integrated. This is not Cox-Reid fitting, learned perturbation response, calibrated biological intervals or a joint RNA composition.")
    }
    /// Apply a frozen, source-bound calibration to a distinct donor's control
    /// counts. Condition identity is explicit; treated response is not inferred.
    public static func posterior(counts: [UInt64],libraryCounts: [UInt64],featureID: String,
        queryDonorID: String,conditionID: String,model: VivoCountObservationCalibrationModel) throws -> VivoCountObservationPosterior {
        guard vivoOmicsID(queryDonorID),!model.trainingDonorIDs.contains(queryDonorID),conditionID==model.conditionID,
              model.method=="donor-stratified-cell-nb2-moments-gamma-new-donor-rate-v1",
              let f=model.features.first(where: { $0.featureID==featureID }),
              f.status=="availableConditionalMomentCalibration",let phi=f.cellDispersion,
              let a=f.gammaPriorShape,let b=f.gammaPriorRatePerCPM else {
            throw VivoOmicsError.invalid("count calibration query identity, overlapping donor or unavailable parameters")
        }
        return try VivoCountObservation.posterior(counts: counts,libraryCounts: libraryCounts,cellDispersion: phi,gammaPriorShape: a,gammaPriorRatePerCPM: b)
    }
}
