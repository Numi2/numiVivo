import Foundation

public enum VivoOmicsNBSupportOutcome: String, Codable, Sendable {
    case ready, insufficientReplication, rankDeficientDesign, rankDeficientPositiveSupport
}
public struct VivoOmicsNBSupportResolution: Codable, Sendable, Equatable {
    public let method: String
    public var outcome: VivoOmicsNBSupportOutcome
    /// Indices address the contrast's original design.observations, not cells.
    public let retainedObservationIndices: [Int]
    public let excludedObservationIndices: [Int]
    public let excludedDonorIDs: [String]
    public let activeDonorIDs: [String]
    public var columnNames: [String]?
    public var rows: [[Double]]?
    public var contrast: [Double]?
    public var residualDegreesOfFreedom: Int?
    public let meanNormalizedCount: Double
    public let qualification: String
}

enum VivoOmicsNBSupport {
    /// Resolve complete zero-total pairs only. Partial pairs and nonzero low
    /// counts are never removed. A nil resolution preserves the original design.
    static func resolve(counts: [UInt64],design: VivoOmicsDesignMatrix,request: VivoOmicsExpressionContrast) throws -> VivoOmicsNBSupportResolution? {
        guard request.design == .pairedDonors, counts.count==design.observations.count,
              counts.count==design.sizeFactorValues.count,design.observations.allSatisfy({ $0.donorID != nil }) else { throw VivoOmicsError.invalid("NB active-donor source dimensions/design") }
        let donors=Set(design.observations.compactMap(\.donorID)).sorted()
        var excludedDonors: [String]=[],excludedRows: [Int]=[]
        for donor in donors {
            let rows=design.observations.indices.filter { design.observations[$0].donorID==donor }
            guard rows.count==2,Set(rows.map { design.observations[$0].condition })==[request.controlCondition,request.treatmentCondition] else {
                throw VivoOmicsError.invalid("NB active-donor policy requires complete paired observations")
            }
            if rows.allSatisfy({ counts[$0]==0 }) { excludedDonors.append(donor);excludedRows+=rows }
        }
        guard !excludedDonors.isEmpty else { return nil }
        let excluded=Set(excludedRows),retained=design.observations.indices.filter { !excluded.contains($0) }
        let active=donors.filter { !excludedDonors.contains($0) }
        let mean=retained.reduce(0.0) { $0+Double(counts[$1])/design.sizeFactorValues[$1]/Double(retained.count) }
        var result=VivoOmicsNBSupportResolution(method: "complete-zero-total-paired-donor-boundary-v1",outcome: .insufficientReplication,
            retainedObservationIndices: retained,excludedObservationIndices: excludedRows.sorted(),excludedDonorIDs: excludedDonors,activeDonorIDs: active,
            meanNormalizedCount: mean,
            qualification: "Gene-specific active-donor adjusted profile; zero-total nuisance donor rates lie at the likelihood boundary. Original counts/offsets retained, full-support prior cohort unchanged. No selection-adjusted FDR or cross-donor generalization qualification.")
        guard active.count>=request.minimumReplicatesPerCondition else { return result }
        do {
            let built=try VivoPseudobulkDifferentialExpression.makeDesign(observations: retained.map { design.observations[$0] },request: request)
            result.columnNames=built.columnNames;result.rows=built.rows;result.contrast=built.contrast
            result.residualDegreesOfFreedom=retained.count-built.columnNames.count
            result.outcome=VivoOmicsNegativeBinomial.positiveSupportIsRankDeficient(counts: retained.map { counts[$0] },design: built.rows) ? .rankDeficientPositiveSupport : .ready
        } catch is CancellationError { throw CancellationError() }
        catch { result.outcome = .rankDeficientDesign }
        return result
    }
}
