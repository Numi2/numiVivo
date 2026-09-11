import Foundation

public struct VivoDurationExposure: Codable, Sendable, Equatable {
    public let condition: String
    public let hours: Double
    public init(condition: String, hours: Double) { self.condition = condition; self.hours = hours }
    private enum CodingKeys: String, CodingKey { case condition, hours }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["condition", "hours"])
        let c = try decoder.container(keyedBy: CodingKeys.self)
        condition = try c.decode(String.self, forKey: .condition); hours = try c.decode(Double.self, forKey: .hours)
    }
}

public struct VivoDurationPlan: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let mapping: VivoH5ADImportPlan
    public let featureNamespace: String
    public let perturbationID: String
    public let controlCondition: String
    public let exposures: [VivoDurationExposure]
    public let responseFeatureIDs: [String]?
    public let provenance: String
    public init(mapping: VivoH5ADImportPlan, featureNamespace: String, perturbationID: String,
                controlCondition: String, exposures: [VivoDurationExposure], responseFeatureIDs: [String]? = nil, provenance: String) {
        schemaVersion = 1; self.mapping = mapping; self.featureNamespace = featureNamespace; self.perturbationID = perturbationID
        self.controlCondition = controlCondition; self.exposures = exposures; self.responseFeatureIDs = responseFeatureIDs; self.provenance = provenance
    }
    private enum CodingKeys: String, CodingKey { case schemaVersion, mapping, featureNamespace, perturbationID, controlCondition, exposures, responseFeatureIDs, provenance }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["schemaVersion", "mapping", "featureNamespace", "perturbationID", "controlCondition", "exposures", "responseFeatureIDs", "provenance"])
        let c = try decoder.container(keyedBy: CodingKeys.self)
        schemaVersion = try c.decode(Int.self, forKey: .schemaVersion); mapping = try c.decode(VivoH5ADImportPlan.self, forKey: .mapping)
        featureNamespace = try c.decode(String.self, forKey: .featureNamespace); perturbationID = try c.decode(String.self, forKey: .perturbationID)
        controlCondition = try c.decode(String.self, forKey: .controlCondition); exposures = try c.decode([VivoDurationExposure].self, forKey: .exposures)
        responseFeatureIDs = try c.decodeIfPresent([String].self, forKey: .responseFeatureIDs); provenance = try c.decode(String.self, forKey: .provenance)
    }
    func validate() throws {
        guard schemaVersion == 1, [featureNamespace, perturbationID, controlCondition].allSatisfy(vivoOmicsID),
              !provenance.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty, provenance.utf8.count <= 16_384,
              (2...63).contains(exposures.count), exposures.allSatisfy({ vivoOmicsID($0.condition) && $0.condition != controlCondition && $0.hours.isFinite && $0.hours > 0 && $0.hours <= 87_600 }),
              Set(exposures.map(\.condition)).count == exposures.count, Set(exposures.map(\.hours)).count == exposures.count,
              mapping.samples.allSatisfy({ $0.condition == controlCondition || exposures.map(\.condition).contains($0.condition) }) else {
            throw VivoOmicsError.invalid("duration plan identity, provenance, exposure or sample conditions")
        }
        if let ids = responseFeatureIDs {
            guard !ids.isEmpty, ids.count <= 100_000, Set(ids).count == ids.count, ids.allSatisfy(vivoOmicsID) else { throw VivoOmicsError.invalid("duration response panel") }
        }
    }
}

public struct VivoDurationQueryPlan: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let mapping: VivoH5ADImportPlan
    public let featureNamespace: String
    public let perturbationID: String
    public let hours: [Double]
    public init(mapping: VivoH5ADImportPlan, featureNamespace: String, perturbationID: String, hours: [Double]) {
        schemaVersion = 1; self.mapping = mapping; self.featureNamespace = featureNamespace; self.perturbationID = perturbationID; self.hours = hours
    }
    private enum CodingKeys: String, CodingKey { case schemaVersion, mapping, featureNamespace, perturbationID, hours }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["schemaVersion", "mapping", "featureNamespace", "perturbationID", "hours"])
        let c = try decoder.container(keyedBy: CodingKeys.self)
        schemaVersion = try c.decode(Int.self, forKey: .schemaVersion); mapping = try c.decode(VivoH5ADImportPlan.self, forKey: .mapping)
        featureNamespace = try c.decode(String.self, forKey: .featureNamespace); perturbationID = try c.decode(String.self, forKey: .perturbationID); hours = try c.decode([Double].self, forKey: .hours)
    }
    func validate() throws {
        guard schemaVersion == 1, vivoOmicsID(featureNamespace), vivoOmicsID(perturbationID), (1...64).contains(hours.count),
              hours.allSatisfy({ $0.isFinite && $0 >= 0 && $0 <= 87_600 }), Set(hours).count == hours.count else { throw VivoOmicsError.invalid("duration query identity or hours") }
    }
}

public struct VivoDurationCurve: Codable, Sendable, Equatable {
    public let donorID: String
    public let hours: [Double]
    public let responses: [[Double]]
}
public struct VivoDurationModel: Codable, Sendable, Equatable {
    public let method: String
    public let plan: VivoDurationPlan
    public let source: VivoFingerprint
    public let featureIDs: [String]
    public let organism: String
    public let cellGroup: String?
    public let curves: [VivoDurationCurve]
    public let selectedFeatureIndices: [Int]
    public let contextCenters: [Double]
    public let contextScales: [Double]
    public let contexts: [[Double]]
    public let ridgeInverse: [[Double]]
    public let maximumSolveResidual: Double
    public let maximumSupportedHours: Double
    public let qualification: String
}
public struct VivoDurationPrediction: Codable, Sendable, Equatable {
    public let hours: Double
    public let prediction: VivoPerturbationPrediction
}
public struct VivoDurationReport: Codable, Sendable, Equatable {
    public let method: String
    public let referenceSource: VivoFingerprint
    public let featureIDs: [String]
    public let predictions: [VivoDurationPrediction]
    public let qualification: String
}

/// Time enters the existing population RNA-response boundary, not a new biological simulator.
public enum VivoDurationPerturbation {
    static func interpolate(_ curve: VivoDurationCurve, hours: Double) throws -> [Double] {
        guard hours.isFinite, hours >= 0, let last = curve.hours.last, hours <= last else { throw VivoOmicsError.invalid("duration extrapolation is unsupported") }
        if let exact = curve.hours.firstIndex(of: hours) { return curve.responses[exact] }
        guard let upper = curve.hours.firstIndex(where: { $0 > hours }), upper > 0 else { throw VivoOmicsError.invalid("duration curve bracket") }
        let lower = upper - 1, fraction = (log1p(hours) - log1p(curve.hours[lower])) / (log1p(curve.hours[upper]) - log1p(curve.hours[lower]))
        return zip(curve.responses[lower], curve.responses[upper]).map { $0 + fraction * ($1 - $0) }
    }
    static func model(_ bulk: VivoPseudobulkCounts, plan: VivoDurationPlan, source: VivoFingerprint) throws -> VivoDurationModel {
        try plan.validate()
        let groups = bulk.groups, features = plan.responseFeatureIDs ?? bulk.featureIDs, m = features.count
        guard !features.isEmpty, Set(bulk.featureIDs).count == bulk.featureIDs.count, Set(features).isSubset(of: Set(bulk.featureIDs)), bulk.countUnit == plan.mapping.countUnit,
              Set(groups.map(\.organism)).count == 1, Set(groups.map(\.cellGroup)).count == 1, groups.allSatisfy({ $0.donorID != nil }) else { throw VivoOmicsError.invalid("duration training axes, units or population") }
        let donors = Set(groups.compactMap(\.donorID)).sorted(), n = donors.count
        guard (2...64).contains(n), n <= 100_000_000 / m / n else { throw VivoOmicsError.limit("duration fit work") }
        let durations = Dictionary(uniqueKeysWithValues: plan.exposures.map { ($0.condition, $0.hours) })
        guard groups.allSatisfy({ $0.condition == plan.controlCondition || durations[$0.condition] != nil }),
              Set(plan.exposures.map(\.condition)).isSubset(of: Set(groups.map(\.condition))) else { throw VivoOmicsError.invalid("duration training condition coverage") }
        let (counts, logs, _) = try VivoPerturbation.dense(bulk)
        let lookup = Dictionary(uniqueKeysWithValues: bulk.featureIDs.enumerated().map { ($0.element, $0.offset) }), order = features.map { lookup[$0]! }
        var controls: [Int] = [], donorRows: [[Int]] = [], curves: [VivoDurationCurve] = []
        for donor in donors {
            let rows = groups.indices.filter { groups[$0].donorID == donor }, control = rows.filter { groups[$0].condition == plan.controlCondition }
            let treated = rows.filter { groups[$0].condition != plan.controlCondition }.sorted { durations[groups[$0].condition]! < durations[groups[$1].condition]! }
            guard control.count == 1, treated.count >= 2, Set(treated.map { groups[$0].condition }).count == treated.count,
                  rows.allSatisfy({ groups[$0].biologicalReplicateID == groups[control[0]].biologicalReplicateID }) else { throw VivoOmicsError.invalid("duration donor needs one control and at least two unique matched exposures") }
            controls.append(control[0]); donorRows.append(rows)
            let values = [[Double](repeating: 0, count: m)] + treated.map { row in order.map { logs[row][$0] - logs[control[0]][$0] } }
            curves.append(.init(donorID: donor, hours: [0] + treated.map { durations[groups[$0].condition]! }, responses: values))
        }
        var selected: [Int] = []
        for (j, column) in order.enumerated() {
            if j % 256 == 0 { try Task.checkCancellation() }
            var total: UInt64 = 0, expressedDonors = 0
            for rows in donorRows {
                var donorTotal: UInt64 = 0
                for row in rows { donorTotal = try vivoOmicsSum(donorTotal, counts[row][column]) }
                total = try vivoOmicsSum(total, donorTotal); if donorTotal > 0 { expressedDonors += 1 }
            }
            if total >= 10 && expressedDonors >= 2 { selected.append(j) }
        }
        guard !selected.isEmpty else { throw VivoOmicsError.invalid("duration context has no training-expressed genes") }
        var centers = [Double](repeating: 0, count: selected.count), scales = centers, contexts = Array(repeating: centers, count: n)
        for (k, j) in selected.enumerated() {
            let column = order[j]
            for row in controls { centers[k] += logs[row][column] }; centers[k] /= Double(n)
            for row in controls { scales[k] += pow(logs[row][column] - centers[k], 2) }; scales[k] = sqrt(scales[k] / Double(n))
            if controls.allSatisfy({ logs[$0][column] == logs[controls[0]][column] }) { centers[k] = logs[controls[0]][column]; scales[k] = 1 }
            else if scales[k] == 0 { scales[k] = 1 }
            for i in 0..<n { contexts[i][k] = (logs[controls[i]][column] - centers[k]) / scales[k] / sqrt(Double(selected.count)) }
        }
        var system = [Double](repeating: 0, count: n*n), inverse = Array(repeating: [Double](repeating: 0, count: n), count: n)
        for i in 0..<n { for j in 0..<n { system[i*n+j] = VivoSingleCellReduction.dot(contexts[i], contexts[j]) + (i == j ? 1 : 0) } }
        let eigen = try VivoSingleCellReduction.symmetricEigen(system, n: n)
        guard eigen.values.allSatisfy({ $0.isFinite && $0 > 0 }) else { throw VivoOmicsError.invalid("duration ridge system") }
        for i in 0..<n { for j in 0..<n { for k in 0..<n { inverse[i][j] += eigen.vectors[i*n+k] * eigen.vectors[j*n+k] / eigen.values[k] } } }
        var residual = 0.0
        for i in 0..<n { for j in 0..<n { var value = 0.0; for k in 0..<n { value += system[i*n+k] * inverse[k][j] }; residual = max(residual, abs(value - (i == j ? 1 : 0))) } }
        guard residual.isFinite, residual <= 1e-9 else { throw VivoOmicsError.invalid("duration ridge inverse residual") }
        return .init(method: "donor-curves-log1p-hours-linear-context-ridge-alpha1-v1", plan: plan, source: source, featureIDs: features,
                     organism: groups[0].organism, cellGroup: groups[0].cellGroup, curves: curves, selectedFeatureIndices: selected,
                     contextCenters: centers, contextScales: scales, contexts: contexts, ridgeInverse: inverse, maximumSolveResidual: residual,
                     maximumSupportedHours: curves.map { $0.hours.last! }.min()!,
                     qualification: "Training-donor empirical RNA curves; zero response at observed control, interpolation linear in log1p(hours). No extrapolation, kinetic mechanism, dose model or causal qualification. Unequal observation schedules remain explicit.")
    }
    static func evaluate(_ bulk: VivoPseudobulkCounts, plan: VivoDurationQueryPlan, model: VivoDurationModel) throws -> VivoDurationReport {
        try plan.validate()
        guard plan.featureNamespace == model.plan.featureNamespace, plan.perturbationID == model.plan.perturbationID,
              plan.mapping.countUnit == model.plan.mapping.countUnit, bulk.countUnit == model.plan.mapping.countUnit,
              plan.mapping.samples.allSatisfy({ $0.condition == model.plan.controlCondition }),
              plan.hours.allSatisfy({ $0 <= model.maximumSupportedHours }), Set(bulk.featureIDs).count == bulk.featureIDs.count,
              Set(model.featureIDs).isSubset(of: Set(bulk.featureIDs)), model.plan.responseFeatureIDs != nil || Set(model.featureIDs) == Set(bulk.featureIDs),
              bulk.groups.allSatisfy({ $0.donorID != nil && $0.condition == model.plan.controlCondition && $0.organism == model.organism && $0.cellGroup == model.cellGroup }),
              Set(bulk.groups.compactMap(\.donorID)).count == bulk.groups.count,
              Set(bulk.groups.compactMap(\.donorID)).isDisjoint(with: Set(model.curves.map(\.donorID))) else { throw VivoOmicsError.invalid("duration query control, identity, donor, panel or supported time mismatch") }
        let m = model.featureIDs.count, n = model.curves.count
        guard bulk.groups.count <= 2_000_000 / m / plan.hours.count,
              bulk.groups.count <= 100_000_000 / m / n / plan.hours.count else { throw VivoOmicsError.limit("duration prediction work or output budget") }
        let (_, logs, totals) = try VivoPerturbation.dense(bulk), lookup = Dictionary(uniqueKeysWithValues: bulk.featureIDs.enumerated().map { ($0.element, $0.offset) })
        let order = model.featureIDs.map { lookup[$0]! }, critical = try VivoOmicsLinearStatistics.studentCriticalValue(degreesOfFreedom: Double(n-1), coverage: 0.95)
        var results: [VivoDurationPrediction] = []
        for row in bulk.groups.indices {
            let control = order.map { logs[row][$0] }, selected = model.selectedFeatureIndices
            let context = selected.indices.map { (control[selected[$0]] - model.contextCenters[$0]) / model.contextScales[$0] / sqrt(Double(selected.count)) }
            let similarities = model.contexts.map { VivoSingleCellReduction.dot(context, $0) }
            let weights = (0..<n).map { j in (0..<n).reduce(0.0) { $0 + similarities[$1] * model.ridgeInverse[$1][j] } }
            for hours in plan.hours {
                try Task.checkCancellation()
                let responses = try model.curves.map { try interpolate($0, hours: hours) }
                var mean = [Double](repeating: 0, count: m), median = mean, ridge = mean
                var lower = [Double?](repeating: nil, count: m), upper = lower, treatedLower = lower, treatedUpper = lower, unavailable: [Int] = []
                for gene in 0..<m {
                    if gene % 256 == 0 { try Task.checkCancellation() }
                    let values = responses.map { $0[gene] }, sorted = values.sorted()
                    mean[gene] = values.reduce(0,+) / Double(n); median[gene] = n % 2 == 0 ? (sorted[n/2-1]+sorted[n/2])/2 : sorted[n/2]
                    ridge[gene] = mean[gene] + (0..<n).reduce(0.0) { $0 + weights[$1] * (values[$1] - mean[gene]) }
                    let variance = values.reduce(0.0) { $0 + pow($1 - mean[gene], 2) } / Double(n-1)
                    if values.allSatisfy({ $0 == values[0] }) || variance == 0 { unavailable.append(gene); continue }
                    guard variance.isFinite, variance > 0 else { throw VivoOmicsError.invalid("duration response variance") }
                    let half = critical * sqrt(variance * (1 + 1 / Double(n)))
                    lower[gene] = mean[gene] - half; upper[gene] = mean[gene] + half
                    treatedLower[gene] = max(0, control[gene] + lower[gene]!); treatedUpper[gene] = max(0, control[gene] + upper[gene]!)
                }
                var estimates: [VivoPerturbationEstimate] = []
                for (name, response) in [("noChange", [Double](repeating: 0, count: m)), ("durationMean", mean), ("durationMedian", median), ("durationContextRidge", ridge)] {
                    let treated = zip(control, response).map { max(0, $0 + $1) }, applied = zip(treated, control).map(-), implied = treated.reduce(0.0) { $0 + expm1($1) }
                    guard response.allSatisfy(\.isFinite), treated.allSatisfy(\.isFinite), implied.isFinite else { throw VivoOmicsError.invalid("nonfinite duration prediction") }
                    estimates.append(.init(baseline: name, unclippedResponse: response, predictedTreated: treated, predictedResponse: applied, impliedCPMSum: implied))
                }
                let interval = VivoDonorResponsePredictiveInterval(method: "normal-independent-interpolated-donor-curve-predictive-t-v1", nominalCoverage: 0.95, trainingDonors: n, degreesOfFreedom: n-1, studentCriticalValue: critical, unclippedResponseLower: lower, unclippedResponseUpper: upper, predictedTreatedLower: treatedLower, predictedTreatedUpper: treatedUpper, unavailableFeatureIndices: unavailable)
                results.append(.init(hours: hours, prediction: .init(group: bulk.groups[row], libraryCounts: totals[row], control: control, estimates: estimates, meanResponsePredictiveInterval: interval)))
            }
        }
        return .init(method: model.method, referenceSource: model.source, featureIDs: model.featureIDs, predictions: results,
                     qualification: "Control-only held-out donors. Equal donor weighting after within-donor log-time interpolation; ridge uses training-control contexts only. Full-source RNA normalization precedes panel projection. Nominal95percent pointwise future-donor intervals assume independent normal interpolated responses; interpolation error and control uncertainty are not modeled. No calibrated, simultaneous, mechanistic, dose-transfer or phenotype claim.")
    }
}
