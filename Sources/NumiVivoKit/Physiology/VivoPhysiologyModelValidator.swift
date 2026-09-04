import Foundation

/// Checks the executable prepared tables, not just their declared dimensions.
/// Codable initializers do not run the source compiler's validation methods.
public enum VivoPhysiologyModelValidator {
    public static func validate(_ model: PreparedVivoPhysiologyModel) throws {
        func require(_ condition: Bool, _ message: String) throws {
            guard condition else { throw VivoRuntimeError.invalidConfiguration("prepared physiology: " + message) }
        }
        let pairs = UInt64(model.analytes.count).multipliedReportingOverflow(by: UInt64(model.compartments.count))
        let elements = UInt64(model.pairCount) * UInt64(model.environmentCount)
        try require(model.schema == PreparedVivoPhysiologyModel.schema && model.environmentCount > 0 && model.pairCount > 0,
                    "invalid schema or empty shape")
        try require(!pairs.overflow && pairs.partialValue == UInt64(model.pairCount) &&
                    model.pairCount <= 1_048_576 && elements <= 268_435_456 && elements <= UInt64(Int.max),
                    "pair/product dimensions exceed executable bounds")
        try require(model.incidenceOffsets.count == Int(model.pairCount) + 1 &&
                    model.clearances.count == Int(model.pairCount) && UInt64(model.initialState.count) == elements &&
                    model.incidence.count <= 16_777_216 && model.doses.count <= 16_777_216,
                    "table shape or bounded record count mismatch")
        try require(model.minimumTimeStepSeconds.isFinite && model.preferredTimeStepSeconds.isFinite &&
                    model.maximumTimeStepSeconds.isFinite && model.linearStabilityLimitSeconds.isFinite &&
                    model.minimumTimeStepSeconds > 0 && Float(model.minimumTimeStepSeconds) > 0 &&
                    model.preferredTimeStepSeconds >= model.minimumTimeStepSeconds &&
                    model.maximumTimeStepSeconds >= model.preferredTimeStepSeconds &&
                    model.maximumTimeStepSeconds <= Double(Float.greatestFiniteMagnitude) && model.linearStabilityLimitSeconds > 0,
                    "invalid time-step or stability interval")
        var compartmentNames = Set<String>(), analyteNames = Set<String>(), doseNames = Set<String>()
        for (index, compartment) in model.compartments.enumerated() {
            try require(compartment.index == UInt32(index) && !compartment.identifier.isEmpty &&
                        compartmentNames.insert(compartment.identifier).inserted &&
                        compartment.volumeLitres.isFinite && compartment.volumeLitres > 0,
                        "invalid compartment index, identity or volume")
        }
        for (index, analyte) in model.analytes.enumerated() {
            try require(analyte.index == UInt32(index) && !analyte.identifier.isEmpty &&
                        analyteNames.insert(analyte.identifier).inserted && !analyte.concentrationUnit.isEmpty &&
                        analyte.minimum.isFinite && analyte.maximum.isFinite && analyte.minimum >= 0 &&
                        analyte.minimum <= analyte.maximum && analyte.minimum <= Double(Float.greatestFiniteMagnitude),
                        "invalid analyte index, identity or bounds")
        }
        try require(model.incidenceOffsets.first == 0 && model.incidenceOffsets.last == UInt32(model.incidence.count),
                    "CSR endpoints do not span the incidence table")
        for pair in 0..<Int(model.pairCount) {
            let begin = model.incidenceOffsets[pair], end = model.incidenceOffsets[pair + 1]
            try require(begin <= end && end <= UInt32(model.incidence.count), "CSR offsets are nonmonotone or out of bounds")
            for entry in model.incidence[Int(begin)..<Int(end)] {
                try require(entry.sourcePairIndex < model.pairCount && entry.coefficientPerSecond.isFinite &&
                            entry.flags == 0 && entry.reserved == 0, "invalid incidence source, coefficient or extension")
            }
            let clearance = model.clearances[pair]
            try require(clearance.firstOrderRate.isFinite && clearance.firstOrderRate >= 0 &&
                        clearance.maximumRate.isFinite && clearance.maximumRate >= 0 &&
                        clearance.halfSaturation.isFinite && clearance.halfSaturation >= 0 && clearance.flags & ~1 == 0,
                        "invalid clearance parameter or flags")
            if clearance.flags & 1 != 0 { try require(clearance.halfSaturation > 0, "saturable clearance requires a positive half-saturation") }
            let analyte = model.analytes[pair / model.compartments.count]
            let minimum = Float(analyte.minimum)
            let maximum = Float(min(analyte.maximum, Double(Float.greatestFiniteMagnitude)))
            let start = pair * Int(model.environmentCount)
            for value in model.initialState[start..<(start + Int(model.environmentCount))] {
                try require(value.isFinite && value >= minimum && value <= maximum, "initial state violates its FP32 analyte bounds")
            }
        }
        var lastTime = -Double.infinity
        var expandedDoses: UInt64 = 0
        for dose in model.doses {
            try require(!dose.identifier.isEmpty && doseNames.insert(dose.identifier).inserted &&
                        dose.timeSeconds.isFinite && dose.endTimeSeconds.isFinite &&
                        dose.timeSeconds >= 0 && dose.endTimeSeconds >= dose.timeSeconds && dose.timeSeconds >= lastTime &&
                        dose.pairIndex < model.pairCount && dose.value.isFinite && dose.value >= 0 &&
                        !dose.environments.isEmpty && Set(dose.environments).count == dose.environments.count &&
                        dose.environments.allSatisfy({ $0 < model.environmentCount }),
                        "invalid, unsorted or duplicated dose record")
            if dose.kind == .concentrationInfusion { try require(dose.endTimeSeconds > dose.timeSeconds, "infusion duration must be positive") }
            lastTime = dose.timeSeconds
            expandedDoses += UInt64(dose.environments.count)
            try require(expandedDoses <= 16_777_216, "expanded dose environment capacity exceeded")
        }
    }
}
