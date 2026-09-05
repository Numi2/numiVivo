import Foundation

public struct VivoAmberRestart: Sendable, Equatable {
    public let title: String
    public let atomCount: UInt32
    public let timePS: Double?
    public let positionsNM: [VivoVector3D]
    /// Preserved raw AMBER restart velocity numbers. Wave B will own the exact
    /// conversion/execution semantics; structure import never silently treats them
    /// as nm/ps or discards their existence.
    public let rawVelocityValues: [Double]?
    public let periodicCell: VivoPeriodicCell?

    public init(data: Data) throws {
        guard let text = String(data: data, encoding: .utf8) else {
            throw VivoArtifactValidationError.invalid("AMBER restart is not UTF-8/ASCII text")
        }
        let lines = text.split(separator: "\n", omittingEmptySubsequences: false).map(String.init)
        guard lines.count >= 2 else { throw VivoArtifactValidationError.invalid("AMBER restart is truncated") }
        title = lines[0].trimmingCharacters(in: .whitespacesAndNewlines)
        let header = lines[1].split(whereSeparator: \.isWhitespace)
        guard let first = header.first, let count = UInt32(first), count > 0 else {
            throw VivoArtifactValidationError.invalid("AMBER restart has invalid atom count")
        }
        atomCount = count
        if header.count > 1 {
            guard let value = Double(header[1]), value.isFinite else {
                throw VivoArtifactValidationError.invalid("AMBER restart has invalid time")
            }
            timePS = value
        } else { timePS = nil }

        var values: [Double] = []
        for line in lines.dropFirst(2) {
            // rst7 uses fixed width 12 for coordinates, but accepting whitespace
            // here also supports files written with expanded exponent fields.
            let tokens = line.split(whereSeparator: \.isWhitespace)
            for token in tokens {
                let normalized = token.replacingOccurrences(of: "D", with: "E").replacingOccurrences(of: "d", with: "e")
                guard let value = Double(normalized), value.isFinite else {
                    throw VivoArtifactValidationError.invalid("AMBER restart contains invalid numeric field")
                }
                values.append(value)
            }
        }
        let coordinateCount = Int(count) * 3
        guard values.count >= coordinateCount else {
            throw VivoArtifactValidationError.invalid("AMBER restart has fewer than 3*NATOM coordinates")
        }
        var positions: [VivoVector3D] = []
        positions.reserveCapacity(Int(count))
        for atom in 0..<Int(count) {
            let base = atom * 3
            positions.append(.init(values[base] * 0.1, values[base + 1] * 0.1, values[base + 2] * 0.1))
        }
        positionsNM = positions

        let remainder = Array(values.dropFirst(coordinateCount))
        if remainder.isEmpty {
            rawVelocityValues = nil; periodicCell = nil
        } else if remainder.count == 6 {
            rawVelocityValues = nil
            periodicCell = try Self.cell(remainder)
        } else if remainder.count == coordinateCount {
            rawVelocityValues = remainder; periodicCell = nil
        } else if remainder.count == coordinateCount + 6 {
            rawVelocityValues = Array(remainder.prefix(coordinateCount))
            periodicCell = try Self.cell(Array(remainder.suffix(6)))
        } else {
            throw VivoArtifactValidationError.incompatible(
                "AMBER restart trailing payload is neither velocities nor six-value box data"
            )
        }
    }

    private static func cell(_ values: [Double]) throws -> VivoPeriodicCell {
        guard values.count == 6 else { throw VivoArtifactValidationError.invalid("AMBER restart box needs six values") }
        return try VivoPDB.cellFromParameters(aA: values[0], bA: values[1], cA: values[2],
                                              alpha: values[3], beta: values[4], gamma: values[5])
    }
}
