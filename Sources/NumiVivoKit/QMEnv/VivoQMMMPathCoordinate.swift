import Foundation

public struct VivoQMMMPathNode: Codable, Sendable, Equatable {
    public var identifier: String
    /// One value per component coordinate, in nm.
    public var valuesNM: [Double]
    public init(identifier: String, valuesNM: [Double]) {
        self.identifier = identifier
        self.valuesNM = valuesNM
    }
}

/// A smooth path collective variable built from the existing mapped scalar QM/MM
/// coordinates. Dimension scales make the squared path metric dimensionless.
public struct VivoQMMMPathCoordinate: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/qmmm-path-coordinate/v1"
    public var schema: String
    public var identifier: String
    public var components: [VivoQMMMReactionCoordinate]
    public var scalesNM: [Double]
    public var nodes: [VivoQMMMPathNode]
    /// Positive smoothing strength multiplying dimensionless squared distance.
    public var lambda: Double

    public init(identifier: String,
                components: [VivoQMMMReactionCoordinate],
                scalesNM: [Double],
                nodes: [VivoQMMMPathNode],
                lambda: Double) {
        schema = Self.schema
        self.identifier = identifier
        self.components = components
        self.scalesNM = scalesNM
        self.nodes = nodes
        self.lambda = lambda
    }

    public func validate() throws {
        guard schema == Self.schema,
              !identifier.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty,
              identifier.utf8.count <= 512,
              components.count >= 2, components.count <= 64,
              Set(components.map(\.identifier)).count == components.count,
              scalesNM.count == components.count,
              scalesNM.allSatisfy({ $0.isFinite && $0 > 1e-12 && $0 <= 1000 }),
              nodes.count >= 2, nodes.count <= 100_000,
              Set(nodes.map(\.identifier)).count == nodes.count,
              lambda.isFinite, lambda > 0, lambda <= 1e8 else {
            throw VivoChemistryError.invalid("QM/MM PathCV identity, dimensions, scales, nodes or smoothing")
        }
        for component in components { try component.validate() }
        for node in nodes {
            guard !node.identifier.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty,
                  node.valuesNM.count == components.count,
                  node.valuesNM.allSatisfy(\.isFinite) else {
                throw VivoChemistryError.invalid("QM/MM PathCV node")
            }
        }
        for pair in zip(nodes, nodes.dropFirst()) {
            let distance = zip(pair.0.valuesNM, pair.1.valuesNM).enumerated().reduce(0.0) { partial, item in
                let d = (item.element.1 - pair.0.valuesNM[item.offset]) / scalesNM[item.offset]
                return partial + d * d
            }
            guard distance.isFinite, distance > 1e-16 else {
                throw VivoChemistryError.invalid("QM/MM PathCV adjacent nodes must be distinct")
            }
        }
    }
}

public struct VivoQMMMPathCoordinateEvaluation: Sendable, Equatable {
    /// Smooth progress from 0 at the first node to 1 at the last node.
    public let progress: Double
    /// Smooth dimensionless distance-like free coordinate: -log(sum w_i)/lambda.
    public let distanceFromPath: Double
    public let componentValuesNM: [Double]
    public let progressGradients: [UInt32: VivoVector3D]
    public let distanceGradients: [UInt32: VivoVector3D]
    public let normalizedNodeWeights: [Double]
}

public struct VivoQMMMResolvedPathCoordinate: Sendable {
    public let source: VivoQMMMPathCoordinate
    public let components: [VivoQMMMResolvedCoordinate]

    public init(source: VivoQMMMPathCoordinate, system: VivoClassicalSystem) throws {
        try source.validate()
        self.source = source
        self.components = try source.components.map { try VivoQMMMResolvedCoordinate(source: $0, system: system) }
    }

    public func evaluate(_ geometry: VivoMDCandidateGeometry) throws -> VivoQMMMPathCoordinateEvaluation {
        let scalar = try components.map { try $0.evaluate(geometry) }
        let q = scalar.map(\.valueNM)
        let dimensions = q.count
        let nodeCount = source.nodes.count

        var squaredDistances = [Double](repeating: 0, count: nodeCount)
        for i in 0..<nodeCount {
            var value = 0.0
            for d in 0..<dimensions {
                let delta = (q[d] - source.nodes[i].valuesNM[d]) / source.scalesNM[d]
                value += delta * delta
            }
            guard value.isFinite else { throw VivoChemistryError.convergence("PathCV distance overflow") }
            squaredDistances[i] = value
        }
        // Stable softmax over -lambda*D_i.
        let logWeights = squaredDistances.map { -source.lambda * $0 }
        guard let maximum = logWeights.max(), maximum.isFinite else {
            throw VivoChemistryError.convergence("PathCV weight overflow")
        }
        let shifted = logWeights.map { exp($0 - maximum) }
        let sum = shifted.reduce(0, +)
        guard sum.isFinite, sum > 0 else { throw VivoChemistryError.convergence("PathCV zero weight") }
        let weights = shifted.map { $0 / sum }
        let logWeightSum = maximum + log(sum)
        let progress = weights.enumerated().reduce(0.0) { partial, item in
            let t = Double(item.offset) / Double(nodeCount - 1)
            return partial + t * item.element
        }
        let z = -logWeightSum / source.lambda
        guard progress.isFinite, z.isFinite else { throw VivoChemistryError.convergence("PathCV value overflow") }

        var dsDq = [Double](repeating: 0, count: dimensions)
        var dzDq = [Double](repeating: 0, count: dimensions)
        for d in 0..<dimensions {
            let inverseScale2 = 1 / (source.scalesNM[d] * source.scalesNM[d])
            var sDerivative = 0.0
            var zDerivative = 0.0
            for i in 0..<nodeCount {
                let delta = q[d] - source.nodes[i].valuesNM[d]
                let dDistanceDq = 2 * delta * inverseScale2
                let dLogWeightDq = -source.lambda * dDistanceDq
                let t = Double(i) / Double(nodeCount - 1)
                sDerivative += weights[i] * (t - progress) * dLogWeightDq
                zDerivative += weights[i] * dDistanceDq
            }
            dsDq[d] = sDerivative
            dzDq[d] = zDerivative
        }

        var progressGradients: [UInt32: VivoVector3D] = [:]
        var distanceGradients: [UInt32: VivoVector3D] = [:]
        func add(_ map: inout [UInt32: VivoVector3D], particle: UInt32, value: VivoVector3D) {
            map[particle] = (map[particle] ?? .zero) + value
        }
        for d in 0..<dimensions {
            for (particle, gradient) in scalar[d].gradients {
                add(&progressGradients, particle: particle, value: gradient * dsDq[d])
                add(&distanceGradients, particle: particle, value: gradient * dzDq[d])
            }
        }
        guard progressGradients.values.allSatisfy(\.isFinite),
              distanceGradients.values.allSatisfy(\.isFinite) else {
            throw VivoChemistryError.convergence("PathCV gradient overflow")
        }
        return .init(progress: progress,
                     distanceFromPath: z,
                     componentValuesNM: q,
                     progressGradients: progressGradients,
                     distanceGradients: distanceGradients,
                     normalizedNodeWeights: weights)
    }
}
