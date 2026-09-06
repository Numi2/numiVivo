import Foundation

/// Portable, user-authored composition of existing native operations. A recipe
/// is data: operation identifiers select an allowlisted registry, never a shell.
public struct VivoWorkflowRecipe: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/workflow-recipe/v1"
    public var schema: String
    public var identifier: String
    public var artifacts: [VivoWorkflowArtifact]
    public var nodes: [VivoWorkflowRecipeNode]
    public var outputs: [VivoWorkflowExport]
    public var policy: VivoWorkflowSchedulingPolicy
    public init(identifier: String, artifacts: [VivoWorkflowArtifact], nodes: [VivoWorkflowRecipeNode],
                outputs: [VivoWorkflowExport], policy: VivoWorkflowSchedulingPolicy = .init()) {
        schema = Self.schema; self.identifier = identifier; self.artifacts = artifacts
        self.nodes = nodes; self.outputs = outputs; self.policy = policy
    }
    /// Ordering independent identity. Array order within scientific payloads,
    /// including atom order, weights and time series, is never canonicalized away.
    public func fingerprint() throws -> VivoFingerprint {
        var copy = self
        copy.artifacts.sort { $0.identifier < $1.identifier }
        copy.nodes.sort { $0.identifier < $1.identifier }
        copy.outputs.sort { $0.name < $1.name }
        return try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(copy))
    }
}
public struct VivoWorkflowArtifact: Codable, Sendable, Equatable {
    public var identifier: String
    public var source: VivoWorkflowArtifactSource
    public init(identifier: String, source: VivoWorkflowArtifactSource) {
        self.identifier = identifier; self.source = source
    }
}
public enum VivoWorkflowArtifactSource: Codable, Sendable, Equatable {
    case json(kind: String, payload: VivoJSONValue)
    case stored(kind: String, fingerprint: VivoFingerprint)
    public var kind: String {
        switch self { case .json(let kind, _), .stored(let kind, _): return kind }
    }
}
public enum VivoWorkflowBinding: Codable, Sendable, Equatable {
    case artifact(identifier: String)
    case output(node: String, port: String)
}
public struct VivoWorkflowRecipeNode: Codable, Sendable, Equatable {
    public var identifier: String
    public var operation: String
    public var version: String
    public var inputs: [String: VivoWorkflowBinding]
    public var configuration: VivoJSONValue
    public var resources: VivoChemistryResourceContract
    public init(identifier: String, operation: String, version: String = "1",
                inputs: [String: VivoWorkflowBinding], configuration: VivoJSONValue = .object([:]),
                resources: VivoChemistryResourceContract = .init()) {
        self.identifier = identifier; self.operation = operation; self.version = version
        self.inputs = inputs; self.configuration = configuration; self.resources = resources
    }
    /// Conservative admission reservation, not a measured process/GPU RSS. A
    /// solver must still enforce its own numerical and allocation contracts.
    public func reservationBytes() throws -> Int {
        try resources.budget.validate()
        guard resources.maximumInputBytes > 0, resources.maximumOutputBytes > 0 else {
            throw VivoChemistryError.invalid("workflow node input/output byte limits must be positive")
        }
        let (first, a) = resources.budget.maximumBytes.addingReportingOverflow(resources.maximumInputBytes)
        let (total, b) = first.addingReportingOverflow(resources.maximumOutputBytes)
        guard !a, !b else { throw VivoChemistryError.resourceLimit("workflow memory reservation overflow") }
        return total
    }
}
public struct VivoWorkflowExport: Codable, Sendable, Equatable {
    public var name: String
    public var node: String
    public var port: String
    public init(name: String, node: String, port: String) { self.name = name; self.node = node; self.port = port }
}
public struct VivoWorkflowSchedulingPolicy: Codable, Sendable, Equatable {
    public var maximumConcurrentTasks: Int
    public var maximumConcurrentMetalTasks: Int
    public var maximumReservedBytes: Int
    public var maximumNodes: Int
    public var maximumInlineBytes: Int
    public init(maximumConcurrentTasks: Int = 2, maximumConcurrentMetalTasks: Int = 1,
                maximumReservedBytes: Int = 2 * 1024 * 1024 * 1024, maximumNodes: Int = 256,
                maximumInlineBytes: Int = 128 * 1024 * 1024) {
        self.maximumConcurrentTasks = maximumConcurrentTasks; self.maximumConcurrentMetalTasks = maximumConcurrentMetalTasks
        self.maximumReservedBytes = maximumReservedBytes; self.maximumNodes = maximumNodes
        self.maximumInlineBytes = maximumInlineBytes
    }
    public func validate() throws {
        guard (1...64).contains(maximumConcurrentTasks), (0...maximumConcurrentTasks).contains(maximumConcurrentMetalTasks),
              (1...4096).contains(maximumNodes), maximumReservedBytes > 0, maximumInlineBytes > 0,
              maximumInlineBytes <= maximumReservedBytes else {
            throw VivoChemistryError.invalid("workflow scheduling bounds")
        }
    }
}
public struct VivoWorkflowOperationDescription: Codable, Sendable, Equatable {
    public let identifier: String
    public let version: String
    public let numericalBackend: String
    public let inputs: [VivoChemistryTaskOutput]
    public let outputs: [VivoChemistryTaskOutput]
    public let summary: String
    public let validationScope: String
}
public struct VivoWorkflowDefinition: Sendable {
    public let operation: VivoChemistryOperation
    public let inputKinds: [String: String]
    public let summary: String
    public let validationScope: String
    public let validateConfiguration: @Sendable (VivoJSONValue) throws -> Void
    public init(operation: VivoChemistryOperation, inputKinds: [String: String], summary: String,
                validationScope: String, validateConfiguration: @escaping @Sendable (VivoJSONValue) throws -> Void) {
        self.operation = operation; self.inputKinds = inputKinds; self.summary = summary
        self.validationScope = validationScope; self.validateConfiguration = validateConfiguration
    }
    public var description: VivoWorkflowOperationDescription {
        .init(identifier: operation.identifier, version: operation.version, numericalBackend: operation.numericalBackend,
              inputs: inputKinds.keys.sorted().map { .init(name: $0, kind: inputKinds[$0]!) },
              outputs: operation.outputs, summary: summary, validationScope: validationScope)
    }
}
public struct VivoWorkflowRegistry: Sendable {
    public let implementationFingerprint: VivoFingerprint
    private let definitions: [String: VivoWorkflowDefinition]
    public init(implementationFingerprint: VivoFingerprint, definitions: [VivoWorkflowDefinition]) throws {
        guard !definitions.isEmpty, definitions.count <= 4096 else { throw VivoChemistryError.invalid("empty or oversized operation registry") }
        var values: [String: VivoWorkflowDefinition] = [:]
        for definition in definitions {
            let operation = definition.operation
            guard values[operation.identifier] == nil, operation.implementationFingerprint == implementationFingerprint,
                  VivoWorkflowPlanner.validName(operation.identifier), VivoWorkflowPlanner.validName(operation.version),
                  !operation.numericalBackend.isEmpty, !operation.outputs.isEmpty,
                  Set(operation.outputs.map(\.name)).count == operation.outputs.count,
                  operation.outputs.allSatisfy({ VivoWorkflowPlanner.validName($0.name) && VivoWorkflowPlanner.validName($0.kind) }),
                  definition.inputKinds.allSatisfy({ VivoWorkflowPlanner.validName($0.key) && VivoWorkflowPlanner.validName($0.value) }) else {
                throw VivoChemistryError.invalid("duplicate, malformed or differently built workflow definition")
            }
            values[operation.identifier] = definition
        }
        self.implementationFingerprint = implementationFingerprint; self.definitions = values
    }
    public var catalog: [VivoWorkflowOperationDescription] {
        definitions.keys.sorted().map { definitions[$0]!.description }
    }
    public func definition(_ identifier: String) throws -> VivoWorkflowDefinition {
        guard let value = definitions[identifier] else { throw VivoChemistryError.unsupported("unregistered native operation: \(identifier)") }
        return value
    }
}
public struct VivoWorkflowPlannedNode: Codable, Sendable, Equatable {
    public let identifier: String
    public let operation: String
    public let dependencies: [String]
    public let reservedBytes: Int
    public let numericalBackend: String
}
public struct VivoWorkflowPlan: Codable, Sendable, Equatable {
    public let schema: String
    public let recipeFingerprint: VivoFingerprint
    public let implementationFingerprint: VivoFingerprint
    public let nodes: [VivoWorkflowPlannedNode]
    public let topologicalOrder: [String]
    public let inlineBytes: Int
    public let declaredMaximumConcurrentTasks: Int
    public let declaredMaximumReservedBytes: Int
    public let scope: String
}
public enum VivoWorkflowPlanner {
    static func validName(_ name: String) -> Bool {
        !name.isEmpty && name.utf8.count <= 240 && name.unicodeScalars.allSatisfy {
            (48...57).contains($0.value) || (65...90).contains($0.value) || (97...122).contains($0.value)
            || [45,46,95].contains($0.value)
        }
    }
    /// Complete graph/type/configuration/resource preflight. No artifact-store
    /// writes or numerical operation executions take place during planning.
    public static func compile(_ recipe: VivoWorkflowRecipe, registry: VivoWorkflowRegistry) throws -> VivoWorkflowPlan {
        try recipe.policy.validate()
        guard recipe.schema == VivoWorkflowRecipe.schema, validName(recipe.identifier),
              !recipe.nodes.isEmpty, recipe.nodes.count <= recipe.policy.maximumNodes,
              recipe.artifacts.count <= 8192, !recipe.outputs.isEmpty, recipe.outputs.count <= 8192,
              Set(recipe.nodes.map(\.identifier)).count == recipe.nodes.count,
              Set(recipe.artifacts.map(\.identifier)).count == recipe.artifacts.count,
              Set(recipe.outputs.map(\.name)).count == recipe.outputs.count else {
            throw VivoChemistryError.invalid("workflow schema, identity, dimensions or duplicate entries")
        }
        let nodes = Dictionary(uniqueKeysWithValues: recipe.nodes.map { ($0.identifier, $0) })
        var artifactKinds: [String: String] = [:], inlineBytes = 0
        for artifact in recipe.artifacts {
            guard validName(artifact.identifier), validName(artifact.source.kind) else { throw VivoChemistryError.invalid("workflow artifact identity or kind") }
            artifactKinds[artifact.identifier] = artifact.source.kind
            if case .json(_, let payload) = artifact.source {
                let size = try VivoCanonicalJSON.encode(payload).count
                guard size <= recipe.policy.maximumInlineBytes - inlineBytes else { throw VivoChemistryError.resourceLimit("workflow inline input byte budget") }
                inlineBytes += size
            }
        }
        var dependencies: [String: Set<String>] = [:], planned: [VivoWorkflowPlannedNode] = []
        for node in recipe.nodes.sorted(by: { $0.identifier < $1.identifier }) {
            let definition = try registry.definition(node.operation), operation = definition.operation
            let reservation = try node.reservationBytes()
            guard validName(node.identifier), node.version == operation.version,
                  node.resources.numericalBackend == operation.numericalBackend,
                  Set(node.inputs.keys) == Set(definition.inputKinds.keys),
                  reservation <= recipe.policy.maximumReservedBytes,
                  !operation.numericalBackend.hasPrefix("metal") || recipe.policy.maximumConcurrentMetalTasks > 0 else {
                throw VivoChemistryError.invalid("workflow node \(node.identifier) version, ports, backend or admission reservation")
            }
            try definition.validateConfiguration(node.configuration)
            guard try VivoCanonicalJSON.encode(node.configuration).count <= node.resources.maximumInputBytes else {
                throw VivoChemistryError.resourceLimit("workflow configuration byte budget")
            }
            var deps = Set<String>()
            for (port, binding) in node.inputs {
                let actual: String
                switch binding {
                case .artifact(let identifier):
                    guard let kind = artifactKinds[identifier] else { throw VivoChemistryError.invalid("workflow input references absent artifact \(identifier)") }
                    actual = kind
                case .output(let identifier, let output):
                    guard identifier != node.identifier, let producer = nodes[identifier] else { throw VivoChemistryError.invalid("workflow unknown/self edge") }
                    let upstream = try registry.definition(producer.operation).operation
                    guard let value = upstream.outputs.first(where: { $0.name == output }) else { throw VivoChemistryError.invalid("workflow absent producer port \(identifier).\(output)") }
                    actual = value.kind; deps.insert(identifier)
                }
                guard actual == definition.inputKinds[port] else { throw VivoChemistryError.invalid("workflow incompatible input kind at \(node.identifier).\(port)") }
            }
            dependencies[node.identifier] = deps
            planned.append(.init(identifier: node.identifier, operation: node.operation, dependencies: deps.sorted(),
                                 reservedBytes: reservation, numericalBackend: operation.numericalBackend))
        }
        var reached = Set<String>(), order: [String] = []
        while order.count < nodes.count {
            let ready = nodes.keys.filter { !reached.contains($0) && dependencies[$0]!.isSubset(of: reached) }.sorted()
            guard !ready.isEmpty else { throw VivoChemistryError.invalid("cyclic workflow; no tasks executed") }
            order.append(contentsOf: ready); reached.formUnion(ready)
        }
        for output in recipe.outputs {
            guard validName(output.name), let node = nodes[output.node],
                  try registry.definition(node.operation).operation.outputs.contains(where: { $0.name == output.port }) else {
                throw VivoChemistryError.invalid("workflow export identity or producer port")
            }
        }
        return .init(schema: "numivivo.org/workflow-plan/v1", recipeFingerprint: try recipe.fingerprint(),
            implementationFingerprint: registry.implementationFingerprint, nodes: planned, topologicalOrder: order,
            inlineBytes: inlineBytes, declaredMaximumConcurrentTasks: recipe.policy.maximumConcurrentTasks,
            declaredMaximumReservedBytes: recipe.policy.maximumReservedBytes,
            scope: "static graph, typed ports, registered operations, configuration and admission checks; external artifact bytes and hardware availability checked at execution; not scientific qualification")
    }
}
