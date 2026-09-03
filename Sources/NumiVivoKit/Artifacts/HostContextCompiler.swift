import Foundation

public struct VivoContextDiagnostic: Codable, Sendable, Equatable, Hashable {
    public enum Severity: String, Codable, Sendable {
        case note
        case warning
    }

    public let severity: Severity
    public let code: String
    public let path: String
    public let message: String

    public init(severity: Severity, code: String, path: String, message: String) {
        self.severity = severity
        self.code = code
        self.path = path
        self.message = message
    }
}

public struct PreparedVivoHostContext: Sendable {
    public let contextFingerprint: VivoFingerprint
    public let programFingerprint: VivoFingerprint
    public let environmentCount: UInt32
    public let laneCount: UInt32
    public let parameterValues: [Float]
    public let transport: [VivoSpeciesTransportABI]
    public let initialCoupling: [VivoCouplingUpdate]
    public let diagnostics: [VivoContextDiagnostic]

    public init(
        contextFingerprint: VivoFingerprint,
        programFingerprint: VivoFingerprint,
        environmentCount: UInt32,
        laneCount: UInt32,
        parameterValues: [Float],
        transport: [VivoSpeciesTransportABI],
        initialCoupling: [VivoCouplingUpdate],
        diagnostics: [VivoContextDiagnostic]
    ) {
        self.contextFingerprint = contextFingerprint
        self.programFingerprint = programFingerprint
        self.environmentCount = environmentCount
        self.laneCount = laneCount
        self.parameterValues = parameterValues
        self.transport = transport
        self.initialCoupling = initialCoupling
        self.diagnostics = diagnostics
    }
}

public struct VivoHostContextCompiler: Sendable {
    public struct Options: Sendable, Codable, Equatable {
        public var strictUnknownHostChannels: Bool
        public var requireSpatialGeometryAgreement: Bool
        public var maximumGeneratedCouplingUpdates: Int

        public init(
            strictUnknownHostChannels: Bool = true,
            requireSpatialGeometryAgreement: Bool = true,
            maximumGeneratedCouplingUpdates: Int = 16_777_216
        ) {
            self.strictUnknownHostChannels = strictUnknownHostChannels
            self.requireSpatialGeometryAgreement = requireSpatialGeometryAgreement
            self.maximumGeneratedCouplingUpdates = maximumGeneratedCouplingUpdates
        }
    }

    private let units: VivoUnitSystem

    public init(units: VivoUnitSystem = .standard) {
        self.units = units
    }

    public func compile(
        _ context: VivoHostContextPack,
        for pack: VivoProgramPack,
        configuration: VivoRuntimeConfiguration,
        options: Options = .init()
    ) throws -> PreparedVivoHostContext {
        let environmentCount = Int(configuration.environmentCount)
        try context.validate(environmentCount: environmentCount)
        try configuration.validate(for: pack)
        guard options.maximumGeneratedCouplingUpdates > 0 else {
            throw VivoArtifactValidationError.invalid("maximumGeneratedCouplingUpdates must be positive")
        }
        try validateGeometry(context.tissue.geometry, configuration: configuration, options: options)

        let contextData = try VivoCanonicalJSON.encode(context)
        let contextFingerprint = try VivoCanonicalJSON.fingerprint(contextData)
        let parameters = try pack.parameterMetadata()
        let species = try pack.speciesMetadata()
        let parameterIndex = try uniqueIndex(parameters.map(\.identifier), label: "ProgramPack parameter")
        let speciesIndex = try uniqueIndex(species.map(\.identifier), label: "ProgramPack species")

        var parameterValues = try pack.parameterValues(environmentCount: environmentCount)
        var diagnostics: [VivoContextDiagnostic] = []
        try applyParameterOverrides(
            context.parameterOverrides,
            metadata: parameters,
            indices: parameterIndex,
            environmentCount: environmentCount,
            values: &parameterValues
        )

        var transport = [VivoSpeciesTransportABI](repeating: .init(), count: species.count)
        try applyTransport(
            context.transport,
            speciesIndices: speciesIndex,
            table: &transport,
            diagnostics: &diagnostics
        )

        var coupling: [VivoCouplingUpdate] = []
        coupling.reserveCapacity(min(options.maximumGeneratedCouplingUpdates, species.count * max(environmentCount, 1)))
        try appendInitialSignals(
            context.initialSignals,
            speciesMetadata: species,
            speciesIndices: speciesIndex,
            configuration: configuration,
            maximumUpdates: options.maximumGeneratedCouplingUpdates,
            destination: &coupling
        )
        try appendHostChannels(
            context.hostChannels,
            speciesMetadata: species,
            speciesIndices: speciesIndex,
            configuration: configuration,
            strict: options.strictUnknownHostChannels,
            maximumUpdates: options.maximumGeneratedCouplingUpdates,
            diagnostics: &diagnostics,
            destination: &coupling
        )

        if context.transport.isEmpty,
           configuration.fidelity.rawValue >= VivoFidelity.spatial.rawValue {
            diagnostics.append(.init(
                severity: .warning,
                code: "NVCTX001",
                path: "$.transport",
                message: "Spatial runtime has no species transport definitions; diffusion, permeability, and extracellular decay default to zero."
            ))
        }
        if context.parameterOverrides.isEmpty {
            diagnostics.append(.init(
                severity: .note,
                code: "NVCTX002",
                path: "$.parameterOverrides",
                message: "No host-specific parameter overrides were supplied; ProgramPack defaults are retained."
            ))
        }

        return PreparedVivoHostContext(
            contextFingerprint: contextFingerprint,
            programFingerprint: pack.header.contentFingerprint,
            environmentCount: configuration.environmentCount,
            laneCount: configuration.laneCount,
            parameterValues: parameterValues,
            transport: transport,
            initialCoupling: coupling,
            diagnostics: diagnostics
        )
    }

    private func applyParameterOverrides(
        _ overrides: [VivoParameterOverride],
        metadata: [VivoProgramPack.ParameterMetadata],
        indices: [String: Int],
        environmentCount: Int,
        values: inout [Float]
    ) throws {
        for override in overrides {
            guard let index = indices[override.parameter] else {
                throw VivoArtifactValidationError.unresolved(
                    "parameter override '\(override.parameter)' is absent from the ProgramPack"
                )
            }
            let parameter = metadata[index]
            for environment in 0..<environmentCount {
                let source = override.values.count == 1 ? override.values[0] : override.values[environment]
                let converted = try units.convert(source, to: parameter.unit)
                guard converted >= parameter.minimum, converted <= parameter.maximum else {
                    throw VivoArtifactValidationError.invalid(
                        "parameter \(parameter.identifier) value \(converted) \(parameter.unit) is outside [\(parameter.minimum), \(parameter.maximum)]"
                    )
                }
                let floatValue = try checkedFloat(
                    converted,
                    label: "parameterOverrides.\(parameter.identifier)[\(environment)]"
                )
                values[index * environmentCount + environment] = floatValue
            }
        }
    }

    private func applyTransport(
        _ definitions: [VivoSpeciesTransport],
        speciesIndices: [String: Int],
        table: inout [VivoSpeciesTransportABI],
        diagnostics: inout [VivoContextDiagnostic]
    ) throws {
        for definition in definitions {
            guard let index = speciesIndices[definition.species] else {
                throw VivoArtifactValidationError.unresolved(
                    "transport species '\(definition.species)' is absent from the ProgramPack"
                )
            }
            let diffusion = try checkedFloat(
                units.convert(definition.diffusion, to: "m2/s"),
                label: "transport.\(definition.species).diffusion"
            )
            let permeability: Float
            if let source = definition.membranePermeability {
                permeability = try checkedFloat(
                    units.convert(source, to: "m/s"),
                    label: "transport.\(definition.species).membranePermeability"
                )
            } else {
                permeability = 0
            }
            let decay: Float
            if let source = definition.extracellularDecayRate {
                decay = try checkedFloat(
                    units.convert(source, to: "1/s"),
                    label: "transport.\(definition.species).extracellularDecayRate"
                )
            } else {
                decay = 0
            }
            guard diffusion >= 0, permeability >= 0, decay >= 0 else {
                throw VivoArtifactValidationError.invalid(
                    "transport coefficients for \(definition.species) must be non-negative"
                )
            }
            table[index] = VivoSpeciesTransportABI(
                diffusion: diffusion,
                membranePermeability: permeability,
                decayRate: decay,
                flags: 0
            )
            if definition.evidence.isEmpty {
                diagnostics.append(.init(
                    severity: .warning,
                    code: "NVCTX003",
                    path: "$.transport.\(definition.species)",
                    message: "Transport coefficients have no attached evidence records."
                ))
            }
        }
    }

    private func appendInitialSignals(
        _ signals: [VivoInitialSignal],
        speciesMetadata: [VivoProgramPack.SpeciesMetadata],
        speciesIndices: [String: Int],
        configuration: VivoRuntimeConfiguration,
        maximumUpdates: Int,
        destination: inout [VivoCouplingUpdate]
    ) throws {
        let environmentCount = Int(configuration.environmentCount)
        let voxelCount = Int(configuration.voxelCount)
        for signal in signals {
            guard let index = speciesIndices[signal.species] else {
                throw VivoArtifactValidationError.unresolved(
                    "initial signal species '\(signal.species)' is absent from the ProgramPack"
                )
            }
            let metadata = speciesMetadata[index]
            guard metadata.isExternallyOwned || metadata.isInput else {
                throw VivoArtifactValidationError.incompatible(
                    "initial signal \(signal.species) targets an internally owned species"
                )
            }
            guard signal.laneValues.count == 1 || signal.laneValues.count == environmentCount else {
                throw VivoArtifactValidationError.invalid(
                    "initial signal \(signal.species) requires one value or one value per environment"
                )
            }
            for environment in 0..<environmentCount {
                let source = signal.laneValues.count == 1 ? signal.laneValues[0] : signal.laneValues[environment]
                let converted = try units.convert(source, to: metadata.unit)
                let value = try checkedSpeciesValue(converted, metadata: metadata, path: "initialSignals.\(signal.species)")
                try appendAcrossVoxels(
                    speciesIndex: UInt32(index),
                    environment: environment,
                    voxelCount: voxelCount,
                    mode: signal.mode,
                    value: value,
                    maximumUpdates: maximumUpdates,
                    destination: &destination
                )
            }
        }
    }

    private func appendHostChannels(
        _ channels: [String: VivoQuantity],
        speciesMetadata: [VivoProgramPack.SpeciesMetadata],
        speciesIndices: [String: Int],
        configuration: VivoRuntimeConfiguration,
        strict: Bool,
        maximumUpdates: Int,
        diagnostics: inout [VivoContextDiagnostic],
        destination: inout [VivoCouplingUpdate]
    ) throws {
        let voxelCount = Int(configuration.voxelCount)
        let environmentCount = Int(configuration.environmentCount)
        for name in channels.keys.sorted() {
            guard let quantity = channels[name] else { continue }
            guard let index = speciesIndices[name] else {
                if strict {
                    throw VivoArtifactValidationError.unresolved(
                        "host channel '\(name)' is absent from the ProgramPack"
                    )
                }
                diagnostics.append(.init(
                    severity: .warning,
                    code: "NVCTX004",
                    path: "$.hostChannels.\(name)",
                    message: "Host channel was not consumed because no ProgramPack species has this identifier."
                ))
                continue
            }
            let metadata = speciesMetadata[index]
            guard metadata.isExternallyOwned || metadata.isInput else {
                throw VivoArtifactValidationError.incompatible(
                    "host channel \(name) targets an internally owned species"
                )
            }
            let converted = try units.convert(quantity, to: metadata.unit)
            let value = try checkedSpeciesValue(converted, metadata: metadata, path: "hostChannels.\(name)")
            for environment in 0..<environmentCount {
                try appendAcrossVoxels(
                    speciesIndex: UInt32(index),
                    environment: environment,
                    voxelCount: voxelCount,
                    mode: .replace,
                    value: value,
                    maximumUpdates: maximumUpdates,
                    destination: &destination
                )
            }
        }
    }

    private func appendAcrossVoxels(
        speciesIndex: UInt32,
        environment: Int,
        voxelCount: Int,
        mode: VivoCouplingMode,
        value: Float,
        maximumUpdates: Int,
        destination: inout [VivoCouplingUpdate]
    ) throws {
        guard destination.count <= maximumUpdates - voxelCount else {
            throw VivoArtifactValidationError.invalid(
                "compiled host context exceeds maximumGeneratedCouplingUpdates (\(maximumUpdates))"
            )
        }
        let laneBase = try checkedProduct(environment, voxelCount, label: "lane base")
        for voxel in 0..<voxelCount {
            let lane = laneBase + voxel
            guard lane <= Int(UInt32.max) else {
                throw VivoArtifactValidationError.invalid("compiled lane index exceeds UInt32")
            }
            destination.append(.init(
                speciesIndex: speciesIndex,
                laneIndex: UInt32(lane),
                mode: mode,
                value: value
            ))
        }
    }

    private func checkedSpeciesValue(
        _ value: Double,
        metadata: VivoProgramPack.SpeciesMetadata,
        path: String
    ) throws -> Float {
        guard value >= Double(metadata.minimum), value <= Double(metadata.maximum) else {
            throw VivoArtifactValidationError.invalid(
                "\(path) value \(value) \(metadata.unit) is outside [\(metadata.minimum), \(metadata.maximum)]"
            )
        }
        let converted = try checkedFloat(value, label: path)
        if metadata.isCountValued, abs(converted.rounded() - converted) > 1e-4 {
            throw VivoArtifactValidationError.invalid("\(path) targets a count-valued species but is not integral")
        }
        return metadata.isCountValued ? converted.rounded() : converted
    }

    private func validateGeometry(
        _ geometry: VivoTissueContext.Geometry,
        configuration: VivoRuntimeConfiguration,
        options: Options
    ) throws {
        guard options.requireSpatialGeometryAgreement else { return }
        switch configuration.fidelity {
        case .logic, .deterministic, .stochastic:
            guard geometry.kind == .wellMixed || configuration.spatialGrid == nil else {
                throw VivoArtifactValidationError.incompatible(
                    "non-spatial runtime cannot consume a voxel or external geometry"
                )
            }
        case .spatial, .tissue:
            guard let runtimeGrid = configuration.spatialGrid else {
                throw VivoArtifactValidationError.incompatible("spatial runtime has no grid")
            }
            if geometry.kind == .voxelGrid {
                guard geometry.dimensions == [runtimeGrid.width, runtimeGrid.height, runtimeGrid.depth] else {
                    throw VivoArtifactValidationError.incompatible(
                        "host-context voxel dimensions do not match runtime configuration"
                    )
                }
                let expectedUnits = ["m", "m", "m"]
                let runtimeSpacing = [runtimeGrid.spacingX, runtimeGrid.spacingY, runtimeGrid.spacingZ]
                for axis in 0..<3 {
                    let meters = try units.convert(geometry.spacing[axis], to: expectedUnits[axis])
                    let expected = Double(runtimeSpacing[axis])
                    let tolerance = max(abs(expected) * 1e-5, 1e-12)
                    guard abs(meters - expected) <= tolerance else {
                        throw VivoArtifactValidationError.incompatible(
                            "host-context voxel spacing axis \(axis) does not match runtime configuration"
                        )
                    }
                }
            } else if geometry.kind != .externalNumanXDomain && geometry.kind != .externalNumiTissueDomain && geometry.kind != .mesh {
                throw VivoArtifactValidationError.incompatible(
                    "F3/F4 runtime requires voxelGrid, mesh, or an external simulation domain"
                )
            }
        }
    }

    private func uniqueIndex(_ identifiers: [String], label: String) throws -> [String: Int] {
        var result: [String: Int] = [:]
        result.reserveCapacity(identifiers.count)
        for (index, identifier) in identifiers.enumerated() {
            guard result.updateValue(index, forKey: identifier) == nil else {
                throw VivoArtifactValidationError.invalid("duplicate \(label) identifier '\(identifier)'")
            }
        }
        return result
    }

    private func checkedFloat(_ value: Double, label: String) throws -> Float {
        guard value.isFinite, abs(value) <= Double(Float.greatestFiniteMagnitude) else {
            throw VivoArtifactValidationError.invalid("\(label) cannot be represented as finite FP32")
        }
        return Float(value)
    }

    private func checkedProduct(_ lhs: Int, _ rhs: Int, label: String) throws -> Int {
        let result = lhs.multipliedReportingOverflow(by: rhs)
        guard !result.overflow else {
            throw VivoArtifactValidationError.invalid("\(label) integer overflow")
        }
        return result.partialValue
    }
}
