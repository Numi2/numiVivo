import Foundation

public struct VivoInterval: Codable, Sendable, Equatable, Hashable {
    public var minimum: Double
    public var maximum: Double

    public init(minimum: Double, maximum: Double) {
        self.minimum = minimum
        self.maximum = maximum
    }

    public func validate(label: String) throws {
        guard minimum.isFinite, maximum.isFinite, minimum <= maximum else {
            throw VivoArtifactValidationError.invalid("\(label) interval must be finite and ordered")
        }
    }
}

public struct VivoQuantity: Codable, Sendable, Equatable, Hashable {
    public var value: Double
    public var unit: String
    public var plausibleRange: VivoInterval?
    public var evidence: VivoEvidenceReference?

    public init(
        value: Double,
        unit: String,
        plausibleRange: VivoInterval? = nil,
        evidence: VivoEvidenceReference? = nil
    ) {
        self.value = value
        self.unit = unit
        self.plausibleRange = plausibleRange
        self.evidence = evidence
    }

    public func validate(label: String, nonnegative: Bool = false) throws {
        guard value.isFinite, !unit.isEmpty else {
            throw VivoArtifactValidationError.invalid("\(label) must contain a finite value and explicit unit")
        }
        if nonnegative, value < 0 {
            throw VivoArtifactValidationError.invalid("\(label) cannot be negative")
        }
        try plausibleRange?.validate(label: label)
        if let plausibleRange,
           !(plausibleRange.minimum...plausibleRange.maximum).contains(value) {
            throw VivoArtifactValidationError.invalid("\(label) value lies outside its declared plausible range")
        }
    }
}

public struct VivoOrganismContext: Codable, Sendable, Equatable {
    public var taxonomyIdentifier: String
    public var speciesName: String
    public var strainOrPopulation: String?
    public var developmentalStage: String?
    public var age: VivoQuantity?
    public var geneticBackground: [String: String]
    public var physiologicalState: [String: VivoQuantity]
    public var annotations: [String: String]

    public init(
        taxonomyIdentifier: String,
        speciesName: String,
        strainOrPopulation: String? = nil,
        developmentalStage: String? = nil,
        age: VivoQuantity? = nil,
        geneticBackground: [String: String] = [:],
        physiologicalState: [String: VivoQuantity] = [:],
        annotations: [String: String] = [:]
    ) {
        self.taxonomyIdentifier = taxonomyIdentifier
        self.speciesName = speciesName
        self.strainOrPopulation = strainOrPopulation
        self.developmentalStage = developmentalStage
        self.age = age
        self.geneticBackground = geneticBackground
        self.physiologicalState = physiologicalState
        self.annotations = annotations
    }

    public func validate() throws {
        guard !taxonomyIdentifier.isEmpty, !speciesName.isEmpty else {
            throw VivoArtifactValidationError.invalid("organism taxonomyIdentifier and speciesName are required")
        }
        try age?.validate(label: "organism.age", nonnegative: true)
        for (key, value) in physiologicalState {
            try value.validate(label: "organism.physiologicalState.\(key)")
        }
    }
}

public struct VivoTissueContext: Codable, Sendable, Equatable {
    public struct Geometry: Codable, Sendable, Equatable {
        public enum Kind: String, Codable, Sendable {
            case wellMixed
            case voxelGrid
            case mesh
            case externalNumanXDomain
            case externalNumiTissueDomain
        }

        public var kind: Kind
        public var dimensions: [UInt32]
        public var spacing: [VivoQuantity]
        public var externalDomainIdentifier: String?

        public init(
            kind: Kind,
            dimensions: [UInt32] = [],
            spacing: [VivoQuantity] = [],
            externalDomainIdentifier: String? = nil
        ) {
            self.kind = kind
            self.dimensions = dimensions
            self.spacing = spacing
            self.externalDomainIdentifier = externalDomainIdentifier
        }

        public func validate() throws {
            switch kind {
            case .wellMixed:
                break
            case .voxelGrid:
                guard dimensions.count == 3, dimensions.allSatisfy({ $0 > 0 }), spacing.count == 3 else {
                    throw VivoArtifactValidationError.invalid("voxelGrid geometry requires three positive dimensions and three spacing values")
                }
                for (index, value) in spacing.enumerated() {
                    try value.validate(label: "tissue.geometry.spacing[\(index)]", nonnegative: true)
                    guard value.value > 0 else {
                        throw VivoArtifactValidationError.invalid("voxel spacing must be positive")
                    }
                }
            case .mesh, .externalNumanXDomain, .externalNumiTissueDomain:
                guard let externalDomainIdentifier, !externalDomainIdentifier.isEmpty else {
                    throw VivoArtifactValidationError.invalid("external or mesh geometry requires an externalDomainIdentifier")
                }
            }
        }
    }

    public var tissueIdentifier: String
    public var ontologyTerm: String?
    public var compartment: String
    public var geometry: Geometry
    public var temperature: VivoQuantity
    public var pH: VivoQuantity
    public var oxygen: VivoQuantity?
    public var perfusion: VivoQuantity?
    public var extracellularVolumeFraction: VivoQuantity?
    public var extracellularMatrix: [String: VivoQuantity]
    public var fields: [String: VivoQuantity]

    public init(
        tissueIdentifier: String,
        ontologyTerm: String? = nil,
        compartment: String,
        geometry: Geometry,
        temperature: VivoQuantity,
        pH: VivoQuantity,
        oxygen: VivoQuantity? = nil,
        perfusion: VivoQuantity? = nil,
        extracellularVolumeFraction: VivoQuantity? = nil,
        extracellularMatrix: [String: VivoQuantity] = [:],
        fields: [String: VivoQuantity] = [:]
    ) {
        self.tissueIdentifier = tissueIdentifier
        self.ontologyTerm = ontologyTerm
        self.compartment = compartment
        self.geometry = geometry
        self.temperature = temperature
        self.pH = pH
        self.oxygen = oxygen
        self.perfusion = perfusion
        self.extracellularVolumeFraction = extracellularVolumeFraction
        self.extracellularMatrix = extracellularMatrix
        self.fields = fields
    }

    public func validate() throws {
        guard !tissueIdentifier.isEmpty, !compartment.isEmpty else {
            throw VivoArtifactValidationError.invalid("tissueIdentifier and compartment are required")
        }
        try geometry.validate()
        try temperature.validate(label: "tissue.temperature")
        try pH.validate(label: "tissue.pH")
        try oxygen?.validate(label: "tissue.oxygen", nonnegative: true)
        try perfusion?.validate(label: "tissue.perfusion", nonnegative: true)
        try extracellularVolumeFraction?.validate(label: "tissue.extracellularVolumeFraction", nonnegative: true)
        if let fraction = extracellularVolumeFraction, fraction.value > 1, fraction.unit == "fraction" {
            throw VivoArtifactValidationError.invalid("extracellular volume fraction cannot exceed one")
        }
        for (key, value) in extracellularMatrix {
            try value.validate(label: "tissue.extracellularMatrix.\(key)")
        }
        for (key, value) in fields {
            try value.validate(label: "tissue.fields.\(key)")
        }
    }
}

public struct VivoCellPopulationContext: Codable, Sendable, Equatable {
    public var identifier: String
    public var cellType: String
    public var ontologyTerm: String?
    public var compartment: String
    public var abundance: VivoQuantity
    public var meanCellVolume: VivoQuantity?
    public var genotype: [String: String]
    public var baselineState: [String: VivoQuantity]
    public var receptorAbundance: [String: VivoQuantity]
    public var fractions: [String: Double]
    public var evidence: [VivoEvidenceReference]

    public init(
        identifier: String,
        cellType: String,
        ontologyTerm: String? = nil,
        compartment: String,
        abundance: VivoQuantity,
        meanCellVolume: VivoQuantity? = nil,
        genotype: [String: String] = [:],
        baselineState: [String: VivoQuantity] = [:],
        receptorAbundance: [String: VivoQuantity] = [:],
        fractions: [String: Double] = [:],
        evidence: [VivoEvidenceReference] = []
    ) {
        self.identifier = identifier
        self.cellType = cellType
        self.ontologyTerm = ontologyTerm
        self.compartment = compartment
        self.abundance = abundance
        self.meanCellVolume = meanCellVolume
        self.genotype = genotype
        self.baselineState = baselineState
        self.receptorAbundance = receptorAbundance
        self.fractions = fractions
        self.evidence = evidence
    }

    public func validate() throws {
        guard !identifier.isEmpty, !cellType.isEmpty, !compartment.isEmpty else {
            throw VivoArtifactValidationError.invalid("cell population identifier, cellType, and compartment are required")
        }
        try abundance.validate(label: "cellPopulations.\(identifier).abundance", nonnegative: true)
        try meanCellVolume?.validate(label: "cellPopulations.\(identifier).meanCellVolume", nonnegative: true)
        for (key, value) in baselineState {
            try value.validate(label: "cellPopulations.\(identifier).baselineState.\(key)")
        }
        for (key, value) in receptorAbundance {
            try value.validate(label: "cellPopulations.\(identifier).receptorAbundance.\(key)", nonnegative: true)
        }
        guard fractions.values.allSatisfy({ $0.isFinite && $0 >= 0 && $0 <= 1 }) else {
            throw VivoArtifactValidationError.invalid("cell population fractions must lie in [0,1]")
        }
    }
}

public struct VivoDeliveryContext: Codable, Sendable, Equatable {
    public enum Mode: String, Codable, Sendable, CaseIterable {
        case retrievableImplant
        case encapsulatedCells
        case exVivoEngineeredCells
        case transientRNA
        case nonIntegratingVector
        case integratingVector
        case permanentGenomeEdit
        case engineeredMicrobe
    }

    public enum Locality: String, Codable, Sendable {
        case local
        case regional
        case systemic
    }

    public var mode: Mode
    public var locality: Locality
    public var route: String
    public var administeredAmount: VivoQuantity
    public var distributionVolume: VivoQuantity?
    public var clearanceHalfLife: VivoQuantity?
    public var retrievalSupported: Bool
    public var externalShutdownSupported: Bool
    public var physicalContainment: [String]
    public var molecularContainment: [String]
    public var assumptions: [String]

    public init(
        mode: Mode,
        locality: Locality,
        route: String,
        administeredAmount: VivoQuantity,
        distributionVolume: VivoQuantity? = nil,
        clearanceHalfLife: VivoQuantity? = nil,
        retrievalSupported: Bool,
        externalShutdownSupported: Bool,
        physicalContainment: [String] = [],
        molecularContainment: [String] = [],
        assumptions: [String] = []
    ) {
        self.mode = mode
        self.locality = locality
        self.route = route
        self.administeredAmount = administeredAmount
        self.distributionVolume = distributionVolume
        self.clearanceHalfLife = clearanceHalfLife
        self.retrievalSupported = retrievalSupported
        self.externalShutdownSupported = externalShutdownSupported
        self.physicalContainment = physicalContainment
        self.molecularContainment = molecularContainment
        self.assumptions = assumptions
    }

    public var isPersistent: Bool {
        switch mode {
        case .integratingVector, .permanentGenomeEdit, .engineeredMicrobe:
            true
        default:
            false
        }
    }

    public func validate() throws {
        guard !route.isEmpty else {
            throw VivoArtifactValidationError.invalid("delivery route is required")
        }
        try administeredAmount.validate(label: "delivery.administeredAmount", nonnegative: true)
        try distributionVolume?.validate(label: "delivery.distributionVolume", nonnegative: true)
        try clearanceHalfLife?.validate(label: "delivery.clearanceHalfLife", nonnegative: true)
        if isPersistent && !externalShutdownSupported && molecularContainment.isEmpty {
            throw VivoArtifactValidationError.invalid(
                "persistent delivery requires external shutdown support or an explicit molecular containment mechanism"
            )
        }
        if mode == .retrievableImplant && !retrievalSupported {
            throw VivoArtifactValidationError.invalid("retrievableImplant mode must declare retrievalSupported")
        }
    }
}

public struct VivoParameterOverride: Codable, Sendable, Equatable {
    public var parameter: String
    public var values: [VivoQuantity]
    public var evidence: VivoEvidenceReference

    public init(parameter: String, values: [VivoQuantity], evidence: VivoEvidenceReference) {
        self.parameter = parameter
        self.values = values
        self.evidence = evidence
    }

    public func validate(environmentCount: Int) throws {
        guard !parameter.isEmpty else {
            throw VivoArtifactValidationError.invalid("parameter override name cannot be empty")
        }
        guard values.count == 1 || values.count == environmentCount else {
            throw VivoArtifactValidationError.invalid(
                "parameter override \(parameter) must provide one value or one value per environment"
            )
        }
        for (index, value) in values.enumerated() {
            try value.validate(label: "parameterOverrides.\(parameter)[\(index)]")
        }
    }
}

public struct VivoSpeciesTransport: Codable, Sendable, Equatable {
    public var species: String
    public var diffusion: VivoQuantity
    public var membranePermeability: VivoQuantity?
    public var extracellularDecayRate: VivoQuantity?
    public var evidence: [VivoEvidenceReference]

    public init(
        species: String,
        diffusion: VivoQuantity,
        membranePermeability: VivoQuantity? = nil,
        extracellularDecayRate: VivoQuantity? = nil,
        evidence: [VivoEvidenceReference] = []
    ) {
        self.species = species
        self.diffusion = diffusion
        self.membranePermeability = membranePermeability
        self.extracellularDecayRate = extracellularDecayRate
        self.evidence = evidence
    }

    public func validate() throws {
        guard !species.isEmpty else {
            throw VivoArtifactValidationError.invalid("species transport identifier cannot be empty")
        }
        try diffusion.validate(label: "transport.\(species).diffusion", nonnegative: true)
        try membranePermeability?.validate(label: "transport.\(species).membranePermeability", nonnegative: true)
        try extracellularDecayRate?.validate(label: "transport.\(species).extracellularDecayRate", nonnegative: true)
    }
}

public struct VivoInitialSignal: Codable, Sendable, Equatable {
    public var species: String
    public var laneValues: [VivoQuantity]
    public var mode: VivoCouplingMode
    public var evidence: VivoEvidenceReference?

    public init(
        species: String,
        laneValues: [VivoQuantity],
        mode: VivoCouplingMode = .replace,
        evidence: VivoEvidenceReference? = nil
    ) {
        self.species = species
        self.laneValues = laneValues
        self.mode = mode
        self.evidence = evidence
    }
}

public struct VivoHostContextPack: Codable, Sendable, Equatable {
    public var organism: VivoOrganismContext
    public var tissue: VivoTissueContext
    public var cellPopulations: [VivoCellPopulationContext]
    public var delivery: VivoDeliveryContext
    public var parameterOverrides: [VivoParameterOverride]
    public var transport: [VivoSpeciesTransport]
    public var initialSignals: [VivoInitialSignal]
    public var hostChannels: [String: VivoQuantity]
    public var uncertaintyTags: [String: VivoInterval]

    public init(
        organism: VivoOrganismContext,
        tissue: VivoTissueContext,
        cellPopulations: [VivoCellPopulationContext],
        delivery: VivoDeliveryContext,
        parameterOverrides: [VivoParameterOverride] = [],
        transport: [VivoSpeciesTransport] = [],
        initialSignals: [VivoInitialSignal] = [],
        hostChannels: [String: VivoQuantity] = [:],
        uncertaintyTags: [String: VivoInterval] = [:]
    ) {
        self.organism = organism
        self.tissue = tissue
        self.cellPopulations = cellPopulations
        self.delivery = delivery
        self.parameterOverrides = parameterOverrides
        self.transport = transport
        self.initialSignals = initialSignals
        self.hostChannels = hostChannels
        self.uncertaintyTags = uncertaintyTags
    }

    public func validate(environmentCount: Int) throws {
        guard environmentCount > 0 else {
            throw VivoArtifactValidationError.invalid("environmentCount must be positive")
        }
        try organism.validate()
        try tissue.validate()
        try delivery.validate()
        guard !cellPopulations.isEmpty else {
            throw VivoArtifactValidationError.invalid("host context requires at least one cell population")
        }
        var cellIdentifiers = Set<String>()
        for population in cellPopulations {
            try population.validate()
            guard cellIdentifiers.insert(population.identifier).inserted else {
                throw VivoArtifactValidationError.invalid("duplicate cell population \(population.identifier)")
            }
        }
        var parameterNames = Set<String>()
        for override in parameterOverrides {
            try override.validate(environmentCount: environmentCount)
            guard parameterNames.insert(override.parameter).inserted else {
                throw VivoArtifactValidationError.invalid("duplicate parameter override \(override.parameter)")
            }
        }
        var transportSpecies = Set<String>()
        for value in transport {
            try value.validate()
            guard transportSpecies.insert(value.species).inserted else {
                throw VivoArtifactValidationError.invalid("duplicate transport definition for \(value.species)")
            }
        }
        for signal in initialSignals {
            guard !signal.species.isEmpty,
                  signal.laneValues.count == 1 || signal.laneValues.count == environmentCount else {
                throw VivoArtifactValidationError.invalid(
                    "initial signal \(signal.species) must provide one value or one value per environment"
                )
            }
            for (index, value) in signal.laneValues.enumerated() {
                try value.validate(label: "initialSignals.\(signal.species)[\(index)]")
            }
        }
        for (key, value) in hostChannels {
            try value.validate(label: "hostChannels.\(key)")
        }
        for (key, value) in uncertaintyTags {
            try value.validate(label: "uncertaintyTags.\(key)")
        }
    }
}
