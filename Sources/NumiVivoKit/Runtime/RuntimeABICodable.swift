import Foundation

extension VivoSpeciesTransportABI: Codable, Equatable {
    private enum CodingKeys: String, CodingKey {
        case diffusion
        case membranePermeability
        case decayRate
        case flags
    }

    public static func == (lhs: VivoSpeciesTransportABI, rhs: VivoSpeciesTransportABI) -> Bool {
        lhs.diffusion == rhs.diffusion &&
        lhs.membranePermeability == rhs.membranePermeability &&
        lhs.decayRate == rhs.decayRate &&
        lhs.flags == rhs.flags
    }

    public init(from decoder: Decoder) throws {
        let container = try decoder.container(keyedBy: CodingKeys.self)
        self.init(
            diffusion: try container.decode(Float.self, forKey: .diffusion),
            membranePermeability: try container.decode(Float.self, forKey: .membranePermeability),
            decayRate: try container.decode(Float.self, forKey: .decayRate),
            flags: try container.decode(UInt32.self, forKey: .flags)
        )
    }

    public func encode(to encoder: Encoder) throws {
        var container = encoder.container(keyedBy: CodingKeys.self)
        try container.encode(diffusion, forKey: .diffusion)
        try container.encode(membranePermeability, forKey: .membranePermeability)
        try container.encode(decayRate, forKey: .decayRate)
        try container.encode(flags, forKey: .flags)
    }
}

extension VivoCouplingUpdateABI: Codable, Equatable {
    private enum CodingKeys: String, CodingKey {
        case speciesIndex
        case laneIndex
        case mode
        case value
    }

    public static func == (lhs: VivoCouplingUpdateABI, rhs: VivoCouplingUpdateABI) -> Bool {
        lhs.speciesIndex == rhs.speciesIndex &&
        lhs.laneIndex == rhs.laneIndex &&
        lhs.mode == rhs.mode &&
        lhs.value == rhs.value
    }

    public init(from decoder: Decoder) throws {
        let container = try decoder.container(keyedBy: CodingKeys.self)
        self.init(
            speciesIndex: try container.decode(UInt32.self, forKey: .speciesIndex),
            laneIndex: try container.decode(UInt32.self, forKey: .laneIndex),
            mode: try container.decode(UInt32.self, forKey: .mode),
            value: try container.decode(Float.self, forKey: .value)
        )
    }

    public func encode(to encoder: Encoder) throws {
        var container = encoder.container(keyedBy: CodingKeys.self)
        try container.encode(speciesIndex, forKey: .speciesIndex)
        try container.encode(laneIndex, forKey: .laneIndex)
        try container.encode(mode, forKey: .mode)
        try container.encode(value, forKey: .value)
    }
}

extension VivoPublicationRequestABI: Codable, Equatable {
    private enum CodingKeys: String, CodingKey {
        case speciesIndex
        case laneIndex
        case outputIndex
        case flags
    }

    public static func == (lhs: VivoPublicationRequestABI, rhs: VivoPublicationRequestABI) -> Bool {
        lhs.speciesIndex == rhs.speciesIndex &&
        lhs.laneIndex == rhs.laneIndex &&
        lhs.outputIndex == rhs.outputIndex &&
        lhs.flags == rhs.flags
    }

    public init(from decoder: Decoder) throws {
        let container = try decoder.container(keyedBy: CodingKeys.self)
        self.init(
            speciesIndex: try container.decode(UInt32.self, forKey: .speciesIndex),
            laneIndex: try container.decode(UInt32.self, forKey: .laneIndex),
            outputIndex: try container.decode(UInt32.self, forKey: .outputIndex),
            flags: try container.decode(UInt32.self, forKey: .flags)
        )
    }

    public func encode(to encoder: Encoder) throws {
        var container = encoder.container(keyedBy: CodingKeys.self)
        try container.encode(speciesIndex, forKey: .speciesIndex)
        try container.encode(laneIndex, forKey: .laneIndex)
        try container.encode(outputIndex, forKey: .outputIndex)
        try container.encode(flags, forKey: .flags)
    }
}
